"""
Manages material properties and preparation for OpenTrons sessions.
"""

import logging
import traceback
import math
import pandas as pd
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from rdkit import Chem
from rdkit.Chem import Descriptors

from ....models import AddAction, Plate, Well, CompoundOrder, SolventPrep, Product
from ....chem_utils import canon_smiles, strip_salts, are_equivalent_structures
from ....db_utils import check_previous_reaction_products

logger = logging.getLogger(__name__)


class MaterialCalculator:
    """Handles calculations related to material properties."""

    def recalculate_amount_for_salt_form(
        self, original_smiles, actual_smiles, original_volume
    ):
        """
        Recalculates the amount needed when switching to a different salt form.

        Parameters
        ----------
        original_smiles: str
            Original SMILES string
        actual_smiles: str
            Actual SMILES string found in custom plate
        original_volume: float
            Original volume calculated

        Returns
        -------
        new_volume: float
            Recalculated volume based on molecular weight ratio
        """
        try:
            # Get molecular weights
            mol_original = Chem.MolFromSmiles(original_smiles)
            mol_actual = Chem.MolFromSmiles(actual_smiles)

            if not mol_original or not mol_actual:
                logger.warning(
                    f"Could not parse SMILES for recalculation: {original_smiles} or {actual_smiles}"
                )
                return original_volume

            mw_original = Descriptors.MolWt(mol_original)
            mw_actual = Descriptors.MolWt(mol_actual)

            # Calculate volume ratio based on molecular weight
            # The ratio is determined by the number of moles needed, which is constant
            # If MW increases, volume must increase proportionally
            volume_ratio = mw_actual / mw_original
            new_volume = original_volume * volume_ratio

            logger.info(
                f"Recalculated volume for salt form. "
                f"Original MW: {mw_original:.2f}, Actual MW: {mw_actual:.2f}, "
                f"Ratio: {volume_ratio:.3f}, "
                f"Original volume: {original_volume:.2f}µL, New volume: {new_volume:.2f}µL"
            )

            return new_volume

        except Exception as e:
            logger.error(f"Error recalculating amount for salt form: {str(e)}")
            return original_volume


class SaltMatchingService:
    """
    Handles salt matching and material identification.
    """

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def find_matching_materials(self, canonical_smiles, potential_matches, volume):
        """Finds matching materials considering salt variations."""
        try:
            if not potential_matches:
                self.logger.debug("No potential matches found")
                return (False, [], None, volume, None)  # Added None for actual_smiles

            # Clean the query SMILES
            cleaned_smiles = strip_salts(canonical_smiles)

            matching_wells = []
            containing_plate = None
            total_available_volume = 0
            actual_smiles = None  # Track the actual SMILES found
            
            # Process the potential matches
            for well in potential_matches:
                # Canonicalize and strip salts from the well's SMILES
                well_canonical = canon_smiles(well.smiles)

               # Use our comprehensive structure comparison
                if are_equivalent_structures(cleaned_smiles, well_canonical):
                    matching_wells.append(well)
                    if containing_plate is None and hasattr(well, "plate_id"):
                        containing_plate = well.plate_id
                    if well.volume is not None:
                        total_available_volume += well.volume

                    # Store the actual SMILES if not set yet (prioritize first match)
                    if actual_smiles is None:
                        actual_smiles = well.smiles
                        self.logger.info(
                            f"Using SMILES from matched well: {actual_smiles}"
                        )
            # Calculate remaining volume needed
            remaining_volume_needed = max(0, volume - total_available_volume)

            # Check if we have enough volume
            if remaining_volume_needed <= 0:
                self.logger.info(
                    f"Found enough material for structure: {cleaned_smiles} "
                    f"Required: {volume}µL, Available: {total_available_volume}µL"
                )
                return (True, matching_wells, containing_plate, 0, actual_smiles)
            elif matching_wells:
                self.logger.info(
                    f"Found material matching structure {cleaned_smiles} but insufficient volume. "
                    f"Required: {volume}µL, Available: {total_available_volume}µL, "
                    f"Still need: {remaining_volume_needed}µL"
                )
                return (
                    False,
                    matching_wells,
                    containing_plate,
                    remaining_volume_needed,
                    actual_smiles,
                )
            else:
                self.logger.info(
                    f"No existing material found matching structure: {cleaned_smiles}"
                )
                return (False, [], None, volume, None)

        except Exception as e:
            self.logger.error(f"Error in salt matching: {str(e)}")
            self.logger.error(traceback.format_exc())
            return (False, [], None, volume, None)


class MaterialManager:
    """
    Manages material properties, concentrations, and preparation.
    """

    def __init__(self, session):
        """
        Initialize with a reference to the parent session.

        Parameters
        ----------
        session: BaseSession
            The parent session object
        """
        self.session = session
        self.calculator = MaterialCalculator()
        self.salt_matcher = SaltMatchingService()

    def get_add_actions_material_dataframe(self, product_exists: bool) -> pd.DataFrame:
        """
        Aggregates all add actions materials and sums up volume requirements using solvent type and
        concentration. Checks existing materials and only requests additional volume if needed.

        Parameters
        ----------
        product_exists: bool
            Set to true to check if the add action material needed is a product from
            a previous reaction

        Returns
        -------
        materials_df: pd.DataFrame
            The add action material as dataframe grouping materials by SMILES, concentration
            and solvent, with volumes adjusted based on existing materials
        """
        try:
            # Get add actions dataframe from the data manager
            if hasattr(self.session, "addactionsdf"):
                addactionsdf = self.session.addactionsdf
            else:
                # If add actions DataFrame doesn't exist yet, create it
                addactionqueryset = self.session.data_manager.get_add_action_query_set(
                    reaction_ids=self.session.reaction_ids,
                    actionsession_ids=self.session.actionsession_ids,
                )
                addactionsdf = self.session.data_manager.get_add_actions_dataframe(
                    addactionqueryset=addactionqueryset
                )

            if addactionsdf.empty:
                return pd.DataFrame()

            # Select only add actions that are not products from previous reactions
            if product_exists:
                condition = addactionsdf.apply(
                    lambda row: check_previous_reaction_products(
                        reaction_id=row["reaction_id_id"],
                        smiles=row["smiles"],
                    ),
                    axis=1,
                )
                filtered_df = addactionsdf[~condition]
            else:
                filtered_df = addactionsdf

            if filtered_df.empty:
                return pd.DataFrame()

            # Begin processing the data
            addactionsdf["uniquesolution"] = addactionsdf.apply(
                self.combine_strings, axis=1
            )

            # Group by the unique solution identifier
            materials_df = filtered_df.groupby(["uniquesolution"]).agg(
                {
                    # "reaction_id_id": "first",
                    "smiles": "first",
                    "volume": "sum",
                    "solvent": "first",
                    "concentration": "first",
                    "molecularweight": "first",
                }
            )

            # Add this line to reset the index, making uniquesolution a regular column
            materials_df = materials_df.reset_index()

            logger.debug(f"Grouped materials dataframe: {materials_df.head()}")

            # Canonicalize SMILES for consistent comparison
            materials_df["SMILES"] = materials_df["smiles"].apply(canon_smiles)

            # Check existing materials and get remaining volume needed
            adjusted_materials = []
            updated_add_actions = {}  # Track which add actions were updated

            for _, row in materials_df.iterrows():
                logger.debug(f"Processing row: {row}")
                # Add 10% safety margin to volume
                total_volume_needed = row["volume"] * 1.1

                # Check if this material already exists
                (
                    exists,
                    matching_wells,
                    plate,
                    remaining_volume,
                    actual_smiles,
                ) = self.check_starting_material_exists(
                    smiles=row["smiles"],
                    volume=total_volume_needed,
                    concentration=row["concentration"],
                    solvent=row["solvent"],
                )

                # Handle salt form matches - update add actions if needed
                if exists and actual_smiles and actual_smiles != row["smiles"]:
                    # Find which add action(s) used this material
                    add_actions_to_update = addactionsdf[
                        (addactionsdf["uniquesolution"] == row["uniquesolution"])
                    ]

                    for idx, add_action_row in add_actions_to_update.iterrows():
                        add_action_id = add_action_row["id"]
                        logger.info(
                            f"Updating add action {add_action_id} for salt form match: "
                            f"{row['smiles']} -> {actual_smiles}"
                        )

                        # Only update if we haven't already processed this add action
                        if add_action_id not in updated_add_actions:
                            try:
                                add_action = AddAction.objects.get(id=add_action_id)

                                # Calculate new volume based on molecular weight differences
                                new_volume = (
                                    self.calculator.recalculate_amount_for_salt_form(
                                        row["smiles"], actual_smiles, add_action.volume
                                    )
                                )

                                # Update the add action
                                add_action.smiles = actual_smiles
                                add_action.volume = new_volume
                                add_action.save()

                                logger.info(
                                    f"Updated add action {add_action_id}: "
                                    f"SMILES from {row['smiles']} to {actual_smiles}, "
                                    f"Amount from {add_action_row['volume']} to {new_volume}µL"
                                )

                                # Track that we've updated this add action
                                updated_add_actions[add_action_id] = True

                            except Exception as e:
                                logger.error(
                                    f"Error updating add action {add_action_id}: {str(e)}"
                                )

                # Only include materials that need additional volume
                if not exists or remaining_volume > 0:
                    # Copy the row and update the volume to what's additionally needed
                    adjusted_row = row.copy()
                    if exists:
                        # Only request the additional volume needed
                        adjusted_row["volume"] = remaining_volume

                    adjusted_materials.append(adjusted_row)

            if not adjusted_materials:
                return pd.DataFrame()

            logger.debug(
                f"Adjusted materials after checking existing ones: {adjusted_materials}"
            )
            # Create final dataframe with adjusted volumes
            result_df = pd.DataFrame(adjusted_materials)

            # Calculate mass needed for each material
            result_df["mass_mg"] = result_df.apply(self.calc_mass, axis=1)

            # Add molecular weight column calculated from the SMILES
            result_df["molecular_weight"] = result_df["smiles"].apply(
                lambda s: Descriptors.MolWt(Chem.MolFromSmiles(strip_salts(s)))
            )

            # Use the actual required volume (convert μL to mL for consistency)
            result_df["stock_volume_ml"] = (
                result_df["volume"] / 1000.0
            )  # Convert μL to mL
            result_df["stock_concentration_M"] = result_df.apply(
                lambda row: (row["mass_mg"] / row["molecular_weight"])
                / row["stock_volume_ml"]
                if row["stock_volume_ml"] > 0
                else 0,
                axis=1,
            )

            # Add safety margin to mass (10%)
            result_df["mass_mg_with_margin"] = result_df["mass_mg"] * 1.1

            # Sort by solvent and volume
            result_df = result_df.sort_values(["solvent", "volume"], ascending=False)

            return result_df

        except Exception as e:
            logger.exception(f"Error in get_add_actions_material_dataframe: {str(e)}")
            return pd.DataFrame()

    def get_product_smiles(self, reaction_ids: list) -> list:
        """
        Get product smiles of reactions.

        Parameters
        ----------
        reaction_ids: list
            The reactions to get product smiles for

        Returns
        -------
        product_smiles: list
            The list of product smiles
        """
        try:
            product_smiles = Product.objects.filter(
                reaction_id__in=reaction_ids
            ).values_list("smiles", flat=True)

            if product_smiles:
                # Canonicalize the SMILES strings for consistent comparison
                canonicalized_smiles = [
                    canon_smiles(smiles) for smiles in product_smiles
                ]
                return list(set(canonicalized_smiles))  # Remove duplicates
            else:
                return []

        except Exception as e:
            logger.error(f"Error getting product SMILES: {str(e)}")
            return []

    def check_starting_material_exists(self, smiles, volume, concentration, solvent):
        """Coordinates the process of finding materials with salt handling"""
        canonical_smiles = canon_smiles(smiles)
        # Query DB for potential matches
        potential_matches = self._get_potential_matches(concentration, solvent)
        # Delegate the complex salt matching to specialized service
        (
            exists,
            matching_wells,
            containing_plate,
            remaining_volume,
            actual_smiles,
        ) = self.salt_matcher.find_matching_materials(
            canonical_smiles, potential_matches, volume
        )

        # Log when a salt form is found
        if actual_smiles and actual_smiles != smiles:
            logger.info(
                f"Found salt form match. Original: {smiles}, Found: {actual_smiles}"
            )

        return exists, matching_wells, containing_plate, remaining_volume, actual_smiles

    def _get_potential_matches(self, concentration, solvent):
        """
        Retrieves potential matches from the database based on concentration and solvent.

        Parameters
        ----------
        concentration: float
            The concentration needed
        solvent: str
            The solvent type needed

        Returns
        -------
        QuerySet
            Potential matches from the database
        """
        try:
            # Get custom starting material plates
            custom_compound_orders = CompoundOrder.objects.filter(
                otsession_id__otbatchprotocol_id=self.session.otbatchprotocolobj,
                iscustomSMplate=True,
            )

            if not custom_compound_orders.exists():
                logger.debug("No custom starting material compound orders found")
                return []

            # Get session IDs from the compound orders
            custom_session_ids = custom_compound_orders.values_list(
                "otsession_id", flat=True
            )

            # Get only custom starting material plates
            plates = Plate.objects.filter(
                otbatchprotocol_id=self.session.otbatchprotocolobj,
                role="startingmaterial",
                otsession_id__in=custom_session_ids,
            )

            if not plates.exists():
                logger.debug("No custom starting material plates found")
                return []

            # Retrieve wells with matching concentration and solvent
            potential_matches = []
            for plate in plates:
                wells = Well.objects.filter(
                    plate_id=plate.id,
                    concentration=concentration,
                    solvent=solvent,
                )
                potential_matches.extend(wells)

            return potential_matches

        except Exception as e:
            logger.exception(f"Error retrieving potential matches: {str(e)}")
            return []

    def get_max_well_volume(self, plateobj: Plate) -> float:
        """
        Get max well volume of a well plate.

        Parameters
        ----------
        plateobj: Plate
            The plate to get the max well volume of

        Returns
        -------
        maxwellvolume: float
            The maximum well volume of a well plate
        """
        maxwellvolume = plateobj.maxwellvolume
        return maxwellvolume

    def get_dead_volume(self, maxwellvolume: float) -> float:
        """
        Calculates the dead volume (5%) of a well.

        Parameters
        ----------
        maxwellvolume: float
            The well's maximum volume

        Returns
        -------
        deadvolume: float
            The dead volume of the well
        """
        deadvolume = maxwellvolume * 0.05
        return deadvolume

    def calc_mass(self, row) -> float:
        """
        Calculates the mass of material (mg) from the
        concentration (mol/L) and volume (uL) needed.

        Parameters
        ----------
        row: DataFrame row
            The row from the dataframe containing the
            concentration and volume information

        Returns
        -------
        mass_mg: float
            The mass of the material needed
        """
        try:
            # Convert units: concentration (mol/L) * volume (µL) * 1e-6 (L/µL) = moles
            volume_field = "amount-uL" if "amount-uL" in row else "volume"
            mols = row["concentration"] * row[volume_field] * 1e-6

            # Get SMILES string
            smiles_field = "SMILES" if "SMILES" in row else "smiles"
            smiles = row[smiles_field]

            # Strip salts before calculating molecular weight
            cleaned_smiles = strip_salts(smiles)
            mol = Chem.MolFromSmiles(cleaned_smiles)

            if mol:
                # Calculate molecular weight and mass
                mw = Descriptors.MolWt(mol)
                mass_mg = mols * mw * 1e3  # Convert to mg
                return round(mass_mg, 3)  # Increased precision
            else:
                logger.warning(f"Could not calculate mass for SMILES: {smiles}")
                return 0.0
        except Exception as e:
            logger.error(f"Error calculating mass: {str(e)}")
            return 0.0

    def create_solvent_prep_model(self, solvent_df: pd.DataFrame):
        """
        Creates a Django solvent prep object - a solvent prep file.

        Parameters
        ----------
        solvent_df: DataFrame
            The solvent dataframe containing plate information,
            well indices, solvents, and volumes required

        Returns
        -------
        solvent_prep_obj: SolventPrep
            The created solvent preparation object
        """
        try:
            if solvent_df.empty:
                logger.warning(
                    "Empty solvent dataframe, not creating solvent prep model"
                )
                return None

            solvent_prep_obj = SolventPrep()
            solvent_prep_obj.otsession_id = self.session.otsessionobj

            # Convert to CSV and save
            csv_data = solvent_df.to_csv(encoding="utf-8", index=False)
            file_name = (
                f"{self.session.actionsessiontype}-session-solventplate-"
                f"for-batch-{self.session.batchobj.batchtag}-"
                f"reactionstep-{self.session.reactionstep}-"
                f"sessionid-{str(self.session.otsessionobj.id)}.csv"
            )

            order_csv = default_storage.save(
                f"solventprep/{file_name}",
                ContentFile(csv_data),
            )

            solvent_prep_obj.solventprepcsv = order_csv
            solvent_prep_obj.save()

            logger.info(f"Created solvent preparation model: {solvent_prep_obj.id}")
            return solvent_prep_obj

        except Exception as e:
            logger.error(f"Error creating solvent prep model: {str(e)}")
            return None

    def combine_strings(self, row):
        """
        Combine SMILES, solvent and concentration into a unique identifier string.

        Parameters
        ----------
        row: DataFrame row
            The row from the dataframe containing smiles, solvent, and concentration

        Returns
        -------
        combined_string: str
            The combined string identifier
        """
        try:
            return (
                str(row["smiles"])
                + "-"
                + str(row["solvent"])
                + "-"
                + str(row["concentration"])
            )
        except Exception as e:
            logger.error(f"Error combining strings: {str(e)}")
            return "unknown"

    def get_number_vials(self, max_volume_vial: float, volume_material: float) -> int:
        """
        Gets the total number of vials needed to prepare a starter plate.

        Parameters
        ----------
        max_volume_vial: float
            The maximum volume of a vial
        volume_material: float
            The volume of the material that needs to be stored in a vial

        Returns
        -------
        no_vials_needed: int
            The number of vials needed to store the material
        """
        # Add logging to understand inputs
        logger.debug(
            f"Calculating vials needed for volume {volume_material}µL in vial of {max_volume_vial}µL"
        )

        # Special handling for microvolumes (< 5 µL)
        is_microvolume = volume_material < 5

        # For tiny volumes, just use 1 vial regardless
        if is_microvolume:
            logger.debug(f"Microvolume detected ({volume_material}µL), using 1 vial")
            return 1

        # Regular calculation for normal volumes
        if max_volume_vial > volume_material:
            logger.debug(
                f"Volume {volume_material}µL fits in one vial of {max_volume_vial}µL"
            )
            return 1
        else:
            # Calculate dead volume - use a fixed minimum to avoid issues with very small wells
            dead_volume = self.get_dead_volume(max_volume_vial)

            # Calculate effective capacity per vial
            effective_capacity = max_volume_vial - dead_volume

            # Calculate needed vials
            no_vials_needed_raw = volume_material / effective_capacity
            no_vials_needed = math.ceil(no_vials_needed_raw)

            logger.debug(
                f"Volume {volume_material}µL requires {no_vials_needed} vials "
                f"(effective capacity: {effective_capacity}µL per vial, "
                f"dead volume: {dead_volume}µL)"
            )

            return no_vials_needed

    def prepare_materials_for_reaction(self):
        """
        Prepare all materials needed for a reaction session.
        Creates solvent preparation models and compound orders as needed.
        """
        try:
            # Get materials that need to be prepared
            materials_df = self.get_add_actions_material_dataframe(product_exists=True)

            if not materials_df.empty:
                # Create solvent preparation models
                self.create_solvent_prep_model(solvent_df=materials_df)

                # Create compound order if needed
                self.session.data_manager.create_compound_order_model(
                    orderdf=materials_df, is_custom_starter_plate=False
                )

            return materials_df

        except Exception as e:
            logger.error(f"Error preparing materials for reaction: {str(e)}")
            return pd.DataFrame()

    def prepare_materials_for_workup(self):
        """
        Prepare all materials needed for a workup session.
        """
        try:
            # Get materials for workup
            materials_df = self.get_add_actions_material_dataframe(product_exists=False)

            if not materials_df.empty:
                # Create solvent preparation models
                self.create_solvent_prep_model(solvent_df=materials_df)

                # Create compound order if needed
                self.session.data_manager.create_compound_order_model(
                    orderdf=materials_df, is_custom_starter_plate=False
                )

            return materials_df

        except Exception as e:
            logger.error(f"Error preparing materials for workup: {str(e)}")
            return pd.DataFrame()
