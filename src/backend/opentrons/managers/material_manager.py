"""
Manages material properties and preparation for OpenTrons sessions.
"""

import logging
import math
import pandas as pd
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from rdkit import Chem
from rdkit.Chem import Descriptors

from ...models import (
    Plate, Well, CompoundOrder, SolventPrep, 
    Product
)
from ...utils import (
    canonSmiles, stripSalts, getInchiKey, getChemicalName,
    checkPreviousReactionProducts, getReaction
)

logger = logging.getLogger(__name__)

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
            if hasattr(self.session, 'addactionsdf'):
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
                
            # Begin processing the data
            # Select only add actions that are not products from previous reactions
            # if product_exists is True, or all add actions if False
            if product_exists:
                condition = addactionsdf["smiles"].apply(
                    lambda x: checkPreviousReactionProducts(
                        smiles=x, reaction_ids=self.session.reaction_ids
                    )
                )
                filtered_df = addactionsdf[~condition]
            else:
                filtered_df = addactionsdf
                
            if filtered_df.empty:
                return pd.DataFrame()
            
            # Group by SMILES, concentration, and solvent, summing the volumes
            materials_df = filtered_df.groupby(["smiles", "concentration", "solvent"]).agg(
                {
                    "reaction_id_id": "first",
                    "volume": "sum",
                    "molecularweight": "first",
                }
            ).reset_index()
            
            # Create a unique identifier for each material
            materials_df["uniquesolution"] = materials_df.apply(
                self.combine_strings, axis=1
            )
            
            # Canonicalize SMILES for consistent comparison
            materials_df["SMILES"] = materials_df["smiles"].apply(canonSmiles)
            
            # Check existing materials and get remaining volume needed
            adjusted_materials = []
            for _, row in materials_df.iterrows():
                # Add 10% safety margin to volume
                total_volume_needed = row["volume"] * 1.1
                
                # Check if this material already exists
                exists, matching_wells, plate, remaining_volume = self.check_starting_material_exists(
                    smiles=row["smiles"],
                    volume=total_volume_needed,
                    concentration=row["concentration"],
                    solvent=row["solvent"],
                )
                
                # Only include materials that need additional volume
                if not exists or remaining_volume > 0:
                    # Copy the row and update the volume to what's additionally needed
                    adjusted_row = row.copy()
                    if exists:
                        # Only request the additional volume needed
                        adjusted_row["volume"] = remaining_volume
                        logger.info(
                            f"Found partial volume for {row['smiles']}. "
                            f"Total needed: {total_volume_needed}µL, "
                            f"Already available: {total_volume_needed - remaining_volume}µL, "
                            f"Still need: {remaining_volume}µL"
                        )
                    
                    adjusted_materials.append(adjusted_row)
            
            if not adjusted_materials:
                return pd.DataFrame()
                
            # Create final dataframe with adjusted volumes
            result_df = pd.DataFrame(adjusted_materials)
            result_df = result_df.sort_values(["solvent", "volume"], ascending=False)
            
            return result_df
            
        except Exception as e:
            logger.error(f"Error in get_add_actions_material_dataframe: {str(e)}")
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
                canonicalized_smiles = [canonSmiles(smiles) for smiles in product_smiles]
                return list(set(canonicalized_smiles))  # Remove duplicates
            else:
                return []
                
        except Exception as e:
            logger.error(f"Error getting product SMILES: {str(e)}")
            return []
    
    def check_starting_material_exists(
        self, smiles: str, volume: float, concentration: float, solvent: str
    ) -> tuple:
        """
        Checks if starting material exists with enough total volume across wells in current OT batch protocol.
        Only considers custom starter plates loaded from CSV files, not auto-generated plates from previous sessions.
        
        Parameters
        ----------
        smiles: str
            The SMILES string of the starting material to check
        volume: float
            The total volume needed in microliters
        concentration: float
            The concentration in moles per liter
        solvent: str
            The solvent used
            
        Returns
        -------
        tuple
            (exists: bool, matching_wells: list[Well], plate: Plate, remaining_volume_needed: float)
            Returns if material exists with enough volume, list of matching wells, plate containing wells,
            and remaining volume needed (0 if all volume is available)
        """
        try:
            # Canonicalize the SMILES for consistent comparison
            canonical_smiles = canonSmiles(smiles)
            
            # Get custom starting material plates
            custom_compound_orders = CompoundOrder.objects.filter(
                otsession_id__otbatchprotocol_id=self.session.otbatchprotocolobj,
                iscustomSMplate=True,
            )
            
            # Get session IDs from the compound orders
            custom_session_ids = custom_compound_orders.values_list(
                "otsession_id", flat=True
            )
            
            # Get only custom starting material plates
            plates = Plate.objects.filter(
                otbatchprotocol_id=self.session.otbatchprotocolobj,
                type="startingmaterial",
                otsession_id__in=custom_session_ids,
            )
            
            if not plates.exists():
                return (False, [], None, volume)
            
            # Track total volume and matching wells across all plates
            total_available_volume = 0
            all_matching_wells = []
            containing_plate = None
            
            # Check each plate for wells containing the desired material
            for plate in plates:
                matching_wells = Well.objects.filter(
                    plate_id=plate.id,
                    smiles=canonical_smiles,
                    concentration=concentration,
                    solvent=solvent,
                )
                
                if matching_wells.exists():
                    containing_plate = plate
                    logger.debug(f"Found {matching_wells.count()} matching wells in plate {plate.name}")
                    
                    for well in matching_wells:
                        if well.volume is not None:
                            all_matching_wells.append(well)
                            total_available_volume += well.volume
            
            # Calculate remaining volume needed
            remaining_volume_needed = max(0, volume - total_available_volume)
            
            # Check if we have enough volume
            if remaining_volume_needed <= 0:
                logger.info(
                    f"Found enough material in custom plates: {canonical_smiles} "
                    f"Required: {volume}µL, Available: {total_available_volume}µL"
                )
                return (True, all_matching_wells, containing_plate, 0)
            else:
                if all_matching_wells:
                    logger.info(
                        f"Found material {canonical_smiles} but insufficient volume. "
                        f"Required: {volume}µL, Available: {total_available_volume}µL, "
                        f"Still need: {remaining_volume_needed}µL"
                    )
                    return (False, all_matching_wells, containing_plate, remaining_volume_needed)
                else:
                    logger.info(f"No existing material found for: {canonical_smiles}")
                    return (False, [], None, volume)
                    
        except Exception as e:
            logger.error(f"Error checking if starting material exists: {str(e)}")
            return (False, [], None, volume)
    
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
        concentration (mol/L) and volume (ul) needed.
        
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
            cleaned_smiles = stripSalts(smiles)
            mol = Chem.MolFromSmiles(cleaned_smiles)
            
            if mol:
                # Calculate molecular weight and mass
                mw = Descriptors.MolWt(mol)
                mass_mg = mols * mw * 1e3  # Convert to mg
                return round(mass_mg, 2)
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
                logger.warning("Empty solvent dataframe, not creating solvent prep model")
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
        if max_volume_vial > volume_material:
            no_vials_needed = 1
        else:
            volumes_to_add = []
            dead_volume = self.session.plate_manager.get_dead_volume(max_well_volume=max_volume_vial)
            no_vials_needed_ratio = volume_material / (max_volume_vial - dead_volume)
            frac, whole = math.modf(no_vials_needed_ratio)
            volumes_to_add = [max_volume_vial for _ in range(int(whole))]
            volumes_to_add.append(frac * max_volume_vial + dead_volume)
            no_vials_needed = len(volumes_to_add)
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