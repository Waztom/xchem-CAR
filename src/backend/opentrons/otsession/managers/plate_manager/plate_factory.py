import logging
import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

from .....models import Plate, AddAction
from .....utils import sanitize_for_python_var
from ....labwareavailable import labware_plates

logger = logging.getLogger(__name__)


class PlateFactory:
    """Factory for creating different types of plates."""

    def __init__(self, session):
        """Initialize with session reference."""
        self.session = session

    def create_plate_model(self, platetype: str, platename: str, labwaretype: str):
        """
        Creates Plate model if Deck index is available.

        Parameters
        ----------
        platetype: str
            The type of plate (e.g., reaction, startingmaterial)
        platename: str
            The name of the plate
        labwaretype: str
            The labware type (e.g., plateone_96_wellplate_2500ul)

        Returns
        -------
        plate_obj: Plate or None
            The created plate object, or None if no deck slots available
        """
        index_slot = self.session.deck_manager.check_deck_slot_available()
        if index_slot:
            plate_index = index_slot
            number_wells_in_column = labware_plates[labwaretype]["no_wells_in_column"]
            max_well_volume = labware_plates[labwaretype]["volume_well"]
            number_wells = labware_plates[labwaretype]["no_wells"]
            number_columns = labware_plates[labwaretype]["no_columns"]
            sanitized_platename = sanitize_for_python_var(platename)

            plate_obj = Plate()
            plate_obj.otbatchprotocol_id = self.session.otbatchprotocolobj
            plate_obj.otsession_id = self.session.otsessionobj
            plate_obj.deck_id = self.session.deckobj
            plate_obj.labware = labwaretype
            plate_obj.index = plate_index
            plate_obj.type = platetype
            plate_obj.maxwellvolume = max_well_volume
            plate_obj.numberwells = number_wells
            plate_obj.numberwellsincolumn = number_wells_in_column
            plate_obj.numbercolumns = number_columns
            plate_obj.save()

            plate_obj.name = f"Reaction_step_{self.session.reactionstep}_{sanitized_platename}_{plate_obj.id}"
            plate_obj.save()

            return plate_obj
        else:
            logger.warning("CreatePlateModel - No more deck slots available")
            return None

    def create_plate_by_reaction_class_recipe(self, reaction_queryset, platetype: str):
        """
        Creates plates by reaction class and recipe.

        Parameters
        ----------
        reaction_queryset: QuerySet[Reaction]
            The reactions to create plates for
        platetype: str
            The type of plate to create

        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Get grouped reactions by class and recipe
        grouped_reaction_querysets = (
            self.session.data_manager.get_grouped_reaction_by_class_recipe(
                reactionqueryset=reaction_queryset
            )
        )

        # Store created plates
        plate_objs = []

        # Create a plate for each group
        for group_key, group_reactions in grouped_reaction_querysets.items():
            # Extract class and recipe from group key
            reaction_class, recipe = group_key.split("-", 1)

            # Generate plate name
            plate_name = f"{platetype}-{reaction_class}-{recipe}"

            # Get volumes for the reaction group
            reaction_volumes = self.session.data_manager.get_rounded_reaction_volumes(
                reactionqueryset=group_reactions
            )

            # Determine best labware type
            labware_type = self.session.labware_selector.get_plate_type(
                platetype="reaction", volumes=reaction_volumes, temperature=25
            )

            # Create the plate
            plate_obj = self.create_plate_model(
                platetype=platetype, platename=plate_name, labwaretype=labware_type
            )

            if plate_obj:
                plate_objs.append(plate_obj)

                # Create a column for the reaction class-recipe
                column_index = (
                    self.session.column_manager.get_plate_current_column_index(
                        plate_obj=plate_obj
                    )
                )
                if column_index is not False:
                    column_obj = self.session.column_manager.create_column_model(
                        plate_obj=plate_obj,
                        columnindex=column_index,
                        columntype=platetype,
                        reactionclass=reaction_class,
                    )

                    # Update column index
                    self.session.column_manager.update_plate_column_index_available(
                        plate_obj=plate_obj, columnindexupdate=column_index + 1
                    )

                    # Add wells for each reaction in this group
                    for reaction in group_reactions:
                        # Get product SMILES for this reaction
                        product_smiles = (
                            self.session.material_manager.get_product_smiles(
                                [reaction.id]
                            )[0]
                        )

                        # Log product information
                        logger.debug(
                            f"Creating well for reaction {reaction.id} with product SMILES: "
                            f"{product_smiles}"
                        )

                        # Check if well index is available
                        well_index = (
                            self.session.well_manager.get_plate_well_index_available(
                                plate_obj=plate_obj
                            )
                        )
                        if well_index is False:
                            # No wells available - create a new plate with the same parameters
                            logger.info(
                                f"Plate {plate_obj.name} is full. Creating another plate."
                            )

                            # Create new plate with same parameters but incremented name
                            new_plate_name = f"{plate_name}-continued"
                            plate_obj = self.create_plate_model(
                                platetype=platetype,
                                platename=new_plate_name,
                                labwaretype=labware_type,
                            )

                            if not plate_obj:
                                logger.error(
                                    "Failed to create additional plate - no more deck slots?"
                                )
                                break

                            # Add the new plate to our result list
                            plate_objs.append(plate_obj)

                            # Create a new column in the new plate
                            column_index = self.session.column_manager.get_plate_current_column_index(
                                plate_obj=plate_obj
                            )
                            if column_index is False:
                                logger.error("Failed to get column index for new plate")
                                break

                            column_obj = (
                                self.session.column_manager.create_column_model(
                                    plate_obj=plate_obj,
                                    columnindex=column_index,
                                    columntype=platetype,
                                    reactionclass=reaction_class,
                                )
                            )

                            # Update column index
                            self.session.column_manager.update_plate_column_index_available(
                                plate_obj=plate_obj, columnindexupdate=column_index + 1
                            )

                            # Try to get well index again on new plate
                            well_index = self.session.well_manager.get_plate_well_index_available(
                                plate_obj=plate_obj
                            )
                            if well_index is False:
                                logger.error("New plate has no available wells")
                                break

                        # Create well with product information
                        well_obj = self.session.well_manager.create_well_model(
                            plate_obj=plate_obj,
                            welltype=platetype,
                            wellindex=well_index,
                            reactionobj=reaction,
                            columnobj=column_obj,
                            smiles=product_smiles,  # Store product SMILES
                            reactantfornextstep=True,  # Mark as available for next step
                        )

                        # Update well index
                        self.session.well_manager.update_plate_well_index(
                            plate_obj=plate_obj, wellindexupdate=well_index + 1
                        )
            else:
                logger.warning(f"Could not create plate for {group_key}")

        return plate_objs

    def create_reaction_plate(self, platetype: str = "reaction"):
        """
        Creates a reaction plate.

        Parameters
        ----------
        platetype: str
            The type of plate to create, default is "reaction"

        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Get grouped reactions by temperature
        grouped_reaction_querysets = (
            self.session.data_manager.get_grouped_temperature_reactions(
                reactionqueryset=self.session.groupreactionqueryset
            )
        )

        # Store created plates
        plate_objs = []

        # Create reaction plates based on temperature groups
        for temp, reaction_queryset in grouped_reaction_querysets.items():
            # Create plates by reaction class and recipe
            temp_plates = self.create_plate_by_reaction_class_recipe(
                reaction_queryset=reaction_queryset, platetype=platetype
            )

            plate_objs.extend(temp_plates)

        return plate_objs

    def create_workup_plate(self, platetype: str):
        """
        Creates a workup plate.

        Parameters
        ----------
        platetype: str
            The type of workup plate (e.g., "workup1", "workup2")

        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Check if we need this plate type for the current actions
        action_session_queryset = (
            self.session.plate_query_service.get_action_session_by_plate_type(
                platetype=platetype
            )
        )

        if not action_session_queryset.exists():
            logger.info(f"No action sessions require {platetype} plate")
            return []

        # Get volumes for determining plate type
        add_action_queryset = self.session.data_manager.get_add_action_query_set(
            reaction_ids=self.session.reaction_ids,
            actionsession_ids=self.session.actionsession_ids,
        )

        rounded_volumes = self.session.data_manager.get_rounded_add_action_volumes(
            addactionqueryset=add_action_queryset
        )

        # Determine best labware type
        labware_type = self.session.labware_selector.get_plate_type(
            platetype="reaction",  # Use reaction type for labware compatibility
            volumes=rounded_volumes,
            temperature=25,
        )

        # Create the plate
        plate_obj = self.create_plate_model(
            platetype=platetype,
            platename=f"{platetype}-plate",
            labwaretype=labware_type,
        )

        if plate_obj:
            return [plate_obj]
        else:
            logger.warning(f"Could not create {platetype} plate")
            return []

    def create_analyse_plate(self):
        """
        Creates an analysis plate.

        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Get volumes for determining plate type
        extract_action_queryset = (
            self.session.data_manager.get_extract_action_query_set(
                reaction_ids=self.session.reaction_ids,
                actionsession_ids=self.session.actionsession_ids,
            )
        )

        rounded_volumes = self.session.data_manager.get_rounded_extract_action_volumes(
            extractactionqueryset=extract_action_queryset
        )

        # Determine best labware type - analysis plates often use smaller volumes
        labware_type = self.session.labware_selector.get_plate_type(
            platetype="analyse", volumes=rounded_volumes, temperature=25
        )

        # Create the plate
        plate_obj = self.create_plate_model(
            platetype="analyse", platename="analyse-plate", labwaretype=labware_type
        )

        if plate_obj:
            return [plate_obj]
        else:
            logger.warning("Could not create analyse plate")
            return []

    def create_solvent_plate(self, materials_df: pd.DataFrame):
        """
        Creates a solvent plate. Used to create solvent plates
        for workup/cleanup solvents.

        Parameters
        ----------
        materials_df: DataFrame
            Dataframe containing material information

        Returns
        -------
        plate_obj: Plate
            The created solvent plate
        """
        if materials_df.empty:
            logger.info("No materials provided for solvent plate")
            return None

        # Group materials by solvent
        solvent_groups = materials_df.groupby("solvent")

        solvent_dicts_list = []

        for solvent_group, group_df in solvent_groups:
            # Calculate total volume needed for this solvent
            total_volume = group_df["volume"].sum() * 1.1  # Add 10% safety margin

            # Create or use existing plate based on volume
            labware_type = self.session.labware_selector.get_plate_type(
                platetype="startingmaterial", volumes=[total_volume], temperature=25
            )

            # Create a solvent plate (or reuse existing one)
            plate_obj = None
            existing_solvent_plates = Plate.objects.filter(
                otsession_id=self.session.otsessionobj, type="solvent"
            )

            if existing_solvent_plates.exists():
                for p in existing_solvent_plates:
                    if p.labware == labware_type:
                        plate_obj = p
                        break

            if not plate_obj:
                plate_obj = self.create_plate_model(
                    platetype="solvent",
                    platename=f"solvent-{solvent_group}",
                    labwaretype=labware_type,
                )

            if not plate_obj:
                logger.warning(f"Could not create solvent plate for {solvent_group}")
                continue

            # Check for available well
            index_well_available = (
                self.session.well_manager.get_plate_well_index_available(
                    plate_obj=plate_obj
                )
            )
            if index_well_available is False:
                logger.warning(f"No wells available on solvent plate {plate_obj.name}")
                continue

            # Create well for this solvent
            well_obj = self.session.well_manager.create_well_model(
                plate_obj=plate_obj,
                welltype="solvent",
                wellindex=index_well_available,
                volume=total_volume,
                solvent=solvent_group,
            )

            # Update well index
            self.session.well_manager.update_plate_well_index(
                plate_obj=plate_obj, wellindexupdate=index_well_available + 1
            )

            # Add to solvent prep list
            solvent_dicts_list.append(
                {
                    "name": plate_obj.name,
                    "labware": plate_obj.labware,
                    "well-index": well_obj.index,
                    "well-name": well_obj.name,
                    "solvent": solvent_group,
                    "amount-uL": total_volume,
                }
            )

        # Create solvent prep model if we have solvents
        if solvent_dicts_list:
            solvent_df = pd.DataFrame(solvent_dicts_list)
            self.session.material_manager.create_solvent_prep_model(
                solvent_df=solvent_df
            )

        return plate_obj

    def create_starting_material_plates_from_csv(self, csv_path: str):
        """
        Creates starting material plates from a CSV file.

        Parameters
        ----------
        csv_path: str
            Path to the CSV file containing starting material information

        Returns
        -------
        plate_obj: Plate
            The created starting material plate
        """
        try:
            # Read the CSV file
            materials_df = pd.read_csv(csv_path)

            if materials_df.empty:
                logger.warning("CSV file is empty")
                return None

            # Make sure required columns exist
            required_columns = [
                "plate-ID",
                "labware-type",
                "well-index",
                "SMILES",
                "amount-uL",
            ]
            missing_columns = [
                col for col in required_columns if col not in materials_df.columns
            ]

            if missing_columns:
                logger.error(
                    f"CSV is missing required columns: {', '.join(missing_columns)}"
                )
                return None

            # Get volumes for plate sizing
            volumes = materials_df["amount-uL"].tolist()

            # Determine best labware type
            labware_type = self.session.labware_selector.get_plate_type(
                platetype="startingmaterial", volumes=volumes, temperature=25
            )

            # Create the plate
            plate_obj = self.create_plate_model(
                platetype="startingmaterial",
                platename="custom-starting-materials",
                labwaretype=labware_type,
            )

            if not plate_obj:
                logger.warning("Could not create custom starting material plate")
                return None

            # Add each material to the plate
            for _, row in materials_df.iterrows():
                # Check for available well
                index_well_available = (
                    self.session.well_manager.get_plate_well_index_available(
                        plate_obj=plate_obj
                    )
                )
                if index_well_available is False:
                    logger.warning(
                        f"No wells available on custom starting plate {plate_obj.name}"
                    )
                    break

                # Create well for this material
                well_obj = self.session.well_manager.create_well_model(
                    plate_obj=plate_obj,
                    welltype="startingmaterial",
                    wellindex=index_well_available,
                    volume=row["amount-uL"],
                    smiles=row["SMILES"],
                    concentration=row["concentration"],
                    solvent=row["solvent"],
                )

                # Update well index
                self.session.well_manager.update_plate_well_index(
                    plate_obj=plate_obj, wellindexupdate=index_well_available + 1
                )

            # Create compound order record
            self.session.data_manager.create_compound_order_model(
                order_df=materials_df, is_custom_starter_plate=True
            )

            return plate_obj

        except Exception as e:
            logger.error(f"Error creating starting material plates from CSV: {str(e)}")
            return None

    def create_reaction_starting_plate(self):
        """
        Creates a starting material plate for reaction materials, respecting the maximum
        fill percentage defined in labware configuration.

        Returns
        -------
        plate_obj: Plate
            The created starting material plate, or None if no new materials needed
        """
        # Get materials that need to be prepared
        materials_df = self.session.material_manager.get_add_actions_material_dataframe(
            product_exists=True
        )

        logger.debug(f"Materials needed for starting plate: {materials_df}")
        if materials_df.empty:
            logger.info("No materials needed for starting plate")
            return None

        # Get total volumes for plate sizing
        volumes = materials_df["volume"].tolist()

        # Determine best labware type
        labware_type = self.session.labware_selector.get_plate_type(
            platetype="startingmaterial", volumes=volumes, temperature=25
        )

        # Create the plate
        plate_obj = self.create_plate_model(
            platetype="startingmaterial",
            platename="reaction-starting-materials",
            labwaretype=labware_type,
        )

        if not plate_obj:
            logger.warning("Could not create starting material plate")
            return None

        # Get the maximum well volume with headspace consideration
        effective_max_volume = self.session.well_manager.get_max_well_volume(plate_obj)
        raw_max_volume = plate_obj.maxwellvolume

        # Log the effective max volume with max fill percentage applied
        logger.info(
            f"Using effective max volume of {effective_max_volume}µL for plate with raw max volume {raw_max_volume}µL"
        )

        # Calculate dead volume for this plate type
        dead_volume = self.session.labware_selector.get_dead_volume(raw_max_volume)
        logger.debug(f"Using dead volume of {dead_volume}µL")

        # Sort materials to optimize well usage (by similar properties)
        sorted_materials_df = materials_df.sort_values(
            by=["smiles", "concentration", "solvent"]
        )

        # Tracking for compound order
        order_dicts_list = []

        # Process each material - ALWAYS use a new well for each material
        for _, row in sorted_materials_df.iterrows():
            # Add 5% extra volume as error margin
            extra_error_volume = row["volume"] * 0.05
            total_volume_needed = row["volume"] + extra_error_volume + dead_volume

            # Get next available well index
            well_idx = self.session.well_manager.get_plate_well_index_available(
                plate_obj
            )

            # If no well is available, create a new plate
            if well_idx is False:
                plate_obj = self.create_plate_model(
                    platetype="startingmaterial",
                    platename="reaction-starting-materials",
                    labwaretype=labware_type,
                )
                if not plate_obj:
                    logger.warning(
                        "Could not create additional starting material plate"
                    )
                    break

                well_idx = self.session.well_manager.get_plate_well_index_available(
                    plate_obj
                )
                effective_max_volume = self.session.well_manager.get_max_well_volume(
                    plate_obj
                )

            # Ensure we don't exceed maximum volume limit
            if total_volume_needed > effective_max_volume:
                logger.warning(
                    f"Material volume {total_volume_needed}µL exceeds effective max volume {effective_max_volume}µL; "
                    f"Material will be split across multiple wells"
                )

                # Split into multiple wells if needed
                remaining_volume = total_volume_needed
                while remaining_volume > 0:
                    volume_this_well = min(remaining_volume, effective_max_volume)

                    # Create well
                    well_obj = self.session.well_manager.create_well_model(
                        plate_obj=plate_obj,
                        welltype="startingmaterial",
                        wellindex=well_idx,
                        volume=volume_this_well,
                        smiles=row["smiles"],
                        concentration=row["concentration"],
                        solvent=row["solvent"],
                    )

                    if not well_obj:
                        logger.warning(
                            f"Could not create well for material {row['smiles']}"
                        )
                        break

                    # Add to order list
                    order_dict = {
                        "SMILES": row["smiles"],
                        "plate-name": plate_obj.name,
                        "labware": plate_obj.labware,
                        "well-index": well_obj.index,
                        "well-name": well_obj.name,
                        "concentration": row["concentration"],
                        "solvent": row["solvent"],
                        "molecularweight": row.get("molecularweight"),
                        "amount-uL": round(volume_this_well, 2),
                    }
                    order_dicts_list.append(order_dict)

                    # Update remaining volume
                    remaining_volume -= volume_this_well

                    # CRITICAL FIX: Always update the well index after creating a well,
                    # regardless of whether there's more volume to process
                    self.session.well_manager.update_plate_well_index(
                        plate_obj=plate_obj, wellindexupdate=well_idx + 1
                    )

                    # Only get a new well index if there's more volume to process
                    if remaining_volume > 0:
                        well_idx = (
                            self.session.well_manager.get_plate_well_index_available(
                                plate_obj
                            )
                        )
                        if well_idx is False:
                            plate_obj = self.create_plate_model(
                                platetype="startingmaterial",
                                platename="reaction-starting-materials",
                                labwaretype=labware_type,
                            )
                            if not plate_obj:
                                logger.warning(
                                    "Could not create additional starting material plate"
                                )
                                break

                            well_idx = self.session.well_manager.get_plate_well_index_available(
                                plate_obj
                            )
            else:
                # Create a single well for this material
                well_obj = self.session.well_manager.create_well_model(
                    plate_obj=plate_obj,
                    welltype="startingmaterial",
                    wellindex=well_idx,
                    volume=total_volume_needed,
                    smiles=row["smiles"],
                    concentration=row["concentration"],
                    solvent=row["solvent"],
                )

                if not well_obj:
                    logger.warning(
                        f"Could not create well for material {row['smiles']}"
                    )
                    continue

                # Add to order list
                order_dict = {
                    "SMILES": row["smiles"],
                    "plate-name": plate_obj.name,
                    "labware": plate_obj.labware,
                    "well-index": well_obj.index,
                    "well-name": well_obj.name,
                    "concentration": row["concentration"],
                    "solvent": row["solvent"],
                    "molecularweight": row.get("molecularweight"),
                    "amount-uL": round(total_volume_needed, 2),
                }
                order_dicts_list.append(order_dict)

            # Update plate's next available well index after creating this well
            self.session.well_manager.update_plate_well_index(
                plate_obj=plate_obj, wellindexupdate=well_idx + 1
            )

        # Create compound order from well assignments
        if order_dicts_list:
            order_df = pd.DataFrame(order_dicts_list)

            # Calculate additional fields
            if (
                "molecularweight" not in order_df.columns
                or order_df["molecularweight"].isna().any()
            ):
                order_df["molecularweight"] = order_df["SMILES"].apply(
                    lambda smiles: Descriptors.MolWt(Chem.MolFromSmiles(smiles))
                )

            # Calculate mass
            order_df["mass-mg"] = order_df.apply(
                self.session.material_manager.calc_mass, axis=1
            )

            # Add inchikey and compound name if available
            if hasattr(self.session, "data_manager") and hasattr(
                self.session.data_manager, "get_inchi_key"
            ):
                order_df["inchikey"] = order_df["SMILES"].apply(
                    self.session.data_manager.get_inchi_key
                )
                order_df["compound-name"] = order_df["inchikey"].apply(
                    self.session.data_manager.get_chemical_name
                )

            # Create the compound order
            self.session.data_manager.create_compound_order_model(
                order_df=order_df, is_custom_starter_plate=False
            )

        logger.info(
            f"Created starting material plate with {len(order_dicts_list)} material wells"
        )
        return plate_obj

    def create_plates_by_temperature(
        self, grouped_reaction_temperature_querysets, platetype="reaction"
    ):
        """
        Creates reaction plates for each temperature group using appropriate labware types.

        This function processes each temperature group, determines the appropriate labware type
        based on volumes and temperature, and creates plates for each group.

        Parameters
        ----------
        grouped_reaction_temperature_querysets: dict
            Dictionary mapping temperature to QuerySets of reactions at that temperature
        platetype: str
            Type of plate to create (e.g. "reaction", "workup1", etc.)

        Returns
        -------
        created_plates: list
            List of created plate objects
        """
        logger.info(
            f"Creating {platetype} plates for {len(grouped_reaction_temperature_querysets)} temperature groups"
        )
        created_plates = []

        for (
            temp,
            reaction_temp_queryset,
        ) in grouped_reaction_temperature_querysets.items():
            # Get volumes for this temperature group
            volumes = self.session.data_manager.get_rounded_reaction_volumes(
                reaction_temp_queryset
            )

            # Get the temperature from the first reaction in group
            reaction_temperature = temp

            # Determine appropriate labware based on volumes and temperature
            labware_platetype = self.session.labware_selector.get_plate_type(
                platetype=platetype,
                volumes=volumes,
                temperature=reaction_temperature,
            )

            logger.info(
                f"Creating {platetype} plate for temperature {reaction_temperature}°C using {labware_platetype}"
            )

            # Create the plate using the existing class-recipe function
            plate = self.create_plate_by_reaction_class_recipe(
                reaction_queryset=reaction_temp_queryset,
                platetype=platetype,
            )

            if plate:
                created_plates.extend(plate if isinstance(plate, list) else [plate])

        logger.info(f"Created {len(created_plates)} {platetype} plates by temperature")
        return created_plates

    def update_plate_deck_ot_session_ids(self, plate_queryset):
        """
        Updates plates to link to the current Deck and OT session.

        Parameters
        ----------
        plate_queryset: QuerySet[Plate]
            The plates to update
        """
        for plate_obj in plate_queryset:
            index_slot = self.session.deck_manager.check_deck_slot_available()
            if index_slot:
                well_queryset = self.session.well_manager.get_plate_wells(
                    plate_obj=plate_obj
                )
                column_queryset = self.session.column_manager.get_plate_columns(
                    plate_obj=plate_obj
                )
                previous_type = plate_obj.type
                plate_obj.deck_id = self.session.deckobj
                plate_obj.otsession_id = self.session.otsessionobj
                if previous_type == "spefilter":
                    plate_obj.labware = "plateone_96_wellplate_2500ul"
                plate_obj.index = index_slot
                plate_obj.save()
                self.session.column_manager.update_column_ot_session_ids(
                    column_queryset=column_queryset, plate_obj=plate_obj
                )
                self.session.well_manager.update_well_ot_session_ids(
                    well_queryset=well_queryset, plate_obj=plate_obj
                )
            else:
                logger.warning("cloneInputPlate - No more deck slots available")

        # Make sure all plates are linked to the current session and deck
        for plate_obj in plate_queryset:
            plate_obj.deck_id = self.session.deckobj
            plate_obj.otsession_id = self.session.otsessionobj
            plate_obj.save()

        return plate_queryset
