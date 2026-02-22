import logging
import math
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

from .....models import Plate, Well
from .....conversions import sanitize_for_python_var
from .....recipe_utils import unparse_plate_type
from ....labwareavailable import labware_plates
from ..multichannel_analyzer import MultichannelAnalyzer

logger = logging.getLogger(__name__)


class PlateFactory:
    """Factory for creating different types of plates."""

    def __init__(self, session):
        """Initialize with session reference."""
        self.session = session

    def create_plate_model(
        self,
        platename: str,
        labwaretype: str,
        *,
        role: str,
        role_index: int = 1,
    ):
        """
        Creates Plate model if Deck index is available.

        Parameters
        ----------
        platename: str
            The name of the plate
        labwaretype: str
            The labware type (e.g., plateone_96_wellplate_2500ul)
        role: str
            The plate role (e.g., "reaction", "workup").
        role_index: int, optional
            The plate role index (default 1).

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
            plate_obj.role = role
            plate_obj.role_index = role_index
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

   
    def create_workup_plate(
        self,
        *,
        role: str = "workup",
        role_index: int = 1,
    ):
        """
        Creates a workup plate.

        Parameters
        ----------
        role: str
            The plate role (default "workup").
        role_index: int
            The workup stage index (1, 2, or 3). Default 1.

        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Generate platetype string for display/logging
        platetype_str = unparse_plate_type(role, role_index)

        # Check if we need this plate type for the current actions
        action_session_queryset = (
            self.session.plate_query_service.get_action_session_by_plate_role(
                role=role,
                role_index=role_index,
            )
        )

        if not action_session_queryset.exists():
            logger.info(f"No action sessions require {platetype_str} plate")
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
            role="reaction",  # Use reaction type for labware compatibility
            volumes=rounded_volumes,
            temperature=25,
        )

        # Create the plate
        plate_obj = self.create_plate_model(
            role=role,
            role_index=role_index,
            platename=f"{platetype_str}-plate",
            labwaretype=labware_type,
        )

        if plate_obj:
            return [plate_obj]
        else:
            logger.warning(f"Could not create {platetype_str} plate")
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
            role="analyse", volumes=rounded_volumes, temperature=25
        )

        # Create the plate using role/role_index
        plate_obj = self.create_plate_model(
            role="analyse",
            role_index=1,
            platename="analyse-plate",
            labwaretype=labware_type,
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
                role="startingmaterial", volumes=[total_volume], temperature=25
            )

            # Create a solvent plate (or reuse existing one)
            plate_obj = None
            existing_solvent_plates = Plate.objects.filter(
                otsession_id=self.session.otsessionobj, role="solvent"
            )

            if existing_solvent_plates.exists():
                for p in existing_solvent_plates:
                    if p.labware == labware_type:
                        plate_obj = p
                        break

            if not plate_obj:
                plate_obj = self.create_plate_model(
                    role="solvent",
                    role_index=1,
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
                role="solvent",
                role_index=1,
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
        Creates starting material plates from a CSV file, handling multiple plate IDs
        and different labware types. Supports both numeric indices and well names.

        Parameters
        ----------
        csv_path: str
            Path to the CSV file containing starting material information

        Returns
        -------
        plate_objects: dict
            Dictionary mapping plate IDs to created plate objects
        """
        try:
            # Read the CSV file
            materials_df = pd.read_csv(csv_path)

            if materials_df.empty:
                logger.warning("CSV file is empty")
                return {}

            # Make sure required columns exist
            required_columns = [
                "plate-ID", 
                "labware-type", 
                "well-index", 
                "SMILES", 
                "amount-uL",
            ]
            missing_columns = [col for col in required_columns if col not in materials_df.columns]

            if missing_columns:
                logger.error(f"CSV is missing required columns: {', '.join(missing_columns)}")
                return {}
                
            # Group materials by plate ID
            plate_groups = materials_df.groupby("plate-ID")
            
            # Store created plates by ID
            created_plates = {}
            
            # Process each plate group
            for plate_id, plate_df in plate_groups:
                # Get labware type for this plate group
                labware_type = plate_df["labware-type"].iloc[0]
                
                # Check for consistent labware type within the plate
                if not plate_df["labware-type"].eq(labware_type).all():
                    logger.warning(
                        f"Inconsistent labware types for plate {plate_id}. "
                        f"Using '{labware_type}' from first row."
                    )
                    
                # Create the plate
                plate_obj = self.create_plate_model(
                    role="startingmaterial",
                    role_index=1,
                    platename=f"custom-starting-materials-{plate_id}",
                    labwaretype=labware_type,
                )
                
                if not plate_obj:
                    logger.warning(f"Could not create custom starting plate for ID {plate_id}")
                    continue
                    
                # Store the created plate
                created_plates[plate_id] = plate_obj
                
                # Track highest well index for this plate
                highest_well_idx = 0
                
                # Create wells for this plate group
                for _, row in plate_df.iterrows():
                    well_idx = row["well-index"]
                    
                    # Track highest index for plate's next available well
                    highest_well_idx = max(highest_well_idx, well_idx)
                    
                    # Create well for this material
                    well_obj = self.session.well_manager.create_well_model(
                        plate_obj=plate_obj,
                        role="startingmaterial",
                        role_index=1,
                        wellindex=well_idx,
                        volume=row["amount-uL"],
                        smiles=row["SMILES"],
                        concentration=row["concentration"] if "concentration" in row else None,
                        solvent=row["solvent"] if "solvent" in row else None,
                    )
                    
                    if not well_obj:
                        logger.warning(f"Plate {plate_id}: Could not create well at index {well_idx}")
                        continue
                
                # Update plate's next available well index
                self.session.well_manager.update_plate_well_index(
                    plate_obj=plate_obj, 
                    wellindexupdate=highest_well_idx + 1
                )
                
                # Create compound order for this plate
                self.session.data_manager.create_compound_order_model(
                    order_df=plate_df, 
                    is_custom_starter_plate=True
                )
                
                logger.info(f"Created custom starter plate {plate_id} with {len(plate_df)} wells")
                
            return created_plates
            
        except Exception as e:
            logger.error(f"Error creating starting material plates from CSV: {str(e)}")
            import traceback
            logger.error(traceback.format_exc())
            return {}
        
    def create_reaction_starting_plate(self):
        """
        Creates a starting material plate for reaction materials, respecting the maximum
        fill percentage defined in labware configuration.

        When ``self.session.use_multichannel`` is True, delegates to
        :meth:`_create_multichannel_starting_plate` which lays out
        reagents in full-column groups for 8-channel pipette transfers.

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

        # Branch: multichannel-aware layout or standard sequential layout
        if getattr(self.session, "use_multichannel", False):
            return self._create_multichannel_starting_plate(materials_df)

        return self._create_sequential_starting_plate(materials_df)

    # ------------------------------------------------------------------
    # Multichannel-aware starter plate layout
    # ------------------------------------------------------------------

    def _create_multichannel_starting_plate(self, materials_df):
        """Create a starting material plate with multichannel-aware column layout.

        Multichannel-eligible reagents (used in ≥ ``wells_per_column`` reactions
        with the same volume) are placed in full, uniform columns so that an
        8-channel pipette can aspirate from the entire column at once.  Each
        source well in a multichannel column holds enough volume for *N*
        transfers, where *N = ceil(reaction_count / wells_per_column)*.

        Remaining reagents (too few reactions, or leftover fractions) revert
        to the standard sequential approach: one well per reagent with the
        summed total volume.

        Parameters
        ----------
        materials_df : pd.DataFrame
            Output of ``MaterialManager.get_add_actions_material_dataframe()``.

        Returns
        -------
        plate_obj : Plate or None
        """
        # --- 1. Build per-transfer materials list from raw add actions ---
        addactionsdf = self.session.addactionsdf
        if addactionsdf is None or addactionsdf.empty:
            logger.warning("No add actions dataframe available for multichannel analysis")
            return self._create_sequential_starting_plate(materials_df)

        # Ensure uniquesolution column exists
        if "uniquesolution" not in addactionsdf.columns:
            addactionsdf["uniquesolution"] = addactionsdf.apply(
                self.session.material_manager.combine_strings, axis=1
            )

        # Group raw add actions by (uniquesolution, volume) to get
        # per-transfer volume and reaction counts
        raw_groups = (
            addactionsdf.groupby(["uniquesolution", "volume"])
            .agg(
                smiles=("smiles", "first"),
                solvent=("solvent", "first"),
                concentration=("concentration", "first"),
                reaction_count=("reaction_id_id", "nunique"),
            )
            .reset_index()
        )

        # Only include materials that still need wells (present in materials_df)
        materials_unique_solutions = set(materials_df["uniquesolution"].tolist())

        materials_for_analysis = []
        for _, row in raw_groups.iterrows():
            if row["uniquesolution"] not in materials_unique_solutions:
                continue
            materials_for_analysis.append(
                {
                    "smiles": row["smiles"],
                    "volume": row["volume"],
                    "reaction_count": int(row["reaction_count"]),
                    "solvent": row["solvent"] if pd.notna(row["solvent"]) else None,
                    "concentration": (
                        float(row["concentration"])
                        if pd.notna(row["concentration"])
                        else None
                    ),
                }
            )

        if not materials_for_analysis:
            logger.info("No materials available for multichannel analysis, "
                        "falling back to sequential layout")
            return self._create_sequential_starting_plate(materials_df)

        # --- 2. Determine plate labware ---
        volumes = materials_df["volume"].tolist()
        labware_type = self.session.labware_selector.get_plate_type(
            role="startingmaterial", volumes=volumes, temperature=25
        )

        plate_obj = self.create_plate_model(
            role="startingmaterial",
            role_index=1,
            platename="reaction-starting-materials-mc",
            labwaretype=labware_type,
        )
        if not plate_obj:
            logger.warning("Could not create starting material plate")
            return None

        wells_per_column = plate_obj.numberwellsincolumn
        num_columns = plate_obj.numbercolumns
        effective_max_volume = self.session.well_manager.get_max_well_volume(plate_obj)
        raw_max_volume = plate_obj.maxwellvolume
        dead_volume = self.session.labware_selector.get_dead_volume(raw_max_volume)

        # --- 3. Run MultichannelAnalyzer ---
        analyzer = MultichannelAnalyzer(
            wells_per_column=wells_per_column,
            num_columns=num_columns,
        )
        wells_per_group = analyzer.wells_per_group
        sub_columns_per_column = analyzer.sub_columns_per_column

        mc_materials, sc_materials = analyzer.prioritize_materials(
            materials_for_analysis
        )

        logger.info(
            f"Multichannel analysis: {len(mc_materials)} MC groups, "
            f"{len(sc_materials)} SC groups "
            f"(wells_per_group={wells_per_group}, "
            f"sub_columns_per_column={sub_columns_per_column})"
        )

        order_dicts_list = []
        next_column = 0       # Next available physical column
        next_sub_column = 0   # Next available sub-column within that column

        total_mc_slots = num_columns * sub_columns_per_column

        # --- 4. Phase A: Place multichannel materials as sub-column groups ---
        sc_leftover = []  # MC leftovers that fall back to sequential
        mc_slots_used = 0

        for mat in mc_materials:
            # Volume per well in the source sub-column.
            # Each channel serves ceil(reaction_count / wells_per_group)
            # transfers from its well.
            transfers_per_channel = math.ceil(
                mat.reaction_count / wells_per_group
            )
            per_well_volume = (
                mat.volume * transfers_per_channel * 1.15 + dead_volume
            )

            # If volume exceeds well capacity, we need multiple source groups
            # (or fall back to sequential for the excess).
            if per_well_volume > effective_max_volume:
                # Max transfers one well can support
                max_transfers = max(
                    1,
                    int((effective_max_volume - dead_volume) / (mat.volume * 1.15)),
                )
                groups_from_volume = math.ceil(
                    mat.reaction_count / (wells_per_group * max_transfers)
                )
                per_well_volume_capped = (
                    mat.volume * max_transfers * 1.15 + dead_volume
                )
                logger.info(
                    f"MC material {mat.smiles}: volume per well {per_well_volume:.1f}µL "
                    f"exceeds max {effective_max_volume:.1f}µL; "
                    f"using {groups_from_volume} groups with {max_transfers} transfers/channel"
                )
            else:
                groups_from_volume = 1
                per_well_volume_capped = per_well_volume

            actual_groups = max(groups_from_volume, 1)

            for _ in range(actual_groups):
                if mc_slots_used >= total_mc_slots:
                    logger.warning(
                        f"Plate full at column {next_column}; "
                        f"remaining MC reactions for {mat.smiles} → sequential"
                    )
                    # Remaining goes to single-channel
                    sc_leftover.append(mat)
                    break

                sub_col_well_indices = analyzer.get_sub_column_well_indices(
                    next_column, next_sub_column
                )

                for well_idx in sub_col_well_indices:
                    well_obj = self.session.well_manager.create_well_model(
                        plate_obj=plate_obj,
                        role="startingmaterial",
                        role_index=1,
                        wellindex=well_idx,
                        volume=round(per_well_volume_capped, 2),
                        smiles=mat.smiles,
                        concentration=mat.concentration,
                        solvent=mat.solvent,
                        transfer_type="multichannel",
                    )
                    if well_obj:
                        order_dicts_list.append(
                            {
                                "SMILES": mat.smiles,
                                "plate-name": plate_obj.name,
                                "labware": plate_obj.labware,
                                "well-index": well_obj.index,
                                "well-name": well_obj.name,
                                "concentration": mat.concentration,
                                "solvent": mat.solvent,
                                "molecularweight": None,
                                "amount-uL": round(per_well_volume_capped, 2),
                            }
                        )

                mc_slots_used += 1
                next_sub_column += 1
                if next_sub_column >= sub_columns_per_column:
                    # Finished all sub-columns for this physical column.
                    # Update plate indices past this physical column.
                    next_well_after_column = (next_column + 1) * wells_per_column
                    self.session.well_manager.update_plate_well_index(
                        plate_obj=plate_obj,
                        wellindexupdate=next_well_after_column,
                    )
                    self.session.column_manager.update_plate_column_index_available(
                        plate_obj=plate_obj,
                        columnindexupdate=next_column + 1,
                    )
                    next_sub_column = 0
                    next_column += 1

            # Leftovers (reactions that don't fill a full column) go to SC
            if mat.leftover_wells > 0:
                sc_leftover.append(
                    type(mat)(
                        identity_key=mat.identity_key,
                        smiles=mat.smiles,
                        volume=mat.volume,
                        concentration=mat.concentration,
                        solvent=mat.solvent,
                        reaction_count=mat.leftover_wells,
                        columns_needed=0,
                        leftover_wells=0,
                    )
                )

        # --- 4.5. Phase A.5: Heterogeneous MC sub-columns ---
        # For destination sub-columns where each reaction needs a
        # DIFFERENT material (e.g. 8 unique amines), create a MC source
        # sub-column with a different reagent at each position, matching
        # the destination well layout.  This lets the 8-channel pipette
        # transfer 8 different materials in one go.
        phase_a5_volume_served = {}  # identity_key → total volume served by MC
        phase_a5_solutions = set()
        try:
            phase_a5_volume_served, phase_a5_solutions = (
                self._create_heterogeneous_mc_sub_columns(
                    plate_obj=plate_obj,
                    addactionsdf=addactionsdf,
                    analyzer=analyzer,
                    mc_materials=mc_materials,
                    wells_per_column=wells_per_column,
                    wells_per_group=wells_per_group,
                    sub_columns_per_column=sub_columns_per_column,
                    effective_max_volume=effective_max_volume,
                    dead_volume=dead_volume,
                    order_dicts_list=order_dicts_list,
                    next_column=next_column,
                    next_sub_column=next_sub_column,
                    mc_slots_used=mc_slots_used,
                    total_mc_slots=total_mc_slots,
                )
            )
        except Exception as e:
            logger.warning(
                f"Phase A.5 heterogeneous MC failed, continuing: {e}"
            )

        # --- 5. Phase B: Place single-channel materials sequentially ---
        # Merge SC leftover + original SC materials.
        # For sequential placement, use one well per unique reagent with
        # the TOTAL volume (same as the existing approach).
        # Retrieve the *adjusted* total per unique solution from materials_df.
        sc_unique_solutions = set()

        # Add MC-leftover unique solutions (their total volume is a fraction)
        for mat in sc_leftover:
            sc_unique_solutions.add(
                analyzer.make_identity_key(mat.smiles, mat.solvent, mat.concentration)
            )
        # Add genuinely single-channel unique solutions
        for mat in sc_materials:
            sc_unique_solutions.add(
                analyzer.make_identity_key(mat.smiles, mat.solvent, mat.concentration)
            )

        # Process SC materials from materials_df (which has correct adjusted total volumes)
        for _, row in materials_df.iterrows():
            row_key = analyzer.make_identity_key(
                row["smiles"],
                row["solvent"] if pd.notna(row.get("solvent")) else None,
                float(row["concentration"]) if pd.notna(row.get("concentration")) else None,
            )
            if row_key not in sc_unique_solutions:
                continue

            # Subtract volume already served by Phase A.5 heterogeneous MC
            raw_volume = row["volume"]
            mc_served = phase_a5_volume_served.get(row_key, 0)
            adjusted_volume = raw_volume - mc_served
            if adjusted_volume <= 0:
                logger.info(
                    f"Material {row_key} fully served by heterogeneous MC "
                    f"— skipping SC well"
                )
                continue

            extra_error_volume = adjusted_volume * 0.15
            total_volume_needed = adjusted_volume + extra_error_volume + dead_volume

            well_idx = self.session.well_manager.get_plate_well_index_available(
                plate_obj
            )
            if well_idx is False:
                plate_obj = self.create_plate_model(
                    role="startingmaterial",
                    role_index=1,
                    platename="reaction-starting-materials-mc",
                    labwaretype=labware_type,
                )
                if not plate_obj:
                    logger.warning("Could not create additional starting material plate")
                    break
                well_idx = self.session.well_manager.get_plate_well_index_available(
                    plate_obj
                )
                effective_max_volume = self.session.well_manager.get_max_well_volume(
                    plate_obj
                )

            # Handle volume overflow — split into multiple wells
            if total_volume_needed > effective_max_volume:
                remaining_volume = total_volume_needed
                while remaining_volume > 0:
                    volume_this_well = min(remaining_volume, effective_max_volume)
                    well_obj = self.session.well_manager.create_well_model(
                        plate_obj=plate_obj,
                        role="startingmaterial",
                        role_index=1,
                        wellindex=well_idx,
                        volume=round(volume_this_well, 2),
                        smiles=row["smiles"],
                        concentration=row["concentration"],
                        solvent=row["solvent"],
                        transfer_type="single",
                    )
                    if well_obj:
                        order_dicts_list.append(
                            {
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
                        )
                    remaining_volume -= volume_this_well
                    self.session.well_manager.update_plate_well_index(
                        plate_obj=plate_obj, wellindexupdate=well_idx + 1
                    )
                    if remaining_volume > 0:
                        well_idx = (
                            self.session.well_manager.get_plate_well_index_available(
                                plate_obj
                            )
                        )
                        if well_idx is False:
                            plate_obj = self.create_plate_model(
                                role="startingmaterial",
                                role_index=1,
                                platename="reaction-starting-materials-mc",
                                labwaretype=labware_type,
                            )
                            if not plate_obj:
                                break
                            well_idx = (
                                self.session.well_manager.get_plate_well_index_available(
                                    plate_obj
                                )
                            )
            else:
                well_obj = self.session.well_manager.create_well_model(
                    plate_obj=plate_obj,
                    role="startingmaterial",
                    role_index=1,
                    wellindex=well_idx,
                    volume=round(total_volume_needed, 2),
                    smiles=row["smiles"],
                    concentration=row["concentration"],
                    solvent=row["solvent"],
                    transfer_type="single",
                )
                if well_obj:
                    order_dicts_list.append(
                        {
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
                    )
                self.session.well_manager.update_plate_well_index(
                    plate_obj=plate_obj, wellindexupdate=well_idx + 1
                )

        # --- 6. Create compound order ---
        self._create_compound_order(order_dicts_list)

        mc_well_count = mc_slots_used * wells_per_group
        sc_well_count = len(order_dicts_list) - mc_well_count
        logger.info(
            f"Created multichannel starting plate: {mc_slots_used} MC groups "
            f"({mc_well_count} wells), {sc_well_count} SC wells"
        )
        return plate_obj

    # ------------------------------------------------------------------
    # Heterogeneous MC sub-column creation (Phase A.5)
    # ------------------------------------------------------------------

    def _create_heterogeneous_mc_sub_columns(
        self,
        plate_obj,
        addactionsdf,
        analyzer,
        mc_materials,
        wells_per_column,
        wells_per_group,
        sub_columns_per_column,
        effective_max_volume,
        dead_volume,
        order_dicts_list,
        next_column,
        next_sub_column,
        mc_slots_used,
        total_mc_slots,
    ):
        """Create heterogeneous MC sub-columns for groups of different
        reagents that share a destination sub-column.

        For each destination sub-column of exactly ``wells_per_group``
        reactions, if the materials at a given add-action step are NOT
        all the same (i.e. not covered by a Phase A homogeneous MC
        sub-column), and all transfer volumes match, a MC source
        sub-column is created with a different reagent at each position
        — matching the position of the corresponding destination well.

        This enables the 8-channel pipette to transfer 8 different
        materials (e.g. unique amines) in a single operation.

        Returns
        -------
        tuple
            ``(volume_served, solutions_set)`` —
            ``volume_served`` is a dict mapping identity_key → total
            volume served by heterogeneous MC, and ``solutions_set``
            collects the identity keys of materials placed.
        """
        volume_served = {}  # identity_key → total volume served
        solutions_set = set()

        # Build reaction → dest well position mapping from reaction plate
        rxn_wells = Well.objects.filter(
            otsession_id=self.session.otsession_id,
            role="reaction",
            reaction_id__isnull=False,
        ).select_related("plate_id")

        if not rxn_wells.exists():
            logger.info("No reaction wells found for Phase A.5")
            return volume_served, solutions_set

        rxn_to_dest = {}  # reaction_id → (plate_id, col, sub_col, pos)
        for w in rxn_wells:
            wpc = w.plate_id.numberwellsincolumn
            sc_per_col = max(1, wpc // wells_per_group)
            col = w.index // wpc
            pos_in_col = w.index % wpc
            if sc_per_col > 1:
                sc = pos_in_col % sc_per_col
                pisc = pos_in_col // sc_per_col
            else:
                sc = 0
                pisc = pos_in_col
            rxn_to_dest[w.reaction_id_id] = (w.plate_id_id, col, sc, pisc)

        # Identity keys already covered by Phase A homogeneous MC
        phase_a_keys = {
            analyzer.make_identity_key(
                m.smiles, m.solvent, m.concentration
            )
            for m in mc_materials
        }

        # Group add actions by (action_number, dest_plate, col, sub_col)
        from collections import defaultdict
        dest_step_groups = defaultdict(dict)

        for _, row in addactionsdf.iterrows():
            rxn_id = row["reaction_id_id"]
            if rxn_id not in rxn_to_dest:
                continue
            _pid, col, sc, pisc = rxn_to_dest[rxn_id]
            action_num = row["number"]
            sol = row["solvent"] if pd.notna(row.get("solvent")) else None
            conc = (
                float(row["concentration"])
                if pd.notna(row.get("concentration"))
                else None
            )
            group_key = (action_num, _pid, col, sc)
            dest_step_groups[group_key][pisc] = {
                "smiles": row["smiles"],
                "volume": row["volume"],
                "solvent": sol,
                "concentration": conc,
            }

        # Identify heterogeneous MC candidates
        candidates = []
        for (action_num, _pid, _col, _sc), pos_map in sorted(
            dest_step_groups.items()
        ):
            if len(pos_map) < wells_per_group:
                continue  # partial sub-column

            # Volume uniformity check
            vols = {d["volume"] for d in pos_map.values()}
            if len(vols) > 1:
                continue

            # Skip if all materials are identical (Phase A handles those)
            materials = {
                (d["smiles"], d["solvent"], d["concentration"])
                for d in pos_map.values()
            }
            if len(materials) == 1:
                mat_key = analyzer.make_identity_key(
                    *next(iter(materials))
                )
                if mat_key in phase_a_keys:
                    continue

            # Skip if ALL materials are already covered by Phase A
            all_in_a = all(
                analyzer.make_identity_key(
                    d["smiles"], d["solvent"], d["concentration"]
                )
                in phase_a_keys
                for d in pos_map.values()
            )
            if all_in_a:
                continue

            vol = next(iter(vols))
            candidates.append(
                {"pos_map": pos_map, "volume": vol}
            )

        if not candidates:
            logger.info("No heterogeneous MC candidates found")
            return volume_served, solutions_set

        logger.info(
            f"Phase A.5: {len(candidates)} heterogeneous MC candidate(s)"
        )

        for cand in candidates:
            if mc_slots_used >= total_mc_slots:
                logger.warning("Plate full — stopping heterogeneous MC")
                break

            pos_map = cand["pos_map"]
            vol = cand["volume"]
            per_well_volume = vol * 1.15 + dead_volume

            if per_well_volume > effective_max_volume:
                logger.info(
                    "Heterogeneous MC well volume "
                    f"{per_well_volume:.1f} µL exceeds max "
                    f"{effective_max_volume:.1f} µL — skipping"
                )
                continue

            sub_col_well_indices = analyzer.get_sub_column_well_indices(
                next_column, next_sub_column
            )

            for pos in range(wells_per_group):
                if pos not in pos_map:
                    continue

                mat_data = pos_map[pos]
                well_idx = sub_col_well_indices[pos]

                well_obj = self.session.well_manager.create_well_model(
                    plate_obj=plate_obj,
                    role="startingmaterial",
                    role_index=1,
                    wellindex=well_idx,
                    volume=round(per_well_volume, 2),
                    smiles=mat_data["smiles"],
                    concentration=mat_data["concentration"],
                    solvent=mat_data["solvent"],
                    transfer_type="multichannel",
                )
                if well_obj:
                    order_dicts_list.append(
                        {
                            "SMILES": mat_data["smiles"],
                            "plate-name": plate_obj.name,
                            "labware": plate_obj.labware,
                            "well-index": well_obj.index,
                            "well-name": well_obj.name,
                            "concentration": mat_data["concentration"],
                            "solvent": mat_data["solvent"],
                            "molecularweight": None,
                            "amount-uL": round(per_well_volume, 2),
                        }
                    )

                    # Track served volume
                    ik = analyzer.make_identity_key(
                        mat_data["smiles"],
                        mat_data["solvent"],
                        mat_data["concentration"],
                    )
                    volume_served[ik] = volume_served.get(ik, 0) + vol
                    solutions_set.add(ik)

            mc_slots_used += 1
            next_sub_column += 1
            if next_sub_column >= sub_columns_per_column:
                next_well_after_column = (
                    (next_column + 1) * wells_per_column
                )
                self.session.well_manager.update_plate_well_index(
                    plate_obj=plate_obj,
                    wellindexupdate=next_well_after_column,
                )
                self.session.column_manager.update_plate_column_index_available(
                    plate_obj=plate_obj,
                    columnindexupdate=next_column + 1,
                )
                next_sub_column = 0
                next_column += 1

        logger.info(
            f"Phase A.5 complete: {len(solutions_set)} material(s) placed "
            f"in heterogeneous MC sub-columns"
        )
        return volume_served, solutions_set

    # ------------------------------------------------------------------
    # Standard sequential starter plate layout (original logic)
    # ------------------------------------------------------------------

    def _create_sequential_starting_plate(self, materials_df):
        """Standard sequential starting plate layout (one well per reagent).

        Parameters
        ----------
        materials_df : pd.DataFrame
            Output of ``MaterialManager.get_add_actions_material_dataframe()``.

        Returns
        -------
        plate_obj : Plate or None
        """

        # Get total volumes for plate sizing
        volumes = materials_df["volume"].tolist()

        # Determine best labware type
        labware_type = self.session.labware_selector.get_plate_type(
            role="startingmaterial", volumes=volumes, temperature=25
        )

        # Create the plate
        plate_obj = self.create_plate_model(
            role="startingmaterial",
            role_index=1,
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
            # Add 15% extra volume as error margin
            extra_error_volume = row["volume"] * 0.15
            total_volume_needed = row["volume"] + extra_error_volume + dead_volume

            # Get next available well index
            well_idx = self.session.well_manager.get_plate_well_index_available(
                plate_obj
            )

            # If no well is available, create a new plate
            if well_idx is False:
                plate_obj = self.create_plate_model(
                    role="startingmaterial",
                    role_index=1,
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
                        role="startingmaterial",
                        role_index=1,
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
                                role="startingmaterial",
                                role_index=1,
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
                    role="startingmaterial",
                    role_index=1,
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
        self._create_compound_order(order_dicts_list)

        logger.info(
            f"Created starting material plate with {len(order_dicts_list)} material wells"
        )
        return plate_obj

    # ------------------------------------------------------------------
    # Shared compound-order helper
    # ------------------------------------------------------------------

    def _create_compound_order(self, order_dicts_list):
        """Build and persist a CompoundOrder from a list of order dicts.

        Parameters
        ----------
        order_dicts_list : list[dict]
            Each dict has at minimum: ``SMILES``, ``plate-name``,
            ``labware``, ``well-index``, ``well-name``, ``concentration``,
            ``solvent``, ``molecularweight``, ``amount-uL``.
        """
        if not order_dicts_list:
            return

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

    def create_plates_by_temperature(
        self,
        grouped_reaction_temperature_querysets,
        *,
        role: str,
        role_index: int = 1,
    ):
        """
        Creates reaction plates for each temperature group using appropriate labware types.

        This function processes each temperature group, determines the appropriate labware type
        based on volumes and temperature, and creates plates for each group.

        Parameters
        ----------
        grouped_reaction_temperature_querysets: dict
            Dictionary mapping temperature to QuerySets of reactions at that temperature
        role: str
            The plate role (e.g., "reaction", "workup").
        role_index: int
            The plate role index (default 1).

        Returns
        -------
        created_plates: list
            List of created plate objects
        """
        # Generate platetype string for logging/naming
        platetype_str = unparse_plate_type(role, role_index)

        logger.info(
            f"Creating {platetype_str} plates for {len(grouped_reaction_temperature_querysets)} temperature groups"
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
                role=role,
                volumes=volumes,
                temperature=reaction_temperature,
            )

            logger.info(
                f"Creating {platetype_str} plate for temperature {reaction_temperature}°C using {labware_platetype}"
            )

            # CREATE PLATE DIRECTLY BY TEMPERATURE - NOT by class/recipe
            plate_name = f"{platetype_str}-{reaction_temperature}C"
            
            # Create the plate
            plate_obj = self.create_plate_model(
                role=role,
                role_index=role_index,
                platename=plate_name,
                labwaretype=labware_platetype,
            )

            if plate_obj:
                created_plates.append(plate_obj)

                # NOW GROUP REACTIONS BY CLASS WITHIN THIS TEMPERATURE PLATE
                # Group reactions by class for column organization
                reaction_class_groups = {}
                for reaction in reaction_temp_queryset:
                    reaction_class = getattr(reaction, 'reactionclass', 'unknown')
                    if reaction_class not in reaction_class_groups:
                        reaction_class_groups[reaction_class] = []
                    reaction_class_groups[reaction_class].append(reaction)

                # Create columns for each reaction class within the temperature plate
                for reaction_class, class_reactions in reaction_class_groups.items():
                    column_index = (
                        self.session.column_manager.get_plate_current_column_index(
                            plate_obj=plate_obj
                        )
                    )
                    if column_index is not False:
                        column_obj = self.session.column_manager.create_column_model(
                            plate_obj=plate_obj,
                            columnindex=column_index,
                            role=role,
                            role_index=role_index,
                            reactionclass=reaction_class,
                        )

                        # Update column index
                        self.session.column_manager.update_plate_column_index_available(
                            plate_obj=plate_obj, columnindexupdate=column_index + 1
                        )

                        # Add wells for each reaction in this class
                        for reaction in class_reactions:
                            # Get product SMILES for this reaction
                            product_smiles = (
                                self.session.material_manager.get_product_smiles(
                                    [reaction.id]
                                )[0]
                            )

                            logger.debug(
                                f"Creating well for reaction {reaction.id} (class: {reaction_class}, "
                                f"temp: {reaction_temperature}°C) with product SMILES: {product_smiles}"
                            )

                            # Check if well index is available
                            well_index = (
                                self.session.well_manager.get_plate_well_index_available(
                                    plate_obj=plate_obj
                                )
                            )
                            
                            if well_index is False:
                                # Create a new plate for the same temperature
                                new_plate_name = f"{plate_name}-continued"
                                plate_obj = self.create_plate_model(
                                    role=role,
                                    role_index=role_index,
                                    platename=new_plate_name,
                                    labwaretype=labware_platetype,
                                )
                                
                                if not plate_obj:
                                    logger.error("Failed to create additional plate")
                                    break
                                    
                                created_plates.append(plate_obj)
                                
                                well_index = self.session.well_manager.get_plate_well_index_available(
                                    plate_obj=plate_obj
                                )
                                if well_index is False:
                                    logger.error("New plate has no available wells")
                                    break

                            # Create well directly (you may want to add column management here)
                            well_obj = self.session.well_manager.create_well_model(
                                plate_obj=plate_obj,
                                role=role,
                                role_index=role_index,
                                wellindex=well_index,
                                reactionobj=reaction,
                                smiles=product_smiles,
                                reactantfornextstep=True,
                            )

                            # Update well index
                            self.session.well_manager.update_plate_well_index(
                                plate_obj=plate_obj, wellindexupdate=well_index + 1
                            )

        logger.info(f"Created {len(created_plates)} {platetype_str} plates by temperature")
        # return created_plates

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
                previous_type = plate_obj.role
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
