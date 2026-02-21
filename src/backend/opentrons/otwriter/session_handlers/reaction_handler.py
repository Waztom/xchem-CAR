"""
Reaction Session Handler for OpenTrons protocols.

This module provides methods for handling reaction sessions.
"""

import logging
import math
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from django.db.models import QuerySet

from backend.models import AddAction, MixAction, RecipeAddAction, RecipeMixAction
from backend.db_utils import get_reaction, get_previous_reaction_query_sets, get_reaction_query_set
from backend.recipe_utils import get_session_recipe_actions

from .base_handler import SessionHandler

logger = logging.getLogger(__name__)


# Number of tips on an 8-channel multi-channel pipette.
_MC_CHANNELS = 8


class ReactionSessionHandler(SessionHandler):
    """
    Handles reaction session processing.

    This class implements the logic for generating reaction session commands.
    """

    def process_session(self, actionsession_queryset: QuerySet) -> None:
        """
        Process the reaction session(s).

        When a multichannel pipette is configured and multichannel source
        wells exist, processing switches to **step-wise** mode.  Instead of
        completing all add actions per-reaction before moving to the next
        reaction, it iterates through the recipe steps *in order* and, for
        each step, decides whether to use multichannel (MC) or
        single-channel (SC) transfers:

        * MC is used for add steps where **all reactions share the same
          reagent** (SMILES + solvent + concentration) at the same volume,
          AND the destination sub-column is fully occupied (8 wells for an
          8-channel pipette).
        * SC is used for all other cases (unique reagents per reaction,
          partial sub-columns, or no MC source wells).

        This preserves add-action ordering while maximising multichannel
        efficiency.

        Parameters
        ----------
        actionsession_queryset : QuerySet
            QuerySet of reaction action sessions to process
        """
        # Use base helper to log session start
        self.log_session_start(actionsession_queryset, "reaction")

        # Get session information
        session_number = self.get_session_number(actionsession_queryset)
        reaction_step = self.script_generator.reactionstep
        logger.info(
            f"Processing reaction step {reaction_step}, session {session_number}"
        )

        # If this is a later reaction step, add dilution processing
        if reaction_step > 1:
            logger.info(
                f"This is reaction step {reaction_step}, adding dilution processing"
            )
            self.process_dilution_step(session_number)

            # Add a pause after dilution
            logger.info("Adding pause after dilution")
            self.add_command(
                self.command_generator.pause_protocol(
                    message="Addition of dilution solvent complete. Confirm dilution complete to restart protocol."
                )
            )

        # --- Decide between step-wise (MC-aware) or per-reaction (SC-only) ---
        mc_pipette_name = self.script_generator.mc_pipettename
        use_step_wise = False
        mc_wells_qs = None

        if mc_pipette_name:
            mc_wells_qs = self.query_service.get_multichannel_source_wells()
            if mc_wells_qs and mc_wells_qs.exists():
                use_step_wise = True
                logger.info(
                    f"MC pipette '{mc_pipette_name}' and {mc_wells_qs.count()} "
                    f"MC source wells detected — using step-wise processing"
                )

        if use_step_wise:
            self._process_session_step_wise(
                actionsession_queryset, session_number, mc_wells_qs
            )
        else:
            self._process_session_per_reaction(
                actionsession_queryset, session_number
            )

        self.log_session_end("reaction")

    # ------------------------------------------------------------------
    # Processing mode A: Original per-reaction flow (single-channel only)
    # ------------------------------------------------------------------

    def _process_session_per_reaction(
        self, actionsession_queryset: QuerySet, session_number: int
    ) -> None:
        """Process reactions one at a time with all actions per-reaction.

        This is the original (non-MC) flow: for each reaction, execute
        actions 1, 2, 3 … sequentially before moving to the next reaction.
        """
        action_count = actionsession_queryset.count()
        logger.info(
            f"Per-reaction processing: {action_count} reaction action session(s)"
        )

        for i, actionsession_obj in enumerate(actionsession_queryset):
            reaction_id = actionsession_obj.reaction_id.id
            logger.info(f"Processing reaction {reaction_id} ({i+1}/{action_count})")
            reaction_obj = get_reaction(reaction_id=reaction_id)

            # --- Human-readable header for this reaction ---
            self.add_command(
                self.annotation_generator.reaction_header(
                    reaction_obj, session_label="Reaction"
                )
            )

            self.process_reaction_actions(
                actionsession_obj, reaction_obj, session_number,
            )

    # ------------------------------------------------------------------
    # Processing mode B: Step-wise flow (multichannel-aware)
    # ------------------------------------------------------------------

    def _process_session_step_wise(
        self,
        actionsession_queryset: QuerySet,
        session_number: int,
        mc_wells_qs: QuerySet,
    ) -> None:
        """Process add actions step-by-step across all reactions.

        Instead of completing every action for reaction 1, then reaction 2,
        etc., this method processes one recipe **step** at a time across all
        reactions.  For each step it decides:

        * **MC-eligible** — all reactions share the same reagent for this
          step and MC source wells exist → full sub-columns use the
          8-channel pipette, partial sub-columns fall back to SC.
        * **SC-only** — the reagent varies per reaction → every reaction
          is transferred individually.

        This preserves the add-action ordering (step 1 before step 2 before
        step 3) while maximising multichannel efficiency.
        """
        mc_source_map = self._build_mc_source_map(mc_wells_qs)

        # Get the canonical recipe action list from the first reaction.
        # All reactions in the session are expected to share the same recipe.
        first_as = actionsession_queryset.first()
        first_rxn = get_reaction(reaction_id=first_as.reaction_id.id)
        reaction_actions = get_session_recipe_actions(
            reaction_class=first_rxn.reactionclass,
            name=first_rxn.recipe,
            session_type="reaction",
            session_number=session_number,
            molecular_context=(
                "intramolecular" if first_rxn.intramolecular else "intermolecular"
            ),
        )

        total_actions = len(reaction_actions)
        logger.info(
            f"Step-wise processing: {total_actions} recipe steps, "
            f"{actionsession_queryset.count()} reactions"
        )

        for step_index, reaction_action in enumerate(reaction_actions):
            if isinstance(reaction_action, RecipeAddAction):
                action_number = reaction_action.action_number

                # Fetch all AddActions for this step number across reactions
                all_add_actions = list(
                    AddAction.objects.filter(
                        actionsession_id__in=actionsession_queryset,
                        number=action_number,
                    )
                    .select_related("reaction_id", "actionsession_id")
                    .order_by("reaction_id__id")
                )

                if not all_add_actions:
                    logger.warning(
                        f"No AddActions for step {action_number} – skipping"
                    )
                    continue

                mc_key = self._check_step_homogeneous(all_add_actions)

                if mc_key and mc_key in mc_source_map:
                    self._process_step_with_mc(
                        all_add_actions,
                        mc_key,
                        mc_source_map,
                        reaction_action,
                        step_index,
                        total_actions,
                        session_number,
                    )
                else:
                    self._process_step_single_channel(
                        all_add_actions,
                        reaction_action,
                        step_index,
                        reaction_actions,
                        session_number,
                    )

            elif isinstance(reaction_action, RecipeMixAction):
                action_number = reaction_action.action_number
                for actionsession_obj in actionsession_queryset:
                    reaction_obj = get_reaction(
                        reaction_id=actionsession_obj.reaction_id.id
                    )
                    self.process_mix_action(
                        action_number,
                        actionsession_obj,
                        reaction_obj,
                        reaction_obj.id,
                    )
            else:
                logger.warning(
                    f"Unknown action type at step {step_index}: "
                    f"{type(reaction_action).__name__}"
                )

    # ------------------------------------------------------------------
    # Step-wise helper: build MC source well lookup
    # ------------------------------------------------------------------

    def _build_mc_source_map(
        self, mc_wells_qs: QuerySet
    ) -> Dict[Tuple, List[dict]]:
        """Build a lookup of multichannel source well groups.

        Groups MC-tagged source wells by their material identity
        ``(smiles, solvent, concentration)`` and source sub-column position.

        Returns
        -------
        dict
            ``{(smiles, solvent, concentration): [
                {'plate': Plate, 'col': int, 'sub_col': int, 'wells': [Well, …]},
                …
            ]}``
        """
        # Group wells by (material_key, plate_id, col, sub_col)
        temp: Dict[tuple, list] = defaultdict(list)
        for well in mc_wells_qs:
            plate = well.plate_id
            wpc = plate.numberwellsincolumn
            sub_cols = max(1, wpc // _MC_CHANNELS)

            col = well.index // wpc
            pos = well.index % wpc
            sub_col = pos % sub_cols if sub_cols > 1 else 0

            mat_key = (well.smiles, well.solvent, well.concentration)
            group_key = (mat_key, plate.id, col, sub_col)
            temp[group_key].append(well)

        mc_map: Dict[Tuple, List[dict]] = defaultdict(list)
        for (mat_key, _plate_id, col, sub_col), wells in temp.items():
            mc_map[mat_key].append(
                {"plate": wells[0].plate_id, "col": col, "sub_col": sub_col, "wells": wells}
            )

        logger.info(
            f"MC source map built: {len(mc_map)} material(s), "
            f"{sum(len(v) for v in mc_map.values())} sub-column group(s)"
        )
        return dict(mc_map)

    # ------------------------------------------------------------------
    # Step-wise helper: check if a step has homogeneous material
    # ------------------------------------------------------------------

    @staticmethod
    def _check_step_homogeneous(
        all_add_actions: list,
    ) -> Optional[Tuple[str, Optional[str], Optional[float]]]:
        """Return the material key if every AddAction shares the same reagent.

        Parameters
        ----------
        all_add_actions : list[AddAction]
            All AddActions for one recipe step, across all reactions.

        Returns
        -------
        tuple or None
            ``(smiles, solvent, concentration)`` when all actions share
            the same material; ``None`` otherwise.
        """
        first = all_add_actions[0]
        ref = (first.smiles, first.solvent, first.concentration)
        for aa in all_add_actions[1:]:
            if (aa.smiles, aa.solvent, aa.concentration) != ref:
                return None
        return ref

    # ------------------------------------------------------------------
    # Step-wise helper: MC-eligible step (full sub-cols MC, partials SC)
    # ------------------------------------------------------------------

    def _process_step_with_mc(
        self,
        all_add_actions: list,
        mc_key: Tuple,
        mc_source_map: Dict,
        reaction_action: RecipeAddAction,
        step_index: int,
        total_actions: int,
        session_number: int,
    ) -> None:
        """Process an MC-eligible add step.

        Full destination sub-columns (exactly ``_MC_CHANNELS`` reactions)
        are transferred with the multichannel pipette.  Partial
        sub-columns fall back to single-channel.

        Parameters
        ----------
        all_add_actions : list[AddAction]
            All AddActions for this step.
        mc_key : tuple
            ``(smiles, solvent, concentration)`` identifying the shared material.
        mc_source_map : dict
            Lookup built by :meth:`_build_mc_source_map`.
        reaction_action : RecipeAddAction
            The recipe-level action template for this step.
        step_index : int
            0-based index of this step within the recipe.
        total_actions : int
            Total recipe steps.
        session_number : int
            Current session number.
        """
        mc_pipette_name = self.script_generator.mc_pipettename
        smiles, solvent, concentration = mc_key
        volume = all_add_actions[0].volume

        mc_sources = mc_source_map[mc_key]

        # --- Step banner ---
        self.add_command(
            self.annotation_generator.multichannel_step_header(
                step_index + 1, total_actions, smiles, solvent, volume
            )
        )

        # Group destination wells by (plate_id, physical_column, sub_column)
        dest_groups: Dict[Tuple, List[Tuple]] = defaultdict(list)
        sc_fallback: List = []

        for aa in all_add_actions:
            rxn_id = aa.reaction_id.id
            try:
                to_well = self.well_finder.find_reaction_well(
                    reaction_id=rxn_id,
                    role=reaction_action.to_plate_role,
                    role_index=reaction_action.to_plate_role_index,
                )
            except Exception:
                logger.warning(
                    f"No reaction well for reaction {rxn_id} – SC fallback"
                )
                sc_fallback.append(aa)
                continue

            to_plate = to_well.plate_id
            wpc = to_plate.numberwellsincolumn
            sub_cols = max(1, wpc // _MC_CHANNELS)

            dest_col = to_well.index // wpc
            pos = to_well.index % wpc
            dest_sub_col = pos % sub_cols if sub_cols > 1 else 0

            dest_groups[(to_plate.id, dest_col, dest_sub_col)].append(
                (aa, to_well)
            )

        # Process each destination sub-column group
        for (dest_plate_id, dest_col, dest_sub_col), group in sorted(
            dest_groups.items()
        ):
            dest_plate = self.query_service.get_plate_by_id(plateid=dest_plate_id)

            if len(group) < _MC_CHANNELS:
                # Partial sub-column → SC fallback
                logger.info(
                    f"Partial sub-column: {dest_plate.name} col {dest_col} "
                    f"sub-col {dest_sub_col} ({len(group)}/{_MC_CHANNELS}) → SC"
                )
                sc_fallback.extend([aa for aa, _w in group])
                continue

            # Find a matching MC source sub-column
            src_group = self._find_mc_source_for_sub_col(
                mc_sources, dest_sub_col
            )
            if src_group is None:
                logger.warning(
                    f"No MC source sub-column {dest_sub_col} for {smiles[:30]}… – SC"
                )
                sc_fallback.extend([aa for aa, _w in group])
                continue

            src_plate = src_group["plate"]
            src_col = src_group["col"]
            src_sub_col = src_group["sub_col"]

            logger.info(
                f"MC transfer step {step_index + 1}: "
                f"{src_plate.name} col {src_col} sub-col {src_sub_col} → "
                f"{dest_plate.name} col {dest_col} sub-col {dest_sub_col} "
                f"({volume:.1f} µL, {len(group)} reactions)"
            )

            # Pick up MC tip
            self.add_command(self.command_generator.pick_up_tip(suffix="MC"))

            # Column-to-column transfer with correct sub-column indices
            self.add_command(
                self.command_generator.transfer_fluid_multi(
                    aspirateplatename=src_plate.name,
                    dispenseplatename=dest_plate.name,
                    aspiratecolumnindex=src_col,
                    dispensecolumnindex=dest_col,
                    transvolume=volume,
                    pipette_name_override=mc_pipette_name,
                    aspirate_sub_column_index=src_sub_col,
                    dispense_sub_column_index=dest_sub_col,
                )
            )

            # Record each transfer in the ledger
            for aa, to_well in group:
                rxn_id = aa.reaction_id.id
                src_well_index = self._mc_source_well_for_dest(
                    src_plate, src_col, to_well, dest_plate
                )
                src_well = self.query_service.get_well_by_plate_and_index(
                    plate_id=src_plate, well_index=src_well_index,
                )

                self.script_generator.transfer_ledger.record(
                    action_type="add",
                    source_plate_name=src_plate.name,
                    source_plate_role="startingmaterial",
                    source_well_index=src_well_index,
                    source_well_name=(
                        (getattr(src_well, "name", "") or "") if src_well else ""
                    ),
                    dest_plate_name=dest_plate.name,
                    dest_plate_role="reaction",
                    dest_well_index=to_well.index,
                    dest_well_name=getattr(to_well, "name", "") or "",
                    volume=volume,
                    smiles=smiles,
                    solvent=solvent,
                    reaction_id=rxn_id,
                    reaction_class=getattr(
                        aa.reaction_id, "reactionclass", None
                    ),
                    recipe=getattr(aa.reaction_id, "recipe", None),
                    transfer_mode="multichannel",
                )

                self.volume_manager.update_well_reactant_status(to_well, True)
                if src_well:
                    self.volume_manager.update_well_volume(
                        wellobj=src_well, transfervolume=volume
                    )

            # Drop MC tip
            self.add_command(self.command_generator.drop_tip(suffix="MC"))

        # ---- SC fallback for partial sub-columns ----
        if sc_fallback:
            logger.info(
                f"Processing {len(sc_fallback)} SC fallback transfers "
                f"for step {step_index + 1}"
            )
            self._process_step_single_channel(
                sc_fallback,
                reaction_action,
                step_index,
                # Fetch canonical reaction_actions for tip logic
                self._get_canonical_reaction_actions(
                    sc_fallback[0], session_number
                ),
                session_number,
            )

    # ------------------------------------------------------------------
    # Step-wise helper: SC-only step
    # ------------------------------------------------------------------

    def _process_step_single_channel(
        self,
        add_actions: list,
        reaction_action: RecipeAddAction,
        step_index: int,
        reaction_actions: list,
        session_number: int,
    ) -> None:
        """Process an add step with single-channel transfers for each reaction.

        Parameters
        ----------
        add_actions : list[AddAction]
            AddActions for this step (may be a subset for MC partial fallback).
        reaction_action : RecipeAddAction
            The recipe-level action template.
        step_index : int
            0-based index of this step within the recipe.
        reaction_actions : list
            Full recipe action list (used for next-is-mix tip logic).
        session_number : int
            Current session number.
        """
        # Emit a SC step banner when called directly (not as MC fallback)
        total_actions = len(reaction_actions)

        for aa in add_actions:
            rxn_id = aa.reaction_id.id
            actionsession_obj = aa.actionsession_id
            reaction_obj = aa.reaction_id

            # Reaction header + add action via the existing SC method
            self.add_command(
                self.annotation_generator.reaction_header(
                    reaction_obj, session_label="Reaction"
                )
            )
            self.process_add_action(
                reaction_action,
                aa.number,
                step_index,
                reaction_actions,
                actionsession_obj,
                reaction_obj,
                rxn_id,
            )

    # ------------------------------------------------------------------
    # Step-wise helper: find a MC source group matching a sub-column
    # ------------------------------------------------------------------

    @staticmethod
    def _find_mc_source_for_sub_col(
        mc_sources: List[dict], target_sub_col: int
    ) -> Optional[dict]:
        """Return the first MC source group whose sub-column matches."""
        for sg in mc_sources:
            if sg["sub_col"] == target_sub_col:
                return sg
        return None

    # ------------------------------------------------------------------
    # Step-wise helper: map dest well → source well for MC ledger
    # ------------------------------------------------------------------

    @staticmethod
    def _mc_source_well_for_dest(
        src_plate, src_col: int, dest_well, dest_plate
    ) -> int:
        """Compute the source well index that is paired with a dest well.

        In a multichannel transfer the tip at physical position *N*
        aspirates from source row *N* and dispenses to destination row
        *N*.  For same-type plates this is a direct position mapping;
        for cross-type (e.g. 96 → 384) the row pitch ratio is applied.
        """
        src_wpc = src_plate.numberwellsincolumn
        dest_wpc = dest_plate.numberwellsincolumn
        dest_pos = dest_well.index % dest_wpc

        if src_wpc == dest_wpc:
            src_pos = dest_pos
        elif dest_wpc > src_wpc:
            # e.g. 96-well source → 384-well dest
            ratio = dest_wpc // src_wpc
            src_pos = dest_pos // ratio
        else:
            ratio = src_wpc // dest_wpc
            src_pos = dest_pos * ratio

        return src_col * src_wpc + src_pos

    # ------------------------------------------------------------------
    # Step-wise helper: canonical reaction actions for a given AddAction
    # ------------------------------------------------------------------

    def _get_canonical_reaction_actions(self, add_action, session_number):
        """Return the recipe action list for the reaction owning *add_action*."""
        rxn = add_action.reaction_id
        return get_session_recipe_actions(
            reaction_class=rxn.reactionclass,
            name=rxn.recipe,
            session_type="reaction",
            session_number=session_number,
            molecular_context=(
                "intramolecular" if rxn.intramolecular else "intermolecular"
            ),
        )

    def process_dilution_step(self, session_number: int) -> None:
        """
        Process dilution for a reaction step.

        Parameters
        ----------
        session_number : int
            The session number
        """
        logger.info(f"Starting dilution processing for session {session_number}")

        # --- Human-readable dilution header ---
        self.add_command(
            self.annotation_generator.dilution_header()
        )

        try:
            # Get reaction data
            reaction_queryset = get_reaction_query_set(
                reaction_ids=self.script_generator.reaction_ids
            )
            reaction_count = reaction_queryset.count()
            logger.info(f"Found {reaction_count} reactions for dilution processing")

            grouped_reaction_querysets = (
                self.query_service.get_grouped_reaction_by_class_recipe(
                    reactionqueryset=reaction_queryset
                )
            )

            group_count = len(grouped_reaction_querysets)
            logger.info(
                f"Reactions grouped into {group_count} sets based on class/recipe"
            )

            # Process each group of reactions
            for i, group_reaction_class_queryset in enumerate(
                grouped_reaction_querysets
            ):
                group_size = group_reaction_class_queryset.count()
                logger.info(
                    f"Processing reaction group {i+1}/{group_count} with {group_size} reactions"
                )

                logger.info("Picking up tip for dilution transfers")
                self.add_command(self.command_generator.pick_up_tip())

                action_session_type = "reaction"
                reaction_classes = group_reaction_class_queryset.values_list(
                    "reactionclass", flat=True
                )

                if not all(reaction_classes):
                    logger.error("Inconsistent reaction classes in grouped queryset")
                    raise Exception(
                        "Expected one type of reaction class in grouped class queryset"
                    )

                reaction_class = reaction_classes[0]
                logger.info(f"Processing dilution for reaction class: {reaction_class}")

                recipe_types = group_reaction_class_queryset.values_list(
                    "recipe", flat=True
                )
                recipe_count = len(set(recipe_types))
                logger.info(f"Found {recipe_count} unique recipe type(s)")

                if not all(recipe_types):
                    logger.info(
                        "Multiple recipe types detected, processing individually"
                    )
                    self.process_multiple_recipe_types(
                        group_reaction_class_queryset,
                        recipe_types,
                        reaction_class,
                        action_session_type,
                        session_number,
                    )
                else:
                    # All reactions have the same recipe
                    recipe_type = recipe_types[0]
                    logger.info(f"All reactions share recipe type: {recipe_type}")
                    self.process_single_recipe_type(
                        group_reaction_class_queryset,
                        reaction_class,
                        recipe_type,
                        action_session_type,
                        session_number,
                    )

                logger.info("Dropping tip after dilution transfers")
                self.add_command(self.command_generator.drop_tip())

            logger.info("Dilution processing completed")

        except Exception as e:
            logger.error(f"Error processing dilution step: {str(e)}")
            raise

    def process_multiple_recipe_types(
        self,
        group_reaction_class_queryset,
        recipe_types,
        reaction_class,
        action_session_type,
        session_number,
    ):
        """Process multiple recipe types for dilution step."""
        unique_recipe_types = set(recipe_types)
        logger.info(
            f"Processing {len(unique_recipe_types)} different recipe types for dilution"
        )

        for i, recipe_type in enumerate(unique_recipe_types):
            if not recipe_type:
                logger.warning(f"Empty recipe type found, skipping")
                continue

            logger.info(
                f"Processing recipe type {i+1}/{len(unique_recipe_types)}: {recipe_type}"
            )

            reaction_subset = group_reaction_class_queryset.filter(recipe=recipe_type)
            reaction_count = reaction_subset.count()
            logger.info(f"Found {reaction_count} reactions with recipe {recipe_type}")

            try:
                intramolecular = reaction_subset.values_list(
                    "intramolecular", flat=True
                ).distinct()[0]

                reaction_action_search = (
                    "intramolecular" if intramolecular else "intermolecular"
                )
                logger.info(
                    f"Using {reaction_action_search} actions for recipe {recipe_type}"
                )

                # Get actions from Recipe DB
                try:
                    reaction_actions = get_session_recipe_actions(
                        reaction_class=reaction_class,
                        name=recipe_type,
                        session_type=action_session_type,
                        session_number=session_number,
                        molecular_context=reaction_action_search,
                    )

                    action_count = len(reaction_actions)
                    logger.info(
                        f"Found {action_count} actions for {reaction_action_search} recipe {recipe_type}"
                    )

                    # Get add actions that have a material SMARTS
                    reaction_add_actions = [
                        action
                        for action in reaction_actions
                        if isinstance(action, RecipeAddAction) and action.material_smarts
                    ]

                    add_action_count = len(reaction_add_actions)
                    logger.info(
                        f"Found {add_action_count} add actions with material SMARTS"
                    )

                    self.process_add_actions_for_dilution(
                        reaction_subset, reaction_add_actions, action_session_type
                    )
                except (KeyError, IndexError) as e:
                    logger.error(
                        f"Error accessing recipe data for {recipe_type}: {str(e)}"
                    )
            except IndexError:
                logger.error(f"No intramolecular value found for recipe {recipe_type}")

    def process_single_recipe_type(
        self,
        group_reaction_class_queryset,
        reaction_class,
        recipe_type,
        action_session_type,
        session_number,
    ):
        """Process a single recipe type for dilution step."""
        try:
            intramolecular = group_reaction_class_queryset.values_list(
                "intramolecular", flat=True
            ).distinct()[0]

            reaction_action_search = (
                "intramolecular" if intramolecular else "intermolecular"
            )
            logger.info(
                f"Using {reaction_action_search} actions for recipe {recipe_type}"
            )

            # Get actions from Recipe DB
            try:
                reaction_actions = get_session_recipe_actions(
                    reaction_class=reaction_class,
                    name=recipe_type,
                    session_type=action_session_type,
                    session_number=session_number,
                    molecular_context=reaction_action_search,
                )

                action_count = len(reaction_actions)
                logger.info(
                    f"Found {action_count} actions for {reaction_action_search} recipe {recipe_type}"
                )

                # Get add actions that have a material SMARTS
                reaction_add_actions = [
                    action
                    for action in reaction_actions
                    if isinstance(action, RecipeAddAction) and action.material_smarts
                ]

                add_action_count = len(reaction_add_actions)
                logger.info(
                    f"Found {add_action_count} add actions with material SMARTS"
                )

                self.process_add_actions_for_dilution(
                    group_reaction_class_queryset,
                    reaction_add_actions,
                    action_session_type,
                )
            except (KeyError, IndexError) as e:
                logger.error(f"Error accessing recipe data for {recipe_type}: {str(e)}")
        except IndexError:
            logger.error(f"No intramolecular value found for recipe {recipe_type}")

    def process_add_actions_for_dilution(
        self, group_reaction_class_queryset, reaction_add_actions, action_session_type
    ):
        """Process add actions for dilution step."""
        reaction_step = self.script_generator.reactionstep
        add_action_count = len(reaction_add_actions)
        logger.info(f"Processing {add_action_count} add actions for dilution")

        for i, reaction_add_action in enumerate(reaction_add_actions):
            action_number = reaction_add_action.action_number
            logger.info(
                f"Processing add action {i+1}/{add_action_count}, number {action_number}"
            )

            # Get add actions for this step
            add_action_queryset = self.query_service.get_add_action_query_set(
                reaction_ids=group_reaction_class_queryset,
                actionsessiontype=action_session_type,
                actionnumber=action_number,
            )

            db_action_count = add_action_queryset.count()
            logger.info(
                f"Found {db_action_count} database add actions for action number {action_number}"
            )

            for j, add_action_obj in enumerate(add_action_queryset):
                reaction_obj = add_action_obj.reaction_id
                reaction_id = reaction_obj.id
                smiles = add_action_obj.smiles

                logger.info(
                    f"Processing add action {j+1}/{db_action_count} for reaction {reaction_id}"
                )
                logger.info(f"SMILES: {smiles[:20]}...")

                # Check if there's a previous reaction for this SMILES
                previous_reaction_queryset = get_previous_reaction_query_sets(
                    reaction_id=reaction_obj.id, smiles=smiles
                )

                prev_reaction_count = (
                    len(previous_reaction_queryset) if previous_reaction_queryset else 0
                )

                if previous_reaction_queryset:
                    logger.info(
                        f"Found {prev_reaction_count} previous reaction(s) for this material"
                    )
                    solvent = add_action_obj.solvent
                    transfer_volume = add_action_obj.volume
                    concentration = add_action_obj.concentration

                    logger.info(
                        f"Dilution parameters: solvent={solvent}, volume={transfer_volume} µL, concentration={concentration}"
                    )

                    # Find starting material wells
                    logger.info("Finding source wells for material")
                    from_well_info = self.well_finder.find_starting_material_wells(
                        reaction_step_no=reaction_step,
                        reaction_id=reaction_obj.id,
                        smiles=smiles,
                        solvent=solvent,
                        concentration=concentration,
                        transfer_volume=transfer_volume,
                    )

                    from_well_count = len(from_well_info)
                    if from_well_count == 0:
                        logger.warning(
                            f"No source wells found for reaction {reaction_id}"
                        )
                        continue

                    logger.info(f"Found {from_well_count} source well(s) for material")

                    # Process each well
                    for k, well_info in enumerate(from_well_info):
                        previous_reaction_objs = well_info[0]
                        from_well_obj = well_info[1]
                        transfer_volume = well_info[2]

                        from_well_index = from_well_obj.index
                        from_plate_id = from_well_obj.plate_id.id
                        logger.info(
                            f"Processing source well {k+1}/{from_well_count}: plate {from_plate_id}, well {from_well_index}"
                        )

                        if previous_reaction_objs:
                            # Find solvent well for dilution
                            logger.info(
                                f"Finding solvent wells for {solvent}, volume={transfer_volume} µL"
                            )
                            from_solvent_well_info = (
                                self.well_finder.find_solvent_wells(
                                    solvent=solvent,
                                    transfer_volume=transfer_volume,
                                )
                            )

                            solvent_well_count = len(from_solvent_well_info)
                            if solvent_well_count == 0:
                                logger.warning(f"No solvent wells found for {solvent}")
                                continue

                            logger.info(f"Found {solvent_well_count} solvent well(s)")

                            # Process each solvent well
                            for m, solvent_well_info in enumerate(
                                from_solvent_well_info
                            ):
                                from_solvent_well_obj = solvent_well_info[0]
                                current_transfer_volume = solvent_well_info[1]
                                to_well_obj = from_well_obj

                                from_solvent_well_index = from_solvent_well_obj.index
                                from_solvent_plate_id = (
                                    from_solvent_well_obj.plate_id.id
                                )
                                to_well_index = to_well_obj.index
                                to_plate_id = to_well_obj.plate_id.id

                                logger.info(
                                    f"Processing solvent well {m+1}/{solvent_well_count}: plate {from_solvent_plate_id}, well {from_solvent_well_index}"
                                )
                                logger.info(
                                    f"Transferring {current_transfer_volume} µL to plate {to_plate_id}, well {to_well_index}"
                                )

                                # Get plate objects
                                from_plate_obj = self.query_service.get_plate_by_id(
                                    plateid=from_solvent_well_obj.plate_id.id
                                )
                                to_plate_obj = self.query_service.get_plate_by_id(
                                    plateid=to_well_obj.plate_id.id
                                )

                                # Generate transfer command
                                aspirate_plate_name = from_plate_obj.name
                                dispense_plate_name = to_plate_obj.name
                                aspirate_well_index = from_solvent_well_obj.index
                                dispense_well_index = to_well_obj.index

                                logger.info(
                                    f"Adding dilution transfer command: {aspirate_plate_name}:{aspirate_well_index} → {dispense_plate_name}:{dispense_well_index}"
                                )

                                # --- Human-readable dilution transfer comment ---
                                self.add_command(
                                    self.annotation_generator.dilution_transfer(
                                        volume=current_transfer_volume,
                                        solvent=solvent,
                                        aspirateplatename=aspirate_plate_name,
                                        aspiratewellindex=aspirate_well_index,
                                        dispenseplatename=dispense_plate_name,
                                        dispensewellindex=dispense_well_index,
                                    )
                                )

                                self.add_command(
                                    self.command_generator.transfer_fluid_single(
                                        aspirateplatename=aspirate_plate_name,
                                        dispenseplatename=dispense_plate_name,
                                        aspiratewellindex=aspirate_well_index,
                                        dispensewellindex=dispense_well_index,
                                        transvolume=current_transfer_volume,
                                        transfertype="dilution",
                                    )
                                )

                                # Record dilution in the ledger
                                self.script_generator.transfer_ledger.record(
                                    action_type="dilution",
                                    source_plate_name=aspirate_plate_name,
                                    source_plate_role="solvent",
                                    source_well_index=aspirate_well_index,
                                    source_well_name=getattr(from_solvent_well_obj, "name", "") or "",
                                    dest_plate_name=dispense_plate_name,
                                    dest_plate_role=getattr(to_well_obj, "role", "") or "",
                                    dest_well_index=dispense_well_index,
                                    dest_well_name=getattr(to_well_obj, "name", "") or "",
                                    volume=current_transfer_volume,
                                    smiles=smiles,
                                    solvent=solvent,
                                    reaction_id=reaction_id,
                                )
                else:
                    logger.info(f"No previous reactions found for material")

    def process_reaction_actions(
        self, actionsession_obj, reaction_obj, session_number,
    ):
        """
        Process actions for a single reaction.

        Used by the per-reaction (non-MC) processing mode.

        Parameters
        ----------
        actionsession_obj : ActionSession
            The action session object
        reaction_obj : Reaction
            The reaction object
        session_number : int
            The session number
        """
        reaction_id = reaction_obj.id
        reaction_class = reaction_obj.reactionclass
        recipe_type = reaction_obj.recipe
        intramolecular = reaction_obj.intramolecular

        logger.info(
            f"Processing reaction {reaction_id} (class: {reaction_class}, recipe: {recipe_type})"
        )

        # Determine which set of actions to use based on intramolecular property
        reaction_action_search = (
            "intramolecular" if intramolecular else "intermolecular"
        )
        logger.info(f"Using {reaction_action_search} actions")

        try:
            # Get actions from Recipe DB
            reaction_actions = get_session_recipe_actions(
                reaction_class=reaction_class,
                name=recipe_type,
                session_type="reaction",
                session_number=session_number,
                molecular_context=reaction_action_search,
            )

            action_count = len(reaction_actions)
            logger.info(f"Found {action_count} actions for reaction")

            # Process each action
            for index, reaction_action in enumerate(reaction_actions):
                if isinstance(reaction_action, RecipeAddAction):
                    action_number = reaction_action.action_number

                    logger.info(
                        f"Processing action {index+1}/{action_count}: add {action_number}"
                    )
                    self.process_add_action(
                        reaction_action,
                        action_number,
                        index,
                        reaction_actions,
                        actionsession_obj,
                        reaction_obj,
                        reaction_id,
                    )

                elif isinstance(reaction_action, RecipeMixAction):
                    action_number = reaction_action.action_number
                    logger.info(
                        f"Processing action {index+1}/{action_count}: mix {action_number}"
                    )
                    self.process_mix_action(
                        action_number, actionsession_obj, reaction_obj, reaction_id
                    )
                else:
                    logger.warning(f"Unknown action type: {type(reaction_action).__name__}")
        except (KeyError, IndexError) as e:
            logger.error(
                f"Error accessing recipe data for reaction {reaction_id}: {str(e)}"
            )

    def process_add_action(
        self,
        reaction_action,
        action_number,
        index,
        reaction_actions,
        actionsession_obj,
        reaction_obj,
        reaction_id,
    ):
        """Process an add action for reaction session."""
        try:
            logger.info(
                f"Processing add action {action_number} for reaction {reaction_id}"
            )

            to_plate_role = reaction_action.to_plate_role
            to_plate_role_index = reaction_action.to_plate_role_index
            from_plate_role = reaction_action.from_plate_role
            from_plate_role_index = reaction_action.from_plate_role_index
            logger.info(f"Transfer plates: {from_plate_role}{from_plate_role_index} → {to_plate_role}{to_plate_role_index}")

            # Get the add action object
            add_action_obj = AddAction.objects.get(
                actionsession_id=actionsession_obj,
                reaction_id=reaction_obj,
                number=action_number,
            )

            smiles = add_action_obj.smiles
            solvent = add_action_obj.solvent
            transfer_volume = add_action_obj.volume
            concentration = add_action_obj.concentration

            # --- Human-readable per-action comment ---
            total_actions = len(reaction_actions)
            desc = self.annotation_generator.add_action_description(
                volume=transfer_volume,
                smiles=smiles,
                solvent=solvent,
                from_plate_role=from_plate_role,
                to_plate_role=to_plate_role,
            )
            self.add_command(
                self.annotation_generator.action_summary("Add", index + 1, total_actions, desc)
            )

            logger.info(
                f"Transfer parameters: SMILES={smiles[:20]}..., volume={transfer_volume} µL"
            )
            if solvent:
                logger.info(f"Solvent: {solvent}, concentration: {concentration}")

            # Find starting material wells
            logger.info("Finding source wells for material")
            from_well_info = self.well_finder.find_starting_material_wells(
                reaction_step_no=self.script_generator.reactionstep,
                reaction_id=reaction_id,
                smiles=smiles,
                solvent=solvent,
                concentration=concentration,
                transfer_volume=transfer_volume,
            )

            well_count = len(from_well_info)
            if well_count == 0:
                logger.warning(f"No source wells found for reaction {reaction_id}")
                return

            logger.info(f"Found {well_count} source well(s) for material")

            for i, well_info in enumerate(from_well_info):
                from_well_obj = well_info[1]
                current_transfer_volume = well_info[2]

                from_well_index = from_well_obj.index
                from_plate_id = from_well_obj.plate_id.id
                logger.info(
                    f"Processing source well {i+1}/{well_count}: plate {from_plate_id}, well {from_well_index}"
                )

                # Get destination well
                logger.info(
                    f"Finding destination well for reaction {reaction_id}, role={to_plate_role}, index={to_plate_role_index}"
                )
                to_well_obj = self.well_finder.find_reaction_well(
                    reaction_id=reaction_id,
                    role=to_plate_role,
                    role_index=to_plate_role_index,
                )

                to_well_index = to_well_obj.index
                to_plate_id = to_well_obj.plate_id.id
                logger.info(
                    f"Found destination well: plate {to_plate_id}, well {to_well_index}"
                )

                # Get plate objects
                from_plate_obj = self.query_service.get_plate_by_id(
                    plateid=from_well_obj.plate_id.id
                )
                to_plate_obj = self.query_service.get_plate_by_id(
                    plateid=to_well_obj.plate_id.id
                )

                aspirate_plate_name = from_plate_obj.name
                dispense_plate_name = to_plate_obj.name
                aspirate_well_index = from_well_obj.index
                dispense_well_index = to_well_obj.index

                # Add commands for transfer
                logger.info("Picking up tip for transfer")
                self.add_command(self.command_generator.pick_up_tip())

                logger.info(
                    f"Adding transfer command: {aspirate_plate_name}:{aspirate_well_index} → {dispense_plate_name}:{dispense_well_index} ({current_transfer_volume} µL)"
                )
                self.add_command(
                    self.command_generator.transfer_fluid_single(
                        aspirateplatename=aspirate_plate_name,
                        dispenseplatename=dispense_plate_name,
                        aspiratewellindex=aspirate_well_index,
                        dispensewellindex=dispense_well_index,
                        transvolume=current_transfer_volume,
                    )
                )

                # Record transfer in the ledger
                self.script_generator.transfer_ledger.record(
                    action_type="add",
                    source_plate_name=aspirate_plate_name,
                    source_plate_role=from_plate_role,
                    source_well_index=aspirate_well_index,
                    source_well_name=getattr(from_well_obj, "name", "") or "",
                    dest_plate_name=dispense_plate_name,
                    dest_plate_role=to_plate_role,
                    dest_well_index=dispense_well_index,
                    dest_well_name=getattr(to_well_obj, "name", "") or "",
                    volume=current_transfer_volume,
                    smiles=smiles,
                    solvent=solvent,
                    reaction_id=reaction_id,
                    reaction_class=reaction_obj.reactionclass,
                    recipe=reaction_obj.recipe,
                )

                # Drop tip if this isn't the action before a mix step
                next_action_is_mix = (
                    index + 1 < len(reaction_actions)
                ) and isinstance(reaction_actions[index + 1], RecipeMixAction)
                if not next_action_is_mix:
                    logger.info("No mix action follows, dropping tip")
                    self.add_command(self.command_generator.drop_tip())
                else:
                    logger.info("Mix action follows, keeping tip")

                # Update well status
                logger.info(f"Updating well statuses for reaction {reaction_id}")
                self.volume_manager.update_well_reactant_status(to_well_obj, True)
                self.volume_manager.update_well_reactant_status(from_well_obj, False)
        except Exception as e:
            logger.error(f"Error processing add action {action_number}: {str(e)}")

    def process_mix_action(
        self, action_number, actionsession_obj, reaction_obj, reaction_id
    ):
        """Process a mix action for reaction session."""
        try:
            logger.info(
                f"Processing mix action {action_number} for reaction {reaction_id}"
            )

            # Get mix action object
            mix_action_obj = MixAction.objects.get(
                actionsession_id=actionsession_obj,
                reaction_id=reaction_obj,
                number=action_number,
            )

            plate_role = mix_action_obj.plate_role
            plate_role_index = mix_action_obj.plate_role_index
            repetitions = mix_action_obj.repetitions

            logger.info(
                f"Mix parameters: role={plate_role}, index={plate_role_index}, repetitions={repetitions}"
            )

            # Find the well to mix
            logger.info(f"Finding well to mix for reaction {reaction_id}")
            mix_well_obj = self.well_finder.find_reaction_well(
                reaction_id=reaction_id, role=plate_role, role_index=plate_role_index
            )

            mix_well_index = mix_well_obj.index
            mix_plate_id = mix_well_obj.plate_id.id
            logger.info(
                f"Found well to mix: plate {mix_plate_id}, well {mix_well_index}"
            )

            mix_plate_obj = self.query_service.get_plate_by_id(
                plateid=mix_well_obj.plate_id.id
            )
            mix_plate_name = mix_plate_obj.name

            # --- Human-readable mix comment ---
            self.add_command(
                self.annotation_generator.reaction_mix(
                    repetitions=repetitions,
                    plate_role=plate_role,
                    plate_name=mix_plate_name,
                    well_index=mix_well_index,
                )
            )

            # Add command for mixing
            logger.info(
                f"Adding mix command: {mix_plate_name}:{mix_well_index} ({repetitions} repetitions)"
            )
            self.add_command(
                self.command_generator.mix_well(
                    wellindex=mix_well_index,
                    nomixes=repetitions,
                    plate=mix_plate_name,
                )
            )

            # Always drop tip after mixing
            logger.info("Dropping tip after mixing")
            self.add_command(self.command_generator.drop_tip())
        except Exception as e:
            logger.error(f"Error processing mix action {action_number}: {str(e)}")
