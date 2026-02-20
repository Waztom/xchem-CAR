"""
Reaction Session Handler for OpenTrons protocols.

This module provides methods for handling reaction sessions.
"""

import logging
from django.db.models import QuerySet

from backend.models import AddAction, MixAction, RecipeAddAction, RecipeMixAction
from backend.db_utils import get_reaction, get_previous_reaction_query_sets, get_reaction_query_set
from backend.recipe_utils import get_session_recipe_actions

from .base_handler import SessionHandler

logger = logging.getLogger(__name__)


class ReactionSessionHandler(SessionHandler):
    """
    Handles reaction session processing.

    This class implements the logic for generating reaction session commands.
    """

    def process_session(self, actionsession_queryset: QuerySet) -> None:
        """
        Process the reaction session(s).

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

        # Process each action session
        action_count = actionsession_queryset.count()
        logger.info(f"Processing {action_count} individual reaction action session(s)")

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
                actionsession_obj, reaction_obj, session_number
            )

        self.log_session_end("reaction")

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

    def process_reaction_actions(self, actionsession_obj, reaction_obj, session_number):
        """
        Process actions for a single reaction.

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
