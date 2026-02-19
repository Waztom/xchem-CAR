"""
Workup Session Handler for OpenTrons protocols.

This module provides methods for handling workup sessions which typically
involve extract and mix operations.
"""

import logging
from django.db.models import QuerySet

from backend.models import ExtractAction, MixAction, RecipeExtractAction, RecipeMixAction
from backend.db_utils import get_reaction
from backend.recipe_utils import get_session_recipe_actions

from .base_handler import SessionHandler

logger = logging.getLogger(__name__)


class WorkupSessionHandler(SessionHandler):
    """
    Handles workup session processing.

    This class implements the logic for generating workup session commands
    like extraction and mixing.
    """

    def process_session(self, actionsession_queryset: QuerySet) -> None:
        """
        Process the workup session(s).

        Parameters
        ----------
        actionsession_queryset : QuerySet
            QuerySet of workup action sessions to process
        """
        # Use base helper to log session start
        self.log_session_start(actionsession_queryset, "workup")

        # Add intro comment
        self.add_command("\n\t# Processing workup actions")

        # Get session information
        session_number = self.get_session_number(actionsession_queryset)
        reaction_step = self.script_generator.reactionstep
        logger.info(f"Processing workup step {reaction_step}, session {session_number}")

        # Process each action session
        action_count = actionsession_queryset.count()
        logger.info(f"Processing {action_count} individual workup action session(s)")

        for i, actionsession_obj in enumerate(actionsession_queryset):
            reaction_id = actionsession_obj.reaction_id.id
            logger.info(
                f"Processing workup for reaction {reaction_id} ({i+1}/{action_count})"
            )
            reaction_obj = get_reaction(reaction_id=reaction_id)
            self.process_workup_actions(actionsession_obj, reaction_obj, session_number)

        self.log_session_end("workup")

    def process_workup_actions(self, actionsession_obj, reaction_obj, session_number):
        """
        Process workup actions for a reaction.

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
            f"Processing workup for reaction {reaction_id} (class: {reaction_class}, recipe: {recipe_type})"
        )

        # Determine which set of actions to use based on intramolecular property
        reaction_action_search = (
            "intramolecular" if intramolecular else "intermolecular"
        )
        logger.info(f"Using {reaction_action_search} workup actions")

        # Get actions from Recipe DB
        try:
            workup_actions = get_session_recipe_actions(
                reaction_class=reaction_class,
                name=recipe_type,
                session_type="workup",
                session_number=session_number,
                molecular_context=reaction_action_search,
            )

            if not workup_actions:
                logger.warning(
                    f"No workup actions defined for reaction {reaction_id}, recipe {recipe_type}"
                )
                self.add_command(
                    f"\n\t# No workup actions found for reaction {reaction_id}"
                )
                return

            action_count = len(workup_actions)
            logger.info(
                f"Found {action_count} workup actions for reaction {reaction_id}"
            )

            # Process each workup action
            for index, workup_action in enumerate(workup_actions):
                if isinstance(workup_action, RecipeExtractAction):
                    action_number = workup_action.action_number
                    logger.info(
                        f"Processing workup action {index+1}/{len(workup_actions)}: extract {action_number}"
                    )
                    self.process_extract_action(
                        workup_action,
                        action_number,
                        index,
                        workup_actions,
                        actionsession_obj,
                        reaction_obj,
                        reaction_id,
                    )

                elif isinstance(workup_action, RecipeMixAction):
                    action_number = workup_action.action_number
                    logger.info(
                        f"Processing workup action {index+1}/{len(workup_actions)}: mix {action_number}"
                    )
                    self.process_mix_action(
                        action_number, actionsession_obj, reaction_obj, reaction_id
                    )
                else:
                    logger.warning(f"Unknown workup action type: {type(workup_action).__name__}")

        except (KeyError, IndexError) as e:
            logger.error(
                f"Error accessing workup recipe data for reaction {reaction_id}: {str(e)}"
            )
            self.add_command(
                f"\n\t# Error finding workup actions for reaction {reaction_id}: {str(e)}"
            )

    def process_extract_action(
        self,
        workup_action,
        action_number,
        index,
        workup_actions,
        actionsession_obj,
        reaction_obj,
        reaction_id,
    ):
        """
        Process an extract action for workup session.

        Parameters
        ----------
        workup_action : Dict
            The workup action data
        action_number : int
            The action number
        index : int
            The index of the action in the action list
        workup_actions : List[Dict]
            All workup actions for the session
        actionsession_obj : ActionSession
            The action session object
        reaction_obj : Reaction
            The reaction object
        reaction_id : int
            The reaction ID
        """
        logger.info(
            f"Processing extract action {action_number} for reaction {reaction_id}"
        )

        try:
            # Get extract action object
            extract_action_obj = ExtractAction.objects.get(
                actionsession_id=actionsession_obj,
                reaction_id=reaction_obj,
                number=action_number,
            )

            from_plate_role = extract_action_obj.from_plate_role
            from_plate_role_index = extract_action_obj.from_plate_role_index
            to_plate_role = extract_action_obj.to_plate_role
            to_plate_role_index = extract_action_obj.to_plate_role_index
            extract_layer = extract_action_obj.layer
            bottom_layer_volume = extract_action_obj.bottomlayervolume

            logger.info(
                f"Extract parameters: from={from_plate_role}{from_plate_role_index}, to={to_plate_role}{to_plate_role_index}, "
                + f"layer={extract_layer}, bottom_volume={bottom_layer_volume} µL"
            )

            # Find source well
            logger.info(
                f"Finding source well for reaction {reaction_id}, role={from_plate_role}, index={from_plate_role_index}"
            )
            from_well_obj = self.well_finder.find_reaction_well(
                reaction_id=reaction_id,
                role=from_plate_role,
                role_index=from_plate_role_index,
            )

            from_well_index = from_well_obj.index
            from_plate_id = from_well_obj.plate_id.id
            logger.info(
                f"Found source well: plate {from_plate_id}, well {from_well_index}"
            )

            # Find destination well
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

            # Calculate aspirate height for extraction
            logger.info(
                f"Calculating aspirate height for bottom layer volume {bottom_layer_volume} µL"
            )
            aspirate_height = self.volume_manager.calculate_aspirate_height(
                bottomlayervolume=bottom_layer_volume,
                labware=from_plate_obj.labware,
            )
            logger.info(f"Calculated aspirate height: {aspirate_height} mm")

            # Get volume to transfer
            transfer_volume = extract_action_obj.volume
            logger.info(f"Extract volume: {transfer_volume} µL")

            # Add extraction commands
            logger.info("Picking up tip for extraction")
            self.add_command(self.command_generator.pick_up_tip())

            aspirate_plate_name = from_plate_obj.name
            dispense_plate_name = to_plate_obj.name
            aspirate_well_index = from_well_obj.index
            dispense_well_index = to_well_obj.index

            logger.info(
                f"Adding extraction command: {aspirate_plate_name}:{aspirate_well_index} → "
                + f"{dispense_plate_name}:{dispense_well_index} ({transfer_volume} µL, layer: {extract_layer})"
            )

            self.add_command(
                self.command_generator.transfer_fluid_single(
                    aspirateplatename=aspirate_plate_name,
                    dispenseplatename=dispense_plate_name,
                    aspiratewellindex=aspirate_well_index,
                    dispensewellindex=dispense_well_index,
                    transvolume=transfer_volume,
                    aspirateheight=aspirate_height,
                    transfertype=f"extraction-{extract_layer}",
                )
            )

            # Drop tip if this isn't the action before a mix step
            next_action_is_mix = (index + 1 < len(workup_actions)) and isinstance(
                workup_actions[index + 1], RecipeMixAction
            )
            if not next_action_is_mix:
                logger.info("No mix action follows, dropping tip")
                self.add_command(self.command_generator.drop_tip())
            else:
                logger.info("Mix action follows, keeping tip")

            # Update well status
            logger.info(f"Updating well statuses for reaction {reaction_id}")
            self.volume_manager.update_well_reactant_status(to_well_obj, True)
            self.volume_manager.update_well_reactant_status(from_well_obj, False)

        except ExtractAction.DoesNotExist:
            logger.error(
                f"Extract action not found for reaction {reaction_id}, action number {action_number}"
            )
            self.add_command(
                f"\n\t# Error: Extract action not found for reaction {reaction_id}"
            )

        except Exception as e:
            logger.error(f"Error processing extract action: {str(e)}")
            self.add_command(f"\n\t# Error processing extract action: {str(e)}")

    def process_mix_action(
        self, action_number, actionsession_obj, reaction_obj, reaction_id
    ):
        """
        Process a mix action for workup session.

        Parameters
        ----------
        action_number : int
            The action number
        actionsession_obj : ActionSession
            The action session object
        reaction_obj : Reaction
            The reaction object
        reaction_id : int
            The reaction ID
        """
        logger.info(f"Processing mix action {action_number} for reaction {reaction_id}")

        try:
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
            logger.info(
                f"Finding well to mix for reaction {reaction_id}, role={plate_role}, index={plate_role_index}"
            )
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

        except MixAction.DoesNotExist:
            logger.error(
                f"Mix action not found for reaction {reaction_id}, action number {action_number}"
            )
            self.add_command(
                f"\n\t# Error: Mix action not found for reaction {reaction_id}"
            )

        except Exception as e:
            logger.error(f"Error processing mix action: {str(e)}")
            self.add_command(f"\n\t# Error processing mix action: {str(e)}")
