"""
Analysis Session Handler for OpenTrons protocols.

This module provides methods for handling analysis sessions which typically
involve transferring samples to analysis plates.
"""

import logging
from typing import List, Dict, Any, Optional, Union, Tuple
from django.db.models import QuerySet, Q

from backend.models import ActionSession, AddAction
from backend.utils import getReaction, getReactionQuerySet
from backend.recipebuilder.encodedrecipes import encoded_recipes

from .base_handler import SessionHandler

logger = logging.getLogger(__name__)


class AnalysisSessionHandler(SessionHandler):
    """
    Handles analysis session processing.

    This class implements the logic for generating analysis session commands
    like sample transfers to analysis plates.
    """

    def process_session(self, actionsession_queryset: QuerySet) -> None:
        """
        Process the analysis session(s).

        Parameters
        ----------
        actionsession_queryset : QuerySet
            QuerySet of analysis action sessions to process
        """
        # Use base helper to log session start
        self.log_session_start(actionsession_queryset, "analysis")

        # Add intro comment
        self.add_command("\n\t# Processing analysis actions")

        # Get session information
        session_number = self.get_session_number(actionsession_queryset)
        reaction_step = self.script_generator.reactionstep
        logger.info(
            f"Processing analysis step {reaction_step}, session {session_number}"
        )

        # Process each action session
        action_count = actionsession_queryset.count()
        logger.info(f"Processing {action_count} individual analysis action session(s)")

        for i, actionsession_obj in enumerate(actionsession_queryset):
            reaction_id = actionsession_obj.reaction_id.id
            logger.info(
                f"Processing analysis for reaction {reaction_id} ({i+1}/{action_count})"
            )
            reaction_obj = getReaction(reaction_id=reaction_id)
            self.process_analysis_actions(
                actionsession_obj, reaction_obj, session_number
            )

        self.log_session_end("analysis")

    def process_analysis_actions(self, actionsession_obj, reaction_obj, session_number):
        """
        Process analysis actions for a reaction.

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
            f"Processing analysis for reaction {reaction_id} (class: {reaction_class}, recipe: {recipe_type})"
        )

        # Determine which set of actions to use based on intramolecular property
        reaction_action_search = (
            "intramolecular" if intramolecular else "intermolecular"
        )
        logger.info(f"Using {reaction_action_search} analysis actions")

        # Get actions from encoded recipes
        try:
            action_sessions = encoded_recipes[reaction_class]["recipes"][recipe_type][
                "actionsessions"
            ]
            analysis_actions = None

            try:
                analysis_actions = [
                    actionsession[reaction_action_search]["actions"]
                    for actionsession in action_sessions
                    if actionsession["type"] == "analyse"
                    and actionsession["sessionnumber"] == session_number
                ][0]

                action_count = len(analysis_actions)
                logger.info(
                    f"Found {action_count} analysis actions for reaction {reaction_id}"
                )
            except (IndexError, KeyError) as e:
                logger.warning(
                    f"No analysis actions defined for reaction {reaction_id}, recipe {recipe_type}: {str(e)}"
                )
                self.add_command(
                    f"\n\t# No analysis actions found for reaction {reaction_id}"
                )
                return

            # Process each analysis action
            for index, analysis_action in enumerate(analysis_actions):
                action_type = analysis_action["type"]
                action_number = analysis_action["actionnumber"]

                logger.info(
                    f"Processing analysis action {index+1}/{len(analysis_actions)}: {action_type} {action_number}"
                )

                if action_type == "add":
                    self.process_add_action(
                        analysis_action,
                        action_number,
                        index,
                        analysis_actions,
                        actionsession_obj,
                        reaction_obj,
                        reaction_id,
                    )
                else:
                    logger.warning(f"Unknown analysis action type: {action_type}")

        except (KeyError, IndexError) as e:
            logger.error(
                f"Error accessing analysis recipe data for reaction {reaction_id}: {str(e)}"
            )
            self.add_command(
                f"\n\t# Error finding analysis actions for reaction {reaction_id}: {str(e)}"
            )

    def process_add_action(
        self,
        analysis_action,
        action_number,
        index,
        analysis_actions,
        actionsession_obj,
        reaction_obj,
        reaction_id,
    ):
        """
        Process an add action for analysis session.

        Parameters
        ----------
        analysis_action : Dict
            The analysis action data
        action_number : int
            The action number
        index : int
            The index of the action in the action list
        analysis_actions : List[Dict]
            All analysis actions for the session
        actionsession_obj : ActionSession
            The action session object
        reaction_obj : Reaction
            The reaction object
        reaction_id : int
            The reaction ID
        """
        logger.info(
            f"Processing analysis add action {action_number} for reaction {reaction_id}"
        )

        try:
            # Get add action object
            add_action_obj = AddAction.objects.get(
                actionsession_id=actionsession_obj,
                reaction_id=reaction_obj,
                number=action_number,
            )

            try:
                from_plate_type = analysis_action["content"]["plates"]["fromplatetype"]
                to_plate_type = analysis_action["content"]["plates"]["toplatetype"]
            except (KeyError, TypeError) as e:
                logger.error(
                    f"Error accessing plate types in analysis action: {str(e)}"
                )
                from_plate_type = add_action_obj.fromplatetype
                to_plate_type = add_action_obj.toplatetype
                logger.info(
                    f"Using plate types from database: from={from_plate_type}, to={to_plate_type}"
                )

            logger.info(
                f"Analysis transfer: from={from_plate_type}, to={to_plate_type}"
            )

            # Find source well (reaction/workup well)
            logger.info(
                f"Finding source well for reaction {reaction_id}, type {from_plate_type}"
            )
            from_well_obj = self.well_finder.find_reaction_well(
                reaction_id=reaction_id,
                well_type=from_plate_type,
            )

            from_well_index = from_well_obj.index
            from_plate_id = from_well_obj.plate_id.id
            logger.info(
                f"Found source well: plate {from_plate_id}, well {from_well_index}"
            )

            # Find destination well (analysis well)
            logger.info(
                f"Finding destination well for reaction {reaction_id}, type {to_plate_type}"
            )
            to_well_obj = self.well_finder.find_reaction_well(
                reaction_id=reaction_id,
                well_type=to_plate_type,
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

            # Get transfer volume
            transfer_volume = add_action_obj.volume
            logger.info(f"Analysis transfer volume: {transfer_volume} µL")

            # Add transfer commands
            logger.info("Picking up tip for analysis transfer")
            self.add_command(self.command_generator.pick_up_tip())

            aspirate_plate_name = from_plate_obj.name
            dispense_plate_name = to_plate_obj.name
            aspirate_well_index = from_well_obj.index
            dispense_well_index = to_well_obj.index

            logger.info(
                f"Adding analysis transfer command: {aspirate_plate_name}:{aspirate_well_index} → "
                + f"{dispense_plate_name}:{dispense_well_index} ({transfer_volume} µL)"
            )

            self.add_command(
                self.command_generator.transfer_fluid_single(
                    aspirateplatename=aspirate_plate_name,
                    dispenseplatename=dispense_plate_name,
                    aspiratewellindex=aspirate_well_index,
                    dispensewellindex=dispense_well_index,
                    transvolume=transfer_volume,
                    transfertype="analysis",
                )
            )

            # Add mix command if specified
            if add_action_obj.mix:
                mix_repetitions = 3  # Standard number of mixes
                logger.info(
                    f"Adding mix command in destination well ({mix_repetitions} repetitions)"
                )
                self.add_command(
                    self.command_generator.mix_well(
                        wellindex=dispense_well_index,
                        nomixes=mix_repetitions,
                        plate=dispense_plate_name,
                    )
                )

            # Always drop tip after analysis transfer
            logger.info("Dropping tip after analysis transfer")
            self.add_command(self.command_generator.drop_tip())

            # For analysis transfers, we typically don't update reactant status
            logger.info(
                "Analysis transfer complete (well status not updated for analysis plates)"
            )

        except AddAction.DoesNotExist:
            logger.error(
                f"Add action not found for reaction {reaction_id}, action number {action_number}"
            )
            self.add_command(
                f"\n\t# Error: Add action not found for reaction {reaction_id}"
            )

        except Exception as e:
            logger.error(f"Error processing analysis add action: {str(e)}")
            self.add_command(f"\n\t# Error processing analysis add action: {str(e)}")
