"""
Implements reaction session functionality for OpenTrons protocols.
"""

import logging

from .base_session import BaseSession

logger = logging.getLogger(__name__)


class ReactionSession(BaseSession):
    """
    Session class for executing reaction protocols.
    """

    def execute(self):
        """
        Execute the reaction session protocol.
        """
        logger.info(f"Executing reaction session for reaction step {self.reactionstep}")

        # Setup common resources (inherited from BaseSession)
        self.setup_common_resources()

        # Group reactions by temperature
        self.grouped_reaction_temperature_querysets = (
            self.data_manager.get_grouped_temperature_reactions(
                reactionqueryset=self.groupreactionqueryset
            )
        )

        # Get add actions
        self.addactionqueryset = self.data_manager.get_add_action_query_set(
            reaction_ids=self.reaction_ids, actionsession_ids=self.actionsession_ids
        )

        logger.info(f"Add action queryset: {self.addactionqueryset}")

        # Create dataframe from add actions
        self.addactionsdf = self.data_manager.get_add_actions_dataframe(
            addactionqueryset=self.addactionqueryset
        )

        # Get extract actions
        self.extractactionqueryset = self.data_manager.get_extract_action_query_set(
            reaction_ids=self.reaction_ids, actionsession_ids=self.actionsession_ids
        )

        # Calculate rounded volumes for add and extract actions
        self.roundedaddvolumes = self.data_manager.get_rounded_add_action_volumes(
            addactionqueryset=self.addactionqueryset
        )
        self.roundedextractvolumes = (
            self.data_manager.get_rounded_extract_action_volumes(
                extractactionqueryset=self.extractactionqueryset
            )
        )
        self.roundedvolumes = self.roundedaddvolumes + self.roundedextractvolumes

        # Setup pipettes and tip racks based on volumes
        tiprack_type = self.pipette_manager.get_tip_rack_type(
            rounded_volumes=self.roundedvolumes
        )
        self.tipracktype = tiprack_type
        self.pipettetype = self.pipette_manager.get_pipette_type(
            rounded_volumes=self.roundedvolumes
        )

        # Create tip racks
        self.pipette_manager.create_tip_racks(tiprack_type=tiprack_type)

        # Handle continuation actions (reactions continuing from previous steps)
        continuation_sessions = self.actionsessionqueryset.filter(continuation=True)
        non_continuation_sessions = self.actionsessionqueryset.filter(
            continuation=False
        )

        # Determine which SMILES to search for in existing plates
        if continuation_sessions.exists():
            searchsmiles = self.material_manager.get_product_smiles(
                reaction_ids=self.reaction_ids
            )
            if searchsmiles:
                searchsmiles += list(
                    self.addactionqueryset.values_list("smiles", flat=True)
                )
        else:
            searchsmiles = list(self.addactionqueryset.values_list("smiles", flat=True))
            # Create reaction plate for non-continuation reactions
            self.plate_factory.create_plates_by_temperature(
                grouped_reaction_temperature_querysets=self.grouped_reaction_temperature_querysets,
                role="reaction",
                role_index=1,
            )

        # Get input plates needed based on SMILES
        input_plates = self.plate_query_service.get_input_plates_needed(
            searchsmiles=searchsmiles,
            groupreactionqueryset=self.groupreactionqueryset,
        )
        if input_plates:
            self.plate_factory.update_plate_deck_ot_session_ids(
                plate_queryset=input_plates
            )

        # Create pipette model
        self.pipette_manager.create_pipette_model()

        # Create custom starting material plates if CSV path provided
        if self.customSMcsvpath:
            self.plate_factory.create_starting_material_plates_from_csv(
                csv_path=self.customSMcsvpath
            )

        # Create reaction starting plate for new materials
        self.plate_factory.create_reaction_starting_plate()

        # For steps after the first, create solvent plate for previous products
        if self.reactionstep > 1:
            logger.info(f"Creating solvent plate for step {self.reactionstep}")
            # Get solvent materials for the current step
            self.solventmaterialsdf = (
                self.material_manager.get_add_actions_material_dataframe(
                    product_exists=True
                )
            )
            if (
                self.solventmaterialsdf is not None
                and not self.solventmaterialsdf.empty
            ):
                self.plate_factory.create_solvent_plate(
                    materials_df=self.solventmaterialsdf
                )

        logger.info(
            f"Reaction session execution completed for step {self.reactionstep}"
        )
        return True
