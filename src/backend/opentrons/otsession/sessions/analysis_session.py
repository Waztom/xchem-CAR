"""
Implements analysis session functionality for OpenTrons protocols.
"""

import logging
from django.db.models import Q

from .base_session import BaseSession

logger = logging.getLogger(__name__)


class AnalysisSession(BaseSession):
    """
    Session class for executing analysis protocols.
    """

    def execute(self):
        """
        Execute the analysis session protocol.
        """
        logger.info(f"Executing analysis session for reaction step {self.reactionstep}")

        # Setup common resources (inherited from BaseSession)
        self.setup_common_resources()

        # Get extract actions (analysis primarily uses extractions)
        self.extractactionqueryset = self.data_manager.get_extract_action_query_set(
            reaction_ids=self.reaction_ids, actionsession_ids=self.actionsession_ids
        )

        # Get add actions for diluent/standards if any
        self.addactionqueryset = self.data_manager.get_add_action_query_set(
            reaction_ids=self.reaction_ids, actionsession_ids=self.actionsession_ids
        )

        # Calculate rounded volumes for extract actions
        self.roundedextractvolumes = (
            self.data_manager.get_rounded_extract_action_volumes(
                extractactionqueryset=self.extractactionqueryset
            )
        )

        # Get add action volumes if any
        if self.addactionqueryset.exists():
            self.roundedaddvolumes = self.data_manager.get_rounded_add_action_volumes(
                addactionqueryset=self.addactionqueryset
            )
            self.addactionsdf = self.data_manager.get_add_actions_dataframe(
                addactionqueryset=self.addactionqueryset
            )
        else:
            self.roundedaddvolumes = []

        # Combine volumes for pipette selection
        self.roundedvolumes = self.roundedextractvolumes + self.roundedaddvolumes

        # Select appropriate pipette based on volumes
        tiprack_type = self.pipette_manager.get_tip_rack_type(
            rounded_volumes=self.roundedvolumes
        )
        self.tipracktype = tiprack_type
        self.pipettetype = self.pipette_manager.get_pipette_type(
            rounded_volumes=self.roundedvolumes
        )

        # Create tip racks
        self.pipette_manager.create_tip_racks(tiprack_type=tiprack_type)

        # Find source plates (workup or reaction plates)
        source_plates = self.get_source_plates_for_analysis()
        if source_plates:
            self.plate_factory.update_plate_deck_ot_session_ids(
                plate_queryset=source_plates
            )

        # Create analysis plate
        analysis_plates = self.plate_factory.create_analyse_plate()

        # Create pipette model
        self.pipette_manager.create_pipette_model()

        # Create solvent/standard plate if needed
        if (
            self.addactionqueryset.exists()
            and hasattr(self, "addactionsdf")
            and not self.addactionsdf.empty
        ):
            self.solventmaterialsdf = (
                self.material_manager.get_add_actions_material_dataframe(
                    product_exists=False
                )
            )
            if not self.solventmaterialsdf.empty:
                self.plate_factory.create_solvent_plate(
                    materials_df=self.solventmaterialsdf
                )

        logger.info(
            f"Analysis session execution completed for step {self.reactionstep}"
        )
        return True

    def get_source_plates_for_analysis(self):
        """
        Get source plates that contain material for analysis.

        Returns
        -------
        source_plates: list
            The plates containing material to be analyzed
        """
        # Find all plates from previous steps in this batch protocol
        all_plates = self.plate_query_service.get_all_ot_batch_protocol_plates(
            otbatchprotocol_id=self.otbatchprotocolobj
        )

        source_plates = []

        # Get workup plates sorted by role_index descending (highest index = most processed)
        workup_plates = sorted(
            [p for p in all_plates if p.role == "workup"],
            key=lambda p: p.role_index,
            reverse=True
        )

        if workup_plates:
            # Use the highest indexed workup plate(s)
            highest_index = workup_plates[0].role_index
            source_plates = [p for p in workup_plates if p.role_index == highest_index]
            logger.info(
                f"Using {len(source_plates)} workup plate(s) with index {highest_index} for analysis"
            )
        else:
            # Fall back to reaction plates
            reaction_plates = [p for p in all_plates if p.role == "reaction"]
            if reaction_plates:
                source_plates = reaction_plates
                logger.info(
                    f"Using {len(source_plates)} reaction plate(s) for analysis"
                )

        if not source_plates:
            logger.warning("No suitable source plates found for analysis")

        return source_plates
