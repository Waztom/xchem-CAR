"""
Implements workup session functionality for OpenTrons protocols.
"""

import logging

from .base_session import BaseSession

logger = logging.getLogger(__name__)

class WorkupSession(BaseSession):
    """
    Session class for executing workup protocols.
    """
    
    def execute(self):
        """
        Execute the workup session protocol.
        """
        logger.info(f"Executing workup session for reaction step {self.reactionstep}")
        
        # Setup common resources (inherited from BaseSession)
        self.setup_common_resources()
        
        # Get add and extract actions
        self.addactionqueryset = self.data_manager.get_add_action_query_set(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids
        )
        
        self.extractactionqueryset = self.data_manager.get_extract_action_query_set(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids
        )
        
        # Create dataframe from add actions
        self.addactionsdf = self.data_manager.get_add_actions_dataframe(
            addactionqueryset=self.addactionqueryset
        )
        
        # Calculate rounded volumes for add and extract actions
        self.roundedaddvolumes = self.data_manager.get_rounded_add_action_volumes(
            addactionqueryset=self.addactionqueryset
        )
        self.roundedextractvolumes = self.data_manager.get_rounded_extract_action_volumes(
            extractactionqueryset=self.extractactionqueryset
        )
        self.roundedvolumes = self.roundedaddvolumes + self.roundedextractvolumes
        
        # Setup pipettes and tip racks based on volumes
        tiprack_type = self.pipette_manager.get_tip_rack_type(rounded_volumes=self.roundedvolumes)
        self.tipracktype = tiprack_type
        self.pipettetype = self.pipette_manager.get_pipette_type(rounded_volumes=self.roundedvolumes)
        
        # Create tip racks
        self.pipette_manager.create_tip_racks(tiprack_type=tiprack_type)
        
        # Determine workup stage number
        workup_stage = self.determine_workup_stage()
        logger.info(f"Executing workup stage {workup_stage} for step {self.reactionstep}")
        
        # Find source reaction plates
        source_plates = self.get_source_plates_for_workup()
        if source_plates:
            self.plate_manager.update_plate_deck_ot_session_ids(plate_queryset=source_plates)
        
        # Create workup plates based on stage
        workup_plate_type = f"workup{workup_stage}"
        volumes = self.data_manager.get_rounded_reaction_volumes(reactionqueryset=self.groupreactionqueryset)
        
        # Get appropriate labware type for this workup
        labware_type = self.plate_manager.get_plate_type(
            platetype="workup", 
            volumes=volumes
        )
        
        # Create the workup plate
        workup_plates = self.plate_manager.create_workup_plate(
            reaction_queryset=self.groupreactionqueryset,
            labware_platetype=labware_type,
            platetype=workup_plate_type
        )
        
        # Create pipette model
        self.pipette_manager.create_pipette_model()
        
        # Create solvent plate if needed for this workup
        if self.addactionqueryset.exists():
            self.solventmaterialsdf = self.material_manager.get_add_actions_material_dataframe(
                product_exists=False
            )
            if not self.solventmaterialsdf.empty:
                self.plate_manager.create_solvent_plate(materialsdf=self.solventmaterialsdf)
        
        logger.info(f"Workup session execution completed for step {self.reactionstep}, stage {workup_stage}")
        return True
    
    def determine_workup_stage(self):
        """
        Determine the workup stage number for this session.
        
        Returns
        -------
        workup_stage: int
            The workup stage number (1, 2, or 3)
        """
        # Check to plates this workup targets
        to_plate_types = set()
        
        if self.addactionqueryset.exists():
            to_plate_types.update(
                self.addactionqueryset.values_list('toplatetype', flat=True).distinct()
            )
            
        if self.extractactionqueryset.exists():
            to_plate_types.update(
                self.extractactionqueryset.values_list('toplatetype', flat=True).distinct()
            )
        
        if 'workup1' in to_plate_types:
            return 1
        elif 'workup2' in to_plate_types:
            return 2
        elif 'workup3' in to_plate_types:
            return 3
        else:
            # Default to stage 1 if not explicitly specified
            logger.warning("Could not determine workup stage, defaulting to 1")
            return 1
    
    def get_source_plates_for_workup(self):
        """
        Get source plates that contain material for this workup session.
        
        Returns
        -------
        source_plates: QuerySet
            The plates containing material to be worked up
        """
        # For first workup stage, get reaction plates
        workup_stage = self.determine_workup_stage()
        
        if workup_stage == 1:
            # Source is reaction plates
            plate_type = "reaction"
        else:
            # Source is previous workup stage
            plate_type = f"workup{workup_stage - 1}"
        
        # Find all plates from previous steps in this batch protocol
        all_plates = self.plate_manager.get_all_ot_batch_protocol_plates(
            otbatchprotocol_id=self.otbatchprotocolobj
        )
        
        # Filter for the correct type
        source_plates = [plate for plate in all_plates if plate.type == plate_type]
        
        if not source_plates:
            logger.warning(f"No source plates found of type {plate_type} for workup stage {workup_stage}")
        
        return source_plates