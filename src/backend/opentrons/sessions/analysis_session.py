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
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids
        )
        
        # Get add actions for diluent/standards if any
        self.addactionqueryset = self.data_manager.get_add_action_query_set(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids
        )
        
        # Calculate rounded volumes for extract actions
        self.roundedextractvolumes = self.data_manager.get_rounded_extract_action_volumes(
            extractactionqueryset=self.extractactionqueryset
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
        tiprack_type = self.pipette_manager.get_tip_rack_type(rounded_volumes=self.roundedvolumes)
        self.tipracktype = tiprack_type
        self.pipettetype = self.pipette_manager.get_pipette_type(rounded_volumes=self.roundedvolumes)
        
        # Create tip racks
        self.pipette_manager.create_tip_racks(tiprack_type=tiprack_type)
        
        # Find source plates (workup or reaction plates)
        source_plates = self.get_source_plates_for_analysis()
        if source_plates:
            self.plate_manager.update_plate_deck_ot_session_ids(plate_queryset=source_plates)
        
        # Create analysis plate
        volumes = self.roundedvolumes
        
        # Get appropriate labware type for analysis
        labware_type = self.plate_manager.get_plate_type(
            platetype="analyse", 
            volumes=volumes
        )
        
        # Create the analysis plate
        analysis_plates = self.plate_manager.create_analyse_plate(
            reaction_queryset=self.groupreactionqueryset,
            labware_platetype=labware_type
        )
        
        # Create pipette model
        self.pipette_manager.create_pipette_model()
        
        # Create solvent/standard plate if needed
        if self.addactionqueryset.exists() and hasattr(self, 'addactionsdf') and not self.addactionsdf.empty:
            self.solventmaterialsdf = self.material_manager.get_add_actions_material_dataframe(
                product_exists=False
            )
            if not self.solventmaterialsdf.empty:
                self.plate_manager.create_solvent_plate(materialsdf=self.solventmaterialsdf)
        
        logger.info(f"Analysis session execution completed for step {self.reactionstep}")
        return True
    
    def get_source_plates_for_analysis(self):
        """
        Get source plates that contain material for analysis.
        
        Returns
        -------
        source_plates: list
            The plates containing material to be analyzed
        """
        # Determine source plate types
        # Analysis typically comes from workup3 -> workup2 -> workup1 -> reaction in that order of preference
        source_plate_types = ["workup3", "workup2", "workup1", "reaction"]
        
        # Find all plates from previous steps in this batch protocol
        all_plates = self.plate_manager.get_all_ot_batch_protocol_plates(
            otbatchprotocol_id=self.otbatchprotocolobj
        )
        
        source_plates = []
        
        # Try to find plates of each type in order of preference
        for plate_type in source_plate_types:
            plates_of_type = [plate for plate in all_plates if plate.type == plate_type]
            if plates_of_type:
                # Found plates of this type, use them
                logger.info(f"Using {len(plates_of_type)} plates of type {plate_type} for analysis")
                source_plates.extend(plates_of_type)
                break
        
        if not source_plates:
            logger.warning("No suitable source plates found for analysis")
        
        return source_plates