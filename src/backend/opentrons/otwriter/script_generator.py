"""
Script Generator for OpenTrons protocols.

This module coordinates the generation of OpenTrons protocol scripts
from action sessions.
"""

import logging
from typing import List, Union

from backend.models import (
    OTSession,
)
from backend.db_utils import getReactionQuerySet

from .file_manager import FileManager
from .command_generator import CommandGenerator
from .utils.query_service import QueryService
from .utils.volume_manager import VolumeManager
from .session_handlers.reaction_handler import ReactionSessionHandler
from .session_handlers.workup_handler import WorkupSessionHandler
from .session_handlers.analysis_handler import AnalysisSessionHandler
from .utils.well_finder import WellFinder

logger = logging.getLogger(__name__)


class ScriptGenerator:
    """
    Coordinates the generation of OpenTrons protocol scripts.

    This class is the main entry point for script generation and replaces
    the original OTWrite class with a more modular design.
    """

    def __init__(
        self,
        batchtag: str,
        otsessionobj: OTSession,
        actionsession_ids: List[int],
        apiLevel: str = "2.9",
    ):
        """
        Initialize the script generator.

        Parameters
        ----------
        batchtag : str
            Unique identifier tag for this batch
        otsessionobj : OTSession
            OTSession object containing session information
        actionsession_ids : List[int]
            List of action session IDs to process
        apiLevel : str, optional
            OpenTrons API level to use, by default "2.9"
        """
        logger.info(
            f"Initializing ScriptGenerator for {otsessionobj.sessiontype} session ID {otsessionobj.id}"
        )
        logger.info(
            f"Processing {len(actionsession_ids)} action sessions with batch tag: {batchtag}"
        )

        # Initialize key properties from original OTWrite class
        self.reactionstep = otsessionobj.reactionstep
        self.otsessionobj = otsessionobj
        self.otsession_id = otsessionobj.id
        self.otsessiontype = otsessionobj.sessiontype
        self.batchtag = batchtag
        self.apiLevel = apiLevel
        self.actionsession_ids = actionsession_ids

        # Set up protocol name BEFORE creating FileManager
        self.protocolname = f"{self.otsessiontype}-b{self.batchtag}-r{self.reactionstep}-s{self.otsession_id}"
        logger.info(f"Protocol name: {self.protocolname}")

        logger.info(f"Creating helper components for reaction step {self.reactionstep}")

        logger.info(f"Creating helper components for reaction step {self.reactionstep}")
        # Create helper components
        self.command_generator = CommandGenerator(self)
        self.file_manager = FileManager(self)
        self.query_service = QueryService(self)
        self.volume_manager = VolumeManager(self)
        self.well_finder = WellFinder(self)

        # Initialize session handlers
        self.reaction_handler = ReactionSessionHandler(self)
        self.workup_handler = WorkupSessionHandler(self)
        self.analysis_handler = AnalysisSessionHandler(self)

        # Script content that will be built
        self.content = []

        # Get action session data
        logger.info("Retrieving action sessions and reactions")
        self.actionsessionqueryset = self.query_service.get_action_session_query_set()

        if not self.actionsessionqueryset.exists():
            logger.warning(
                f"No action sessions found for IDs: {self.actionsession_ids}"
            )

        self.reaction_ids = [
            actionsession_obj.reaction_id.id
            for actionsession_obj in self.actionsessionqueryset
        ]

        logger.info(f"Found {len(self.reaction_ids)} reactions to process")
        self.groupreactionqueryset = getReactionQuerySet(reaction_ids=self.reaction_ids)

        # Get required resources
        logger.info("Loading equipment resources (tip racks, pipettes, plates)")
        self.tiprackqueryset = self.query_service.get_tip_racks()
        self.pipetteobj = self.query_service.get_pipette()
        self.platequeryset = self.query_service.get_plates()
        self.pipettename = self.pipetteobj.name

        logger.info(
            f"Using pipette: {self.pipettename} with {self.tiprackqueryset.count()} tip racks"
        )
        logger.info(f"Using {self.platequeryset.count()} plates for protocol")

        # Create file path for script
        self.filepath, self.filename = self.file_manager.create_file_path()
        logger.info(f"Script will be written to: {self.filepath}")

    def generate_script(self) -> str:
        """
        Generate the complete OpenTrons protocol script.

        Returns
        -------
        str
            Path to the generated script file
        """
        logger.info(f"Starting script generation for {self.otsessiontype} session")

        try:
            # Setup script basics
            logger.info("Setting up script basics (imports, metadata, labware)")
            self.setup_script()

            # Write the session actions based on session type
            if self.otsessiontype == "reaction":
                logger.info("Processing reaction session actions")
                self.write_reaction_session()
            elif self.otsessiontype == "workup":
                logger.info("Processing workup session actions")
                self.write_workup_session()
            elif self.otsessiontype == "analyse":
                logger.info("Processing analysis session actions")
                self.write_analyse_session()
            else:
                logger.warning(
                    f"Unknown session type: {self.otsessiontype}, no actions generated"
                )

            # Write the file and create OTScript model
            logger.info("Writing script content to file")
            self.file_manager.write_content(self.content)

            logger.info("Creating OTScript database record")
            script_obj = self.file_manager.create_ot_script_model()
            logger.info(f"Created OTScript record with ID: {script_obj.id}")

            logger.info(f"Script generation completed successfully: {self.filepath}")
            return self.filepath

        except Exception as e:
            logger.error(f"Error generating script: {str(e)}")
            logger.error(
                f"Script generation failed for {self.otsessiontype} session {self.otsession_id}"
            )
            raise

    def setup_script(self):
        """Set up the script with initial imports, metadata, and function definitions."""
        # Add imports and metadata
        logger.info("Adding script imports and metadata")
        self.content.extend(
            self.command_generator.get_script_setup(
                protocolname=self.protocolname, apiLevel=self.apiLevel
            )
        )

        # Add labware
        logger.info(f"Setting up {self.platequeryset.count()} labware items")
        self.content.extend(
            self.command_generator.get_labware_setup(
                [
                    {"name": plateobj.name, "labware": plateobj.labware, "index": idx}
                    for idx, plateobj in zip(
                        range(
                            len(self.tiprackqueryset) + 1,
                            len(self.tiprackqueryset) + len(self.platequeryset) + 1,
                        ),
                        self.platequeryset,
                    )
                ]
            )
        )

        # Add tip racks
        logger.info(f"Setting up {self.tiprackqueryset.count()} tip racks")
        self.content.extend(
            self.command_generator.get_tiprack_setup(
                [
                    {
                        "name": tiprack.name,
                        "labware": tiprack.labware,
                        "index": tiprack.index,
                    }
                    for tiprack in self.tiprackqueryset
                ]
            )
        )

        # Add pipettes
        logger.info(f"Setting up {self.pipettename} pipette")
        self.content.extend(
            self.command_generator.get_pipette_setup(
                pipette_name=self.pipetteobj.name,
                pipette_labware=self.pipetteobj.labware,
                pipette_position=self.pipetteobj.position,
                tiprack_names=[tiprack.name for tiprack in self.tiprackqueryset],
            )
        )

        # Set up tip tracking
        logger.info("Setting up tip tracking system")
        self.content.extend(
            self.command_generator.get_number_tips_available_setup(
                num_tipracks=len(self.tiprackqueryset),
                channel_type=self.pipetteobj.type,
            )
        )

        # Add pick up and drop tip functions
        logger.info("Adding tip handling functions")
        self.content.extend(
            self.command_generator.get_pickup_tip_function(
                pipette_name=self.pipettename
            )
        )

        self.content.extend(
            self.command_generator.get_drop_tip_function(pipette_name=self.pipettename)
        )

        logger.info("Script setup completed")

    def write_reaction_session(self):
        """Write the reaction session actions."""
        if self.actionsessionqueryset.exists():
            action_count = self.actionsessionqueryset.count()
            logger.info(f"Processing {action_count} reaction action sessions")
            self.reaction_handler.process_session(self.actionsessionqueryset)
            logger.info("Reaction session processing completed")
        else:
            logger.warning("No reaction sessions found to process")

    def write_workup_session(self):
        """Write the workup session actions."""
        if self.actionsessionqueryset.exists():
            action_count = self.actionsessionqueryset.count()
            logger.info(f"Processing {action_count} workup action sessions")
            self.workup_handler.process_session(self.actionsessionqueryset)
            logger.info("Workup session processing completed")
        else:
            logger.warning("No workup sessions found to process")

    def write_analyse_session(self):
        """Write the analysis session actions."""
        if self.actionsessionqueryset.exists():
            action_count = self.actionsessionqueryset.count()
            logger.info(f"Processing {action_count} analysis action sessions")
            self.analysis_handler.process_session(self.actionsessionqueryset)
            logger.info("Analysis session processing completed")
        else:
            logger.warning("No analysis sessions found to process")

    def add_command(self, commands: Union[str, List[str]]):
        """
        Add commands to the script content.

        Parameters
        ----------
        commands : Union[str, List[str]]
            Commands to add to script
        """
        if isinstance(commands, str):
            self.content.append(commands)
        elif isinstance(commands, list):
            self.content.extend(commands)
