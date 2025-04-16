"""
Base session class with common functionality for all OpenTrons sessions.
"""

import logging
from abc import ABC, abstractmethod

from ..managers.deck_manager import DeckManager
from ..managers.plate_manager import PlateManager
from ..managers.pipette_manager import PipetteManager
from ..managers.material_manager import MaterialManager
from ..managers.data_manager import DataManager
from ...models import OTSession
from ...utils import getReactionQuerySet

logger = logging.getLogger(__name__)

class BaseSession(ABC):
    """
    Abstract base class for OpenTrons session types.
    Provides common functionality and defines the interface for all session types.
    """
    
    def __init__(self, reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath=None):
        """
        Initialize the session.
        
        Parameters
        ----------
        reactionstep: int
            The reaction step number
        otbatchprotocolobj: OTBatchProtocol
            The batch protocol this session is part of
        actionsessionqueryset: QuerySet[ActionSession]
            The action sessions to execute
        customSMcsvpath: str, optional
            Path to custom starting material CSV file
        """
        self.reactionstep = reactionstep
        self.otbatchprotocolobj = otbatchprotocolobj
        self.actionsessionqueryset = actionsessionqueryset
        self.customSMcsvpath = customSMcsvpath
        
        # Initialize session properties
        self.batchobj = otbatchprotocolobj.batch_id
        self.actionsession_ids = list(actionsessionqueryset.values_list("id", flat=True))
        self.reaction_ids = list(actionsessionqueryset.values_list("reaction_id", flat=True).distinct())
        self.actionsessiontype = actionsessionqueryset.values_list("type", flat=True).distinct()[0]
        
        # Properties that will be set during setup/execution
        self.otsessionobj = None
        self.deckobj = None
        self.groupreactionqueryset = None
        self.tipracktype = None
        self.pipettetype = None
        self.roundedaddvolumes = []
        self.roundedextractvolumes = []
        self.roundedvolumes = []
        self.addactionsdf = None
        self.addactionqueryset = None
        self.extractactionqueryset = None
        
        # Manager instances will be created in setup_common_resources
        self.deck_manager = None
        self.plate_manager = None
        self.pipette_manager = None
        self.material_manager = None
        self.data_manager = None
        
        # Session initialization flag
        self.is_initialized = False
    
    def create_ot_session_model(self):
        """
        Create an OT session model in the database.
        
        Returns
        -------
        otsessionobj: OTSession
            The created session object
        """
        try:
            otsessionobj = OTSession()
            otsessionobj.otbatchprotocol_id = self.otbatchprotocolobj
            otsessionobj.actionsessiontype = self.actionsessiontype
            otsessionobj.reactionstep = self.reactionstep
            otsessionobj.save()
            
            logger.info(f"Created OT session: {otsessionobj.id} for step {self.reactionstep}")
            return otsessionobj
            
        except Exception as e:
            logger.error(f"Error creating OT session model: {str(e)}")
            raise
    
    def get_action_session_type(self):
        """
        Get the type of action session.
        
        Returns
        -------
        action_session_type: str
            The type of action session (reaction, workup, analyse)
        """
        return self.actionsessiontype
    
    def setup_common_resources(self):
        """
        Initialize managers and set up resources needed by all session types.
        """
        try:
            # Create session model first
            self.otsessionobj = self.create_ot_session_model()
            
            # Initialize managers - order matters here due to dependencies
            self.data_manager = DataManager(self)
            self.deck_manager = DeckManager(self)
            self.plate_manager = PlateManager(self)
            self.material_manager = MaterialManager(self)
            self.pipette_manager = PipetteManager(self)
            
            # Create deck
            self.deckobj = self.deck_manager.create_deck_model()
            
            # Get reaction queryset for this session
            self.groupreactionqueryset = getReactionQuerySet(reaction_ids=self.reaction_ids)
            
            self.is_initialized = True
            logger.info(f"Common resources set up for {self.actionsessiontype} session")
            
        except Exception as e:
            logger.error(f"Error setting up common resources: {str(e)}")
            self.cleanup()
            raise
    
    def cleanup(self):
        """
        Clean up any created database entries if an error occurs.
        """
        try:
            if hasattr(self, "otsessionobj") and self.otsessionobj:
                logger.info(f"Cleaning up session {self.otsessionobj.id} due to error")
                self.otsessionobj.delete()  # Will cascade delete related objects
        except Exception as e:
            logger.error(f"Error during cleanup: {str(e)}")
    
    @abstractmethod
    def execute(self):
        """
        Execute the session protocol. Must be implemented by subclasses.
        
        Returns
        -------
        success: bool
            True if execution was successful, False otherwise
        """
        pass