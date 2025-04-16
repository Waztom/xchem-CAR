"""
Compatibility layer for the original CreateOTSession class.
Forwards calls to the new modular implementation.
"""

from .session_orchestrator import SessionOrchestrator

class CreateOTSession:
    """
    Legacy compatibility class that forwards to the new modular implementation.
    """
    
    def __init__(self, reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath=None):
        """
        Initialize using the new SessionOrchestrator.
        """
        # Create session using new orchestrator
        orchestrator = SessionOrchestrator()
        self.session = orchestrator.create_session(
            reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath
        )
        
        # Expose session attributes for backward compatibility
        self.reactionstep = self.session.reactionstep
        self.otbatchprotocolobj = self.session.otbatchprotocolobj
        self.actionsessionqueryset = self.session.actionsessionqueryset
        self.customSMcsvpath = self.session.customSMcsvpath
        self.actionsessionnumber = self.session.actionsessionnumber
        self.actionsession_ids = self.session.actionsession_ids
        self.reaction_ids = self.session.reaction_ids
        self.groupreactionqueryset = self.session.groupreactionqueryset
        self.otsessionqueryset = self.session.otsessionqueryset
        self.batchobj = self.session.batchobj
        self.actionsessiontype = self.session.actionsessiontype
        self.otsessionobj = self.session.otsessionobj
        
        # Map other attributes as needed (roundedvolumes, addactionqueryset, etc.)
        if hasattr(self.session, 'deckobj'):
            self.deckobj = self.session.deckobj
        if hasattr(self.session, 'tipracktype'):
            self.tipracktype = self.session.tipracktype
        if hasattr(self.session, 'pipettetype'):
            self.pipettetype = self.session.pipettetype
        # ...map other attributes as needed
    
    # Forward method calls to the session or appropriate manager
    def getActionSessionType(self):
        return self.session._get_action_session_type()
    
    def createOTSessionModel(self):
        return self.session._create_ot_session_model()
    
    # ... forward other method calls as needed