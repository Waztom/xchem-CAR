"""Orchestrates different OpenTrons session types."""

import logging
from .sessions.reaction_session import ReactionSession
from .sessions.workup_session import WorkupSession
from .sessions.analysis_session import AnalysisSession

logger = logging.getLogger(__name__)

class SessionOrchestrator:
    """
    Orchestrates the creation and execution of different OT session types.
    Acts as the main entry point for session creation.
    """
    
    def create_session(self, reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath=None):
        """
        Factory method to create the appropriate session type based on actionsessionqueryset.
        
        Parameters
        ----------
        reactionstep: int
            The reaction step that the protocol is being created for
        otbatchprotocolobj: OTBatchProtocol
            The Django OTBatchProtocol model object
        actionsessionqueryset: QuerySet[ActionSession]
            The action sessions being executed on the OT
        customSMcsvpath: str (Optional)
            The path to the custom starting material plate csv file
            
        Returns
        -------
        session: BaseSession
            The appropriate session object (Reaction, Workup, or Analysis)
        """
        # Determine session type
        action_session_type = actionsessionqueryset.values_list("type", flat=True).distinct()[0]
        logger.info(f"Creating {action_session_type} session for reaction step {reactionstep}")
        
        # Create and return appropriate session type
        if action_session_type == "reaction":
            session = ReactionSession(reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath)
            session.execute()
            return session
        elif action_session_type == "workup":
            session = WorkupSession(reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath)
            session.execute()
            return session
        elif action_session_type == "analyse":
            session = AnalysisSession(reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath)
            session.execute()
            return session
        else:
            error_msg = f"Unsupported session type: {action_session_type}"
            logger.error(error_msg)
            raise ValueError(error_msg)