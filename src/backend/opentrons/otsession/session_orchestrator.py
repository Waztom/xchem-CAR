"""
SessionOrchestrator handles creation and execution of OpenTrons sessions based on action type.
"""

import logging
from django.db.models import QuerySet

from .sessions.base_session import BaseSession
from .sessions.reaction_session import ReactionSession
from .sessions.workup_session import WorkupSession
from .sessions.analysis_session import AnalysisSession
from ...models import ActionSession, OTBatchProtocol

logger = logging.getLogger(__name__)


class SessionOrchestrator:
    """
    Determines the appropriate session type based on action sessions,
    instantiates it, and manages its execution.
    """

    def __init__(
        self,
        reactionstep: int,
        otbatchprotocolobj: OTBatchProtocol,
        actionsessionqueryset: QuerySet[ActionSession],
        customSMcsvpath: str = None,
    ):
        """
        Initialize the orchestrator.

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

        # Determine the session type
        self.actionsessiontype = self._determine_session_type()

        # Initialize session object for appropriate type
        self.session = self._create_session()

        # Reference to session's OT session object for compatibility with existing code
        self.otsessionobj = (
            self.session.otsessionobj if hasattr(self.session, "otsessionobj") else None
        )

        logger.info(
            f"SessionOrchestrator created for {self.actionsessiontype} session, step {reactionstep}"
        )

    def _determine_session_type(self) -> str:
        """
        Determine the appropriate session type based on action sessions.

        Returns
        -------
        session_type: str
            The type of session to create (reaction, workup, or analysis)
        """
        if not self.actionsessionqueryset.exists():
            logger.error("No action sessions provided")
            raise ValueError("No action sessions provided")

        # Get the session type from the action sessions
        session_types = set(
            self.actionsessionqueryset.values_list("type", flat=True).distinct()
        )

        if len(session_types) > 1:
            logger.warning(
                f"Multiple session types found: {session_types}. Using the first one."
            )

        session_type = list(session_types)[0]
        return session_type

    def _create_session(self) -> BaseSession:
        """
        Create an appropriate session based on the action session type.

        Returns
        -------
        session: BaseSession
            The created session instance
        """
        if self.actionsessiontype == "reaction":
            session = ReactionSession(
                reactionstep=self.reactionstep,
                otbatchprotocolobj=self.otbatchprotocolobj,
                actionsessionqueryset=self.actionsessionqueryset,
                customSMcsvpath=self.customSMcsvpath,
            )
        elif self.actionsessiontype == "workup":
            session = WorkupSession(
                reactionstep=self.reactionstep,
                otbatchprotocolobj=self.otbatchprotocolobj,
                actionsessionqueryset=self.actionsessionqueryset,
                customSMcsvpath=self.customSMcsvpath,
            )
        elif self.actionsessiontype in ["analyse", "analysis"]:
            session = AnalysisSession(
                reactionstep=self.reactionstep,
                otbatchprotocolobj=self.otbatchprotocolobj,
                actionsessionqueryset=self.actionsessionqueryset,
                customSMcsvpath=self.customSMcsvpath,
            )
        else:
            logger.error(f"Unknown session type: {self.actionsessiontype}")
            raise ValueError(f"Unknown session type: {self.actionsessiontype}")

        return session

    def execute(self):
        """
        Execute the session protocol.

        Returns
        -------
        success: bool
            True if execution was successful, False otherwise
        """
        try:
            logger.info(f"Executing {self.actionsessiontype} session via orchestrator")
            result = self.session.execute()

            # Update otsessionobj reference after execution
            self.otsessionobj = (
                self.session.otsessionobj
                if hasattr(self.session, "otsessionobj")
                else None
            )

            if self.otsessionobj is None:
                logger.error("Session execution didn't create an OTSession object")
                raise ValueError("Failed to create OTSession object during execution")

            logger.info(f"Session execution complete with result: {result}")
            return result
        except Exception as e:
            logger.error(f"Error executing session: {str(e)}")
            if hasattr(self.session, "cleanup"):
                self.session.cleanup()
            raise
