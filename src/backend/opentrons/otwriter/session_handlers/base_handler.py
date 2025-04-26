"""
Base Session Handler for OpenTrons protocols.

This module provides a base class for different session handlers.
"""

import logging
from typing import List, Union
from django.db.models import QuerySet

logger = logging.getLogger(__name__)

class SessionHandler:
    """
    Base class for handling different session types.
    
    This class provides common functionality for all session handlers.
    """
    
    def __init__(self, script_generator):
        """
        Initialize the session handler.
        
        Parameters
        ----------
        script_generator : ScriptGenerator
            Parent script generator
        """
        self.script_generator = script_generator
        self.command_generator = script_generator.command_generator
        self.query_service = script_generator.query_service
        self.volume_manager = script_generator.volume_manager
        self.well_finder = script_generator.well_finder
        
        # Log initialization with class name for better tracking
        handler_type = self.__class__.__name__
        session_type = script_generator.otsessiontype
        session_id = script_generator.otsession_id
        logger.info(f"Initialized {handler_type} for {session_type} session ID {session_id}")
        
    def process_session(self, actionsession_queryset: QuerySet) -> None:
        """
        Process the action session(s).
        
        This method should be implemented by subclasses to handle
        specific session types.
        
        Parameters
        ----------
        actionsession_queryset : QuerySet
            QuerySet of action sessions to process
        """
        logger.error(f"{self.__class__.__name__} failed to implement process_session method")
        raise NotImplementedError("Subclasses must implement process_session method")
        
    def add_command(self, commands: Union[str, List[str]]) -> None:
        """
        Add commands to the script content.
        
        Parameters
        ----------
        commands : Union[str, List[str]]
            Commands to add to script
        """
        # For single line command
        if isinstance(commands, str):
            if commands.strip().startswith('#'):
                # Log comments to show section divisions
                logger.info(f"Adding comment section: {commands.strip()}")
            else:
                command_preview = commands.strip()[:50]
                if len(commands) > 50:
                    command_preview += "..."
                logger.info(f"Adding command: {command_preview}")
        # For multi-line commands
        elif isinstance(commands, list):
            command_count = len(commands)
            if command_count > 0:
                first_command = commands[0].strip()[:50]
                if len(commands[0]) > 50:
                    first_command += "..."
                logger.info(f"Adding {command_count} command(s), starting with: {first_command}")
            else:
                logger.warning("Attempted to add empty command list")
                
        self.script_generator.add_command(commands)
        
    def get_session_number(self, actionsession_queryset: QuerySet) -> int:
        """
        Get the session number from the action session queryset.
        
        Parameters
        ----------
        actionsession_queryset : QuerySet
            QuerySet of action sessions
            
        Returns
        -------
        int
            The session number
        """
        if not actionsession_queryset.exists():
            logger.warning("Attempted to get session number from empty queryset")
            return 0
            
        session_numbers = actionsession_queryset.values_list("sessionnumber", flat=True).distinct()
        session_count = len(session_numbers)
        
        if session_count == 0:
            logger.error("No session numbers found in queryset")
            raise ValueError("No session numbers found in action session queryset")
        elif session_count > 1:
            logger.warning(f"Multiple session numbers found: {list(session_numbers)}, using the first one")
        
        session_number = session_numbers[0]
        logger.info(f"Using session number: {session_number}")
        return session_number
        
    def log_session_start(self, actionsession_queryset: QuerySet, session_type: str) -> None:
        """
        Log the start of session processing.
        
        Parameters
        ----------
        actionsession_queryset : QuerySet
            QuerySet of action sessions
        session_type : str
            Type of session being processed
        """
        action_count = actionsession_queryset.count()
        session_number = self.get_session_number(actionsession_queryset)
        
        logger.info(f"Starting processing of {session_type} session {session_number}")
        logger.info(f"Processing {action_count} action(s) for {session_type} session")
        
        # Log reaction IDs being processed
        if action_count > 0:
            reaction_ids = list(actionsession_queryset.values_list('reaction_id', flat=True).distinct())
            logger.info(f"Session involves {len(reaction_ids)} reaction(s): {reaction_ids}")
        
    def log_session_end(self, session_type: str) -> None:
        """
        Log the end of session processing.
        
        Parameters
        ----------
        session_type : str
            Type of session that was processed
        """
        logger.info(f"Completed processing of {session_type} session")