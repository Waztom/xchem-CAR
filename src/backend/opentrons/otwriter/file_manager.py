"""
File Manager for OpenTrons protocol scripts.

This module handles file operations for script generation including
creating file paths, writing content to files, and creating database entries.
"""

import os
import logging
from typing import List, Tuple
from django.core.files.storage import default_storage
from django.conf import settings

from backend.models import OTScript

logger = logging.getLogger(__name__)

class FileManager:
    """
    Manages file operations for OpenTrons protocol scripts.
    """
    
    def __init__(self, script_generator):
        """
        Initialize the file manager.
        
        Parameters
        ----------
        script_generator : ScriptGenerator
            Parent script generator
        """
        self.script_generator = script_generator
        self.batchtag = script_generator.batchtag
        self.otsessionobj = script_generator.otsessionobj
        self.protocolname = script_generator.protocolname
        
        # Create file path
        self.filepath, self.filename = self.create_file_path()
        
        logger.debug(f"File manager initialized with path: {self.filepath}")
    
    def create_file_path(self) -> Tuple[str, str]:
        """
        Create file path for the protocol script.
        
        Returns
        -------
        Tuple[str, str]
            Tuple of (filepath, filename)
        """
        filename = f"{self.protocolname}.txt"
        path = f"tmp/{filename}"
        filepath = str(os.path.join(settings.MEDIA_ROOT, path))
        
        logger.debug(f"Created file path: {filepath}")
        return filepath, filename
    
    def write_content(self, content: List[str]) -> str:
        """
        Write content to protocol script file.
        
        Parameters
        ----------
        content : List[str]
            Lines of content to write
        
        Returns
        -------
        str
            Path to the written file
        """
        try:
            with open(self.filepath, 'w') as f:
                for line in content:
                    f.write(f"{line}\n")
            
            logger.info(f"Protocol script written to: {self.filepath}")
            return self.filepath
            
        except Exception as e:
            logger.error(f"Error writing protocol script: {str(e)}")
            raise
    
    def create_ot_script_model(self) -> OTScript:
        """
        Create OTScript model in database.
        
        Returns
        -------
        OTScript
            Created OTScript object
        """
        try:
            # Ensure file exists and is ready for reading
            with open(self.filepath, "a") as script:
                pass  # Just ensuring file exists
            
            # Create OTScript object
            otscriptobj = OTScript()
            otscriptobj.otsession_id = self.otsessionobj
            
            # Read file and save to storage
            with open(self.filepath, "rb") as otscriptfile:
                otscriptfn = default_storage.save(
                    f"otscripts/{self.filename.strip('.txt')}.py", 
                    otscriptfile
                )
            
            # Save path to database
            otscriptobj.otscript = otscriptfn
            otscriptobj.save()
            
            logger.info(f"OTScript model created: {otscriptobj.id}")
            return otscriptobj
            
        except Exception as e:
            logger.error(f"Error creating OTScript model: {str(e)}")
            raise