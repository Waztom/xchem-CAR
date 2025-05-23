"""
OpenTrons Protocol Writer Package.

This package provides a modular approach to generating OpenTrons protocols
from reaction data in the XChem database.
"""

from .script_generator import ScriptGenerator
from .command_generator import CommandGenerator
from .file_manager import FileManager
from .utils.query_service import QueryService
from .utils.volume_manager import VolumeManager
from .utils.well_finder import WellFinder
from .session_handlers import (
    SessionHandler,
    ReactionSessionHandler,
    WorkupSessionHandler,
    AnalysisSessionHandler,
)

__all__ = [
    "ScriptGenerator",
    "CommandGenerator",
    "FileManager",
    "QueryService",
    "VolumeManager",
    "WellFinder",
    "SessionHandler",
    "ReactionSessionHandler",
    "WorkupSessionHandler",
    "AnalysisSessionHandler",
]

__version__ = "1.0.0"
