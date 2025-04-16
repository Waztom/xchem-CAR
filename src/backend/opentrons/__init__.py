"""
OpenTrons automation package for liquid handling operations.
"""

# Import new classes
from .session_orchestrator import SessionOrchestrator
from .sessions.base_session import BaseSession
from .sessions.reaction_session import ReactionSession
from .sessions.workup_session import WorkupSession
from .sessions.analysis_session import AnalysisSession

from .managers.deck_manager import DeckManager
from .managers.plate_manager import PlateManager
from .managers.pipette_manager import PipetteManager
from .managers.material_manager import MaterialManager
from .managers.data_manager import DataManager

# Import compatibility layer
from .otsession_compat import CreateOTSession

# Export symbols
__all__ = [
    'SessionOrchestrator',
    'BaseSession',
    'ReactionSession',
    'WorkupSession',
    'AnalysisSession',
    'DeckManager',
    'PlateManager',
    'PipetteManager',
    'MaterialManager',
    'DataManager',
    'CreateOTSession',
]