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
from .managers.pipette_manager import PipetteManager
from .managers.material_manager import MaterialManager
from .managers.data_manager import DataManager

# Import refactored plate manager components
from .managers.plate_manager.plate_factory import PlateFactory
from .managers.plate_manager.well_manager import WellManager
from .managers.plate_manager.column_manager import ColumnManager
from .managers.plate_manager.labware_selector import LabwareSelector
from .managers.plate_manager.plate_query_service import PlateQueryService

# Export symbols
__all__ = [
    "SessionOrchestrator",
    "BaseSession",
    "ReactionSession",
    "WorkupSession",
    "AnalysisSession",
    "DeckManager",
    "PlateFactory",
    "WellManager",
    "ColumnManager", 
    "LabwareSelector",
    "PlateQueryService",
    "PipetteManager",
    "MaterialManager",
    "DataManager",
    "CreateOTSession",
]
