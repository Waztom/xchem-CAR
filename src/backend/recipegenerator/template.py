"""
Recipe Template definition and validation.

A RecipeTemplate defines:
- A base recipe to inherit from (reaction class + recipe name)
- Action identifiers mapping semantic names to action locations
- Session identifiers for session reordering
- Variable definitions mapping parameter names to action paths
- Optional ordering presets for common addition sequences
"""

from dataclasses import dataclass, field
from typing import Any, Optional
import copy
import logging

from backend.recipe_utils import (
    recipe_exists,
    get_latest_recipe,
    recipe_to_dict,
)
from .exceptions import (
    TemplateValidationError,
    ActionNotFoundError,
    SessionNotFoundError,
)

logger = logging.getLogger(__name__)


@dataclass
class ActionLocation:
    """Identifies an action's location within a recipe."""
    session: int  # sessionnumber (1-indexed)
    actionnumber: int  # actionnumber within the session
    
    def __post_init__(self):
        if self.session < 1:
            raise ValueError("Session number must be >= 1")
        if self.actionnumber < 1:
            raise ValueError("Action number must be >= 1")


@dataclass
class VariableDefinition:
    """Defines a variable parameter within the template."""
    action_id: str  # References an action in action_ids
    path: str  # Dot-notation path to the value (e.g., "content.material.solvent")
    
    def get_path_parts(self) -> list[str]:
        """Split path into components, handling array indices."""
        parts = []
        for part in self.path.split("."):
            # Handle array indices like "actions[0]"
            if "[" in part:
                name = part[:part.index("[")]
                index = int(part[part.index("[")+1:part.index("]")])
                parts.append(name)
                parts.append(index)
            else:
                parts.append(part)
        return parts


@dataclass
class OrderPreset:
    """Defines a named ordering preset."""
    name: str
    description: Optional[str] = None
    action_orders: dict[str, list[str]] = field(default_factory=dict)
    # Maps session_id to ordered list of action_ids within that session
    session_order: Optional[list[str]] = None
    # Optional reordering of sessions themselves


@dataclass
class RecipeTemplate:
    """
    A recipe template for generating encoded recipe variants.
    
    Attributes:
        name: Template name (e.g., "Amidation-screening")
        base: Base recipe in format "ReactionClass.recipe_name"
        action_ids: Maps semantic action names to their locations
        session_ids: Maps semantic session names to session numbers
        variables: Maps parameter names to their action paths
        order_presets: Named ordering configurations
    """
    name: str
    base: str  # Format: "ReactionClass.recipe_name" or just "ReactionClass" for standard
    action_ids: dict[str, dict[str, int]]  # {"add_acid": {"session": 1, "actionnumber": 1}}
    session_ids: dict[str, int] = field(default_factory=dict)  # {"reaction_1": 1, "stir_1": 2}
    variables: dict[str, dict[str, str]] = field(default_factory=dict)
    order_presets: dict[str, dict] = field(default_factory=dict)
    description: Optional[str] = None
    
    def __post_init__(self):
        """Validate the template on creation."""
        self._validate()
    
    def _validate(self) -> None:
        """Validate template structure and references."""
        # Parse and validate base recipe reference
        self._parse_base()
        
        # Validate action_ids structure
        for action_id, location in self.action_ids.items():
            if "session" not in location or "actionnumber" not in location:
                raise TemplateValidationError(
                    f"Action '{action_id}' must have 'session' and 'actionnumber' keys"
                )
            # Validate action exists in base recipe
            self._validate_action_exists(action_id, location)
        
        # Validate session_ids
        for session_id, sessionnumber in self.session_ids.items():
            self._validate_session_exists(session_id, sessionnumber)
        
        # Validate variables reference valid action_ids
        for var_name, var_def in self.variables.items():
            if "action_id" not in var_def:
                raise TemplateValidationError(
                    f"Variable '{var_name}' must have 'action_id' key"
                )
            if var_def["action_id"] not in self.action_ids:
                raise TemplateValidationError(
                    f"Variable '{var_name}' references unknown action_id '{var_def['action_id']}'"
                )
            if "path" not in var_def:
                raise TemplateValidationError(
                    f"Variable '{var_name}' must have 'path' key"
                )
    
    def _parse_base(self) -> tuple[str, str]:
        """Parse base recipe reference into (reaction_class, recipe_name)."""
        if "." in self.base:
            parts = self.base.split(".", 1)
            reaction_class, recipe_name = parts[0], parts[1]
        else:
            reaction_class = self.base
            recipe_name = "standard"
        
        # Validate reaction class exists in Recipe DB
        if not recipe_exists(reaction_class):
            raise TemplateValidationError(
                f"Reaction class '{reaction_class}' not found in Recipe DB"
            )
        
        # Validate recipe exists
        from backend.models import Recipe
        if not Recipe.objects.filter(
            reaction_class=reaction_class, name=recipe_name
        ).exists():
            raise TemplateValidationError(
                f"Recipe '{recipe_name}' not found in reaction class '{reaction_class}'"
            )
        
        return reaction_class, recipe_name
    
    def get_base_recipe(self) -> dict:
        """Get a deep copy of the base recipe as a clean dict.

        The structure mirrors the DB schema (not the legacy encoded format)::

            {"sessions": [{..."actions": [{...}]}], "estimated_yield": ..., ...}
        """
        reaction_class, recipe_name = self._parse_base()
        return copy.deepcopy(
            recipe_to_dict(reaction_class, recipe_name)
        )
    
    def get_recipe(self):
        """Get the Recipe model instance for this template's base recipe."""
        reaction_class, recipe_name = self._parse_base()
        return get_latest_recipe(reaction_class, recipe_name)
    
    def get_reaction_class(self) -> str:
        """Get the reaction class name."""
        reaction_class, _ = self._parse_base()
        return reaction_class
    
    def get_recipe_name(self) -> str:
        """Get the base recipe name."""
        _, recipe_name = self._parse_base()
        return recipe_name
    
    def _validate_action_exists(self, action_id: str, location: dict) -> None:
        """Validate that an action exists at the specified location in the DB."""
        recipe = self.get_recipe()
        session_num = location["session"]
        action_num = location["actionnumber"]
        
        session = recipe.action_sessions.filter(session_number=session_num).first()
        if session is None:
            raise ActionNotFoundError(action_id, session_num, action_num)
        
        # Check all action types for this action_number
        found = (
            session.add_actions.filter(action_number=action_num).exists()
            or session.stir_actions.filter(action_number=action_num).exists()
            or session.extract_actions.filter(action_number=action_num).exists()
            or session.mix_actions.filter(action_number=action_num).exists()
        )
        if not found:
            raise ActionNotFoundError(action_id, session_num, action_num)
    
    def _validate_session_exists(self, session_id: str, sessionnumber: int) -> None:
        """Validate that a session exists at the specified number in the DB."""
        recipe = self.get_recipe()
        if not recipe.action_sessions.filter(session_number=sessionnumber).exists():
            raise SessionNotFoundError(session_id, sessionnumber)
    
    def get_action_location(self, action_id: str) -> ActionLocation:
        """Get the location of an action by its ID."""
        if action_id not in self.action_ids:
            raise ActionNotFoundError(action_id)
        loc = self.action_ids[action_id]
        return ActionLocation(session=loc["session"], actionnumber=loc["actionnumber"])
    
    def get_variable_definition(self, var_name: str) -> VariableDefinition:
        """Get a variable definition by name."""
        if var_name not in self.variables:
            raise TemplateValidationError(f"Variable '{var_name}' not defined")
        var_def = self.variables[var_name]
        return VariableDefinition(action_id=var_def["action_id"], path=var_def["path"])
    
    def get_actions_in_session(self, session_num: int) -> list[str]:
        """Get all action IDs belonging to a specific session."""
        return [
            action_id for action_id, loc in self.action_ids.items()
            if loc["session"] == session_num
        ]
    
    def to_dict(self) -> dict:
        """Serialize template to dictionary."""
        return {
            "name": self.name,
            "base": self.base,
            "description": self.description,
            "action_ids": self.action_ids,
            "session_ids": self.session_ids,
            "variables": self.variables,
            "order_presets": self.order_presets,
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> "RecipeTemplate":
        """Create template from dictionary."""
        return cls(
            name=data["name"],
            base=data["base"],
            description=data.get("description"),
            action_ids=data.get("action_ids", {}),
            session_ids=data.get("session_ids", {}),
            variables=data.get("variables", {}),
            order_presets=data.get("order_presets", {}),
        )
