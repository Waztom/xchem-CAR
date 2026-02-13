"""
Recipe Generator - generates recipe variants from templates and design matrices.
"""

import copy
import logging
from typing import Any, Union
import pandas as pd

from .template import RecipeTemplate, ActionLocation, VariableDefinition
from .exceptions import (
    RecipeGeneratorError,
    ActionNotFoundError,
    PathNotFoundError,
    OrderingError,
)

logger = logging.getLogger(__name__)


class RecipeGenerator:
    """
    Generates encoded recipe variants from a template and design matrix.
    
    The generator handles:
    - Parameter value substitution at action level
    - Action reordering within sessions
    - Session reordering
    - Named ordering presets
    
    Example:
        template = RecipeTemplate(...)
        generator = RecipeGenerator(template)
        
        # From DataFrame
        recipes = generator.from_dataframe(design_df)
        
        # From CSV file
        recipes = generator.from_csv("design_matrix.csv")
        
        # Single recipe with explicit values
        recipe = generator.generate_single(
            solvent="DMF",
            temperature=60,
            action_order=["add_coupling_agent", "add_acid", "add_base"]
        )
    """
    
    def __init__(self, template: RecipeTemplate):
        """
        Initialize generator with a recipe template.
        
        Args:
            template: A validated RecipeTemplate instance
        """
        self.template = template
    
    def from_dataframe(
        self, 
        design_df: pd.DataFrame,
        recipe_name_column: str = None,
    ) -> list[dict]:
        """
        Generate recipes from a pandas DataFrame design matrix.
        
        Args:
            design_df: DataFrame where columns match variable names from template
            recipe_name_column: Optional column to use as recipe names
            
        Returns:
            List of generated recipe dictionaries
        """
        recipes = []
        
        for idx, row in design_df.iterrows():
            # Convert row to dict, excluding NaN values
            params = {k: v for k, v in row.to_dict().items() if pd.notna(v)}
            
            # Generate recipe name
            if recipe_name_column and recipe_name_column in params:
                recipe_name = str(params.pop(recipe_name_column))
            else:
                recipe_name = f"{self.template.name}_{idx}"
            
            try:
                recipe = self._generate_recipe(params, recipe_name)
                recipes.append(recipe)
            except RecipeGeneratorError as e:
                logger.error(f"Error generating recipe for row {idx}: {e}")
                raise
        
        return recipes
    
    def from_csv(
        self, 
        filepath: str,
        recipe_name_column: str = None,
        **pandas_kwargs
    ) -> list[dict]:
        """
        Generate recipes from a CSV design matrix file.
        
        Args:
            filepath: Path to CSV file
            recipe_name_column: Optional column to use as recipe names
            **pandas_kwargs: Additional arguments passed to pd.read_csv
            
        Returns:
            List of generated recipe dictionaries
        """
        design_df = pd.read_csv(filepath, **pandas_kwargs)
        return self.from_dataframe(design_df, recipe_name_column)
    
    def generate_single(self, recipe_name: str = None, **params) -> dict:
        """
        Generate a single recipe with specified parameter values.
        
        Args:
            recipe_name: Optional name for the recipe
            **params: Variable values (must match template variable names)
            
        Returns:
            Generated recipe dictionary
        """
        if recipe_name is None:
            recipe_name = f"{self.template.name}_single"
        return self._generate_recipe(params, recipe_name)
    
    def _generate_recipe(self, params: dict, recipe_name: str) -> dict:
        """
        Internal method to generate a single recipe.
        
        Args:
            params: Dictionary of parameter values
            recipe_name: Name for the generated recipe
            
        Returns:
            Generated recipe dictionary
        """
        # Start with deep copy of base recipe
        recipe = self.template.get_base_recipe()
        
        # Separate ordering params from value params
        action_order_params = {}
        session_order_param = None
        order_preset_param = None
        value_params = {}
        
        for key, value in params.items():
            if key.endswith("_action_order"):
                # Extract session identifier
                session_key = key.replace("_action_order", "")
                action_order_params[session_key] = value
            elif key == "session_order":
                session_order_param = value
            elif key == "order_preset":
                order_preset_param = value
            else:
                value_params[key] = value
        
        # Apply order preset if specified
        if order_preset_param and order_preset_param != "default":
            self._apply_order_preset(recipe, order_preset_param)
        
        # Apply explicit action orderings (override preset)
        for session_key, order_value in action_order_params.items():
            if order_value and order_value != "default":
                action_order = self._parse_order_value(order_value)
                self._reorder_actions_in_session(recipe, session_key, action_order)
        
        # Apply session ordering
        if session_order_param and session_order_param != "default":
            session_order = self._parse_order_value(session_order_param)
            self._reorder_sessions(recipe, session_order)
        
        # Common metadata columns that should be silently ignored
        metadata_columns = {
            "experiment_id", "recipe_name", "notes", "description", 
            "batch", "plate", "well", "row", "column"
        }
        
        # Apply parameter value substitutions
        for var_name, value in value_params.items():
            if var_name in self.template.variables:
                self._apply_variable(recipe, var_name, value)
            elif var_name not in metadata_columns:
                logger.warning(f"Unknown variable '{var_name}' in parameters, skipping")
        
        return {
            "name": recipe_name,
            "template": self.template.name,
            "base": self.template.base,
            "recipe": recipe,
        }
    
    def _apply_variable(self, recipe: dict, var_name: str, value: Any) -> None:
        """
        Apply a variable value to the recipe.
        
        Args:
            recipe: Recipe dict to modify
            var_name: Variable name from template
            value: Value to set
        """
        var_def = self.template.get_variable_definition(var_name)
        action_loc = self.template.get_action_location(var_def.action_id)
        
        # Find the action
        action = self._find_action(recipe, action_loc)
        if action is None:
            raise ActionNotFoundError(
                var_def.action_id, action_loc.session, action_loc.actionnumber
            )
        
        # Navigate to the path and set value
        self._set_nested_value(action, var_def.path, value, var_def.action_id)
    
    def _find_action(self, recipe: dict, location: ActionLocation) -> dict:
        """Find an action in the recipe by its location."""
        # Find session
        session = None
        for s in recipe.get("actionsessions", []):
            if s.get("sessionnumber") == location.session:
                session = s
                break
        
        if session is None:
            return None
        
        # Find action in session (check all possible action containers)
        for actions_key in ["intermolecular", "Intramolecular", "actions"]:
            if actions_key in session:
                actions_container = session[actions_key]
                # Handle nested "actions" key
                if isinstance(actions_container, dict) and "actions" in actions_container:
                    actions = actions_container["actions"]
                else:
                    actions = actions_container
                
                if isinstance(actions, list):
                    for action in actions:
                        if action.get("actionnumber") == location.actionnumber:
                            return action
        
        return None
    
    def _set_nested_value(
        self, 
        obj: dict, 
        path: str, 
        value: Any, 
        action_id: str
    ) -> None:
        """
        Set a value at a nested path within an object.
        
        Supports dot notation and array indices: "content.material.solvent"
        """
        parts = self._parse_path(path)
        current = obj
        
        try:
            for i, part in enumerate(parts[:-1]):
                if isinstance(part, int):
                    current = current[part]
                else:
                    current = current[part]
            
            # Set the final value
            final_key = parts[-1]
            if isinstance(final_key, int):
                current[final_key] = value
            else:
                current[final_key] = value
                
        except (KeyError, IndexError, TypeError) as e:
            raise PathNotFoundError(path, action_id)
    
    def _parse_path(self, path: str) -> list[Union[str, int]]:
        """Parse a dot-notation path into components."""
        parts = []
        for part in path.split("."):
            if "[" in part:
                # Handle array index like "actions[0]"
                name = part[:part.index("[")]
                index = int(part[part.index("[")+1:part.index("]")])
                parts.append(name)
                parts.append(index)
            else:
                parts.append(part)
        return parts
    
    def _parse_order_value(self, order_value: Any) -> list[str]:
        """Parse an order value into a list of IDs."""
        if isinstance(order_value, list):
            return order_value
        elif isinstance(order_value, str):
            # Handle comma-separated string
            return [x.strip() for x in order_value.split(",")]
        else:
            raise OrderingError(f"Invalid order value type: {type(order_value)}")
    
    def _apply_order_preset(self, recipe: dict, preset_name: str) -> None:
        """Apply a named ordering preset to the recipe."""
        if preset_name not in self.template.order_presets:
            raise OrderingError(f"Unknown order preset: '{preset_name}'")
        
        preset = self.template.order_presets[preset_name]
        
        # Apply action orderings
        if "action_orders" in preset:
            for session_key, action_order in preset["action_orders"].items():
                self._reorder_actions_in_session(recipe, session_key, action_order)
        
        # Apply session ordering
        if "session_order" in preset and preset["session_order"]:
            self._reorder_sessions(recipe, preset["session_order"])
    
    def _reorder_actions_in_session(
        self, 
        recipe: dict, 
        session_identifier: str, 
        new_order: list[str]
    ) -> None:
        """
        Reorder actions within a session.
        
        Args:
            recipe: Recipe to modify
            session_identifier: Either a session_id from template or "session_N" format
            new_order: List of action_ids in desired order
        """
        # Resolve session number
        if session_identifier in self.template.session_ids:
            session_num = self.template.session_ids[session_identifier]
        elif session_identifier.startswith("session_"):
            session_num = int(session_identifier.replace("session_", ""))
        else:
            raise OrderingError(f"Unknown session identifier: '{session_identifier}'")
        
        # Find session
        session = None
        for s in recipe.get("actionsessions", []):
            if s.get("sessionnumber") == session_num:
                session = s
                break
        
        if session is None:
            raise OrderingError(f"Session {session_num} not found in recipe")
        
        # Find the actions container
        actions_key = None
        actions_container = None
        for key in ["intermolecular", "Intramolecular", "actions"]:
            if key in session:
                if isinstance(session[key], dict) and "actions" in session[key]:
                    actions_key = key
                    actions_container = session[key]["actions"]
                elif isinstance(session[key], list):
                    actions_key = key
                    actions_container = session[key]
                if actions_container:
                    break
        
        if actions_container is None:
            raise OrderingError(f"No actions found in session {session_num}")
        
        # Build map of action_id -> action
        action_id_to_action = {}
        for action_id, loc in self.template.action_ids.items():
            if loc["session"] == session_num:
                for action in actions_container:
                    if action.get("actionnumber") == loc["actionnumber"]:
                        action_id_to_action[action_id] = action
                        break
        
        # Validate all IDs in new_order exist
        for action_id in new_order:
            if action_id not in action_id_to_action:
                raise OrderingError(
                    f"Action '{action_id}' not found in session {session_num}"
                )
        
        # Reorder
        reordered = [action_id_to_action[aid] for aid in new_order]
        
        # Include any actions not in new_order at the end (preserve them)
        ordered_actions_set = set(new_order)
        for action_id, action in action_id_to_action.items():
            if action_id not in ordered_actions_set:
                reordered.append(action)
        
        # Renumber actions sequentially
        for i, action in enumerate(reordered, 1):
            action["actionnumber"] = i
        
        # Update recipe
        if actions_key in ["intermolecular", "Intramolecular"]:
            session[actions_key]["actions"] = reordered
        else:
            session["actions"] = reordered
    
    def _reorder_sessions(self, recipe: dict, new_order: list[str]) -> None:
        """
        Reorder sessions in the recipe.
        
        Args:
            recipe: Recipe to modify
            new_order: List of session_ids in desired order
        """
        current_sessions = recipe.get("actionsessions", [])
        
        # Build map of session_id -> session
        session_id_to_session = {}
        for session_id, session_num in self.template.session_ids.items():
            for session in current_sessions:
                if session.get("sessionnumber") == session_num:
                    session_id_to_session[session_id] = session
                    break
        
        # Validate all IDs in new_order exist
        for session_id in new_order:
            if session_id not in session_id_to_session:
                raise OrderingError(f"Session '{session_id}' not found")
        
        # Reorder
        reordered = [session_id_to_session[sid] for sid in new_order]
        
        # Include any sessions not in new_order at the end
        ordered_sessions_set = set(new_order)
        for session_id, session in session_id_to_session.items():
            if session_id not in ordered_sessions_set:
                reordered.append(session)
        
        # Also include any sessions not in session_ids mapping
        mapped_session_nums = set(self.template.session_ids.values())
        for session in current_sessions:
            if session.get("sessionnumber") not in mapped_session_nums:
                if session not in reordered:
                    reordered.append(session)
        
        # Renumber sessions sequentially
        for i, session in enumerate(reordered, 1):
            session["sessionnumber"] = i
        
        recipe["actionsessions"] = reordered
    
    def get_variable_names(self) -> list[str]:
        """Get list of all variable names defined in the template."""
        return list(self.template.variables.keys())
    
    def get_required_columns(self) -> list[str]:
        """
        Get list of column names expected in a design matrix.
        
        Returns variable names plus optional ordering columns.
        """
        columns = list(self.template.variables.keys())
        
        # Add ordering columns based on template
        for session_id in self.template.session_ids:
            columns.append(f"{session_id}_action_order")
        
        if self.template.session_ids:
            columns.append("session_order")
        
        if self.template.order_presets:
            columns.append("order_preset")
        
        return columns
    
    def generate_design_matrix_template(self) -> pd.DataFrame:
        """
        Generate an empty design matrix template DataFrame.
        
        Useful for chemists to fill in with their experimental design.
        """
        columns = ["experiment_id"] + self.get_required_columns()
        return pd.DataFrame(columns=columns)
