"""
Parsers for design matrices and template files.

Supports:
- CSV design matrices
- Excel design matrices
- JSON/YAML template definitions
"""

import json
import logging
from pathlib import Path
from typing import Any, Optional, Union
import pandas as pd

from .template import RecipeTemplate
from .exceptions import DesignMatrixError, TemplateValidationError

logger = logging.getLogger(__name__)


class DesignMatrixParser:
    """
    Parser for design matrix files (CSV, Excel).
    
    Handles validation and type conversion of design matrix data.
    """
    
    # Column type hints for automatic conversion
    NUMERIC_SUFFIXES = ["_equiv", "_temp", "_temperature", "_time", "_duration", "_conc", "_concentration"]
    ORDER_SUFFIXES = ["_order", "_action_order"]
    
    def __init__(self, template: RecipeTemplate = None):
        """
        Initialize parser with optional template for validation.
        
        Args:
            template: Optional RecipeTemplate for column validation
        """
        self.template = template
    
    def parse_csv(
        self, 
        filepath: Union[str, Path],
        validate: bool = True,
        **pandas_kwargs
    ) -> pd.DataFrame:
        """
        Parse a CSV design matrix file.
        
        Args:
            filepath: Path to CSV file
            validate: Whether to validate columns against template
            **pandas_kwargs: Additional arguments for pd.read_csv
            
        Returns:
            Parsed and validated DataFrame
        """
        try:
            df = pd.read_csv(filepath, **pandas_kwargs)
        except Exception as e:
            raise DesignMatrixError(f"Failed to read CSV file: {e}")
        
        return self._process_dataframe(df, validate)
    
    def parse_excel(
        self,
        filepath: Union[str, Path],
        sheet_name: Union[str, int] = 0,
        validate: bool = True,
        **pandas_kwargs
    ) -> pd.DataFrame:
        """
        Parse an Excel design matrix file.
        
        Args:
            filepath: Path to Excel file
            sheet_name: Sheet name or index to read
            validate: Whether to validate columns against template
            **pandas_kwargs: Additional arguments for pd.read_excel
            
        Returns:
            Parsed and validated DataFrame
        """
        try:
            df = pd.read_excel(filepath, sheet_name=sheet_name, **pandas_kwargs)
        except Exception as e:
            raise DesignMatrixError(f"Failed to read Excel file: {e}")
        
        return self._process_dataframe(df, validate)
    
    def parse_dict_list(
        self,
        data: list[dict],
        validate: bool = True
    ) -> pd.DataFrame:
        """
        Parse a list of dictionaries as design matrix.
        
        Args:
            data: List of experiment parameter dictionaries
            validate: Whether to validate columns against template
            
        Returns:
            DataFrame representation
        """
        df = pd.DataFrame(data)
        return self._process_dataframe(df, validate)
    
    def _process_dataframe(self, df: pd.DataFrame, validate: bool) -> pd.DataFrame:
        """Process and optionally validate a DataFrame."""
        # Clean column names
        df.columns = [str(col).strip() for col in df.columns]
        
        # Apply type conversions
        df = self._convert_types(df)
        
        # Validate if template provided and validation requested
        if validate and self.template:
            self._validate_columns(df)
        
        return df
    
    def _convert_types(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply automatic type conversions based on column naming conventions."""
        for col in df.columns:
            # Convert numeric columns
            if any(col.lower().endswith(suffix) for suffix in self.NUMERIC_SUFFIXES):
                df[col] = pd.to_numeric(df[col], errors="coerce")
            
            # Ensure order columns are strings
            if any(col.lower().endswith(suffix) for suffix in self.ORDER_SUFFIXES):
                df[col] = df[col].astype(str).replace("nan", "default")
        
        return df
    
    def _validate_columns(self, df: pd.DataFrame) -> None:
        """Validate that DataFrame columns match template variables."""
        template_vars = set(self.template.variables.keys())
        df_cols = set(df.columns)
        
        # Check for unknown columns (warning only)
        # Exclude common metadata columns
        metadata_cols = {"experiment_id", "recipe_name", "notes", "order_preset", "session_order"}
        order_cols = {col for col in df_cols if col.endswith("_action_order")}
        
        unknown_cols = df_cols - template_vars - metadata_cols - order_cols
        if unknown_cols:
            logger.warning(
                f"Design matrix contains columns not in template: {unknown_cols}"
            )
        
        # Check for missing required columns (warning only, they may be optional)
        missing_vars = template_vars - df_cols
        if missing_vars:
            logger.info(
                f"Template variables not in design matrix (will use base values): {missing_vars}"
            )


class TemplateParser:
    """
    Parser for recipe template definition files.
    
    Supports JSON and YAML formats.
    """
    
    def parse_json(self, filepath: Union[str, Path]) -> RecipeTemplate:
        """
        Parse a JSON template definition file.
        
        Args:
            filepath: Path to JSON file
            
        Returns:
            RecipeTemplate instance
        """
        try:
            with open(filepath, "r") as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise TemplateValidationError(f"Invalid JSON in template file: {e}")
        except Exception as e:
            raise TemplateValidationError(f"Failed to read template file: {e}")
        
        return RecipeTemplate.from_dict(data)
    
    def parse_yaml(self, filepath: Union[str, Path]) -> RecipeTemplate:
        """
        Parse a YAML template definition file.
        
        Args:
            filepath: Path to YAML file
            
        Returns:
            RecipeTemplate instance
        """
        try:
            import yaml
        except ImportError:
            raise TemplateValidationError(
                "PyYAML is required to parse YAML files. Install with: pip install pyyaml"
            )
        
        try:
            with open(filepath, "r") as f:
                data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise TemplateValidationError(f"Invalid YAML in template file: {e}")
        except Exception as e:
            raise TemplateValidationError(f"Failed to read template file: {e}")
        
        return RecipeTemplate.from_dict(data)
    
    def parse(self, filepath: Union[str, Path]) -> RecipeTemplate:
        """
        Parse a template file, auto-detecting format from extension.
        
        Args:
            filepath: Path to template file (.json or .yaml/.yml)
            
        Returns:
            RecipeTemplate instance
        """
        filepath = Path(filepath)
        suffix = filepath.suffix.lower()
        
        if suffix == ".json":
            return self.parse_json(filepath)
        elif suffix in [".yaml", ".yml"]:
            return self.parse_yaml(filepath)
        else:
            raise TemplateValidationError(
                f"Unknown template file format: {suffix}. Use .json or .yaml"
            )


def create_template_from_recipe(
    reaction_class: str,
    recipe_name: str = "standard",
    template_name: str = None,
    auto_map_actions: bool = True,
) -> RecipeTemplate:
    """
    Helper function to create a template from an existing encoded recipe.
    
    Automatically generates action_ids and session_ids mappings.
    
    Args:
        reaction_class: Name of the reaction class
        recipe_name: Name of the recipe (default: "standard")
        template_name: Name for the template (default: "{reaction_class}-template")
        auto_map_actions: Whether to auto-generate action mappings
        
    Returns:
        RecipeTemplate with auto-generated mappings
    """
    from backend.recipe_utils import (
        recipe_exists,
        recipe_to_dict,
    )
    
    if not recipe_exists(reaction_class):
        raise TemplateValidationError(f"Reaction class '{reaction_class}' not found")
    
    from backend.models import Recipe
    if not Recipe.objects.filter(
        reaction_class=reaction_class, name=recipe_name
    ).exists():
        raise TemplateValidationError(
            f"Recipe '{recipe_name}' not found in '{reaction_class}'"
        )
    
    recipe = recipe_to_dict(reaction_class, recipe_name)
    
    action_ids = {}
    session_ids = {}
    variables = {}
    
    if auto_map_actions:
        for session in recipe.get("sessions", []):
            session_num = session.get("session_number")
            session_type = session.get("session_type", "unknown")
            session_id = f"{session_type}_{session_num}"
            session_ids[session_id] = session_num
            
            for action in session.get("actions", []):
                action_num = action.get("action_number")
                action_type = action.get("type", "unknown")
                
                # Generate descriptive action_id
                if action_type == "add":
                    if action.get("material_smarts"):
                        action_id = f"add_reactant_{action_num}"
                    elif action.get("material_smiles"):
                        action_id = f"add_reagent_{action_num}"
                    else:
                        action_id = f"add_{action_num}"
                else:
                    action_id = f"{action_type}_{action_num}"
                
                action_ids[action_id] = {
                    "session": session_num,
                    "actionnumber": action_num
                }
                
                # Auto-generate common variable mappings
                if action_type == "add":
                    if action.get("solvent"):
                        variables[f"{action_id}_solvent"] = {
                            "action_id": action_id,
                            "path": "solvent"
                        }
                    if action.get("concentration"):
                        variables[f"{action_id}_concentration"] = {
                            "action_id": action_id,
                            "path": "concentration"
                        }
                    if action.get("equivalents") is not None:
                        variables[f"{action_id}_equiv"] = {
                            "action_id": action_id,
                            "path": "equivalents"
                        }
                
                elif action_type == "stir":
                    if action.get("temperature") is not None:
                        variables[f"{action_id}_temperature"] = {
                            "action_id": action_id,
                            "path": "temperature"
                        }
                    if action.get("duration") is not None:
                        variables[f"{action_id}_duration"] = {
                            "action_id": action_id,
                            "path": "duration"
                        }
    
    return RecipeTemplate(
        name=template_name or f"{reaction_class}-template",
        base=f"{reaction_class}.{recipe_name}",
        action_ids=action_ids,
        session_ids=session_ids,
        variables=variables,
        description=f"Auto-generated template for {reaction_class} ({recipe_name})"
    )