"""
Recipe Generator Module

This module provides tools for generating encoded recipes from templates
and design matrices, enabling chemists to easily create multiple recipe
variants for experimental designs (full factorial, Plackett-Burman, 
Bayesian optimization, etc.)

Example usage:
    from backend.recipegenerator import RecipeTemplate, RecipeGenerator
    
    # Define template
    template = RecipeTemplate(
        name="Amidation-screening",
        base="Amidation.standard",
        action_ids={
            "add_acid": {"session": 1, "actionnumber": 1},
            "add_coupling_agent": {"session": 1, "actionnumber": 2},
        },
        variables={
            "solvent": {"action_id": "add_acid", "path": "content.material.solvent"},
        }
    )
    
    # Generate recipes from design matrix
    generator = RecipeGenerator(template)
    recipes = generator.from_csv("design_matrix.csv")
"""

from .template import RecipeTemplate
from .generator import RecipeGenerator
from .parsers import DesignMatrixParser
from .exceptions import (
    RecipeGeneratorError,
    TemplateValidationError,
    DesignMatrixError,
    ActionNotFoundError,
    SessionNotFoundError,
)

__all__ = [
    "RecipeTemplate",
    "RecipeGenerator",
    "DesignMatrixParser",
    "RecipeGeneratorError",
    "TemplateValidationError",
    "DesignMatrixError",
    "ActionNotFoundError",
    "SessionNotFoundError",
]
