"""
Unit tests for the recipegenerator module.

Uses Django TestCase for tests that need database access (recipe lookups).
"""

from unittest import TestCase
import json
import tempfile
import os
from pathlib import Path

import pandas as pd

from django.test import TestCase as DjangoTestCase

from backend.recipegenerator import (
    RecipeTemplate,
    RecipeGenerator,
    DesignMatrixParser,
    TemplateValidationError,
    ActionNotFoundError,
)
from backend.recipegenerator.parsers import create_template_from_recipe

# Re-use fixture ingest helper from test_recipe_models
from .test_recipe_models import ingest_recipe_from_json

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


class RecipeTemplateTestCase(DjangoTestCase):
    """Tests for RecipeTemplate class."""
    
    @classmethod
    def setUpTestData(cls):
        """Load amidation fixture into test DB."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)
    
    def setUp(self):
        """Set up test fixtures."""
        # Use Amidation.standard as test base - it has a well-defined structure
        self.valid_template_data = {
            "name": "test-template",
            "base": "Amidation.standard",
            "action_ids": {
                "add_acid": {"session": 1, "actionnumber": 1},
                "add_coupling_agent": {"session": 1, "actionnumber": 2},
                "add_base": {"session": 1, "actionnumber": 3},
            },
            "session_ids": {
                "reaction_1": 1,
                "reaction_2": 2,
                "stir_1": 3,
            },
            "variables": {
                "acid_solvent": {
                    "action_id": "add_acid",
                    "path": "solvent"
                },
                "acid_equiv": {
                    "action_id": "add_acid",
                    "path": "equivalents"
                },
            },
        }
    
    def test_template_creation_valid(self):
        """Test creating a valid template."""
        template = RecipeTemplate.from_dict(self.valid_template_data)
        
        self.assertEqual(template.name, "test-template")
        self.assertEqual(template.base, "Amidation.standard")
        self.assertEqual(len(template.action_ids), 3)
        self.assertEqual(len(template.variables), 2)
    
    def test_template_get_base_recipe(self):
        """Test retrieving the base recipe."""
        template = RecipeTemplate.from_dict(self.valid_template_data)
        recipe = template.get_base_recipe()
        
        self.assertIn("sessions", recipe)
        self.assertIsInstance(recipe["sessions"], list)
    
    def test_template_invalid_reaction_class(self):
        """Test that invalid reaction class raises error."""
        data = self.valid_template_data.copy()
        data["base"] = "NonexistentReaction.standard"
        
        with self.assertRaises(TemplateValidationError):
            RecipeTemplate.from_dict(data)
    
    def test_template_invalid_recipe_name(self):
        """Test that invalid recipe name raises error."""
        data = self.valid_template_data.copy()
        data["base"] = "Amidation.nonexistent_recipe"
        
        with self.assertRaises(TemplateValidationError):
            RecipeTemplate.from_dict(data)
    
    def test_template_invalid_action_reference(self):
        """Test that invalid action location raises error."""
        data = self.valid_template_data.copy()
        data["action_ids"]["invalid_action"] = {"session": 999, "actionnumber": 999}
        
        with self.assertRaises(ActionNotFoundError):
            RecipeTemplate.from_dict(data)
    
    def test_template_variable_references_unknown_action(self):
        """Test that variable referencing unknown action raises error."""
        data = self.valid_template_data.copy()
        data["variables"]["bad_var"] = {
            "action_id": "nonexistent_action",
            "path": "solvent"
        }
        
        with self.assertRaises(TemplateValidationError):
            RecipeTemplate.from_dict(data)
    
    def test_template_serialization(self):
        """Test template serialization round-trip."""
        template = RecipeTemplate.from_dict(self.valid_template_data)
        serialized = template.to_dict()
        restored = RecipeTemplate.from_dict(serialized)
        
        self.assertEqual(template.name, restored.name)
        self.assertEqual(template.base, restored.base)
        self.assertEqual(template.action_ids, restored.action_ids)


class RecipeGeneratorTestCase(DjangoTestCase):
    """Tests for RecipeGenerator class."""
    
    @classmethod
    def setUpTestData(cls):
        """Load amidation fixture into test DB."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)
    
    def setUp(self):
        """Set up test fixtures."""
        self.template = RecipeTemplate(
            name="amidation-test",
            base="Amidation.standard",
            action_ids={
                "add_acid": {"session": 1, "actionnumber": 1},
                "add_coupling_agent": {"session": 1, "actionnumber": 2},
                "add_base": {"session": 1, "actionnumber": 3},
            },
            session_ids={
                "reaction_1": 1,
                "reaction_2": 2,
                "stir_1": 3,
            },
            variables={
                "solvent": {
                    "action_id": "add_acid",
                    "path": "solvent"
                },
                "acid_equiv": {
                    "action_id": "add_acid",
                    "path": "equivalents"
                },
            },
        )
        self.generator = RecipeGenerator(self.template)
    
    def test_generate_single_recipe(self):
        """Test generating a single recipe with parameters."""
        result = self.generator.generate_single(
            recipe_name="test-recipe-1",
            solvent="DMF",
            acid_equiv=1.5
        )
        
        self.assertEqual(result["name"], "test-recipe-1")
        self.assertIn("recipe", result)
        
        # Verify the solvent was changed
        recipe = result["recipe"]
        action = self._find_action(recipe, session=1, action_number=1)
        self.assertEqual(action["solvent"], "DMF")
        self.assertEqual(action["equivalents"], 1.5)
    
    def test_generate_from_dataframe(self):
        """Test generating recipes from a DataFrame."""
        design_df = pd.DataFrame({
            "experiment_id": [1, 2, 3],
            "solvent": ["DMF", "DMA", "NMP"],
            "acid_equiv": [1.0, 1.2, 1.5],
        })
        
        recipes = self.generator.from_dataframe(design_df)
        
        self.assertEqual(len(recipes), 3)
        
        # Check first recipe
        recipe_1 = recipes[0]["recipe"]
        action_1 = self._find_action(recipe_1, session=1, action_number=1)
        self.assertEqual(action_1["solvent"], "DMF")
        
        # Check third recipe
        recipe_3 = recipes[2]["recipe"]
        action_3 = self._find_action(recipe_3, session=1, action_number=1)
        self.assertEqual(action_3["solvent"], "NMP")
    
    def test_generate_from_csv(self):
        """Test generating recipes from a CSV file."""
        # Create temporary CSV
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("experiment_id,solvent,acid_equiv\n")
            f.write("1,DMF,1.0\n")
            f.write("2,DMSO,2.0\n")
            csv_path = f.name
        
        try:
            recipes = self.generator.from_csv(csv_path)
            self.assertEqual(len(recipes), 2)
            
            # Verify values
            recipe_2 = recipes[1]["recipe"]
            action = self._find_action(recipe_2, session=1, action_number=1)
            self.assertEqual(action["solvent"], "DMSO")
            self.assertEqual(action["equivalents"], 2.0)
        finally:
            os.unlink(csv_path)
    
    def test_action_reordering(self):
        """Test reordering actions within a session."""
        result = self.generator.generate_single(
            recipe_name="reorder-test",
            reaction_1_action_order="add_coupling_agent,add_acid,add_base"
        )
        
        recipe = result["recipe"]
        session_1 = next(s for s in recipe["sessions"] if s["session_number"] == 1)
        
        # Get actions in order
        actions = session_1["actions"]
        action_nums = [a["action_number"] for a in actions]
        
        # Actions should be renumbered sequentially
        self.assertEqual(action_nums, [1, 2, 3])
    
    def test_get_variable_names(self):
        """Test retrieving variable names."""
        var_names = self.generator.get_variable_names()
        self.assertIn("solvent", var_names)
        self.assertIn("acid_equiv", var_names)
    
    def test_generate_design_matrix_template(self):
        """Test generating empty design matrix template."""
        df = self.generator.generate_design_matrix_template()
        
        self.assertIn("experiment_id", df.columns)
        self.assertIn("solvent", df.columns)
        self.assertIn("acid_equiv", df.columns)
    
    def _find_action(self, recipe: dict, session: int, action_number: int) -> dict:
        """Helper to find an action in a recipe (new format)."""
        for s in recipe.get("sessions", []):
            if s.get("session_number") == session:
                for action in s.get("actions", []):
                    if action.get("action_number") == action_number:
                        return action
        return None


class DesignMatrixParserTestCase(TestCase):
    """Tests for DesignMatrixParser class."""
    
    def test_parse_csv(self):
        """Test parsing a CSV file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("solvent,acid_equiv,temperature\n")
            f.write("DMF,1.0,25\n")
            f.write("DMA,1.5,40\n")
            csv_path = f.name
        
        try:
            parser = DesignMatrixParser()
            df = parser.parse_csv(csv_path, validate=False)
            
            self.assertEqual(len(df), 2)
            self.assertIn("solvent", df.columns)
        finally:
            os.unlink(csv_path)
    
    def test_type_conversion(self):
        """Test automatic type conversion of columns."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("name_equiv,reaction_1_action_order\n")
            f.write("1.5,add_acid,add_base\n")
            csv_path = f.name
        
        try:
            parser = DesignMatrixParser()
            df = parser.parse_csv(csv_path, validate=False)
            
            # Numeric columns should be converted
            self.assertEqual(df["name_equiv"].dtype, float)
        finally:
            os.unlink(csv_path)


class CreateTemplateFromRecipeTestCase(DjangoTestCase):
    """Tests for create_template_from_recipe helper function."""
    
    @classmethod
    def setUpTestData(cls):
        """Load amidation fixture into test DB."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)
    
    def test_create_template_basic(self):
        """Test creating a template from an existing recipe."""
        template = create_template_from_recipe(
            reaction_class="Amidation",
            recipe_name="standard",
            template_name="auto-amidation"
        )
        
        self.assertEqual(template.name, "auto-amidation")
        self.assertEqual(template.base, "Amidation.standard")
        
        # Should have auto-generated action_ids
        self.assertGreater(len(template.action_ids), 0)
        
        # Should have auto-generated variables
        self.assertGreater(len(template.variables), 0)
    
    def test_create_template_invalid_class(self):
        """Test that invalid reaction class raises error."""
        with self.assertRaises(TemplateValidationError):
            create_template_from_recipe(
                reaction_class="NonexistentReaction",
                recipe_name="standard"
            )
    
    def test_auto_generated_variables(self):
        """Test that auto-generated variables are correct."""
        template = create_template_from_recipe(
            reaction_class="Amidation",
            recipe_name="standard"
        )
        
        # Should have solvent and concentration variables for add actions
        var_names = list(template.variables.keys())
        solvent_vars = [v for v in var_names if "solvent" in v]
        
        self.assertGreater(len(solvent_vars), 0)


class IntegrationTestCase(DjangoTestCase):
    """Integration tests for full workflow."""
    
    @classmethod
    def setUpTestData(cls):
        """Load amidation fixture into test DB."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)
    
    def test_full_workflow(self):
        """Test complete workflow from template creation to recipe generation."""
        # 1. Create template from existing recipe
        template = create_template_from_recipe(
            reaction_class="Amidation",
            recipe_name="standard",
            template_name="integration-test"
        )
        
        # 2. Create generator
        generator = RecipeGenerator(template)
        
        # 3. Generate design matrix template
        design_template = generator.generate_design_matrix_template()
        self.assertIsInstance(design_template, pd.DataFrame)
        
        # 4. Create a simple design matrix
        var_names = generator.get_variable_names()
        if var_names:
            # Use first variable if any exist
            first_var = var_names[0]
            design_df = pd.DataFrame({
                "experiment_id": [1, 2],
                first_var: ["test_value_1", "test_value_2"],
            })
            
            # 5. Generate recipes
            recipes = generator.from_dataframe(design_df)
            
            self.assertEqual(len(recipes), 2)
            self.assertIn("recipe", recipes[0])
            self.assertIn("name", recipes[0])
