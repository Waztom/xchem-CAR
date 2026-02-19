"""Tests for recipe_utils.py - utility functions for querying Recipe DB models.

This module tests:
- Pure functions (no DB required): parse_plate_type, unparse_plate_type
- DB-dependent functions: recipe queries, session actions, plate extraction

Uses the amidation_standard.json and mitsunobu_with_extraction.json fixtures
from test_recipe_models.py.
"""

import json
from pathlib import Path
from unittest import TestCase as UnitTestCase

from django.test import TestCase as DjangoTestCase

from backend.models import (
    Recipe,
    RecipeActionSession,
    RecipeAddAction,
    RecipeStirAction,
    RecipeExtractAction,
    RecipeMixAction,
)
from backend.recipe_utils import (
    parse_plate_type,
    unparse_plate_type,
    get_latest_recipe,
    recipe_exists,
    get_recipe_intramolecular,
    get_recipe_yield,
    get_recipe_smarts,
    get_recipe_stir_temperature,
    get_session_recipe_actions,
    collect_session_actions,
    action_to_dict,
    recipe_to_dict,
    get_session_destination_plates,
    get_session_source_plates,
)

# Re-use fixture ingest helper from test_recipe_models
from .test_recipe_models import ingest_recipe_from_json

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


# ═══════════════════════════════════════════════════════════════════
# PURE FUNCTION TESTS (no DB required)
# ═══════════════════════════════════════════════════════════════════


class ParsePlateTypeTestCase(UnitTestCase):
    """Tests for parse_plate_type() - converts legacy plate strings to (role, index)."""

    def test_startingmaterial_no_index(self):
        """'startingmaterial' → ('startingmaterial', 1)"""
        role, idx = parse_plate_type("startingmaterial")
        self.assertEqual(role, "startingmaterial")
        self.assertEqual(idx, 1)

    def test_starting_material_with_space(self):
        """'starting material' (legacy with space) → ('startingmaterial', 1)"""
        role, idx = parse_plate_type("starting material")
        self.assertEqual(role, "startingmaterial")
        self.assertEqual(idx, 1)

    def test_reaction_no_index(self):
        """'reaction' → ('reaction', 1)"""
        role, idx = parse_plate_type("reaction")
        self.assertEqual(role, "reaction")
        self.assertEqual(idx, 1)

    def test_workup_no_index(self):
        """'workup' (no number) → ('workup', 1)"""
        role, idx = parse_plate_type("workup")
        self.assertEqual(role, "workup")
        self.assertEqual(idx, 1)

    def test_workup1(self):
        """'workup1' → ('workup', 1)"""
        role, idx = parse_plate_type("workup1")
        self.assertEqual(role, "workup")
        self.assertEqual(idx, 1)

    def test_workup2(self):
        """'workup2' → ('workup', 2)"""
        role, idx = parse_plate_type("workup2")
        self.assertEqual(role, "workup")
        self.assertEqual(idx, 2)

    def test_workup3(self):
        """'workup3' → ('workup', 3)"""
        role, idx = parse_plate_type("workup3")
        self.assertEqual(role, "workup")
        self.assertEqual(idx, 3)

    def test_spefilter_no_index(self):
        """'spefilter' → ('spefilter', 1)"""
        role, idx = parse_plate_type("spefilter")
        self.assertEqual(role, "spefilter")
        self.assertEqual(idx, 1)

    def test_lcms_no_index(self):
        """'lcms' → ('lcms', 1)"""
        role, idx = parse_plate_type("lcms")
        self.assertEqual(role, "lcms")
        self.assertEqual(idx, 1)

    def test_xchem_no_index(self):
        """'xchem' → ('xchem', 1)"""
        role, idx = parse_plate_type("xchem")
        self.assertEqual(role, "xchem")
        self.assertEqual(idx, 1)

    def test_nmr_no_index(self):
        """'nmr' → ('nmr', 1)"""
        role, idx = parse_plate_type("nmr")
        self.assertEqual(role, "nmr")
        self.assertEqual(idx, 1)

    def test_solvent_no_index(self):
        """'solvent' → ('solvent', 1)"""
        role, idx = parse_plate_type("solvent")
        self.assertEqual(role, "solvent")
        self.assertEqual(idx, 1)

    def test_analyse_no_index(self):
        """'analyse' → ('analyse', 1)"""
        role, idx = parse_plate_type("analyse")
        self.assertEqual(role, "analyse")
        self.assertEqual(idx, 1)

    def test_case_insensitive_reaction(self):
        """'REACTION' → ('reaction', 1)"""
        role, idx = parse_plate_type("REACTION")
        self.assertEqual(role, "reaction")
        self.assertEqual(idx, 1)

    def test_case_insensitive_workup(self):
        """'Workup2' → ('workup', 2)"""
        role, idx = parse_plate_type("Workup2")
        self.assertEqual(role, "workup")
        self.assertEqual(idx, 2)

    def test_whitespace_trimmed(self):
        """'  reaction  ' → ('reaction', 1)"""
        role, idx = parse_plate_type("  reaction  ")
        self.assertEqual(role, "reaction")
        self.assertEqual(idx, 1)

    def test_unknown_type_passthrough(self):
        """Unknown type falls through with index 1."""
        role, idx = parse_plate_type("custom")
        self.assertEqual(role, "custom")
        self.assertEqual(idx, 1)


class UnparsePlateTypeTestCase(UnitTestCase):
    """Tests for unparse_plate_type() - converts (role, index) to legacy strings."""

    def test_workup_index_1(self):
        """('workup', 1) → 'workup' (number omitted)"""
        self.assertEqual(unparse_plate_type("workup", 1), "workup")

    def test_workup_index_2(self):
        """('workup', 2) → 'workup2'"""
        self.assertEqual(unparse_plate_type("workup", 2), "workup2")

    def test_workup_index_3(self):
        """('workup', 3) → 'workup3'"""
        self.assertEqual(unparse_plate_type("workup", 3), "workup3")

    def test_reaction_index_1(self):
        """('reaction', 1) → 'reaction' (number omitted)"""
        self.assertEqual(unparse_plate_type("reaction", 1), "reaction")

    def test_startingmaterial_index_1(self):
        """('startingmaterial', 1) → 'startingmaterial'"""
        self.assertEqual(unparse_plate_type("startingmaterial", 1), "startingmaterial")

    def test_spefilter_index_1(self):
        """('spefilter', 1) → 'spefilter'"""
        self.assertEqual(unparse_plate_type("spefilter", 1), "spefilter")

    def test_default_index(self):
        """Calling with only role defaults to index 1."""
        self.assertEqual(unparse_plate_type("reaction"), "reaction")


class ParseUnparseRoundTripTestCase(UnitTestCase):
    """Tests that parse_plate_type and unparse_plate_type are inverses."""

    def test_roundtrip_reaction(self):
        """parse → unparse → parse for 'reaction'"""
        role, idx = parse_plate_type("reaction")
        back = unparse_plate_type(role, idx)
        self.assertEqual(back, "reaction")

    def test_roundtrip_workup2(self):
        """parse → unparse → parse for 'workup2'"""
        role, idx = parse_plate_type("workup2")
        back = unparse_plate_type(role, idx)
        self.assertEqual(back, "workup2")

    def test_roundtrip_startingmaterial(self):
        """parse → unparse for 'startingmaterial'"""
        role, idx = parse_plate_type("startingmaterial")
        back = unparse_plate_type(role, idx)
        self.assertEqual(back, "startingmaterial")


# ═══════════════════════════════════════════════════════════════════
# DB-DEPENDENT TESTS (require Django TestCase)
# ═══════════════════════════════════════════════════════════════════


class RecipeLookupTestCase(DjangoTestCase):
    """Tests for recipe lookup functions using the amidation fixture."""

    @classmethod
    def setUpTestData(cls):
        """Load amidation_standard.json into test DB."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)

    def test_get_latest_recipe_found(self):
        """get_latest_recipe returns the recipe."""
        recipe = get_latest_recipe("Amidation", "standard")
        self.assertEqual(recipe.id, self.recipe.id)
        self.assertEqual(recipe.reaction_class, "Amidation")

    def test_get_latest_recipe_not_found(self):
        """get_latest_recipe raises DoesNotExist for unknown class."""
        with self.assertRaises(Recipe.DoesNotExist):
            get_latest_recipe("nonexistent-reaction", "standard")

    def test_get_latest_recipe_specific_version(self):
        """get_latest_recipe can fetch a specific version."""
        recipe = get_latest_recipe("Amidation", "standard", version=1)
        self.assertEqual(recipe.version, 1)

    def test_recipe_exists_true(self):
        """recipe_exists returns True for known reaction class."""
        self.assertTrue(recipe_exists("Amidation"))

    def test_recipe_exists_false(self):
        """recipe_exists returns False for unknown reaction class."""
        self.assertFalse(recipe_exists("nonexistent-reaction"))

    def test_get_recipe_intramolecular(self):
        """get_recipe_intramolecular returns the flag."""
        # Amidation is intramolecular=True per fixture
        result = get_recipe_intramolecular("Amidation", "standard")
        self.assertTrue(result)

    def test_get_recipe_yield(self):
        """get_recipe_yield returns estimated yield as fraction."""
        result = get_recipe_yield("Amidation", "standard")
        # Amidation has estimated_yield=85
        self.assertAlmostEqual(result, 0.85, places=2)

    def test_get_recipe_smarts(self):
        """get_recipe_smarts returns the SMARTS list."""
        smarts = get_recipe_smarts("Amidation", "standard")
        self.assertIsInstance(smarts, list)
        self.assertGreater(len(smarts), 0)
        # Should contain at least one SMARTS pattern
        self.assertIn(">>", smarts[0])

    def test_get_recipe_stir_temperature(self):
        """get_recipe_stir_temperature returns first stir temp."""
        temp = get_recipe_stir_temperature("Amidation", "standard")
        # Amidation session 3 has a stir at 25°C
        self.assertEqual(temp, 25)


class RecipeVersioningTestCase(DjangoTestCase):
    """Tests for recipe versioning with get_latest_recipe."""

    @classmethod
    def setUpTestData(cls):
        """Create multiple versions of a recipe."""
        # Version 1
        Recipe.objects.create(
            reaction_class="test-versioned",
            name="standard",
            version=1,
            estimated_yield=80,
            reaction_smarts=["[A:1]>>[B:1]"],
        )
        # Version 2
        Recipe.objects.create(
            reaction_class="test-versioned",
            name="standard",
            version=2,
            estimated_yield=90,
            reaction_smarts=["[A:1]>>[C:1]"],
        )

    def test_get_latest_returns_highest_version(self):
        """get_latest_recipe returns highest version when none specified."""
        recipe = get_latest_recipe("test-versioned", "standard")
        self.assertEqual(recipe.version, 2)
        self.assertEqual(recipe.estimated_yield, 90)

    def test_get_specific_version(self):
        """Can still retrieve older version explicitly."""
        recipe = get_latest_recipe("test-versioned", "standard", version=1)
        self.assertEqual(recipe.version, 1)
        self.assertEqual(recipe.estimated_yield, 80)


class SessionRecipeActionsTestCase(DjangoTestCase):
    """Tests for get_session_recipe_actions and collect_session_actions."""

    @classmethod
    def setUpTestData(cls):
        """Load amidation_standard.json."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)

    def test_get_reaction_session_actions(self):
        """Session 1 (reaction) has multiple add actions."""
        actions = get_session_recipe_actions(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
        )
        self.assertGreater(len(actions), 0)
        # All should be RecipeAddAction for session 1
        add_actions = [a for a in actions if isinstance(a, RecipeAddAction)]
        self.assertGreater(len(add_actions), 0)

    def test_get_reaction_session_with_context(self):
        """Filtering by molecular_context limits add actions."""
        all_actions = get_session_recipe_actions(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
        )
        inter_actions = get_session_recipe_actions(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
            molecular_context="intermolecular",
        )
        # Intermolecular should have fewer or equal actions
        self.assertLessEqual(len(inter_actions), len(all_actions))

    def test_get_stir_session_actions(self):
        """Session 3 (stir) has stir actions."""
        actions = get_session_recipe_actions(
            reaction_class="Amidation",
            name="standard",
            session_type="stir",
            session_number=3,
        )
        stir_actions = [a for a in actions if isinstance(a, RecipeStirAction)]
        self.assertGreater(len(stir_actions), 0)

    def test_get_nonexistent_session(self):
        """Nonexistent session returns empty list."""
        actions = get_session_recipe_actions(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=99,
        )
        self.assertEqual(actions, [])

    def test_actions_sorted_by_action_number(self):
        """Actions are returned sorted by action_number."""
        actions = get_session_recipe_actions(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
        )
        action_numbers = [a.action_number for a in actions]
        self.assertEqual(action_numbers, sorted(action_numbers))

    def test_collect_session_actions_directly(self):
        """collect_session_actions returns all actions from a session object."""
        session = RecipeActionSession.objects.filter(
            recipe=self.recipe,
            session_type="reaction",
            session_number=1,
        ).first()
        self.assertIsNotNone(session)
        
        actions = collect_session_actions(session)
        self.assertGreater(len(actions), 0)
        # Should include add actions from session 1
        self.assertTrue(any(isinstance(a, RecipeAddAction) for a in actions))

    def test_collect_session_actions_with_molecular_context(self):
        """collect_session_actions filters by molecular_context for reaction sessions."""
        session = RecipeActionSession.objects.filter(
            recipe=self.recipe,
            session_type="reaction",
            session_number=1,
        ).first()
        
        all_actions = collect_session_actions(session)
        inter_actions = collect_session_actions(session, molecular_context="intermolecular")
        
        # Filtered should have fewer or equal actions
        self.assertLessEqual(len(inter_actions), len(all_actions))

    def test_collect_session_actions_stir_session(self):
        """collect_session_actions works on stir sessions."""
        session = RecipeActionSession.objects.filter(
            recipe=self.recipe,
            session_type="stir",
            session_number=3,
        ).first()
        self.assertIsNotNone(session)
        
        actions = collect_session_actions(session)
        self.assertGreater(len(actions), 0)
        self.assertTrue(any(isinstance(a, RecipeStirAction) for a in actions))

    def test_collect_session_actions_sorted(self):
        """collect_session_actions returns actions sorted by action_number."""
        session = RecipeActionSession.objects.filter(
            recipe=self.recipe,
            session_type="reaction",
            session_number=1,
        ).first()
        
        actions = collect_session_actions(session)
        action_numbers = [a.action_number for a in actions]
        self.assertEqual(action_numbers, sorted(action_numbers))


class ActionToDictTestCase(DjangoTestCase):
    """Tests for action_to_dict serialization."""

    @classmethod
    def setUpTestData(cls):
        """Load amidation_standard.json."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)

    def test_add_action_to_dict(self):
        """RecipeAddAction serializes with correct keys."""
        add = RecipeAddAction.objects.filter(session__recipe=self.recipe).first()
        d = action_to_dict(add)
        self.assertEqual(d["type"], "add")
        self.assertIn("action_number", d)
        self.assertIn("from_plate_role", d)
        self.assertIn("to_plate_role", d)
        self.assertIn("equivalents", d)

    def test_stir_action_to_dict(self):
        """RecipeStirAction serializes with correct keys."""
        stir = RecipeStirAction.objects.filter(session__recipe=self.recipe).first()
        d = action_to_dict(stir)
        self.assertEqual(d["type"], "stir")
        self.assertIn("temperature", d)
        self.assertIn("duration", d)
        self.assertIn("plate_role", d)

    def test_mix_action_to_dict(self):
        """RecipeMixAction serializes with correct keys."""
        mix = RecipeMixAction.objects.filter(session__recipe=self.recipe).first()
        d = action_to_dict(mix)
        self.assertEqual(d["type"], "mix")
        self.assertIn("repetitions", d)
        self.assertIn("plate_role", d)


class RecipeToDictTestCase(DjangoTestCase):
    """Tests for recipe_to_dict full serialization."""

    @classmethod
    def setUpTestData(cls):
        """Load amidation_standard.json."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)

    def test_recipe_to_dict_structure(self):
        """recipe_to_dict returns expected top-level keys."""
        d = recipe_to_dict("Amidation", "standard")
        self.assertIn("estimated_yield", d)
        self.assertIn("reaction_smarts", d)
        self.assertIn("sessions", d)

    def test_recipe_to_dict_session_count(self):
        """Session count matches fixture."""
        d = recipe_to_dict("Amidation", "standard")
        fixture_sessions = len(self.amidation_data["action_sessions"])
        self.assertEqual(len(d["sessions"]), fixture_sessions)

    def test_recipe_to_dict_session_has_actions(self):
        """Each session has actions list."""
        d = recipe_to_dict("Amidation", "standard")
        for session in d["sessions"]:
            self.assertIn("actions", session)
            self.assertIsInstance(session["actions"], list)


class DestinationPlatesTestCase(DjangoTestCase):
    """Tests for get_session_destination_plates."""

    @classmethod
    def setUpTestData(cls):
        """Load amidation_standard.json."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)

    def test_reaction_session_destinations(self):
        """Reaction session adds go to reaction plate."""
        plates = get_session_destination_plates(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
        )
        self.assertIn(("reaction", 1), plates)

    def test_stir_session_destinations(self):
        """Stir session targets reaction plate."""
        plates = get_session_destination_plates(
            reaction_class="Amidation",
            name="standard",
            session_type="stir",
            session_number=3,
        )
        self.assertIn(("reaction", 1), plates)

    def test_plates_are_unique(self):
        """Returned plates are unique (no duplicates)."""
        plates = get_session_destination_plates(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
        )
        self.assertEqual(len(plates), len(set(plates)))

    def test_plates_are_sorted(self):
        """Returned plates are sorted."""
        plates = get_session_destination_plates(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
        )
        self.assertEqual(plates, sorted(plates))


class SourcePlatesTestCase(DjangoTestCase):
    """Tests for get_session_source_plates."""

    @classmethod
    def setUpTestData(cls):
        """Load amidation_standard.json."""
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        cls.amidation_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.amidation_data)

    def test_reaction_session_sources(self):
        """Reaction session pulls from startingmaterial/solvent."""
        plates = get_session_source_plates(
            reaction_class="Amidation",
            name="standard",
            session_type="reaction",
            session_number=1,
        )
        # Should include startingmaterial as source
        source_roles = [p[0] for p in plates]
        self.assertIn("startingmaterial", source_roles)

    def test_sources_only_from_add_extract(self):
        """Only add/extract actions contribute source plates."""
        plates = get_session_source_plates(
            reaction_class="Amidation",
            name="standard",
            session_type="stir",
            session_number=3,
        )
        # Stir session has no add/extract, so no source plates
        self.assertEqual(plates, [])


class MitsunobuExtractionTestCase(DjangoTestCase):
    """Tests with mitsunobu fixture for extract actions."""

    @classmethod
    def setUpTestData(cls):
        """Load mitsunobu_with_extraction.json."""
        fixture_path = FIXTURES_DIR / "mitsunobu_with_extraction.json"
        cls.mitsunobu_data = json.loads(fixture_path.read_text())
        cls.recipe = ingest_recipe_from_json(cls.mitsunobu_data)

    def test_extract_action_to_dict(self):
        """RecipeExtractAction serializes correctly."""
        ext = RecipeExtractAction.objects.filter(session__recipe=self.recipe).first()
        d = action_to_dict(ext)
        self.assertEqual(d["type"], "extract")
        self.assertIn("layer", d)
        self.assertIn("volume", d)
        self.assertIn("from_plate_role", d)
        self.assertIn("to_plate_role", d)

    def test_workup_session_has_extracts(self):
        """Workup session contains extract actions."""
        actions = get_session_recipe_actions(
            reaction_class="Mitsunobu",
            name="standard_with_extraction",
            session_type="workup",
            session_number=5,  # Workup session 5 has extracts per fixture
        )
        extract_actions = [a for a in actions if isinstance(a, RecipeExtractAction)]
        self.assertGreater(len(extract_actions), 0)

    def test_workup_destination_plates_multiple(self):
        """Workup with extract can target multiple plates."""
        plates = get_session_destination_plates(
            reaction_class="Mitsunobu",
            name="standard_with_extraction",
            session_type="workup",
            session_number=5,
        )
        # Should have workup plates
        workup_plates = [p for p in plates if p[0] == "workup"]
        self.assertGreater(len(workup_plates), 0)
