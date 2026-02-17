"""Tests for the normalised Recipe DB models.

Uses the Amidation / standard recipe as the reference because it
exercises every feature:
  - intramolecular=True  (both inter- and intramolecular add actions)
  - continuation session  (session 2)
  - stir sessions
  - workup session with SPE filter transfers
  - mix action
"""

import json
from pathlib import Path

from django.test import TestCase

from backend.models import (
    PlateRole,
    Recipe,
    RecipeActionSession,
    RecipeAddAction,
    RecipeStirAction,
    RecipeExtractAction,
    RecipeMixAction,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"


# ── helpers: JSON → DB ingest ──────────────────────────────────────

def ingest_recipe_from_json(data: dict) -> Recipe:
    """Create a Recipe and all child sessions/actions from a JSON dict.

    This is the reference ingest implementation used by the tests.
    The production serializer should produce identical DB state.
    """
    recipe = Recipe.objects.create(
        reaction_class=data["reaction_class"],
        name=data["name"],
        description=data.get("description", ""),
        intramolecular=data.get("intramolecular", False),
        estimated_yield=data.get("estimated_yield"),
        reaction_smarts=data.get("reaction_smarts", []),
        references=data.get("references"),
    )

    for sess_idx, sess_data in enumerate(data["action_sessions"], start=1):
        session = RecipeActionSession.objects.create(
            recipe=recipe,
            session_type=sess_data["session_type"],
            driver=sess_data.get("driver", "robot"),
            continuation=sess_data.get("continuation", False),
        )

        for act_idx, action in enumerate(sess_data.get("actions", []), start=1):
            action_type = action["type"]

            if action_type == "add":
                RecipeAddAction.objects.create(
                    session=session,
                    molecular_context=action.get("molecular_context"),
                    material_smarts=action.get("material_smarts"),
                    material_smiles=action.get("material_smiles"),
                    equivalents=action["amount"],
                    quantity_unit=action.get("quantity_unit", "moleq"),
                    solvent=action.get("solvent"),
                    concentration=action.get("concentration"),
                    density=action.get("density"),
                    from_plate_role=action.get(
                        "from_plate_role", "startingmaterial"
                    ),
                    from_plate_role_index=action.get("from_plate_role_index", 1),
                    to_plate_role=action.get("to_plate_role", "reaction"),
                    to_plate_role_index=action.get("to_plate_role_index", 1),
                )

            elif action_type == "stir":
                RecipeStirAction.objects.create(
                    session=session,
                    temperature=action.get("temperature", 25),
                    temperature_unit=action.get("temperature_unit", "degC"),
                    duration=action["duration"],
                    duration_unit=action.get("duration_unit", "h"),
                    stirring_speed=action.get("stirring_speed", "normal"),
                    plate_role=action.get("plate_role", "reaction"),
                    plate_role_index=action.get("plate_role_index", 1),
                )

            elif action_type == "extract":
                RecipeExtractAction.objects.create(
                    session=session,
                    layer=action.get("layer", "bottom"),
                    volume=action["volume"],
                    bottom_layer_volume=action.get("bottom_layer_volume"),
                    smiles=action.get("smiles"),
                    solvent=action.get("solvent"),
                    concentration=action.get("concentration"),
                    from_plate_role=action.get("from_plate_role", "reaction"),
                    from_plate_role_index=action.get("from_plate_role_index", 1),
                    to_plate_role=action.get("to_plate_role", "workup"),
                    to_plate_role_index=action.get("to_plate_role_index", 1),
                )

            elif action_type == "mix":
                RecipeMixAction.objects.create(
                    session=session,
                    plate_role=action.get("plate_role", "reaction"),
                    plate_role_index=action.get("plate_role_index", 1),
                    repetitions=action["repetitions"],
                )

    return recipe


# ── helpers: DB → JSON export ──────────────────────────────────────

def _omit_defaults(d: dict, defaults: dict) -> dict:
    """Return *d* with entries removed when their value equals the default.

    Also strips any key whose value is ``None``.
    """
    return {
        k: v
        for k, v in d.items()
        if v is not None and v != defaults.get(k, _SENTINEL)
    }


_SENTINEL = object()  # marker so we keep keys with no listed default

# Defaults for each action / session / recipe level.
_ADD_DEFAULTS = {
    "quantity_unit": "moleq",
    "from_plate_role_index": 1,
    "to_plate_role_index": 1,
}
_STIR_DEFAULTS = {
    "temperature": 25,
    "temperature_unit": "degC",
    "duration_unit": "h",
    "stirring_speed": "normal",
    "plate_role_index": 1,
}
_EXTRACT_DEFAULTS = {
    "from_plate_role_index": 1,
    "to_plate_role_index": 1,
}
_MIX_DEFAULTS = {
    "plate_role_index": 1,
}
_SESSION_DEFAULTS = {
    "continuation": False,
}
_RECIPE_DEFAULTS = {
    "description": "",
    "references": None,
}


def export_recipe_to_json(recipe: Recipe) -> dict:
    """Export a Recipe and all child sessions/actions to a JSON-serialisable dict.

    Produces the *lean* canonical format: keys whose values equal the
    documented default are omitted so chemists only see what matters.
    This is the reference export implementation used by the tests.
    """
    sessions_out = []

    for session in recipe.action_sessions.all().order_by("session_number"):
        # Collect every action with its DB action_number for ordering,
        # then strip the number — JSON position is the canonical order.
        tagged: list[tuple[int, str, dict]] = []
        #   (action_number, context_sort_key, action_dict)

        for add in session.add_actions.all().order_by("action_number"):
            action_dict: dict = {
                "type": "add",
            }
            if add.molecular_context is not None:
                action_dict["molecular_context"] = add.molecular_context
            if add.material_smarts is not None:
                action_dict["material_smarts"] = add.material_smarts
            if add.material_smiles is not None:
                action_dict["material_smiles"] = add.material_smiles
            action_dict["amount"] = add.equivalents
            action_dict["quantity_unit"] = add.quantity_unit
            if add.solvent is not None:
                action_dict["solvent"] = add.solvent
            if add.concentration is not None:
                action_dict["concentration"] = add.concentration
            if add.density is not None:
                action_dict["density"] = add.density
            action_dict["from_plate_role"] = add.from_plate_role
            if add.from_plate_role_index != 1:
                action_dict["from_plate_role_index"] = add.from_plate_role_index
            action_dict["to_plate_role"] = add.to_plate_role
            if add.to_plate_role_index != 1:
                action_dict["to_plate_role_index"] = add.to_plate_role_index
            tagged.append((
                add.action_number,
                add.molecular_context or "",  # NULL sorts before "inter…"
                action_dict,
            ))

        for stir in session.stir_actions.all().order_by("action_number"):
            stir_dict: dict = {"type": "stir"}
            # Always include temperature for chemical clarity
            # Normalize to int if it's a whole number for cleaner JSON
            stir_dict["temperature"] = (
                int(stir.temperature)
                if stir.temperature == int(stir.temperature)
                else stir.temperature
            )
            if stir.temperature_unit != "degC":
                stir_dict["temperature_unit"] = stir.temperature_unit
            stir_dict["duration"] = stir.duration
            stir_dict["duration_unit"] = stir.duration_unit
            if stir.stirring_speed != "normal":
                stir_dict["stirring_speed"] = stir.stirring_speed
            stir_dict["plate_role"] = stir.plate_role
            if stir.plate_role_index != 1:
                stir_dict["plate_role_index"] = stir.plate_role_index
            tagged.append((
                stir.action_number,
                "",  # no molecular_context — sorts first
                stir_dict,
            ))

        for extract in session.extract_actions.all().order_by("action_number"):
            extract_dict: dict = {"type": "extract"}
            if extract.layer is not None:
                extract_dict["layer"] = extract.layer
            extract_dict["volume"] = extract.volume
            extract_dict["volume_unit"] = "uL"  # Always uL for extracts
            if extract.bottom_layer_volume is not None:
                extract_dict["bottom_layer_volume"] = extract.bottom_layer_volume
            if extract.smiles is not None:
                extract_dict["smiles"] = extract.smiles
            if extract.solvent is not None:
                extract_dict["solvent"] = extract.solvent
            if extract.concentration is not None:
                extract_dict["concentration"] = extract.concentration
            extract_dict["from_plate_role"] = extract.from_plate_role
            if extract.from_plate_role_index != 1:
                extract_dict["from_plate_role_index"] = extract.from_plate_role_index
            extract_dict["to_plate_role"] = extract.to_plate_role
            if extract.to_plate_role_index != 1:
                extract_dict["to_plate_role_index"] = extract.to_plate_role_index
            tagged.append((
                extract.action_number,
                "",
                extract_dict,
            ))

        for mix in session.mix_actions.all().order_by("action_number"):
            mix_dict: dict = {"type": "mix"}
            mix_dict["repetitions"] = mix.repetitions
            mix_dict["plate_role"] = mix.plate_role
            if mix.plate_role_index != 1:
                mix_dict["plate_role_index"] = mix.plate_role_index
            tagged.append((
                mix.action_number,
                "",
                mix_dict,
            ))

        # Sort by (action_number, molecular_context) so inter always
        # precedes intra at the same position, and mixed types interleave.
        tagged.sort(key=lambda t: (t[0], t[1]))
        actions_out = [entry[2] for entry in tagged]

        sess_dict: dict = {
            "session_type": session.session_type,
            "driver": session.driver,
        }
        if session.continuation:
            sess_dict["continuation"] = True
        sess_dict["actions"] = actions_out
        sessions_out.append(sess_dict)

    result: dict = {
        "reaction_class": recipe.reaction_class,
        "name": recipe.name,
    }
    if recipe.description:
        result["description"] = recipe.description
    result["intramolecular"] = recipe.intramolecular
    if recipe.estimated_yield is not None:
        result["estimated_yield"] = recipe.estimated_yield
    if recipe.reaction_smarts:
        result["reaction_smarts"] = recipe.reaction_smarts
    if recipe.references is not None:
        result["references"] = recipe.references
    result["action_sessions"] = sessions_out
    return result


class AmidationRecipeTestCase(TestCase):
    """Build the full Amidation / standard recipe and verify the ORM
    relationships, defaults, and constraints."""

    @classmethod
    def setUpTestData(cls):
        # ── Recipe ──────────────────────────────────────────────
        cls.recipe = Recipe.objects.create(
            reaction_class="Amidation",
            name="standard",
            version=1,
            intramolecular=True,
            estimated_yield=85,
            reaction_smarts=[
                "[#6:1](=[#8:2])-[#8;H1].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]"
            ],
        )

        # ── Session 1 — reaction (robot) ────────────────────────
        cls.s1 = RecipeActionSession.objects.create(
            recipe=cls.recipe,
            session_type="reaction",
            driver="robot",
        )

        # Intermolecular actions (3 adds)
        RecipeAddAction.objects.create(
            session=cls.s1,
            molecular_context="intermolecular",
            material_smarts="[#6:1](=[#8:2])-[#8;H1]",
            equivalents=1.0,
            solvent="DMA",
            concentration=0.5,
        )
        RecipeAddAction.objects.create(
            session=cls.s1,
            molecular_context="intermolecular",
            material_smiles="CCCP1(=O)OP(=O)(OP(=O)(O1)CCC)CCC",
            equivalents=1.2,
            solvent="DMA",
            concentration=0.5,
        )
        RecipeAddAction.objects.create(
            session=cls.s1,
            molecular_context="intermolecular",
            material_smiles="CCN(C(C)C)C(C)C",
            equivalents=3.5,
            density=0.74,
        )

        # Intramolecular actions (3 adds — different amounts / solvent)
        RecipeAddAction.objects.create(
            session=cls.s1,
            molecular_context="intramolecular",
            material_smarts="[#6:1](=[#8:2])-[#8;H1]",
            equivalents=1.0,
            solvent="DMA",
            concentration=0.5,
        )
        RecipeAddAction.objects.create(
            session=cls.s1,
            molecular_context="intramolecular",
            material_smiles="CCCP1(=O)OP(=O)(OP(=O)(O1)CCC)CCC",
            equivalents=1.2,
            solvent="DMA",
            concentration=0.5,
        )
        RecipeAddAction.objects.create(
            session=cls.s1,
            molecular_context="intramolecular",
            material_smiles="CCN(C(C)C)C(C)C",
            equivalents=3.5,
            solvent="DMA",
            concentration=10,
        )

        # ── Session 2 — reaction continuation (robot) ───────────
        cls.s2 = RecipeActionSession.objects.create(
            recipe=cls.recipe,
            session_type="reaction",
            driver="robot",
            continuation=True,
        )

        RecipeAddAction.objects.create(
            session=cls.s2,
            molecular_context="intermolecular",
            material_smarts="[#7;H3,H2,H1]",
            equivalents=1.1,
            solvent="DMA",
            concentration=0.5,
        )

        # ── Session 3 — stir (human) ───────────────────────────
        cls.s3 = RecipeActionSession.objects.create(
            recipe=cls.recipe,
            session_type="stir",
            driver="human",
        )

        RecipeStirAction.objects.create(
            session=cls.s3,
            temperature=25,
            temperature_unit="degC",
            duration=12,
            duration_unit="h",
        )

        # ── Session 4 — workup: add ACN wash (robot) ───────────
        cls.s4 = RecipeActionSession.objects.create(
            recipe=cls.recipe,
            session_type="workup",
            driver="robot",
        )

        RecipeAddAction.objects.create(
            session=cls.s4,
            material_smiles="CC#N",
            equivalents=200,
            quantity_unit="uL",
            solvent="ACN",
            from_plate_role=PlateRole.SOLVENT,
            to_plate_role=PlateRole.REACTION,
        )

        # ── Session 5 — stir (human) ───────────────────────────
        cls.s5 = RecipeActionSession.objects.create(
            recipe=cls.recipe,
            session_type="stir",
            driver="human",
        )

        RecipeStirAction.objects.create(
            session=cls.s5,
            temperature=25,
            duration=1,
        )

        # ── Session 6 — workup: SPE filter (robot) ─────────────
        cls.s6 = RecipeActionSession.objects.create(
            recipe=cls.recipe,
            session_type="workup",
            driver="robot",
        )

        # Move reaction mix → SPE filter
        RecipeAddAction.objects.create(
            session=cls.s6,
            equivalents=200,
            quantity_unit="uL",
            from_plate_role=PlateRole.REACTION,
            to_plate_role=PlateRole.SPEFILTER,
        )
        # ACN wash → reaction plate
        RecipeAddAction.objects.create(
            session=cls.s6,
            material_smiles="CC#N",
            equivalents=200,
            quantity_unit="uL",
            solvent="ACN",
            from_plate_role=PlateRole.SOLVENT,
            to_plate_role=PlateRole.REACTION,
        )
        # Mix reaction plate
        RecipeMixAction.objects.create(
            session=cls.s6,
            repetitions=3,
        )
        # Move reaction mix → SPE filter (second pass)
        RecipeAddAction.objects.create(
            session=cls.s6,
            equivalents=300,
            quantity_unit="uL",
            from_plate_role=PlateRole.REACTION,
            to_plate_role=PlateRole.SPEFILTER,
        )

    # ── Recipe-level tests ──────────────────────────────────────

    def test_recipe_str(self):
        self.assertEqual(str(self.recipe), "Amidation / standard (v1)")

    def test_recipe_unique_constraint(self):
        """Duplicate (reaction_class, name, version) should fail."""
        from django.db import IntegrityError

        with self.assertRaises(IntegrityError):
            Recipe.objects.create(
                reaction_class="Amidation",
                name="standard",
                version=1,
            )

    def test_recipe_intramolecular_flag(self):
        self.assertTrue(self.recipe.intramolecular)

    def test_recipe_smarts_stored(self):
        self.assertEqual(len(self.recipe.reaction_smarts), 1)
        self.assertIn(">>", self.recipe.reaction_smarts[0])

    def test_recipe_version_default(self):
        r2 = Recipe.objects.create(
            reaction_class="Amidation", name="high-temp"
        )
        self.assertEqual(r2.version, 1)

    def test_recipe_parent_lineage(self):
        child = Recipe.objects.create(
            reaction_class="Amidation",
            name="standard",
            version=2,
            parent=self.recipe,
        )
        self.assertEqual(child.parent, self.recipe)
        self.assertIn(child, self.recipe.derived_versions.all())

    # ── Session-level tests ─────────────────────────────────────

    def test_session_count(self):
        self.assertEqual(self.recipe.action_sessions.count(), 6)

    def test_session_ordering(self):
        sessions = list(self.recipe.action_sessions.all())
        numbers = [s.session_number for s in sessions]
        self.assertEqual(numbers, [1, 2, 3, 4, 5, 6])

    def test_session_types(self):
        types = list(
            self.recipe.action_sessions.values_list("session_type", flat=True)
        )
        self.assertEqual(
            types,
            ["reaction", "reaction", "stir", "workup", "stir", "workup"],
        )

    def test_continuation_session(self):
        self.assertTrue(self.s2.continuation)
        self.assertFalse(self.s1.continuation)

    def test_session_drivers(self):
        self.assertEqual(self.s1.driver, "robot")
        self.assertEqual(self.s3.driver, "human")

    def test_session_auto_numbering(self):
        """session_number is auto-assigned from creation order."""
        self.assertEqual(self.s1.session_number, 1)
        self.assertEqual(self.s2.session_number, 2)
        self.assertEqual(self.s6.session_number, 6)

    def test_action_auto_numbering_session1(self):
        """Session 1 has 6 add actions; each gets a sequential number."""
        numbers = list(
            self.s1.add_actions.values_list("action_number", flat=True)
        )
        self.assertEqual(sorted(numbers), [1, 2, 3, 4, 5, 6])

    def test_action_auto_numbering_cross_type(self):
        """Session 6 interleaves add/mix, numbering spans all types."""
        adds = list(
            self.s6.add_actions.order_by("action_number")
            .values_list("action_number", flat=True)
        )
        mixes = list(
            self.s6.mix_actions.values_list("action_number", flat=True)
        )
        self.assertEqual(adds, [1, 2, 4])
        self.assertEqual(mixes, [3])

    # ── Add-action tests ────────────────────────────────────────

    def test_intermolecular_add_count_session1(self):
        inter = self.s1.add_actions.filter(molecular_context="intermolecular")
        self.assertEqual(inter.count(), 3)

    def test_intramolecular_add_count_session1(self):
        intra = self.s1.add_actions.filter(molecular_context="intramolecular")
        self.assertEqual(intra.count(), 3)

    def test_workup_add_has_null_context(self):
        """Workup add actions have NULL molecular_context (applies to all pathways)."""
        workup_add = self.s4.add_actions.first()
        self.assertIsNone(workup_add.molecular_context)

    def test_add_action_smarts_vs_smiles(self):
        """First action uses SMARTS (reactant match), second uses SMILES (reagent)."""
        actions = self.s1.add_actions.filter(
            molecular_context="intermolecular"
        ).order_by("action_number")
        self.assertIsNotNone(actions[0].material_smarts)
        self.assertIsNone(actions[0].material_smiles)
        self.assertIsNone(actions[1].material_smarts)
        self.assertIsNotNone(actions[1].material_smiles)

    def test_add_action_density_field(self):
        """DIPEA (intermolecular) is neat liquid — has density, no solvent."""
        dipea = self.s1.add_actions.filter(
            molecular_context="intermolecular",
            material_smiles="CCN(C(C)C)C(C)C",
        ).first()
        self.assertAlmostEqual(dipea.density, 0.74)
        self.assertIsNone(dipea.solvent)

    def test_intramolecular_dipea_uses_solvent(self):
        """Intramolecular DIPEA is dissolved in DMA instead of neat."""
        dipea = self.s1.add_actions.filter(
            molecular_context="intramolecular",
            material_smiles="CCN(C(C)C)C(C)C",
        ).first()
        self.assertIsNone(dipea.density)
        self.assertEqual(dipea.solvent, "DMA")
        self.assertAlmostEqual(dipea.concentration, 10)

    def test_continuation_session_adds(self):
        """Session 2 only has one intermolecular add (the amine)."""
        adds = self.s2.add_actions.all()
        self.assertEqual(adds.count(), 1)
        self.assertEqual(adds[0].material_smarts, "[#7;H3,H2,H1]")

    # ── Plate role tests ────────────────────────────────────────

    def test_default_plates_reaction_add(self):
        """Reaction adds default to SM → reaction."""
        add = self.s1.add_actions.filter(
            molecular_context="intermolecular",
            material_smarts="[#6:1](=[#8:2])-[#8;H1]",
        ).first()
        self.assertEqual(add.from_plate_role, PlateRole.STARTINGMATERIAL)
        self.assertEqual(add.to_plate_role, PlateRole.REACTION)
        self.assertEqual(add.from_plate_role_index, 1)
        self.assertEqual(add.to_plate_role_index, 1)

    def test_workup_wash_plates(self):
        """Workup wash: solvent → reaction."""
        wash = self.s4.add_actions.first()
        self.assertEqual(wash.from_plate_role, PlateRole.SOLVENT)
        self.assertEqual(wash.to_plate_role, PlateRole.REACTION)

    def test_spe_filter_transfer(self):
        """SPE filtration: reaction → spefilter."""
        spe_adds = self.s6.add_actions.filter(
            to_plate_role=PlateRole.SPEFILTER
        ).order_by("action_number")
        self.assertEqual(spe_adds.count(), 2)
        for add in spe_adds:
            self.assertEqual(add.from_plate_role, PlateRole.REACTION)
            self.assertEqual(add.to_plate_role, PlateRole.SPEFILTER)

    def test_plate_role_index_defaults(self):
        """All plate indices in this recipe should be 1 (no multi-plate workups)."""
        for add in RecipeAddAction.objects.filter(session__recipe=self.recipe):
            self.assertEqual(add.from_plate_role_index, 1)
            self.assertEqual(add.to_plate_role_index, 1)

    # ── Stir-action tests ──────────────────────────────────────

    def test_stir_session3(self):
        stir = self.s3.stir_actions.first()
        self.assertAlmostEqual(stir.temperature, 25)
        self.assertEqual(stir.temperature_unit, "degC")
        self.assertAlmostEqual(stir.duration, 12)
        self.assertEqual(stir.duration_unit, "h")

    def test_stir_session5(self):
        stir = self.s5.stir_actions.first()
        self.assertAlmostEqual(stir.duration, 1)

    def test_stir_plate_defaults(self):
        stir = self.s3.stir_actions.first()
        self.assertEqual(stir.plate_role, PlateRole.REACTION)
        self.assertEqual(stir.plate_role_index, 1)

    # ── Mix-action tests ────────────────────────────────────────

    def test_mix_action(self):
        mix = self.s6.mix_actions.first()
        self.assertEqual(mix.repetitions, 3)
        self.assertEqual(mix.plate_role, PlateRole.REACTION)

    # ── Cascade delete tests ────────────────────────────────────

    def test_delete_recipe_cascades(self):
        """Deleting a recipe should cascade-delete all sessions and actions."""
        recipe_id = self.recipe.pk
        # Count before
        sessions = RecipeActionSession.objects.filter(recipe_id=recipe_id).count()
        self.assertGreater(sessions, 0)

        # Delete
        Recipe.objects.filter(pk=recipe_id).delete()

        self.assertEqual(
            RecipeActionSession.objects.filter(recipe_id=recipe_id).count(), 0
        )
        self.assertEqual(
            RecipeAddAction.objects.filter(session__recipe_id=recipe_id).count(), 0
        )
        self.assertEqual(
            RecipeStirAction.objects.filter(session__recipe_id=recipe_id).count(), 0
        )
        self.assertEqual(
            RecipeMixAction.objects.filter(session__recipe_id=recipe_id).count(), 0
        )


class RecipeVersioningTestCase(TestCase):
    """Test the version auto-increment pattern and parent lineage."""

    def test_version_increment_pattern(self):
        """Simulates the upload-creates-new-version workflow."""
        v1 = Recipe.objects.create(
            reaction_class="Suzuki coupling",
            name="standard",
            version=1,
        )
        # Simulate uploading a new version of same recipe
        latest = (
            Recipe.objects.filter(
                reaction_class="Suzuki coupling", name="standard"
            )
            .order_by("-version")
            .first()
        )
        v2 = Recipe.objects.create(
            reaction_class="Suzuki coupling",
            name="standard",
            version=latest.version + 1,
            parent=latest,
        )
        self.assertEqual(v2.version, 2)
        self.assertEqual(v2.parent, v1)

    def test_multiple_recipe_names_same_class(self):
        """Different recipe names under the same reaction class are independent."""
        Recipe.objects.create(
            reaction_class="Amidation", name="standard", version=1
        )
        Recipe.objects.create(
            reaction_class="Amidation", name="high-temp", version=1
        )
        self.assertEqual(
            Recipe.objects.filter(reaction_class="Amidation").count(), 2
        )


class RecipeJSONRoundTripTestCase(TestCase):
    """Load amidation_standard.json → DB → export → compare."""

    @classmethod
    def setUpTestData(cls):
        fixture_path = FIXTURES_DIR / "amidation_standard.json"
        with open(fixture_path) as f:
            cls.source_json = json.load(f)
        cls.recipe = ingest_recipe_from_json(cls.source_json)

    # ── ingest verification ─────────────────────────────────────

    def test_recipe_created(self):
        self.assertEqual(self.recipe.reaction_class, "Amidation")
        self.assertEqual(self.recipe.name, "standard")
        self.assertTrue(self.recipe.intramolecular)

    def test_session_count_from_json(self):
        self.assertEqual(self.recipe.action_sessions.count(), 6)

    def test_total_add_actions(self):
        count = RecipeAddAction.objects.filter(
            session__recipe=self.recipe
        ).count()
        expected = sum(
            1
            for s in self.source_json["action_sessions"]
            for a in s["actions"]
            if a["type"] == "add"
        )
        self.assertEqual(count, expected)

    def test_total_stir_actions(self):
        count = RecipeStirAction.objects.filter(
            session__recipe=self.recipe
        ).count()
        expected = sum(
            1
            for s in self.source_json["action_sessions"]
            for a in s["actions"]
            if a["type"] == "stir"
        )
        self.assertEqual(count, expected)

    def test_total_mix_actions(self):
        count = RecipeMixAction.objects.filter(
            session__recipe=self.recipe
        ).count()
        expected = sum(
            1
            for s in self.source_json["action_sessions"]
            for a in s["actions"]
            if a["type"] == "mix"
        )
        self.assertEqual(count, expected)

    def test_total_extract_actions(self):
        """Amidation standard has no extract actions."""
        count = RecipeExtractAction.objects.filter(
            session__recipe=self.recipe
        ).count()
        self.assertEqual(count, 0)

    # ── round-trip: DB → JSON → compare ────────────────────────

    def test_round_trip_recipe_level(self):
        exported = export_recipe_to_json(self.recipe)
        for key in [
            "reaction_class",
            "name",
            "intramolecular",
            "estimated_yield",
            "reaction_smarts",
        ]:
            self.assertEqual(
                exported.get(key),
                self.source_json.get(key),
                f"Mismatch on recipe key '{key}'",
            )
        # Optional keys — only compare when present in the source fixture
        for key in ["description", "references"]:
            if key in self.source_json:
                self.assertEqual(
                    exported.get(key),
                    self.source_json[key],
                    f"Mismatch on recipe key '{key}'",
                )

    def test_round_trip_session_count(self):
        exported = export_recipe_to_json(self.recipe)
        self.assertEqual(
            len(exported["action_sessions"]),
            len(self.source_json["action_sessions"]),
        )

    def test_round_trip_session_metadata(self):
        exported = export_recipe_to_json(self.recipe)
        for idx, (exp_sess, src_sess) in enumerate(
            zip(exported["action_sessions"], self.source_json["action_sessions"]),
            start=1,
        ):
            for key in ["session_type", "driver"]:
                self.assertEqual(
                    exp_sess[key],
                    src_sess[key],
                    f"Session {idx} mismatch on '{key}'",
                )
            # continuation is omitted when False (the default)
            self.assertEqual(
                exp_sess.get("continuation", False),
                src_sess.get("continuation", False),
                f"Session {idx} mismatch on 'continuation'",
            )

    def test_round_trip_action_counts_per_session(self):
        exported = export_recipe_to_json(self.recipe)
        for idx, (exp_sess, src_sess) in enumerate(
            zip(exported["action_sessions"], self.source_json["action_sessions"]),
            start=1,
        ):
            self.assertEqual(
                len(exp_sess["actions"]),
                len(src_sess["actions"]),
                f"Action count mismatch in session {idx}",
            )

    def test_round_trip_full_json_equality(self):
        """The exported JSON should be identical to the source fixture."""
        exported = export_recipe_to_json(self.recipe)

        # Normalise both to sorted JSON strings for comparison
        exported_str = json.dumps(exported, sort_keys=True, indent=2)
        source_str = json.dumps(self.source_json, sort_keys=True, indent=2)

        self.assertEqual(
            exported_str,
            source_str,
            "Full round-trip JSON mismatch — see diff above.",
        )
