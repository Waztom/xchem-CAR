"""
Management command: ingest_recipes

Bulk-load all ~38 recipes from the legacy ``encoded_recipes`` dict
(backend/recipebuilder/encodedrecipes.py) into the normalised Recipe DB
models.

Usage::

    python3 manage.py ingest_recipes              # load all recipes
    python3 manage.py ingest_recipes --dry-run    # report what would be loaded
    python3 manage.py ingest_recipes --class "Amidation"  # one reaction class
    python3 manage.py ingest_recipes --wipe       # delete ALL existing recipes + reload

The converter handles every quirk of the old format:
  - "starting material" (with space) → role=startingmaterial
  - "workup1" / "workup2" → role=workup, index=1/2
  - "Intramolecular" (capital I) and "intramolecular" keys
  - continuation sessions
  - density field inside material dict
  - extract actions with bottomlayerquantity
"""
from __future__ import annotations

import logging

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from backend.models import (
    Recipe,
    RecipeActionSession,
    RecipeAddAction,
    RecipeExtractAction,
    RecipeMixAction,
    RecipeStirAction,
)
from backend.recipe_utils import parse_plate_type
from backend.recipebuilder.encodedrecipes import encoded_recipes

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "Ingest recipes from encodedrecipes.py into Recipe DB models"

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Print what would be created without writing to DB",
        )
        parser.add_argument(
            "--class",
            dest="reaction_class",
            type=str,
            default=None,
            help="Only ingest recipes for this reaction class",
        )
        parser.add_argument(
            "--wipe",
            action="store_true",
            help="Delete ALL existing Recipe rows before ingesting",
        )

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        target_class = options["reaction_class"]
        wipe = options["wipe"]

        if wipe and not dry_run:
            count, _ = Recipe.objects.all().delete()
            self.stdout.write(self.style.WARNING(f"Wiped {count} DB objects"))

        stats = {"recipes": 0, "sessions": 0, "actions": 0, "skipped": 0}

        for reaction_class, class_data in encoded_recipes.items():
            if target_class and reaction_class != target_class:
                continue

            intramolecular_possible = class_data.get("intramolecular", False)
            recipes_dict = class_data.get("recipes", {})

            for recipe_name, recipe_data in recipes_dict.items():
                # Skip if already exists
                if Recipe.objects.filter(
                    reaction_class=reaction_class,
                    name=recipe_name,
                    version=1,
                ).exists():
                    if not wipe:
                        stats["skipped"] += 1
                        if not dry_run:
                            self.stdout.write(
                                f"  SKIP {reaction_class}/{recipe_name} (already exists)"
                            )
                        continue

                if dry_run:
                    sessions = recipe_data.get("actionsessions", [])
                    n_actions = _count_legacy_actions(sessions, intramolecular_possible)
                    self.stdout.write(
                        f"  WOULD CREATE {reaction_class}/{recipe_name} "
                        f"({len(sessions)} sessions, ~{n_actions} actions)"
                    )
                    stats["recipes"] += 1
                    stats["sessions"] += len(sessions)
                    stats["actions"] += n_actions
                    continue

                try:
                    with transaction.atomic():
                        r, n_sess, n_act = _ingest_one_recipe(
                            reaction_class=reaction_class,
                            recipe_name=recipe_name,
                            recipe_data=recipe_data,
                            intramolecular_possible=intramolecular_possible,
                        )
                    stats["recipes"] += 1
                    stats["sessions"] += n_sess
                    stats["actions"] += n_act
                    self.stdout.write(
                        self.style.SUCCESS(
                            f"  ✓ {reaction_class}/{recipe_name} "
                            f"({n_sess} sessions, {n_act} actions)"
                        )
                    )
                except Exception as e:
                    self.stderr.write(
                        self.style.ERROR(
                            f"  ✗ {reaction_class}/{recipe_name}: {e}"
                        )
                    )
                    logger.exception("Failed to ingest %s/%s", reaction_class, recipe_name)

        verb = "Would create" if dry_run else "Created"
        self.stdout.write(
            self.style.SUCCESS(
                f"\n{verb} {stats['recipes']} recipes, "
                f"{stats['sessions']} sessions, "
                f"{stats['actions']} actions. "
                f"Skipped {stats['skipped']}."
            )
        )


# ── Core ingest logic ──────────────────────────────────────────────


def _ingest_one_recipe(
    reaction_class: str,
    recipe_name: str,
    recipe_data: dict,
    intramolecular_possible: bool,
) -> tuple[Recipe, int, int]:
    """Create a :class:`Recipe` and all child sessions/actions from the
    legacy ``encoded_recipes`` dict entry.

    Returns ``(recipe, n_sessions, n_actions)``.
    """
    recipe = Recipe.objects.create(
        reaction_class=reaction_class,
        name=recipe_name,
        version=1,
        intramolecular=intramolecular_possible,
        estimated_yield=recipe_data.get("yield"),
        reaction_smarts=recipe_data.get("reactionSMARTS", []),
        references=recipe_data.get("references"),
    )

    n_sessions = 0
    n_actions = 0

    for sess_data in recipe_data.get("actionsessions", []):
        sess_type = sess_data.get("type", "reaction")
        driver = sess_data.get("driver", "robot")
        continuation = sess_data.get("continuation", False)

        session = RecipeActionSession.objects.create(
            recipe=recipe,
            session_type=sess_type,
            driver=driver,
            continuation=bool(continuation),
        )
        n_sessions += 1

        if sess_type == "reaction":
            # Reaction sessions have inter/intramolecular sub-dicts
            for context_key in ("intermolecular", "intramolecular", "Intramolecular"):
                if context_key not in sess_data:
                    continue
                ctx_data = sess_data[context_key]
                actions_list = ctx_data.get("actions", [])
                mol_context = context_key.lower()  # normalise "Intramolecular" → "intramolecular"
                for action in actions_list:
                    _create_action(session, action, molecular_context=mol_context)
                    n_actions += 1
        else:
            # Non-reaction sessions: actions list at top level
            for action in sess_data.get("actions", []):
                _create_action(session, action, molecular_context=None)
                n_actions += 1

    return recipe, n_sessions, n_actions


def _create_action(
    session: RecipeActionSession,
    action: dict,
    molecular_context: str | None,
) -> None:
    """Dispatch to the correct model creator based on action ``type``."""
    atype = action.get("type")
    if atype == "add":
        _create_add_action(session, action, molecular_context)
    elif atype == "stir":
        _create_stir_action(session, action)
    elif atype == "extract":
        _create_extract_action(session, action)
    elif atype == "mix":
        _create_mix_action(session, action)
    else:
        logger.warning("Unknown action type '%s' in session %s", atype, session)


def _create_add_action(
    session: RecipeActionSession,
    action: dict,
    molecular_context: str | None,
) -> None:
    content = action.get("content", {})
    plates = content.get("plates", {})
    material = content.get("material", {})
    quantity = material.get("quantity", {})

    from_role, from_idx = parse_plate_type(plates.get("fromplatetype", "startingmaterial"))
    to_role, to_idx = parse_plate_type(plates.get("toplatetype", "reaction"))

    RecipeAddAction.objects.create(
        session=session,
        molecular_context=molecular_context,
        material_smarts=material.get("SMARTS"),
        material_smiles=material.get("SMILES"),
        equivalents=quantity.get("value", 0),
        quantity_unit=quantity.get("unit", "moleq"),
        solvent=material.get("solvent"),
        concentration=material.get("concentration"),
        density=material.get("density"),
        from_plate_role=from_role,
        from_plate_role_index=from_idx,
        to_plate_role=to_role,
        to_plate_role_index=to_idx,
    )


def _create_stir_action(session: RecipeActionSession, action: dict) -> None:
    content = action.get("content", {})
    temp = content.get("temperature", {})
    dur = content.get("duration", {})

    plate_type = content.get("platetype", "reaction")
    role, idx = parse_plate_type(plate_type)

    RecipeStirAction.objects.create(
        session=session,
        temperature=temp.get("value", 25),
        temperature_unit=temp.get("unit", "degC"),
        duration=dur.get("value", 1),
        duration_unit=dur.get("unit", "hours"),
        plate_role=role,
        plate_role_index=idx,
    )


def _create_extract_action(session: RecipeActionSession, action: dict) -> None:
    content = action.get("content", {})
    plates = content.get("plates", {})
    material = content.get("material", {})
    quantity = material.get("quantity", {})

    from_role, from_idx = parse_plate_type(plates.get("fromplatetype", "reaction"))
    to_role, to_idx = parse_plate_type(plates.get("toplatetype", "workup"))

    bottom_layer_vol = None
    if "bottomlayerquantity" in material:
        bottom_layer_vol = material["bottomlayerquantity"].get("value")

    RecipeExtractAction.objects.create(
        session=session,
        layer=material.get("layer", "bottom"),
        volume=quantity.get("value", 0),
        bottom_layer_volume=bottom_layer_vol,
        smiles=material.get("extraction_smiles"),
        solvent=material.get("solvent"),
        concentration=material.get("concentration"),
        from_plate_role=from_role,
        from_plate_role_index=from_idx,
        to_plate_role=to_role,
        to_plate_role_index=to_idx,
    )


def _create_mix_action(session: RecipeActionSession, action: dict) -> None:
    content = action.get("content", {})
    reps = content.get("repetitions", {})
    plate_type = content.get("platetype", "reaction")
    role, idx = parse_plate_type(plate_type)

    RecipeMixAction.objects.create(
        session=session,
        plate_role=role,
        plate_role_index=idx,
        repetitions=reps.get("value", 1) if isinstance(reps, dict) else reps,
    )


# ── Helpers ────────────────────────────────────────────────────────


def _count_legacy_actions(sessions: list, intramolecular_possible: bool) -> int:
    """Estimate action count from legacy session list (for dry-run)."""
    total = 0
    for sess in sessions:
        stype = sess.get("type", "")
        if stype == "reaction":
            for key in ("intermolecular", "intramolecular", "Intramolecular"):
                if key in sess:
                    total += len(sess[key].get("actions", []))
        else:
            total += len(sess.get("actions", []))
    return total
