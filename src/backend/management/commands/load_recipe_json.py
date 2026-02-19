"""
Management command: load_recipe_json

Load one or more recipes from JSON files into the Recipe DB.

Usage::

    python3 manage.py load_recipe_json recipes/amidation_standard.json
    python3 manage.py load_recipe_json recipes/*.json
    python3 manage.py load_recipe_json recipes/*.json --dry-run
    python3 manage.py load_recipe_json recipes/*.json --update  # bump version

The JSON format mirrors the DB schema directly.  All ``*_plate_index``
fields are **optional** and default to 1 — the chemist only needs to
specify them for multi-stage workup routing (e.g. workup plate 2).

See ``backend/tests/fixtures/amidation_standard.json`` for the
canonical reference.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

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

logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = "Load recipes from JSON files into the Recipe DB"

    def add_arguments(self, parser):
        parser.add_argument(
            "files",
            nargs="+",
            type=str,
            help="One or more JSON file paths to load",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Validate and report what would be created without writing",
        )
        parser.add_argument(
            "--update",
            action="store_true",
            help="If recipe already exists, create a new version instead of skipping",
        )

    def handle(self, *args, **options):
        dry_run = options["dry_run"]
        update = options["update"]
        files = options["files"]

        stats = {"recipes": 0, "sessions": 0, "actions": 0, "skipped": 0}

        for filepath in files:
            path = Path(filepath)
            if not path.exists():
                self.stderr.write(self.style.ERROR(f"File not found: {filepath}"))
                continue

            try:
                data = json.loads(path.read_text())
            except json.JSONDecodeError as e:
                self.stderr.write(
                    self.style.ERROR(f"Invalid JSON in {filepath}: {e}")
                )
                continue

            try:
                _validate_recipe_json(data, filepath)
            except ValueError as e:
                self.stderr.write(self.style.ERROR(f"{filepath}: {e}"))
                continue

            reaction_class = data["reaction_class"]
            recipe_name = data["name"]

            # Check for existing recipe
            existing = Recipe.objects.filter(
                reaction_class=reaction_class,
                name=recipe_name,
            ).order_by("-version").first()

            if existing and not update:
                stats["skipped"] += 1
                self.stdout.write(
                    f"  SKIP {reaction_class}/{recipe_name} "
                    f"(v{existing.version} exists, use --update to create new version)"
                )
                continue

            version = (existing.version + 1) if existing else 1
            n_sessions = len(data.get("action_sessions", []))
            n_actions = sum(
                len(s.get("actions", []))
                for s in data.get("action_sessions", [])
            )

            if dry_run:
                self.stdout.write(
                    f"  WOULD CREATE {reaction_class}/{recipe_name} v{version} "
                    f"({n_sessions} sessions, {n_actions} actions) from {filepath}"
                )
                stats["recipes"] += 1
                stats["sessions"] += n_sessions
                stats["actions"] += n_actions
                continue

            try:
                with transaction.atomic():
                    recipe, n_sess, n_act = _load_recipe(
                        data=data,
                        version=version,
                        parent=existing,
                    )
                stats["recipes"] += 1
                stats["sessions"] += n_sess
                stats["actions"] += n_act
                self.stdout.write(
                    self.style.SUCCESS(
                        f"  ✓ {reaction_class}/{recipe_name} v{version} "
                        f"({n_sess} sessions, {n_act} actions)"
                    )
                )
            except Exception as e:
                self.stderr.write(
                    self.style.ERROR(
                        f"  ✗ {reaction_class}/{recipe_name}: {e}"
                    )
                )
                logger.exception(
                    "Failed to load %s/%s from %s",
                    reaction_class, recipe_name, filepath,
                )

        verb = "Would create" if dry_run else "Created"
        self.stdout.write(
            self.style.SUCCESS(
                f"\n{verb} {stats['recipes']} recipes, "
                f"{stats['sessions']} sessions, "
                f"{stats['actions']} actions. "
                f"Skipped {stats['skipped']}."
            )
        )


# ── Validation ─────────────────────────────────────────────────────

VALID_ACTION_TYPES = {"add", "stir", "extract", "mix"}
VALID_SESSION_TYPES = {"reaction", "stir", "workup", "analyse"}
VALID_DRIVERS = {"robot", "human"}
VALID_PLATE_ROLES = {
    "reaction", "workup", "spefilter", "lcms", "xchem",
    "nmr", "startingmaterial", "solvent",
}
VALID_QUANTITY_UNITS = {"moleq", "masseq", "uL", "mL", "mg", "g", "M", "uM"}
VALID_MOL_CONTEXTS = {"intermolecular", "intramolecular", None}
VALID_LAYERS = {"top", "bottom"}


def _validate_recipe_json(data: dict, filepath: str) -> None:
    """Raise ``ValueError`` if the JSON structure is invalid."""
    for required in ("reaction_class", "name", "action_sessions"):
        if required not in data:
            raise ValueError(f"Missing required top-level key: '{required}'")

    if not isinstance(data["action_sessions"], list) or not data["action_sessions"]:
        raise ValueError("'action_sessions' must be a non-empty list")

    for i, sess in enumerate(data["action_sessions"]):
        prefix = f"action_sessions[{i}]"

        if "session_type" not in sess:
            raise ValueError(f"{prefix}: missing 'session_type'")
        if sess["session_type"] not in VALID_SESSION_TYPES:
            raise ValueError(
                f"{prefix}: invalid session_type '{sess['session_type']}'"
            )
        if sess.get("driver") and sess["driver"] not in VALID_DRIVERS:
            raise ValueError(f"{prefix}: invalid driver '{sess['driver']}'")

        actions = sess.get("actions", [])
        if not isinstance(actions, list):
            raise ValueError(f"{prefix}: 'actions' must be a list")

        for j, action in enumerate(actions):
            _validate_action(action, f"{prefix}.actions[{j}]")


def _validate_action(action: dict, prefix: str) -> None:
    """Validate a single action dict."""
    atype = action.get("type")
    if atype not in VALID_ACTION_TYPES:
        raise ValueError(f"{prefix}: invalid action type '{atype}'")

    if atype == "add":
        if "amount" not in action:
            raise ValueError(f"{prefix}: add action missing 'amount'")
        if action.get("quantity_unit") and action["quantity_unit"] not in VALID_QUANTITY_UNITS:
            raise ValueError(f"{prefix}: invalid quantity_unit")
        if action.get("molecular_context") not in VALID_MOL_CONTEXTS:
            raise ValueError(f"{prefix}: invalid molecular_context")
        for role_key in ("from_plate_role", "to_plate_role"):
            if action.get(role_key) and action[role_key] not in VALID_PLATE_ROLES:
                raise ValueError(f"{prefix}: invalid {role_key} '{action[role_key]}'")

    elif atype == "stir":
        if "duration" not in action:
            raise ValueError(f"{prefix}: stir action missing 'duration'")

    elif atype == "extract":
        if "volume" not in action:
            raise ValueError(f"{prefix}: extract action missing 'volume'")
        if action.get("layer") and action["layer"] not in VALID_LAYERS:
            raise ValueError(f"{prefix}: invalid layer '{action['layer']}'")

    elif atype == "mix":
        if "repetitions" not in action:
            raise ValueError(f"{prefix}: mix action missing 'repetitions'")


# ── Core load logic ────────────────────────────────────────────────


def _load_recipe(
    data: dict,
    version: int,
    parent: Recipe | None,
) -> tuple[Recipe, int, int]:
    """Create a Recipe and all child sessions/actions from a JSON dict.

    Returns ``(recipe, n_sessions, n_actions)``.
    """
    recipe = Recipe.objects.create(
        reaction_class=data["reaction_class"],
        name=data["name"],
        version=version,
        parent=parent,
        description=data.get("description", ""),
        intramolecular=data.get("intramolecular", False),
        estimated_yield=data.get("estimated_yield"),
        reaction_smarts=data.get("reaction_smarts", []),
        references=data.get("references"),
    )

    n_sessions = 0
    n_actions = 0

    for sess_data in data.get("action_sessions", []):
        session = RecipeActionSession.objects.create(
            recipe=recipe,
            session_type=sess_data["session_type"],
            driver=sess_data.get("driver", "robot"),
            continuation=sess_data.get("continuation", False),
        )
        n_sessions += 1

        for action_data in sess_data.get("actions", []):
            _create_action(session, action_data)
            n_actions += 1

    return recipe, n_sessions, n_actions


def _create_action(session: RecipeActionSession, action: dict) -> None:
    """Dispatch to the correct model creator based on action ``type``."""
    atype = action["type"]
    if atype == "add":
        _create_add_action(session, action)
    elif atype == "stir":
        _create_stir_action(session, action)
    elif atype == "extract":
        _create_extract_action(session, action)
    elif atype == "mix":
        _create_mix_action(session, action)


def _create_add_action(session: RecipeActionSession, action: dict) -> None:
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
        from_plate_role=action.get("from_plate_role", "startingmaterial"),
        from_plate_role_index=action.get("from_plate_role_index", action.get("from_plate_index", 1)),
        to_plate_role=action.get("to_plate_role", "reaction"),
        to_plate_role_index=action.get("to_plate_role_index", action.get("to_plate_index", 1)),
    )


def _create_stir_action(session: RecipeActionSession, action: dict) -> None:
    RecipeStirAction.objects.create(
        session=session,
        temperature=action.get("temperature", 25),
        temperature_unit=action.get("temperature_unit", "degC"),
        duration=action["duration"],
        duration_unit=action.get("duration_unit", "h"),
        stirring_speed=action.get("stirring_speed", "normal"),
        plate_role=action.get("plate_role", "reaction"),
        plate_role_index=action.get("plate_role_index", action.get("plate_index", 1)),
    )


def _create_extract_action(session: RecipeActionSession, action: dict) -> None:
    RecipeExtractAction.objects.create(
        session=session,
        layer=action.get("layer", "bottom"),
        volume=action["volume"],
        bottom_layer_volume=action.get("bottom_layer_volume"),
        smiles=action.get("smiles"),
        solvent=action.get("solvent"),
        concentration=action.get("concentration"),
        from_plate_role=action.get("from_plate_role", "reaction"),
        from_plate_role_index=action.get("from_plate_role_index", action.get("from_plate_index", 1)),
        to_plate_role=action.get("to_plate_role", "workup"),
        to_plate_role_index=action.get("to_plate_role_index", action.get("to_plate_index", 1)),
    )


def _create_mix_action(session: RecipeActionSession, action: dict) -> None:
    RecipeMixAction.objects.create(
        session=session,
        plate_role=action.get("plate_role", "reaction"),
        plate_role_index=action.get("plate_role_index", action.get("plate_index", 1)),
        repetitions=action["repetitions"],
    )
