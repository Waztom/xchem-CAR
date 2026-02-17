"""
Utility functions for querying Recipe DB models.

Replaces direct lookups into ``encoded_recipes`` dict.  Every consumer
that previously read from the dict should use these helpers instead.
"""
from __future__ import annotations

import logging
import re
from typing import Optional

from django.db.models import Q

from .models import (
    PlateRole,
    Recipe,
    RecipeActionSession,
    RecipeAddAction,
    RecipeExtractAction,
    RecipeMixAction,
    RecipeStirAction,
)

logger = logging.getLogger(__name__)

# ── plate-type mapping ─────────────────────────────────────────────
# Old encoded recipes use flat strings like "workup1", "workup2",
# "starting material".  Map them to (PlateRole, index) pairs.

_PLATE_TYPE_RE = re.compile(
    r"^(reaction|workup|spefilter|lcms|xchem|nmr|startingmaterial|"
    r"starting material|solvent|analyse)(\d*)$",
    re.IGNORECASE,
)


def parse_plate_type(raw: str) -> tuple[str, int]:
    """Convert a legacy plate-type string to ``(role, role_index)``.

    This function is primarily used for ingesting legacy recipe data from
    JSON files that use the old string-based plate type format. New code
    should work with ``(role, role_index)`` tuples directly.

    Examples::

        "startingmaterial"   → ("startingmaterial", 1)
        "starting material"  → ("startingmaterial", 1)
        "workup1"            → ("workup", 1)
        "workup2"            → ("workup", 2)
        "reaction"           → ("reaction", 1)
    """
    raw = raw.strip()
    m = _PLATE_TYPE_RE.match(raw)
    if not m:
        # Fall-through: treat as-is with index 1
        return raw.lower().replace(" ", ""), 1
    role = m.group(1).lower().replace(" ", "")
    idx = int(m.group(2)) if m.group(2) else 1
    return role, idx


def unparse_plate_type(role: str, role_index: int = 1) -> str:
    """Convert a ``(role, role_index)`` pair back to a legacy plate-type string.

    This is the inverse of :func:`parse_plate_type`.  It is primarily
    used for logging, display names, and backward-compatibility layers
    that still expect the old string format.

    Examples::

        ("workup", 2)         → "workup2"
        ("reaction", 1)       → "reaction"
        ("startingmaterial", 1) → "startingmaterial"

    DEPRECATED: Prefer passing role/role_index directly to new APIs.
    """
    # For index 1, omit the number (matches legacy convention)
    if role_index == 1:
        return role
    return f"{role}{role_index}"


# ── Recipe lookups ─────────────────────────────────────────────────


def get_latest_recipe(
    reaction_class: str,
    name: str = "standard",
    version: int = None,
) -> Recipe:
    """Return the latest (or a specific) Recipe for a reaction class + name.

    Raises ``Recipe.DoesNotExist`` when not found.
    """
    qs = Recipe.objects.filter(reaction_class=reaction_class, name=name)
    if version is not None:
        return qs.get(version=version)
    return qs.order_by("-version").first() or qs.none().get()  # raises DoesNotExist


def recipe_exists(reaction_class: str) -> bool:
    """Check whether *any* recipe exists for the given reaction class."""
    return Recipe.objects.filter(reaction_class=reaction_class).exists()


def get_recipe_intramolecular(
    reaction_class: str, name: str = "standard"
) -> bool:
    """Return the ``intramolecular`` flag for a recipe."""
    return get_latest_recipe(reaction_class, name).intramolecular


def get_recipe_yield(reaction_class: str, name: str = "standard") -> float:
    """Return estimated yield as a fraction (e.g. 0.85)."""
    recipe = get_latest_recipe(reaction_class, name)
    return (recipe.estimated_yield or 100) / 100


def get_recipe_smarts(
    reaction_class: str, name: str = "standard"
) -> list[str]:
    """Return the ``reaction_smarts`` list."""
    return get_latest_recipe(reaction_class, name).reaction_smarts


def get_recipe_stir_temperature(
    reaction_class: str, name: str = "standard"
) -> Optional[float]:
    """Return the temperature from the first stir session of a recipe."""
    recipe = get_latest_recipe(reaction_class, name)
    stir = (
        RecipeStirAction.objects.filter(session__recipe=recipe)
        .order_by("session__session_number", "action_number")
        .first()
    )
    return stir.temperature if stir else None


# ── Type alias for recipe action model instances ───────────────────

from typing import Union

RecipeAction = Union[
    RecipeAddAction, RecipeStirAction, RecipeExtractAction, RecipeMixAction
]


# ── Session / action queries (returns model instances) ─────────────


def get_session_recipe_actions(
    reaction_class: str,
    name: str,
    session_type: str,
    session_number: int,
    molecular_context: str = None,
) -> list[RecipeAction]:
    """Return a sorted list of recipe action **model instances** for a session.

    Returns a mixed-type list of ``RecipeAddAction``, ``RecipeStirAction``,
    ``RecipeExtractAction``, and ``RecipeMixAction`` instances, sorted by
    ``action_number``.

    For reaction sessions, ``molecular_context`` filters add actions to
    those matching the context (or ``NULL``).  Pass ``None`` for workup /
    stir / analyse sessions.

    Returns an empty list when the session or actions are not found.
    """
    recipe = get_latest_recipe(reaction_class, name)
    session = (
        recipe.action_sessions.filter(
            session_type=session_type,
            session_number=session_number,
        )
        .first()
    )
    if session is None:
        return []

    return collect_session_actions(session, molecular_context)


def collect_session_actions(
    session: RecipeActionSession,
    molecular_context: str = None,
) -> list[RecipeAction]:
    """Collect all action model instances from a session, sorted by
    ``action_number``.

    For reaction sessions with a ``molecular_context``, only add actions
    matching the context (or ``NULL``) are included.  Non-add actions
    (stir, extract, mix) are always included.
    """
    tagged: list[tuple[int, RecipeAction]] = []

    # ── add actions ─────────────────────────────────────────
    add_qs = session.add_actions.all()
    if molecular_context and session.session_type == "reaction":
        add_qs = add_qs.filter(
            Q(molecular_context=molecular_context) | Q(molecular_context__isnull=True)
        )
    for add in add_qs.order_by("action_number"):
        tagged.append((add.action_number, add))

    # ── stir actions ────────────────────────────────────────
    for stir in session.stir_actions.all().order_by("action_number"):
        tagged.append((stir.action_number, stir))

    # ── extract actions ─────────────────────────────────────
    for ext in session.extract_actions.all().order_by("action_number"):
        tagged.append((ext.action_number, ext))

    # ── mix actions ─────────────────────────────────────────
    for mix in session.mix_actions.all().order_by("action_number"):
        tagged.append((mix.action_number, mix))

    tagged.sort(key=lambda t: t[0])
    return [entry[1] for entry in tagged]


# ── Clean dict serialisation (for recipe generator) ────────────────
# The generator needs a mutable dict representation to produce recipe
# variants.  This format mirrors the DB schema — no legacy nesting.


def action_to_dict(action: RecipeAction) -> dict:
    """Serialise a single recipe action model instance to a flat dict."""
    if isinstance(action, RecipeAddAction):
        return {
            "type": "add",
            "action_number": action.action_number,
            "molecular_context": action.molecular_context,
            "material_smarts": action.material_smarts,
            "material_smiles": action.material_smiles,
            "equivalents": action.equivalents,
            "quantity_unit": action.quantity_unit,
            "solvent": action.solvent,
            "concentration": action.concentration,
            "density": action.density,
            "from_plate_role": action.from_plate_role,
            "from_plate_role_index": action.from_plate_role_index,
            "to_plate_role": action.to_plate_role,
            "to_plate_role_index": action.to_plate_role_index,
        }
    if isinstance(action, RecipeStirAction):
        return {
            "type": "stir",
            "action_number": action.action_number,
            "temperature": action.temperature,
            "temperature_unit": action.temperature_unit,
            "duration": action.duration,
            "duration_unit": action.duration_unit,
            "plate_role": action.plate_role,
            "plate_role_index": action.plate_role_index,
        }
    if isinstance(action, RecipeExtractAction):
        return {
            "type": "extract",
            "action_number": action.action_number,
            "layer": action.layer,
            "volume": action.volume,
            "bottom_layer_volume": action.bottom_layer_volume,
            "smiles": action.smiles,
            "solvent": action.solvent,
            "concentration": action.concentration,
            "from_plate_role": action.from_plate_role,
            "from_plate_role_index": action.from_plate_role_index,
            "to_plate_role": action.to_plate_role,
            "to_plate_role_index": action.to_plate_role_index,
        }
    if isinstance(action, RecipeMixAction):
        return {
            "type": "mix",
            "action_number": action.action_number,
            "plate_role": action.plate_role,
            "plate_role_index": action.plate_role_index,
            "repetitions": action.repetitions,
        }
    raise TypeError(f"Unknown recipe action type: {type(action)}")


def recipe_to_dict(reaction_class: str, name: str = "standard") -> dict:
    """Serialise a full recipe from the DB to a clean dict.

    Structure mirrors the DB schema::

        {
            "estimated_yield": int,
            "reaction_smarts": [...],
            "sessions": [
                {
                    "session_number": 1,
                    "session_type": "reaction",
                    "driver": "robot",
                    "continuation": False,
                    "actions": [
                        {"type": "add", "action_number": 1, ...},
                        ...
                    ]
                },
                ...
            ]
        }
    """
    recipe = get_latest_recipe(reaction_class, name)
    sessions = []
    for session in recipe.action_sessions.order_by("session_number"):
        # Get ALL actions (no molecular_context filter) for the full view
        actions = collect_session_actions(session)
        sessions.append({
            "session_number": session.session_number,
            "session_type": session.session_type,
            "driver": session.driver,
            "continuation": session.continuation,
            "actions": [action_to_dict(a) for a in actions],
        })
    return {
        "estimated_yield": recipe.estimated_yield or 100,
        "reaction_smarts": recipe.reaction_smarts,
        "sessions": sessions,
    }


# ── Plate role extraction from recipe actions ──────────────────────


def get_session_destination_plates(
    reaction_class: str,
    name: str,
    session_type: str,
    session_number: int,
    molecular_context: str = None,
) -> list[tuple[str, int]]:
    """Extract unique destination plate (role, index) pairs from a recipe session.

    Scans all actions in the session and collects the ``to_plate_role`` /
    ``to_plate_role_index`` (for add/extract actions) or ``plate_role`` /
    ``plate_role_index`` (for stir/mix actions) values.

    Returns a sorted list of unique ``(role, index)`` tuples representing
    all plates that need to be created for this session.

    Parameters
    ----------
    reaction_class : str
        The reaction class (e.g., "amide-coupling").
    name : str
        Recipe name (e.g., "standard").
    session_type : str
        Session type (e.g., "reaction", "workup", "analyse").
    session_number : int
        The session number within the recipe.
    molecular_context : str, optional
        Filter for intermolecular/intramolecular (reaction sessions only).

    Returns
    -------
    list[tuple[str, int]]
        Sorted list of (plate_role, plate_role_index) pairs.
    """
    actions = get_session_recipe_actions(
        reaction_class=reaction_class,
        name=name,
        session_type=session_type,
        session_number=session_number,
        molecular_context=molecular_context,
    )

    plates: set[tuple[str, int]] = set()

    for action in actions:
        if isinstance(action, RecipeAddAction):
            plates.add((action.to_plate_role, action.to_plate_role_index))
        elif isinstance(action, RecipeExtractAction):
            plates.add((action.to_plate_role, action.to_plate_role_index))
        elif isinstance(action, RecipeStirAction):
            plates.add((action.plate_role, action.plate_role_index))
        elif isinstance(action, RecipeMixAction):
            plates.add((action.plate_role, action.plate_role_index))

    return sorted(plates, key=lambda p: (p[0], p[1]))


def get_session_source_plates(
    reaction_class: str,
    name: str,
    session_type: str,
    session_number: int,
    molecular_context: str = None,
) -> list[tuple[str, int]]:
    """Extract unique source plate (role, index) pairs from a recipe session.

    Scans add/extract actions and collects their ``from_plate_role`` /
    ``from_plate_role_index`` values.

    Returns a sorted list of unique ``(role, index)`` tuples representing
    all source plates referenced by this session's actions.

    Parameters
    ----------
    reaction_class : str
        The reaction class (e.g., "amide-coupling").
    name : str
        Recipe name (e.g., "standard").
    session_type : str
        Session type (e.g., "reaction", "workup", "analyse").
    session_number : int
        The session number within the recipe.
    molecular_context : str, optional
        Filter for intermolecular/intramolecular (reaction sessions only).

    Returns
    -------
    list[tuple[str, int]]
        Sorted list of (plate_role, plate_role_index) pairs.
    """
    actions = get_session_recipe_actions(
        reaction_class=reaction_class,
        name=name,
        session_type=session_type,
        session_number=session_number,
        molecular_context=molecular_context,
    )

    plates: set[tuple[str, int]] = set()

    for action in actions:
        if isinstance(action, RecipeAddAction):
            plates.add((action.from_plate_role, action.from_plate_role_index))
        elif isinstance(action, RecipeExtractAction):
            plates.add((action.from_plate_role, action.from_plate_role_index))

    return sorted(plates, key=lambda p: (p[0], p[1]))
