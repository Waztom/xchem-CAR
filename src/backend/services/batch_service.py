"""Batch lifecycle service.

Encapsulates the business logic for batch/target/method cloning and
reaction status management — previously embedded in api.py ViewSet
methods.

All functions operate on model objects or IDs; none depend on
``request`` or any HTTP/DRF layer.
"""
from __future__ import annotations

from typing import List

from django.db import transaction

from ..models import Batch, Method, Reaction, Target


# ------------------------------------------------------------------
# Cloning
# ------------------------------------------------------------------


def clone_target(target_obj: Target, batch_obj: Batch) -> Target:
    """Deep-clone a Target (and its CatalogEntries) into *batch_obj*.

    The original image file is reused (path copied, not duplicated).

    Parameters
    ----------
    target_obj : Target
        The source target to clone.
    batch_obj : Batch
        The destination batch.

    Returns
    -------
    Target
        The newly saved clone.
    """
    related_catalogentry_queryset = target_obj.catalogentries.all().order_by("id")

    original_image_path = target_obj.image.name if target_obj.image else None

    target_obj.pk = None
    target_obj.batch_id = batch_obj

    if original_image_path:
        target_obj.image.name = original_image_path

    target_obj.save()

    for catalogentry_obj in related_catalogentry_queryset:
        catalogentry_obj.pk = None
        catalogentry_obj.target_id = target_obj
        catalogentry_obj.save()

    return target_obj


def clone_method(method_obj: Method, target_obj: Target) -> None:
    """Deep-clone a Method and its full tree into *target_obj*.

    Clones: Method → Reactions → Products, Reactants → CatalogEntries.
    Image file paths are reused.

    Parameters
    ----------
    method_obj : Method
        The source method.
    target_obj : Target
        The destination target (must already be saved).
    """
    related_reaction_queryset = method_obj.reactions.all().order_by("id")
    method_obj.pk = None
    method_obj.target_id = target_obj
    method_obj.save()

    for reaction_obj in related_reaction_queryset:
        product_obj = reaction_obj.products.all()[0]
        related_reactant_objs = reaction_obj.reactants.all().order_by("id")

        original_reaction_image = (
            reaction_obj.image.name if reaction_obj.image else None
        )
        original_product_image = product_obj.image.name if product_obj.image else None

        reaction_obj.pk = None
        reaction_obj.method_id = method_obj
        if original_reaction_image:
            reaction_obj.image.name = original_reaction_image
        reaction_obj.save()

        product_obj.pk = None
        product_obj.reaction_id = reaction_obj
        if original_product_image:
            product_obj.image.name = original_product_image
        product_obj.save()

        for reactant_obj in related_reactant_objs:
            related_catalogentry_objs = reactant_obj.catalogentries.all().order_by("id")
            reactant_obj.pk = None
            reactant_obj.reaction_id = reaction_obj
            reactant_obj.save()
            for catalog_obj in related_catalogentry_objs:
                catalog_obj.pk = None
                catalog_obj.reactant_id = reactant_obj
                catalog_obj.save()


def clone_batch(method_ids: List[int], batchtag: str) -> Batch:
    """Create a new batch by cloning selected methods.

    Finds the targets owning *method_ids*, creates a new child batch
    under the same project, and deep-clones each target with only the
    selected methods.

    Parameters
    ----------
    method_ids : list[int]
        Primary keys of the :class:`Method` objects to include.
    batchtag : str
        Name for the new batch.

    Returns
    -------
    Batch
        The newly created batch.

    Raises
    ------
    Target.DoesNotExist
        If no targets own the given method IDs.
    """
    target_query_set = (
        Target.objects.filter(methods__id__in=method_ids)
        .distinct()
        .order_by("id")
    )
    batch_obj = target_query_set[0].batch_id
    project_obj = batch_obj.project_id

    with transaction.atomic():
        batch_obj_new = Batch()
        batch_obj_new.project_id = project_obj
        batch_obj_new.batch_id = batch_obj
        batch_obj_new.batchtag = batchtag
        batch_obj_new.save()

        for target_obj in target_query_set:
            method_query_set_to_clone = (
                Method.objects.filter(target_id=target_obj)
                .filter(pk__in=method_ids)
                .order_by("id")
            )
            target_obj_clone = clone_target(
                target_obj=target_obj, batch_obj=batch_obj_new
            )
            for method_obj in method_query_set_to_clone:
                clone_method(method_obj=method_obj, target_obj=target_obj_clone)

    return batch_obj_new


# ------------------------------------------------------------------
# Reaction status management
# ------------------------------------------------------------------


def mark_reactions_failed(reaction_ids: List[int]) -> dict:
    """Bulk-mark reactions as unsuccessful.

    Parameters
    ----------
    reaction_ids : list[int]
        Primary keys of the reactions to mark as failed.

    Returns
    -------
    dict
        ``{"reaction_ids": <list>}`` if rows were updated, or
        ``{"reaction_ids": None}`` if none matched.
    """
    qs = Reaction.objects.filter(id__in=reaction_ids)
    if qs.exists():
        qs.update(success=False)
        return {"reaction_ids": list(reaction_ids)}
    return {"reaction_ids": None}
