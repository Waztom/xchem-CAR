"""Chemistry service.

Centralises SMILES-related look-ups, molecular-weight calculations,
and product-SMILES retrieval that were previously duplicated across
api.py, db_utils.py, and createmodels.py.

All functions are pure queries or computations — no HTTP/DRF
dependency.
"""
from __future__ import annotations

from typing import List

from ..models import Batch, Method, Well


def get_ot_batch_product_smiles(batch_id: int) -> list:
    """Return product SMILES for a batch in well-index (execution) order.

    Traverses Batch → Targets → Methods → Wells (``role="reaction"``),
    and returns SMILES ordered by ``Well.index``.

    Parameters
    ----------
    batch_id : int
        Primary key of the :class:`Batch`.

    Returns
    -------
    list[str]
        Product SMILES in execution order.
    """
    targetqs = Batch.objects.get(id=batch_id).targets.all()
    methodqs = Method.objects.filter(target_id__in=targetqs)
    wellsqs = (
        Well.objects.filter(method_id__in=methodqs, role="reaction")
        .order_by("index")
        .distinct()
    )
    return list(wellsqs.values_list("smiles", flat=True))


def get_batch_product_mws(batch_id: int) -> List[float]:
    """Return molecular weights of batch products.

    Delegates to :func:`backend.chem_utils.get_mws` via
    :func:`backend.db_utils.get_batch_reaction_product_mws`.

    Parameters
    ----------
    batch_id : int
        Primary key of the :class:`Batch`.

    Returns
    -------
    list[float]
        Molecular weights of the batch's reaction products.
    """
    from ..db_utils import get_batch_reaction_product_mws
    return get_batch_reaction_product_mws(batch_id=batch_id)


def get_batch_target_mws(batch_id: int) -> List[float]:
    """Return molecular weights of batch targets.

    Parameters
    ----------
    batch_id : int
        Primary key of the :class:`Batch`.

    Returns
    -------
    list[float]
        Molecular weights of the batch's target SMILES.
    """
    from ..db_utils import get_batch_target_mws as _get_batch_target_mws
    return _get_batch_target_mws(batch_id=batch_id)
