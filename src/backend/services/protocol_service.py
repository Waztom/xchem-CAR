"""Protocol service.

Handles OT project creation orchestration — parsing custom
starting-material files, saving them to temp storage, and dispatching
the ``create_ot_script`` Celery task.

Extracted from ``OTProjectViewSet.create_ot_project`` so the ViewSet
does only request parsing and JSON response assembly.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

from celery.result import AsyncResult

logger = logging.getLogger(__name__)


def save_tmp_file(myfile) -> str:
    """Persist an uploaded file to ``MEDIA_ROOT/tmp/`` and return its path.

    Parameters
    ----------
    myfile
        A Django ``UploadedFile`` (or file-like with ``.name`` and ``.read()``).

    Returns
    -------
    str
        Absolute filesystem path to the saved temporary file.
    """
    import os

    from django.conf import settings
    from django.core.files.base import ContentFile
    from django.core.files.storage import default_storage

    name = myfile.name
    path = default_storage.save("tmp/" + name, ContentFile(myfile.read()))
    return str(os.path.join(settings.MEDIA_ROOT, path))


def initiate_ot_project(
    batch_ids: List[int],
    protocol_name: str,
    custom_files: Optional[Dict[str, "UploadedFile"]] = None,
) -> str:
    """Dispatch the ``create_ot_script`` Celery task.

    Saves any custom starting-material CSVs to temp storage and kicks
    off the asynchronous OT-script generation pipeline.

    Parameters
    ----------
    batch_ids : list[int]
        Batch IDs to include in the protocol.
    protocol_name : str
        Human-readable name for the OT project.
    custom_files : dict or None
        Mapping of ``batch_id`` (str) → ``UploadedFile`` for custom
        starting-material CSVs.  ``None`` when not applicable.

    Returns
    -------
    str
        The Celery task ID for status polling.
    """
    from ..tasks import create_ot_script

    starting_material_paths: Optional[Dict[str, str]] = None
    if custom_files:
        starting_material_paths = {}
        for batch_id_str, uploaded_file in custom_files.items():
            starting_material_paths[batch_id_str] = save_tmp_file(uploaded_file)

    task = create_ot_script.delay(
        batchids=batch_ids,
        protocol_name=protocol_name,
        custom_SM_files=starting_material_paths,
    )
    return task.id


def poll_task_status(task_id: str) -> dict:
    """Check the status of a Celery task and return a JSON-serialisable dict.

    This is a shared helper for the ``get_task_status`` endpoints on
    ``ProjectViewSet``, ``BatchViewSet``, and ``OTProjectViewSet``.

    Parameters
    ----------
    task_id : str
        The Celery task ID.

    Returns
    -------
    dict
        Always contains ``"task_status"``; may also contain
        ``"result"``, ``"traceback"``, or task-specific keys depending
        on the status.
    """
    task = AsyncResult(task_id)

    if task.status == "FAILURE":
        return {"task_status": task.status, "traceback": str(task.traceback)}

    if task.status == "SUCCESS":
        return {"task_status": task.status, "result": task.get()}

    # PENDING / STARTED / RETRY / etc.
    return {"task_status": task.status}
