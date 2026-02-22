"""OT protocol creation tasks: session orchestration, script generation, zip packaging."""
from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING
from zipfile import ZipFile

from celery import shared_task, current_task
from django.conf import settings
from django.db import transaction

from ..createmodels import CreateEncodedActionModels
from ..db_utils import (
    get_action_session_query_set,
    get_action_session_sequence_numbers,
    get_action_session_types,
    get_batch_reactions,
    get_batch_tag,
    get_grouped_action_session_sequences,
    get_grouped_action_session_types,
    get_max_reaction_number,
    get_ot_batch_protocol_query_set,
    get_reactions_to_do,
    group_reactions,
)
from ..models import (
    Batch,
    CompoundOrder,
    OTBatchProtocol,
    OTProject,
    OTScript,
    OTSession,
    SolventPrep,
)
from ..opentrons.otsession import SessionOrchestrator
from ..opentrons.otwriter import ScriptGenerator
from ..recipe_utils import get_recipe_intramolecular

if TYPE_CHECKING:
    from django.db.models import QuerySet
    from ..models import ActionSession

logger = logging.getLogger(__name__)


def _get_custom_sm_csv_path(custom_SM_files, batchid):
    """Return the custom starting-material CSV path for *batchid*, or None."""
    if custom_SM_files and str(batchid) in custom_SM_files:
        return custom_SM_files[str(batchid)]
    return None


def _process_reaction_step_sessions(
    reaction_ids,
    reactionstep,
    otbatchprotocolobj,
    custom_sm_csv_path,
    batchtag,
    use_multichannel=False,
):
    """Shared logic for processing a single reaction step.

    Groups the action sessions by sequence number and type, filters for
    robot-driven sessions, and creates OT sessions for each group.

    Parameters
    ----------
    reaction_ids : list[int]
        Reaction IDs to process for this step.
    reactionstep : int
        1-based reaction step number.
    otbatchprotocolobj : OTBatchProtocol
        The parent batch protocol object.
    custom_sm_csv_path : str or None
        Path to a custom starting-material CSV, if any.
    batchtag : str or None
        Batch tag for session naming.
    use_multichannel : bool, optional
        When True, starter plates are laid out for multichannel
        pipette column transfers where possible (default False).
    """
    actionsessionqueryset = get_action_session_query_set(
        reaction_ids=reaction_ids
    )
    sessionnumbers = get_action_session_sequence_numbers(
        actionsessionqueryset=actionsessionqueryset
    )
    grouped_sequences = get_grouped_action_session_sequences(
        sessionnumbers=sessionnumbers,
        actionsessionqueryset=actionsessionqueryset,
    )
    for group_action_session in grouped_sequences:
        action_session_types = get_action_session_types(
            actionsessionqueryset=group_action_session
        )
        grouped_types = get_grouped_action_session_types(
            actionsessiontypes=action_session_types,
            actionsessionqueryset=group_action_session,
        )
        for action_session_type_qs in grouped_types:
            robot_qs = action_session_type_qs.filter(driver="robot")
            if not robot_qs:
                continue
            create_multiple_ot_sessions(
                reactionstep=reactionstep,
                otbatchprotocolobj=otbatchprotocolobj,
                actionsessionqueryset=robot_qs,
                customSMcsvpath=custom_sm_csv_path,
                batchtag=batchtag,
                use_multichannel=use_multichannel,
            )


@shared_task
def create_ot_script(batchids: list, protocol_name: str, custom_SM_files: dict = None, use_multichannel: bool = False):
    """Create OT scripts and starting plates for a list of batch IDs.

    Parameters
    ----------
    batchids : list
        Batch IDs to create OT scripts for.
    protocol_name : str
        Name of the protocol to create.
    custom_SM_files : dict or None
        Optional custom starting-material files keyed by batch ID (as str)
        to CSV file path.  Format: ``{"batch_id": "path/to/csv"}``
    """
    task_summary = {}

    with transaction.atomic():
        otprojectobj = OTProject()
        projectobj = Batch.objects.get(id=batchids[0]).project_id
        otprojectobj.project_id = projectobj
        otprojectobj.name = protocol_name
        otprojectobj.save()

        for batchid in batchids:
            reactionqueryset = get_batch_reactions(batchid=batchid)
            actionsessionqueryset = get_action_session_query_set(
                reaction_ids=reactionqueryset
            )

            # Ensure action models exist for every reaction in this batch
            if not actionsessionqueryset:
                for reactionobj in reactionqueryset:
                    reactant_pair_smiles = list(
                        reactionobj.reactants.all().values_list(
                            "smiles", flat=True
                        )
                    )
                    intramolecular = (
                        get_recipe_intramolecular(
                            reactionobj.reactionclass, reactionobj.recipe
                        )
                        and len(reactant_pair_smiles) == 1
                    )
                    CreateEncodedActionModels(
                        reaction_class=reactionobj.reactionclass,
                        recipe_name=reactionobj.recipe,
                        intramolecular=intramolecular,
                        target_id=reactionobj.method_id.target_id.id,
                        reaction_id=reactionobj.id,
                        reactant_pair_smiles=reactant_pair_smiles,
                    )

            batchtag = get_batch_tag(batchid=batchid)
            otbatchprotocolqueryset = get_ot_batch_protocol_query_set(
                batch_id=batchid
            )

            # Re-use an existing completed protocol if one exists
            if otbatchprotocolqueryset and otbatchprotocolqueryset[0].zipfile:
                otbatchprotocolobj = otbatchprotocolqueryset[0]
                otbatchprotocolobj.otproject_id = otprojectobj
                otbatchprotocolobj.celery_taskid = current_task.request.id
                otbatchprotocolobj.save()
                task_summary[batchid] = True
                continue

            # Create a fresh batch protocol
            otbatchprotocolobj = OTBatchProtocol()
            otbatchprotocolobj.batch_id = Batch.objects.get(id=batchid)
            otbatchprotocolobj.otproject_id = otprojectobj
            otbatchprotocolobj.celery_taskid = current_task.request.id
            otbatchprotocolobj.save()

            custom_sm_csv_path = _get_custom_sm_csv_path(
                custom_SM_files, batchid
            )
            max_reaction_number = get_max_reaction_number(
                reactionqueryset=reactionqueryset
            )
            grouped_reaction_querysets = group_reactions(
                reactionqueryset=reactionqueryset,
                maxreactionnumber=max_reaction_number,
            )

            for step_index, step_reactions in enumerate(
                grouped_reaction_querysets
            ):
                # For steps beyond the first, filter out reactions that
                # failed QC — if none remain, stop processing further steps.
                if step_index > 0:
                    step_reactions = get_reactions_to_do(
                        groupreactionqueryset=step_reactions
                    )
                    if not step_reactions:
                        break

                reaction_ids = [
                    reaction.id for reaction in step_reactions
                ]
                _process_reaction_step_sessions(
                    reaction_ids=reaction_ids,
                    reactionstep=step_index + 1,
                    otbatchprotocolobj=otbatchprotocolobj,
                    custom_sm_csv_path=custom_sm_csv_path,
                    batchtag=batchtag,
                    use_multichannel=use_multichannel,
                )

            zip_protocol = ZipOTBatchProtocol(
                otbatchprotocolobj=otbatchprotocolobj, batchtag=batchtag
            )
            task_summary[batchid] = not zip_protocol.errors

    # Clean up temporary CSV files outside the transaction
    if custom_SM_files:
        for filepath in custom_SM_files.values():
            try:
                if os.path.exists(filepath):
                    os.remove(filepath)
            except Exception as e:
                logger.warning(f"Error cleaning up file {filepath}: {e}")

    return task_summary, otprojectobj.id


def create_multiple_ot_sessions(
    reactionstep: int,
    otbatchprotocolobj: OTBatchProtocol,
    actionsessionqueryset: "QuerySet[ActionSession]",
    customSMcsvpath: str = None,
    max_reactions_per_session: int = None,
    batchtag: str = None,
    use_multichannel: bool = False,
):
    """Splits a large reaction set into multiple OT sessions to prevent deck overflow.

    For each session, immediately performs OTWrite to ensure database changes from
    one session are available to subsequent sessions.
    """
    # Get unique reaction IDs from action sessions
    reaction_ids = list(
        set(actionsessionqueryset.values_list("reaction_id", flat=True))
    )
    total_reactions = len(reaction_ids)

    logger.info(
        f"Processing {total_reactions} unique reactions across {actionsessionqueryset.count()} action sessions"
    )

    if max_reactions_per_session is None:
        # Start with a smaller default to avoid issues
        max_reactions_per_session = 600 if total_reactions > 600 else total_reactions

    # Calculate initial number of sessions needed
    num_sessions = (
        total_reactions + max_reactions_per_session - 1
    ) // max_reactions_per_session

    logger.info(
        f"Initial plan: Creating {num_sessions} sessions with max {max_reactions_per_session} reactions per session"
    )

    created_sessions = []
    reaction_groups = []

    # Create initial reaction groups
    for i in range(0, total_reactions, max_reactions_per_session):
        batch = reaction_ids[i : i + max_reactions_per_session]
        reaction_groups.append(batch)

    # Try to create sessions for each reaction group
    for group_index, reaction_group in enumerate(reaction_groups):
        try:
            # Find action sessions related to this reaction group
            group_action_sessions = actionsessionqueryset.filter(
                reaction_id__in=reaction_group
            )
            group_action_session_ids = list(
                group_action_sessions.values_list("id", flat=True)
            )

            logger.info(
                f"Attempting session {group_index + 1} with {len(reaction_group)} reactions"
            )

            # Create a session for this group
            orchestrator = SessionOrchestrator(
                reactionstep=reactionstep,
                otbatchprotocolobj=otbatchprotocolobj,
                actionsessionqueryset=group_action_sessions,
                customSMcsvpath=customSMcsvpath,
                use_multichannel=use_multichannel,
            )

            orchestrator.execute()

            # Validate the otsessionobj exists before using it
            if not orchestrator.otsessionobj:
                logger.error(
                    "Orchestrator execution completed but didn't create an OTSession object"
                )
                raise ValueError(
                    "Session execution failed to create necessary database objects"
                )

            # Now it's safe to use orchestrator.otsessionobj.id
            session_batchtag = (
                f"{batchtag}_session_{orchestrator.otsessionobj.id}"
                if batchtag
                else f"session_{orchestrator.otsessionobj.id}"
            )

            logger.info(f"Running ScriptGenerator for session {group_index + 1}")
            try:
                # Create ScriptGenerator instance and generate the script
                script_generator = ScriptGenerator(
                    batchtag=session_batchtag,
                    otsessionobj=orchestrator.otsessionobj,
                    actionsession_ids=group_action_session_ids,
                )

                # Generate the script and get the filepath
                filepath = script_generator.generate_script()

                # Create a result object similar to what OTWrite would have returned
                otwrite_result = {
                    "success": True,
                    "filepath": filepath,
                    "session_id": orchestrator.otsessionobj.id,
                }

                # If we get here, everything worked - store session info
                created_sessions.append(
                    {
                        "session": orchestrator,
                        "action_session_ids": group_action_session_ids,
                        "otwrite_result": otwrite_result,
                    }
                )

                logger.info(f"Successfully created and wrote session {group_index + 1}")

            except Exception as e:
                logger.error(
                    f"Error in script generation for session {group_index + 1}: {str(e)}"
                )
                otwrite_result = {
                    "success": False,
                    "error": str(e),
                    "session_id": orchestrator.otsessionobj.id,
                }
                raise

        except ValueError as ve:
            # Handle deck slot issues by splitting the group
            error_msg = str(ve)
            logger.error(f"ValueError in session {group_index + 1}: {error_msg}")
            if (
                error_msg == "No deck slots available - cannot create more plates"
                and len(reaction_group) > 1
            ):
                logger.warning(f"Deck full error detected: {str(ve)}")
                logger.warning(
                    f"Splitting group of {len(reaction_group)} reactions in half"
                )

                # Recursive call with smaller max_reactions_per_session for just this group
                smaller_max = max(1, len(reaction_group) // 2)
                logger.info(f"Trying with batch size: {smaller_max}")

                sub_action_sessions = actionsessionqueryset.filter(
                    reaction_id__in=reaction_group
                )
                sub_sessions = create_multiple_ot_sessions(
                    reactionstep=reactionstep,
                    otbatchprotocolobj=otbatchprotocolobj,
                    actionsessionqueryset=sub_action_sessions,
                    customSMcsvpath=customSMcsvpath,
                    max_reactions_per_session=smaller_max,
                    batchtag=batchtag,
                    use_multichannel=use_multichannel,
                )
                created_sessions.extend(sub_sessions)
            else:
                # If not a deck slot issue or already at minimum group size, propagate the error
                logger.error(f"ValueError in session {group_index + 1}: {str(ve)}")
                raise

    logger.info(
        f"Created {len(created_sessions)} sessions to handle {total_reactions} reactions"
    )
    return created_sessions


class ZipOTBatchProtocol:
    """Creates a zip archive of compound orders, OT scripts, and solvent
    preps for a batch protocol.

    The zip is built during ``__init__``, written to the
    ``OTBatchProtocol`` model's ``zipfile`` field, and the temp file is
    cleaned up.  Any missing querysets are recorded in ``self.errors``.
    """

    # (Model, file-field attribute name, zip subdirectory)
    _ARTIFACT_TYPES = [
        (SolventPrep, "solventprepcsv", "solventprep"),
        (CompoundOrder, "ordercsv", "compoundorders"),
        (OTScript, "otscript", "otscripts"),
    ]

    def __init__(self, otbatchprotocolobj, batchtag: str):
        self.errors = []
        self.otbatchprotocolobj = otbatchprotocolobj
        self.mediaroot = settings.MEDIA_ROOT
        self.zipfn = f"batch-{batchtag}-protocol.zip"
        self.ziptmpfp = os.path.join(
            settings.MEDIA_ROOT, "tmp", "batchprotocoltmp.zip"
        )
        self._build_zip()

    # ------------------------------------------------------------------
    # Zip construction
    # ------------------------------------------------------------------

    def _build_zip(self):
        """Orchestrate zip creation: query sessions, collect artifacts,
        write to media, then clean up the temp file."""
        ot_sessions = self._get_ot_sessions()
        if not ot_sessions:
            return

        with ZipFile(self.ziptmpfp, "w") as ziparchive:
            for session in ot_sessions:
                for model, file_field, dest_dir in self._ARTIFACT_TYPES:
                    self._add_artifacts_to_zip(
                        ziparchive, session, model, file_field, dest_dir
                    )

            # Add session visualization HTML files (not model-backed)
            self._add_session_visualizations(ziparchive)

        self._write_zip_to_media()
        self._delete_tmp_zip()

    def add_warning(self, function: str, errorwarning: str):
        """Record a warning/error encountered during zip construction."""
        self.errors.append(
            {"function": function, "errorwarning": errorwarning}
        )

    # ------------------------------------------------------------------
    # Queryset helpers
    # ------------------------------------------------------------------

    def _get_ot_sessions(self):
        """Return the OTSession queryset for the batch protocol, or
        ``None`` if no sessions exist (records a warning)."""
        qs = OTSession.objects.filter(
            otbatchprotocol_id=self.otbatchprotocolobj
        )
        if not qs:
            self.add_warning(
                function="_get_ot_sessions",
                errorwarning="No queryset found",
            )
            return None
        return qs

    def _get_related_queryset(self, model, session):
        """Return the queryset for *model* filtered to *session*, or
        ``None`` if empty (records a warning)."""
        qs = model.objects.filter(otsession_id=session)
        if not qs:
            self.add_warning(
                function=f"_get_related_queryset({model.__name__})",
                errorwarning="No queryset found",
            )
            return None
        return qs

    # ------------------------------------------------------------------
    # Zip I/O
    # ------------------------------------------------------------------

    def _add_artifacts_to_zip(
        self, ziparchive, session, model, file_field, dest_dir
    ):
        """Add all artifacts of a given type from *session* to the zip."""
        qs = self._get_related_queryset(model, session)
        if not qs:
            return
        for obj in qs:
            filepath = os.path.join(
                self.mediaroot, getattr(obj, file_field).name
            )
            arcname = os.path.join(dest_dir, os.path.basename(filepath))
            ziparchive.write(filename=filepath, arcname=arcname)

    def _add_session_visualizations(self, ziparchive):
        """Add any session visualization HTML files to the zip archive.

        Visualization files are written to ``MEDIA_ROOT/session_visualizations/``
        by ``SessionVisualizer`` during script generation.  We glob for all
        ``.html`` files in that directory and add them under a
        ``sessionvisualizations/`` subdirectory in the zip.
        """
        viz_dir = os.path.join(self.mediaroot, "session_visualizations")
        if not os.path.isdir(viz_dir):
            return
        for fname in os.listdir(viz_dir):
            if fname.endswith(".html"):
                filepath = os.path.join(viz_dir, fname)
                arcname = os.path.join("sessionvisualizations", fname)
                ziparchive.write(filename=filepath, arcname=arcname)

    def _write_zip_to_media(self):
        """Persist the zip archive to the OTBatchProtocol model."""
        with open(self.ziptmpfp, "rb") as zf:
            self.otbatchprotocolobj.zipfile.save(self.zipfn, zf)
            self.otbatchprotocolobj.save()

    def _delete_tmp_zip(self):
        """Remove the temporary zip file."""
        os.remove(self.ziptmpfp)
