from rest_framework import viewsets
from rest_framework.decorators import action
from django.http import JsonResponse
from django.conf import settings
import os
import json
from celery.result import AsyncResult
import logging

logger = logging.getLogger(__name__)

# from viewer.tasks import check_services
import pandas as pd

from .tasks import (
    validate_file_upload,
    upload_manifold_reaction,
    upload_custom_reaction,
    upload_combi_custom_reaction,
    canonicalize_smiles,
)

# Import standard models
from .models import (
    Project,
    Batch,
    PubChemInfo,
    Target,
    Method,
    Reaction,
    Reactant,
    CatalogEntry,
    Product,
)

# Import action models
from .models import (
    ActionSession,
    AddAction,
    ExtractAction,
    MixAction,
    StirAction,
)

# Import OT Session models
from .models import (
    OTSession,
    Deck,
    Pipette,
    TipRack,
    Plate,
    Column,
    Well,
    OTProject,
    OTBatchProtocol,
    CompoundOrder,
    OTScript,
)

# Import standard serializers
from .serializers import (
    OTProjectSerializer,
    ProjectSerializer,
    ProjectSerializerAll,
    BatchSerializer,
    BatchSerializerAll,
    TargetSerializer,
    TargetSerializerAll,
    MethodSerializer,
    MethodSerializerAll,
    ReactionSerializer,
    ReactionSerializerAll,
    PubChemInfoSerializer,
    ProductSerializer,
    ProductSerializerAll,
    ReactantSerializer,
    ReactantSerializerAll,
    CatalogEntrySerializer,
)

# Import action serializers
from .serializers import (
    ActionSessionSerializer,
    AddActionSerializer,
    ExtractActionSerializer,
    MixActionSerializer,
    StirActionSerializer,
)

# Import OT Session serializers
from .serializers import (
    OTBatchProtocolSerializer,
    OTSessionSerializer,
    DeckSerializer,
    PipetteSerializer,
    TipRackSerializer,
    PlateSerializer,
    ColumnSerializer,
    WellSerializer,
    CompoundOrderSerializer,
    OTScriptSerializer,
)


from .services.batch_service import (
    clone_batch,
    mark_reactions_failed,
)
from .services.chemistry_service import get_ot_batch_product_smiles
from .services.protocol_service import (
    initiate_ot_project,
    save_tmp_file,
)


class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()

    def get_serializer_class(self):
        fetchall = self.request.GET.get("fetchall", None)
        return ProjectSerializerAll if fetchall == "yes" else ProjectSerializer

    @action(methods=["post"], detail=False, url_path="createproject")
    def create_project(self, request, pk=None):
        """Post method to create a new project

        Parameters
        ----------
        request: JSON
            Will have structure:

            {"projectname": str,
             "submittername": str,
             "submitterorganisation": str,
             "proteintarget": str
             "validate_choice": int,
             "API_choice": int,
             "fetch_pubchem": bool,
             "fetch_catalogue": bool,
            }

        """
        # check_services() # Uncomment when in production
        project_info = {}
        project_info["projectname"] = request.data["project_name"]
        project_info["submittername"] = request.data["submitter_name"]
        project_info["submitterorganisation"] = request.data["submitter_organisation"]
        project_info["proteintarget"] = request.data["protein_target"]
        validate_choice = request.data["validate_choice"]
        API_choice = request.data["API_choice"]
        fetch_pubchem = request.data.get("fetch_pubchem", True)
        fetch_catalogue = request.data.get("fetch_catalogue", True)

        csvfile = request.FILES["csv_file"]
        tmp_file = save_tmp_file(csvfile)

        if str(validate_choice) == "0":
            if str(API_choice) == "0":
                task = validate_file_upload.delay(
                    csv_fp=tmp_file, validate_type="retro-API"
                )

            if str(API_choice) == "1":
                task = validate_file_upload.delay(
                    csv_fp=tmp_file, validate_type="custom-chem"
                )

            if str(API_choice) == "2":
                task = validate_file_upload.delay(
                    csv_fp=tmp_file, validate_type="combi-custom-chem"
                )

        if str(validate_choice) == "1":
            if str(API_choice) == "0":
                task = (
                    validate_file_upload.s(
                        csv_fp=tmp_file,
                        validate_type="retro-API",
                        project_info=project_info,
                        validate_only=False,
                    )
                    | upload_manifold_reaction.s(
                        fetch_pubchem=fetch_pubchem,
                    )
                ).apply_async()

            if str(API_choice) == "1":
                task = (
                    validate_file_upload.s(
                        csv_fp=tmp_file,
                        validate_type="custom-chem",
                        project_info=project_info,
                        validate_only=False,
                    )
                    | upload_custom_reaction.s(
                        fetch_pubchem=fetch_pubchem,
                        fetch_catalogue=fetch_catalogue,
                    )
                ).apply_async()

            if str(API_choice) == "2":
                task = (
                    validate_file_upload.s(
                        csv_fp=tmp_file,
                        validate_type="combi-custom-chem",
                        project_info=project_info,
                        validate_only=False,
                    )
                    | upload_combi_custom_reaction.s(
                        fetch_pubchem=fetch_pubchem,
                        fetch_catalogue=fetch_catalogue,
                    )
                ).apply_async()

        data = {"task_id": task.id}
        return JsonResponse(data=data)

    @action(detail=False, methods=["get"], url_path="gettaskstatus")
    def get_task_status(self, request, pk=None):
        """Get method to get the Celery task status"""
        task_id = self.request.GET.get("task_id", None)
        if task_id:
            task = AsyncResult(task_id)
            if task.status == "FAILURE":
                data = {"task_status": task.status, "traceback": str(task.traceback)}
                return JsonResponse(data)

            if task.status == "SUCCESS":
                results = task.get()
                validate_dict = results[0]
                validated = results[1]
                project_info = results[2]

                if not project_info:
                    if validated:
                        data = {"task_status": task.status, "validated": True}
                        return JsonResponse(data)

                    if not validated:
                        errorsummary = json.dumps(validate_dict)
                        data = {
                            "task_status": task.status,
                            "validated": False,
                            "validation_errors": errorsummary,
                        }
                        return JsonResponse(data)

                if project_info:
                    project_id = project_info["project_id"]

                    if validated:
                        data = {
                            "task_status": task.status,
                            "validated": True,
                            "project_id": project_id,
                        }
                        return JsonResponse(data)

                    if not validated:
                        errorsummary = json.dumps(validate_dict)
                        data = {
                            "task_status": task.status,
                            "validated": False,
                            "validation_errors": errorsummary,
                        }

                        return JsonResponse(data)

            if task.status == "PENDING":
                data = {"task_status": task.status}
                return JsonResponse(data)


class BatchViewSet(viewsets.ModelViewSet):
    queryset = Batch.objects.all()
    filterset_fields = ["project_id"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get("fetchall", None)
        return BatchSerializerAll if fetchall == "yes" else BatchSerializer

    def create_batch(
        self, project_obj: Project, batch_node_obj: Batch, batchtag: str
    ) -> Batch:
        """Creates a batch object

        Parameters
        ----------
        project_obj: Project
            The project object that the new batch will be related to
        batch_node_obj: Batch
            The parent batch of the newly created batch
        batchtag: str
            The batch tag for the new batch

        Returns
        -------
        batch_obj: Batch
            The newly created batch object
        """
        batch_obj = Batch()
        batch_obj.project_id = project_obj
        batch_obj.batch_id = batch_node_obj
        batch_obj.batchtag = batchtag
        batch_obj.save()
        return batch_obj

    def create(self, request, **kwargs):
        """Post method to create a new batch by cloning selected methods.

        Parameters
        ----------
        request: JSON
            Will have structure:

            {"methodids": list,
             "batchtag": str,
            }

        The methodids are the methods selected for the new batch
        The batchtag is the name of the new batch to be created
        """
        method_ids = request.data["methodids"]
        batchtag = request.data["batchtag"]
        try:
            batch_obj_new = clone_batch(method_ids=method_ids, batchtag=batchtag)
            serialized_data = BatchSerializer(batch_obj_new).data
            if serialized_data:
                return JsonResponse(data=serialized_data)
            else:
                return JsonResponse(data="Something went wrong")
        except Exception as e:
            return JsonResponse(data={"message": "Something went wrong", "error": e})

    @action(methods=["post"], detail=False, url_path="canonicalizesmiles")
    def canonicalize_smiles(self, request, pk=None):
        """Post method to canonicalise a list or csv file of SMILES"""
        # check_services()
        if request.POST.get("smiles"):
            smiles = request.POST.getlist("smiles")
            task = canonicalize_smiles.delay(smiles=smiles)
            data = {"task_id": task.id}
            return JsonResponse(data=data)
        if len(request.FILES) != 0:
            csvfile = request.FILES["csv_file"]
            tmp_file = save_tmp_file(csvfile)
            task = canonicalize_smiles.delay(csvfile=tmp_file)
            data = {"task_id": task.id}
            return JsonResponse(data=data)

    @action(detail=False, methods=["get"], url_path="gettaskstatus")
    def get_task_status(self, request, pk=None):
        """Get method to check the Celery task status of the the
        SMILES being canonicalised
        """
        task_id = self.request.GET.get("task_id", None)
        if task_id:
            task = AsyncResult(task_id)
            if task.status == "FAILURE":
                data = {"task_status": task.status, "traceback": str(task.traceback)}
                return JsonResponse(data)

            if task.status == "SUCCESS":
                result = task.get()
                validated = result[0]

                if validated:
                    canonicalizedsmiles = result[1]
                    data = {
                        "task_status": task.status,
                        "canonicalizedsmiles": canonicalizedsmiles,
                    }
                    return JsonResponse(data)
                if not validated:
                    error_summary = result[1]
                    data = {"task_status": task.status, "error_summary": error_summary}
                    return JsonResponse(data)

            if task.status == "PENDING":
                data = {"task_status": task.status}
                return JsonResponse(data)

    @action(methods=["post"], detail=False, url_path="updatereactionsuccess")
    def update_reaction_success(self, request, pk=None):
        """Updates reactions to be set to be unsuccessful"""
        if request.POST.get("reaction_ids"):
            reaction_ids = request.POST.getlist("reaction_ids")
        if len(request.FILES) != 0:
            csvfile = request.FILES["csv_file"]
            reaction_ids = pd.read_csv(csvfile)["reaction_id"]
        data = mark_reactions_failed(reaction_ids)
        return JsonResponse(data=data)


class TargetViewSet(viewsets.ModelViewSet):
    queryset = Target.objects.all()
    filterset_fields = ["batch_id"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get("fetchall", None)
        return TargetSerializerAll if fetchall == "yes" else TargetSerializer


class MethodViewSet(viewsets.ModelViewSet):
    queryset = Method.objects.all()
    filterset_fields = ["target_id", "nosteps"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get("fetchall", None)
        return MethodSerializerAll if fetchall == "yes" else MethodSerializer


class ReactionViewSet(viewsets.ModelViewSet):
    queryset = Reaction.objects.all()
    filterset_fields = {"method_id": ["exact"]}

    def get_serializer_class(self):
        fetchall = self.request.GET.get("fetchall", None)
        return ReactionSerializerAll if fetchall == "yes" else ReactionSerializer


class PubChemInfoViewSet(viewsets.ModelViewSet):
    queryset = PubChemInfo.objects.all()
    serializer_class = PubChemInfoSerializer


class ProductViewSet(viewsets.ModelViewSet):
    queryset = Product.objects.all()
    serializer_class = ProductSerializer
    filterset_fields = ["reaction_id"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get("fetchall", None)
        return ProductSerializerAll if fetchall == "yes" else ProductSerializer


class ReactantViewSet(viewsets.ModelViewSet):
    queryset = Reactant.objects.all()
    filterset_fields = ["reaction_id"]

    def get_serializer_class(self):
        fetchall = self.request.GET.get("fetchall", None)
        return ReactantSerializerAll if fetchall == "yes" else ReactantSerializer


class CatalogEntryViewSet(viewsets.ModelViewSet):
    queryset = CatalogEntry.objects.all()
    serializer_class = CatalogEntrySerializer


# Action viewsets
class ActionSessionViewSet(viewsets.ModelViewSet):
    queryset = ActionSession.objects.all()
    serializer_class = ActionSessionSerializer
    filterset_fields = ["reaction_id"]


class AddActionViewSet(viewsets.ModelViewSet):
    queryset = AddAction.objects.all()
    serializer_class = AddActionSerializer
    filterset_fields = ["reaction_id"]


class ExtractActionViewSet(viewsets.ModelViewSet):
    queryset = ExtractAction.objects.all()
    serializer_class = ExtractActionSerializer
    filterset_fields = ["reaction_id"]


class MixActionViewSet(viewsets.ModelViewSet):
    queryset = MixAction.objects.all()
    serializer_class = MixActionSerializer
    filterset_fields = ["reaction_id"]


class StirActionViewSet(viewsets.ModelViewSet):
    queryset = StirAction.objects.all()
    serializer_class = StirActionSerializer
    filterset_fields = ["reaction_id"]


# OT Session viewsets
class OTProjectViewSet(viewsets.ModelViewSet):
    queryset = OTProject.objects.all()
    serializer_class = OTProjectSerializer
    filterset_fields = ["project_id"]

    @action(methods=["post"], detail=False, url_path="createotproject")
    def create_ot_project(self, request, pk=None):
        """Post method to create an OT project.

        Parses the request, extracts batch IDs, protocol name, and any
        custom starting-material CSVs, then delegates to
        ``protocol_service.initiate_ot_project``.
        """
        logger.info("The data is: %s", request.data)
        logger.info("The files are: %s", request.FILES)
        logger.info(
            "The starting  materials request is: %s",
            request.data.get("has_custom_starting_materials"),
        )

        batch_ids = json.loads(request.data["batchids"])
        protocol_name = request.data["protocol_name"]
        has_custom_materials = request.data["has_custom_starting_materials"] == "true"

        logger.info(
            "The OT project has been custom starting materials set to %s",
            has_custom_materials,
        )

        custom_files = None
        if has_custom_materials:
            custom_files = {}
            for batch_id in batch_ids:
                file_key = f"starting_materials_batch_{batch_id}"
                if file_key in request.FILES:
                    custom_files[str(batch_id)] = request.FILES[file_key]

        task_id = initiate_ot_project(
            batch_ids=batch_ids,
            protocol_name=protocol_name,
            custom_files=custom_files,
        )

        return JsonResponse(data={"task_id": task_id})

    @action(detail=False, methods=["get"], url_path="gettaskstatus")
    def get_task_status(self, request, pk=None):
        """Get method for getting the OT project celery task status"""
        task_id = self.request.GET.get("task_id", None)
        if task_id:
            task = AsyncResult(task_id)
            if task.status == "FAILURE":
                data = {"task_status": task.status, "traceback": str(task.traceback)}
                return JsonResponse(data)

            if task.status == "SUCCESS":
                task_summary, otproject_id = task.get()
                data = {
                    "task_status": task.status,
                    "otproject_id": otproject_id,
                    "task_summary": task_summary,
                }
                return JsonResponse(data)

            if task.status == "PENDING":
                data = {"task_status": task.status}
                return JsonResponse(data)


class OTBatchProtocolViewSet(viewsets.ModelViewSet):
    queryset = OTBatchProtocol.objects.all()
    serializer_class = OTBatchProtocolSerializer
    filterset_fields = ["otproject_id", "batch_id", "celery_taskid"]


class OTSessionViewSet(viewsets.ModelViewSet):
    queryset = OTSession.objects.all()
    serializer_class = OTSessionSerializer
    filterset_fields = ["otbatchprotocol_id"]


class DeckViewSet(viewsets.ModelViewSet):
    queryset = Deck.objects.all()
    serializer_class = DeckSerializer
    filterset_fields = ["otsession_id"]


class PipetteViewSet(viewsets.ModelViewSet):
    queryset = Pipette.objects.all()
    serializer_class = PipetteSerializer
    filterset_fields = ["otsession_id"]


class TipRackViewSet(viewsets.ModelViewSet):
    queryset = TipRack.objects.all()
    serializer_class = TipRackSerializer
    filterset_fields = ["otsession_id"]


class PlateViewSet(viewsets.ModelViewSet):
    queryset = Plate.objects.all()
    serializer_class = PlateSerializer
    filterset_fields = ["otbatchprotocol_id"]


class ColumnViewSet(viewsets.ModelViewSet):
    queryset = Column.objects.all()
    serializer_class = ColumnSerializer
    filterset_fields = ["otsession_id"]


class WellViewSet(viewsets.ModelViewSet):
    queryset = Well.objects.all()
    serializer_class = WellSerializer
    filterset_fields = ["otsession_id"]


class CompoundOrderViewSet(viewsets.ModelViewSet):
    queryset = CompoundOrder.objects.all()
    serializer_class = CompoundOrderSerializer
    filterset_fields = ["otsession_id"]


class OTScriptViewSet(viewsets.ModelViewSet):
    queryset = OTScript.objects.all()
    serializer_class = OTScriptSerializer
    filterset_fields = ["otsession_id"]
