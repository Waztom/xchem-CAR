# Import standard models
from .models import Project, MculeQuote, Target, Method, Reaction, Product, AnalyseAction

# Import IBM models
from .models import (
    IBMAddAction,
    IBMCollectLayerAction,
    IBMConcentrateAction,
    IBMDegasAction,
    IBMDrySolidAction,
    IBMDrySolutionAction,
    IBMExtractAction,
    IBMFilterAction,
    IBMMakeSolutionAction,
    IBMPartitionAction,
    IBMpHAction,
    IBMPhaseSeparationAction,
    IBMQuenchAction,
    IBMRefluxAction,
    IBMSetTemperatureAction,
    IBMStirAction,
    IBMStoreAction,
    IBMWaitAction,
    IBMWashAction,
)

# Import OT Session models
from .models import OTSession, Deck, Pipette, TipRack, Plate, Well, CompoundOrder, OTScript

from rest_framework import viewsets, permissions, filters
from rest_framework.decorators import action
from django.http import JsonResponse
from django.core.files.base import ContentFile
from django.db.models import Q

# Import standard serializers
from .serializers import (
    ProjectSerializer,
    MculeQuoteSerializer,
    TargetSerializer,
    MethodSerializer,
    ReactionSerializer,
    ProductSerializer,
    AnalyseActionSerializer,
)

# Import IBM serializers
from .serializers import (
    IBMAddActionSerializer,
    IBMCollectLayerActionSerializer,
    IBMConcentrateActionSerializer,
    IBMDegasActionSerializer,
    IBMDrySolidActionSerializer,
    IBMDrySolutionActionSerializer,
    IBMExtractActionSerializer,
    IBMFilterActionSerializer,
    IBMMakeSolutionActionSerializer,
    IBMPartitionActionSerializer,
    IBMpHActionSerializer,
    IBMPhaseSeparationActionSerializer,
    IBMQuenchActionSerializer,
    IBMRefluxActionSerializer,
    IBMSetTemperatureActionSerializer,
    IBMStirActionSerializer,
    IBMStoreActionSerializer,
    IBMWaitActionSerializer,
    IBMWashActionSerializer,
)

# Import OT Session serializers
from .serializers import (
    OTSessionSerializer,
    DeckSerializer,
    PipetteSerializer,
    TipRackSerializer,
    PlateSerializer,
    WellSerializer,
    CompoundOrderSerializer,
    OTScriptSerializer,
)

from rdkit import Chem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage

from .utils import createSVGString


class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ProjectSerializer


class MculeQuoteViewSet(viewsets.ModelViewSet):
    queryset = MculeQuote.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = MculeQuoteSerializer


class TargetViewSet(viewsets.ModelViewSet):
    queryset = Target.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = TargetSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=project_id__id"]


class MethodViewSet(viewsets.ModelViewSet):
    queryset = Method.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = MethodSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=target_id__id", "nosteps"]


class GroupByStepsViewSet(viewsets.ModelViewSet):
    """ "
    Viewset to filter methods that have a given no of steps in a particular project
    """

    serializer_class = MethodSerializer

    def get_queryset(self):
        projectid = self.request.GET.get("projectid")
        nosteps = self.request.GET.get("nosteps")
        targetqueryset = Target.objects.filter(project_id=projectid).order_by("id")
        targetids = [target.id for target in targetqueryset]
        methodsqueryset = Method.objects.filter(Q(target_id__in=targetids) & Q(nosteps=nosteps))
        return methodsqueryset


class ReactionViewSet(viewsets.ModelViewSet):
    queryset = Reaction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ReactionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=method_id__id"]
    filter_fields = {"successrate": ["gte", "lte"]}


class ProductViewSet(viewsets.ModelViewSet):
    queryset = Product.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ProductSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class AnalyseActionViewSet(viewsets.ModelViewSet):
    queryset = AnalyseAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = AnalyseActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


# IBM viewsets here
class IBMAddActionViewSet(viewsets.ModelViewSet):
    queryset = IBMAddAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMAddActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]

    def get_patch_object(self, pk):
        return IBMAddAction.objects.get(pk=pk)

    # @action(methods=["PATCH"], detail=True)
    def partial_update(self, request, pk):
        addaction = self.get_patch_object(pk)

        if "materialsmiles" in request.data:
            materialsmiles = request.data["materialsmiles"]
            mol = Chem.MolFromSmiles(materialsmiles)
            molecular_weight = Descriptors.ExactMolWt(mol)
            add_svg_string = createSVGString(materialsmiles)
            # Delete previous image
            svg_fp = addaction.materialimage.path
            default_storage.delete(svg_fp)
            # Add new image
            add_svg_fn = default_storage.save(
                "addactionimages/{}.svg".format(materialsmiles),
                ContentFile(add_svg_string),
            )
            addaction.materialsmiles = materialsmiles
            addaction.molecularweight = molecular_weight
            addaction.materialimage = add_svg_fn
            addaction.save()
            serialized_data = IBMAddActionSerializer(addaction).data
            if serialized_data:
                return JsonResponse(data=serialized_data)
            else:
                return JsonResponse(data="wrong parameters")
        else:
            serializer = IBMAddActionSerializer(addaction, data=request.data, partial=True)
            if serializer.is_valid():
                serializer.save()
                return JsonResponse(data=serializer.data)
            else:
                return JsonResponse(data="wrong parameters")


class IBMCollectLayerActionViewSet(viewsets.ModelViewSet):
    queryset = IBMCollectLayerAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMCollectLayerActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMConcentrateActionViewSet(viewsets.ModelViewSet):
    queryset = IBMConcentrateAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMConcentrateActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMDegasActionViewSet(viewsets.ModelViewSet):
    queryset = IBMDegasAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMDegasActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMDrySolidActionViewSet(viewsets.ModelViewSet):
    queryset = IBMDrySolidAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMDrySolidActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMDrySolutionActionViewSet(viewsets.ModelViewSet):
    queryset = IBMDrySolutionAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMDrySolutionActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMExtractActionViewSet(viewsets.ModelViewSet):
    queryset = IBMExtractAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMExtractActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMFilterActionViewSet(viewsets.ModelViewSet):
    queryset = IBMFilterAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMFilterActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMMakeSolutionActionViewSet(viewsets.ModelViewSet):
    queryset = IBMMakeSolutionAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMMakeSolutionActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMPartitionActionViewSet(viewsets.ModelViewSet):
    queryset = IBMPartitionAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMPartitionActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMpHActionViewSet(viewsets.ModelViewSet):
    queryset = IBMpHAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMpHActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMPhaseSeparationActionViewSet(viewsets.ModelViewSet):
    queryset = IBMPhaseSeparationAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMPhaseSeparationActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMQuenchActionViewSet(viewsets.ModelViewSet):
    queryset = IBMQuenchAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMQuenchActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMRefluxActionViewSet(viewsets.ModelViewSet):
    queryset = IBMRefluxAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMRefluxActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMSetTemperatureActionViewSet(viewsets.ModelViewSet):
    queryset = IBMSetTemperatureAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMSetTemperatureActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMStirActionViewSet(viewsets.ModelViewSet):
    queryset = IBMStirAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMStirActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMStoreActionViewSet(viewsets.ModelViewSet):
    queryset = IBMStoreAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMStoreActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMWaitActionViewSet(viewsets.ModelViewSet):
    queryset = IBMWaitAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMWaitActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class IBMWashActionViewSet(viewsets.ModelViewSet):
    queryset = IBMWashAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = IBMWashActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


# OT Session viewsets
class OTSessionViewSet(viewsets.ModelViewSet):
    queryset = OTSession.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = OTSessionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=project_id__id"]


class DeckViewSet(viewsets.ModelViewSet):
    queryset = Deck.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = DeckSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["= otsession_id__id"]


class PipetteViewSet(viewsets.ModelViewSet):
    queryset = Pipette.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = PipetteSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=otsession_id__id"]


class TipRackViewSet(viewsets.ModelViewSet):
    queryset = TipRack.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = TipRackSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=otsession_id__id"]


class PlateViewSet(viewsets.ModelViewSet):
    queryset = Plate.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = PlateSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=otsession_id__id"]


class WellViewSet(viewsets.ModelViewSet):
    queryset = Well.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = WellSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=otsession_id__id"]


class CompoundOrderViewSet(viewsets.ModelViewSet):
    queryset = CompoundOrder.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = CompoundOrderSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=project_id__id"]


class OTScriptViewSet(viewsets.ModelViewSet):
    queryset = OTScript.objects.all()
    permission_classes = [permissions.AllowAny]
    serializer_class = OTScriptSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=otsession_id__id"]
