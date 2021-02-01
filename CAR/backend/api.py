# Import standard models
from .models import Project, Target, Method, Reaction, Product, AnalyseAction

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


from rest_framework import viewsets, permissions, filters

# Import standard serializers
from .serializers import (
    ProjectSerializer,
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


class ProjectViewSet(viewsets.ModelViewSet):
    queryset = Project.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ProjectSerializer


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
    search_fields = ["=target_id__id"]


class ReactionViewSet(viewsets.ModelViewSet):
    queryset = Reaction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ReactionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=method_id__id"]


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
