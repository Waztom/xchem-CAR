from backend.models import (
    Project,
    Target,
    Method,
    Reaction,
    Product,
    AddAction,
    MakeSolutionAction,
    StirAction,
    WashAction,
    DrySolutionAction,
    ConcentrateAction,
    AnalyseAction,
)

from rest_framework import viewsets, permissions, filters
from backend.serializers import (
    ProjectSerializer,
    TargetSerializer,
    MethodSerializer,
    ReactionSerializer,
    ProductSerializer,
    AddActionSerializer,
    MakeSolutionActionSerializer,
    StirActionSerializer,
    WashActionSerializer,
    DrySolutionActionSerializer,
    ConcentrateActionSerializer,
    AnalyseActionSerializer,
)

# Viewsets allow us to create CRUD like usability without explicitly setting
# these up!


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


class AddActionViewSet(viewsets.ModelViewSet):
    queryset = AddAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = AddActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class MakeSolutionActionViewSet(viewsets.ModelViewSet):
    queryset = MakeSolutionAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = MakeSolutionActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class StirActionViewSet(viewsets.ModelViewSet):
    queryset = StirAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = StirActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class WashActionViewSet(viewsets.ModelViewSet):
    queryset = WashAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = WashActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class DrySolutionActionViewSet(viewsets.ModelViewSet):
    queryset = DrySolutionAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = DrySolutionActionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ["=reaction_id__id"]


class ConcentrateActionViewSet(viewsets.ModelViewSet):
    queryset = ConcentrateAction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ConcentrateActionSerializer
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

