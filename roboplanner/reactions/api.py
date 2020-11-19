from reactions.models import Project, Target, Method, Reaction, Reactant
from rest_framework import viewsets, permissions, filters
from reactions.serializers import ProjectSerializer, TargetSerializer, MethodSerializer, ReactionSerializer, ReactantSerializer

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


class MethodViewSet(viewsets.ModelViewSet):
    queryset = Method.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = MethodSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ['target_id__id']


class ReactionViewSet(viewsets.ModelViewSet):
    queryset = Reaction.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ReactionSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ['method_id__id']



class ReactantViewSet(viewsets.ModelViewSet):
    queryset = Reactant.objects.all()
    permission_classes = [
        # Allows all users to access this model - will change when users addded
        permissions.AllowAny
    ]
    serializer_class = ReactantSerializer
    filter_backends = [filters.SearchFilter]
    search_fields = ['reaction_id__id']



