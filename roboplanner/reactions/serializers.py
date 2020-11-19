from rest_framework import serializers
from reactions.models import Project, Target, Method, Reaction, Reactant

class ProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = '__all__'

class TargetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Target
        fields = '__all__'

class MethodSerializer(serializers.ModelSerializer):
    class Meta:
        model = Method
        fields = '__all__'

class ReactionSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reaction
        fields = '__all__'

class ReactantSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reactant
        fields = '__all__'
