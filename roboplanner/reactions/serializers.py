from rest_framework import serializers
from reactions.models import (
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


class ProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = "__all__"


class TargetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Target
        fields = "__all__"


class MethodSerializer(serializers.ModelSerializer):
    class Meta:
        model = Method
        fields = "__all__"


class ReactionSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reaction
        fields = "__all__"


class ProductSerializer(serializers.ModelSerializer):
    class Meta:
        model = Product
        fields = "__all__"


class AddActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = AddAction
        fields = "__all__"


class MakeSolutionActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = MakeSolutionAction
        fields = "__all__"


class StirActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = StirAction
        fields = "__all__"


class WashActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = WashAction
        fields = "__all__"


class DrySolutionActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = DrySolutionAction
        fields = "__all__"


class ConcentrateActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = ConcentrateAction
        fields = "__all__"


class AnalyseActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = AnalyseAction
        fields = "__all__"

