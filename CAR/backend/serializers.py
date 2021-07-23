from rest_framework import serializers

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


class ProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = "__all__"


class MculeQuoteSerializer(serializers.ModelSerializer):
    class Meta:
        model = MculeQuote
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


class AnalyseActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = AnalyseAction
        fields = "__all__"


# IBM models here
class IBMAddActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMAddAction
        fields = "__all__"


class IBMCollectLayerActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMCollectLayerAction
        fields = "__all__"


class IBMConcentrateActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMConcentrateAction
        fields = "__all__"


class IBMDegasActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMDegasAction
        fields = "__all__"


class IBMDrySolidActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMDrySolidAction
        fields = "__all__"


class IBMDrySolutionActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMDrySolutionAction
        fields = "__all__"


class IBMExtractActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMExtractAction
        fields = "__all__"


class IBMFilterActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMFilterAction
        fields = "__all__"


class IBMMakeSolutionActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMMakeSolutionAction
        fields = "__all__"


class IBMPartitionActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMPartitionAction
        fields = "__all__"


class IBMpHActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMpHAction
        fields = "__all__"


class IBMPhaseSeparationActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMPhaseSeparationAction
        fields = "__all__"


class IBMQuenchActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMQuenchAction
        fields = "__all__"


class IBMRefluxActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMRefluxAction
        fields = "__all__"


class IBMSetTemperatureActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMSetTemperatureAction
        fields = "__all__"


class IBMStirActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMStirAction
        fields = "__all__"


class IBMStoreActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMStoreAction
        fields = "__all__"


class IBMWaitActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMWaitAction
        fields = "__all__"


class IBMWashActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = IBMWashAction
        fields = "__all__"
