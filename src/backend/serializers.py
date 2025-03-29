from rest_framework import serializers

# Import standard models
from .models import (
    Project,
    Batch,
    Target,
    Method,
    Reaction,
    PubChemInfo,
    Product,
    Reactant,
    CatalogEntry,
)

# Import action models
from .models import (
    ActionSession,
    AddAction,
    ExtractAction,
    MixAction,
    StirAction,
)

# Import OT session models
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


class CatalogEntrySerializer(serializers.ModelSerializer):
    class Meta:
        model = CatalogEntry
        fields = "__all__"


class PubChemInfoSerializer(serializers.ModelSerializer):
    class Meta:
        model = PubChemInfo
        fields = "__all__"


class ReactantSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reactant
        fields = "__all__"


class ReactantSerializerAll(serializers.ModelSerializer):
    catalogentries = CatalogEntrySerializer(many=True, read_only=True)
    reactantpubcheminfo = serializers.SerializerMethodField()

    class Meta:
        model = Reactant
        fields = "__all__"

    def get_reactantpubcheminfo(self, obj):
        if obj.pubcheminfo_id:
            reactantpubcheminfo = PubChemInfo.objects.get(id=obj.pubcheminfo_id.id)
            return PubChemInfoSerializer(reactantpubcheminfo).data


class ProductSerializer(serializers.ModelSerializer):
    class Meta:
        model = Product
        fields = "__all__"


class ProductSerializerAll(serializers.ModelSerializer):
    productpubcheminfo = serializers.SerializerMethodField()

    class Meta:
        model = Product
        fields = "__all__"

    def get_productpubcheminfo(self, obj):
        if obj.pubcheminfo_id:
            productpubcheminfo = PubChemInfo.objects.get(id=obj.pubcheminfo_id.id)
            return PubChemInfoSerializer(productpubcheminfo).data


class ReactionSerializer(serializers.ModelSerializer):
    class Meta:
        model = Reaction
        fields = "__all__"


class ReactionSerializerAll(serializers.ModelSerializer):
    products = ProductSerializerAll(many=True, read_only=True)
    reactants = ReactantSerializerAll(many=True, read_only=True)

    class Meta:
        model = Reaction
        fields = "__all__"


class MethodSerializer(serializers.ModelSerializer):
    class Meta:
        model = Method
        fields = "__all__"


class MethodSerializerAll(serializers.ModelSerializer):
    reactions = ReactionSerializerAll(many=True, read_only=True)

    class Meta:
        model = Method
        fields = "__all__"


class TargetSerializer(serializers.ModelSerializer):
    class Meta:
        model = Target
        fields = "__all__"


class TargetSerializerAll(serializers.ModelSerializer):
    methods = MethodSerializerAll(many=True, read_only=True)
    catalogentries = CatalogEntrySerializer(many=True, read_only=True)

    class Meta:
        model = Target
        fields = "__all__"


class BatchSerializer(serializers.ModelSerializer):
    class Meta:
        model = Batch
        fields = "__all__"


class BatchSerializerAll(serializers.ModelSerializer):
    targets = TargetSerializerAll(many=True, read_only=True)

    class Meta:
        model = Batch
        fields = "__all__"


class ProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = Project
        fields = "__all__"


class ProjectSerializerAll(serializers.ModelSerializer):
    batches = BatchSerializer(many=True, read_only=True)

    class Meta:
        model = Project
        fields = "__all__"


# Action models here
class ActionSessionSerializer(serializers.ModelSerializer):
    class Meta:
        model = ActionSession
        fields = "__all__"


class AddActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = AddAction
        fields = "__all__"


class ExtractActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = ExtractAction
        fields = "__all__"


class MixActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = MixAction
        fields = "__all__"


class StirActionSerializer(serializers.ModelSerializer):
    class Meta:
        model = StirAction
        fields = "__all__"


# OT Session serializers
class DeckSerializer(serializers.ModelSerializer):
    class Meta:
        model = Deck
        fields = "__all__"


class PipetteSerializer(serializers.ModelSerializer):
    class Meta:
        model = Pipette
        fields = "__all__"


class TipRackSerializer(serializers.ModelSerializer):
    class Meta:
        model = TipRack
        fields = "__all__"


class PlateSerializer(serializers.ModelSerializer):
    class Meta:
        model = Plate
        fields = "__all__"


class ColumnSerializer(serializers.ModelSerializer):
    class Meta:
        model = Column
        fields = "__all__"


class WellSerializer(serializers.ModelSerializer):
    class Meta:
        model = Well
        fields = "__all__"


class CompoundOrderSerializer(serializers.ModelSerializer):
    class Meta:
        model = CompoundOrder
        fields = "__all__"


class OTScriptSerializer(serializers.ModelSerializer):
    class Meta:
        model = OTScript
        fields = "__all__"


class OTSessionSerializer(serializers.ModelSerializer):
    otscripts = OTScriptSerializer(many=True, read_only=True)
    compoundorders = CompoundOrderSerializer(many=True, read_only=True)

    class Meta:
        model = OTSession
        fields = "__all__"


class OTBatchProtocolSerializer(serializers.ModelSerializer):
    class Meta:
        model = OTBatchProtocol
        fields = "__all__"


class OTProjectSerializer(serializers.ModelSerializer):
    class Meta:
        model = OTProject
        fields = "__all__"
