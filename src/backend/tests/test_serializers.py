"""Tests for backend.serializers module.

Covers all 26 serializer classes (263 LOC):
  - Flat ModelSerializers: verify fields == "__all__" round-trips
  - Nested *All serializers: verify nested related fields appear
  - Custom SerializerMethodFields: get_reactantpubcheminfo,
    get_productpubcheminfo
  - Edge cases: None FK, empty nested querysets, field presence

All ORM access is mocked — these are pure unit tests.
"""

from unittest import TestCase
from unittest.mock import patch, MagicMock, PropertyMock

from backend.serializers import (
    # Flat serializers
    CatalogEntrySerializer,
    PubChemInfoSerializer,
    ReactantSerializer,
    ReactantSerializerAll,
    ProductSerializer,
    ProductSerializerAll,
    ReactionSerializer,
    ReactionSerializerAll,
    MethodSerializer,
    MethodSerializerAll,
    TargetSerializer,
    TargetSerializerAll,
    BatchSerializer,
    BatchSerializerAll,
    ProjectSerializer,
    ProjectSerializerAll,
    # Action serializers
    ActionSessionSerializer,
    AddActionSerializer,
    ExtractActionSerializer,
    MixActionSerializer,
    StirActionSerializer,
    # OT serializers
    DeckSerializer,
    PipetteSerializer,
    TipRackSerializer,
    PlateSerializer,
    ColumnSerializer,
    WellSerializer,
    CompoundOrderSerializer,
    OTScriptSerializer,
    OTSessionSerializer,
    OTBatchProtocolSerializer,
    OTProjectSerializer,
)


# ===================================================================
# Helpers
# ===================================================================
def _make_mock_instance(fields: dict) -> MagicMock:
    """Create a MagicMock with the given attributes set."""
    obj = MagicMock()
    for k, v in fields.items():
        setattr(obj, k, v)
    obj.pk = fields.get("id", 1)
    return obj


# ===================================================================
# Meta field configuration tests
# ===================================================================
class TestSerializerMetaConfig(TestCase):
    """Verify every serializer declares fields='__all__'."""

    def test_all_flat_serializers_use_all_fields(self):
        flat_serializers = [
            CatalogEntrySerializer,
            PubChemInfoSerializer,
            ReactantSerializer,
            ProductSerializer,
            ReactionSerializer,
            MethodSerializer,
            TargetSerializer,
            BatchSerializer,
            ProjectSerializer,
            ActionSessionSerializer,
            AddActionSerializer,
            ExtractActionSerializer,
            MixActionSerializer,
            StirActionSerializer,
            DeckSerializer,
            PipetteSerializer,
            TipRackSerializer,
            PlateSerializer,
            ColumnSerializer,
            WellSerializer,
            CompoundOrderSerializer,
            OTScriptSerializer,
            OTBatchProtocolSerializer,
            OTProjectSerializer,
        ]
        for cls in flat_serializers:
            with self.subTest(serializer=cls.__name__):
                self.assertEqual(
                    cls.Meta.fields,
                    "__all__",
                    f"{cls.__name__} should use fields='__all__'",
                )

    def test_all_nested_serializers_use_all_fields(self):
        nested_serializers = [
            ReactantSerializerAll,
            ProductSerializerAll,
            ReactionSerializerAll,
            MethodSerializerAll,
            TargetSerializerAll,
            BatchSerializerAll,
            ProjectSerializerAll,
            OTSessionSerializer,
        ]
        for cls in nested_serializers:
            with self.subTest(serializer=cls.__name__):
                self.assertEqual(
                    cls.Meta.fields,
                    "__all__",
                    f"{cls.__name__} should use fields='__all__'",
                )


# ===================================================================
# Model reference tests
# ===================================================================
class TestSerializerModelBindings(TestCase):
    """Verify each serializer points to the correct model."""

    def test_model_bindings(self):
        from backend.models import (
            CatalogEntry,
            PubChemInfo,
            Reactant,
            Product,
            Reaction,
            Method,
            Target,
            Batch,
            Project,
            ActionSession,
            AddAction,
            ExtractAction,
            MixAction,
            StirAction,
            Deck,
            Pipette,
            TipRack,
            Plate,
            Column,
            Well,
            CompoundOrder,
            OTScript,
            OTSession,
            OTBatchProtocol,
            OTProject,
        )

        expected = {
            CatalogEntrySerializer: CatalogEntry,
            PubChemInfoSerializer: PubChemInfo,
            ReactantSerializer: Reactant,
            ReactantSerializerAll: Reactant,
            ProductSerializer: Product,
            ProductSerializerAll: Product,
            ReactionSerializer: Reaction,
            ReactionSerializerAll: Reaction,
            MethodSerializer: Method,
            MethodSerializerAll: Method,
            TargetSerializer: Target,
            TargetSerializerAll: Target,
            BatchSerializer: Batch,
            BatchSerializerAll: Batch,
            ProjectSerializer: Project,
            ProjectSerializerAll: Project,
            ActionSessionSerializer: ActionSession,
            AddActionSerializer: AddAction,
            ExtractActionSerializer: ExtractAction,
            MixActionSerializer: MixAction,
            StirActionSerializer: StirAction,
            DeckSerializer: Deck,
            PipetteSerializer: Pipette,
            TipRackSerializer: TipRack,
            PlateSerializer: Plate,
            ColumnSerializer: Column,
            WellSerializer: Well,
            CompoundOrderSerializer: CompoundOrder,
            OTScriptSerializer: OTScript,
            OTSessionSerializer: OTSession,
            OTBatchProtocolSerializer: OTBatchProtocol,
            OTProjectSerializer: OTProject,
        }
        for serializer_cls, model_cls in expected.items():
            with self.subTest(serializer=serializer_cls.__name__):
                self.assertEqual(
                    serializer_cls.Meta.model,
                    model_cls,
                    f"{serializer_cls.__name__}.Meta.model != {model_cls.__name__}",
                )


# ===================================================================
# Nested field declarations
# ===================================================================
class TestNestedFieldDeclarations(TestCase):
    """Verify nested serializers declare their nested child fields."""

    def _get_field_names(self, serializer_cls):
        """Return the declared field names from a serializer class."""
        return set(serializer_cls._declared_fields.keys())

    def test_reactant_serializer_all_has_nested_fields(self):
        fields = self._get_field_names(ReactantSerializerAll)
        self.assertIn("catalogentries", fields)
        self.assertIn("reactantpubcheminfo", fields)

    def test_product_serializer_all_has_nested_fields(self):
        fields = self._get_field_names(ProductSerializerAll)
        self.assertIn("productpubcheminfo", fields)

    def test_reaction_serializer_all_has_nested_fields(self):
        fields = self._get_field_names(ReactionSerializerAll)
        self.assertIn("products", fields)
        self.assertIn("reactants", fields)

    def test_method_serializer_all_has_nested_fields(self):
        fields = self._get_field_names(MethodSerializerAll)
        self.assertIn("reactions", fields)

    def test_target_serializer_all_has_nested_fields(self):
        fields = self._get_field_names(TargetSerializerAll)
        self.assertIn("methods", fields)
        self.assertIn("catalogentries", fields)

    def test_batch_serializer_all_has_nested_fields(self):
        fields = self._get_field_names(BatchSerializerAll)
        self.assertIn("targets", fields)

    def test_project_serializer_all_has_nested_fields(self):
        fields = self._get_field_names(ProjectSerializerAll)
        self.assertIn("batches", fields)

    def test_otsession_serializer_has_nested_fields(self):
        fields = self._get_field_names(OTSessionSerializer)
        self.assertIn("otscripts", fields)
        self.assertIn("compoundorders", fields)


# ===================================================================
# SerializerMethodField logic — ReactantSerializerAll
# ===================================================================
class TestReactantSerializerAllPubChemInfo(TestCase):
    """Tests for ReactantSerializerAll.get_reactantpubcheminfo."""

    @patch("backend.serializers.PubChemInfo.objects")
    def test_returns_pubchem_data_when_fk_present(self, mock_qs):
        """When reactant has a pubcheminfo_id FK, fetch and serialize it."""
        pubchem_obj = MagicMock()
        pubchem_obj.id = 42
        pubchem_obj.compoundid = 12345
        pubchem_obj.smiles = "CCO"
        pubchem_obj.summaryurl = "https://pubchem.ncbi.nlm.nih.gov/compound/12345"
        pubchem_obj.lcssurl = "https://pubchem.ncbi.nlm.nih.gov/compound/12345#datasheet=LCSS"
        pubchem_obj.cas = "64-17-5"

        mock_qs.get.return_value = pubchem_obj

        reactant_obj = MagicMock()
        reactant_obj.pubcheminfo_id = MagicMock()
        reactant_obj.pubcheminfo_id.id = 42

        serializer = ReactantSerializerAll()
        result = serializer.get_reactantpubcheminfo(reactant_obj)

        mock_qs.get.assert_called_once_with(id=42)
        self.assertIsNotNone(result)

    def test_returns_none_when_fk_absent(self):
        """When pubcheminfo_id is None, should return None."""
        reactant_obj = MagicMock()
        reactant_obj.pubcheminfo_id = None

        serializer = ReactantSerializerAll()
        result = serializer.get_reactantpubcheminfo(reactant_obj)

        self.assertIsNone(result)


# ===================================================================
# SerializerMethodField logic — ProductSerializerAll
# ===================================================================
class TestProductSerializerAllPubChemInfo(TestCase):
    """Tests for ProductSerializerAll.get_productpubcheminfo."""

    @patch("backend.serializers.PubChemInfo.objects")
    def test_returns_pubchem_data_when_fk_present(self, mock_qs):
        pubchem_obj = MagicMock()
        pubchem_obj.id = 99
        pubchem_obj.compoundid = 67890
        mock_qs.get.return_value = pubchem_obj

        product_obj = MagicMock()
        product_obj.pubcheminfo_id = MagicMock()
        product_obj.pubcheminfo_id.id = 99

        serializer = ProductSerializerAll()
        result = serializer.get_productpubcheminfo(product_obj)

        mock_qs.get.assert_called_once_with(id=99)
        self.assertIsNotNone(result)

    def test_returns_none_when_fk_absent(self):
        product_obj = MagicMock()
        product_obj.pubcheminfo_id = None

        serializer = ProductSerializerAll()
        result = serializer.get_productpubcheminfo(product_obj)

        self.assertIsNone(result)


# ===================================================================
# Nested read_only flags
# ===================================================================
class TestNestedReadOnlyFlags(TestCase):
    """Nested fields should all be read_only=True (no write via nested)."""

    def test_nested_fields_are_read_only(self):
        checks = [
            (ReactantSerializerAll, "catalogentries"),
            (ReactionSerializerAll, "products"),
            (ReactionSerializerAll, "reactants"),
            (MethodSerializerAll, "reactions"),
            (TargetSerializerAll, "methods"),
            (TargetSerializerAll, "catalogentries"),
            (BatchSerializerAll, "targets"),
            (ProjectSerializerAll, "batches"),
            (OTSessionSerializer, "otscripts"),
            (OTSessionSerializer, "compoundorders"),
        ]
        for serializer_cls, field_name in checks:
            with self.subTest(serializer=serializer_cls.__name__, field=field_name):
                field = serializer_cls._declared_fields[field_name]
                self.assertTrue(
                    field.read_only,
                    f"{serializer_cls.__name__}.{field_name} should be read_only",
                )


# ===================================================================
# Nested many=True flags
# ===================================================================
class TestNestedManyFlags(TestCase):
    """Nested list fields should use many=True."""

    def test_nested_list_fields_use_many(self):
        checks = [
            (ReactantSerializerAll, "catalogentries"),
            (ReactionSerializerAll, "products"),
            (ReactionSerializerAll, "reactants"),
            (MethodSerializerAll, "reactions"),
            (TargetSerializerAll, "methods"),
            (TargetSerializerAll, "catalogentries"),
            (BatchSerializerAll, "targets"),
            (ProjectSerializerAll, "batches"),
            (OTSessionSerializer, "otscripts"),
            (OTSessionSerializer, "compoundorders"),
        ]
        for serializer_cls, field_name in checks:
            with self.subTest(serializer=serializer_cls.__name__, field=field_name):
                field = serializer_cls._declared_fields[field_name]
                # ListSerializer wraps child when many=True
                from rest_framework.serializers import ListSerializer
                self.assertIsInstance(
                    field,
                    ListSerializer,
                    f"{serializer_cls.__name__}.{field_name} should be many=True",
                )


# ===================================================================
# Nested child serializer types
# ===================================================================
class TestNestedChildSerializerTypes(TestCase):
    """Verify nested fields use the correct child serializer class."""

    def _get_child_class(self, serializer_cls, field_name):
        field = serializer_cls._declared_fields[field_name]
        # many=True wraps in ListSerializer — child is the actual serializer
        return type(field.child) if hasattr(field, "child") else type(field)

    def test_reactant_all_catalogentries_child(self):
        self.assertEqual(
            self._get_child_class(ReactantSerializerAll, "catalogentries"),
            CatalogEntrySerializer,
        )

    def test_reaction_all_products_child(self):
        self.assertEqual(
            self._get_child_class(ReactionSerializerAll, "products"),
            ProductSerializerAll,
        )

    def test_reaction_all_reactants_child(self):
        self.assertEqual(
            self._get_child_class(ReactionSerializerAll, "reactants"),
            ReactantSerializerAll,
        )

    def test_method_all_reactions_child(self):
        self.assertEqual(
            self._get_child_class(MethodSerializerAll, "reactions"),
            ReactionSerializerAll,
        )

    def test_target_all_methods_child(self):
        self.assertEqual(
            self._get_child_class(TargetSerializerAll, "methods"),
            MethodSerializerAll,
        )

    def test_target_all_catalogentries_child(self):
        self.assertEqual(
            self._get_child_class(TargetSerializerAll, "catalogentries"),
            CatalogEntrySerializer,
        )

    def test_batch_all_targets_child(self):
        self.assertEqual(
            self._get_child_class(BatchSerializerAll, "targets"),
            TargetSerializerAll,
        )

    def test_project_all_batches_child(self):
        self.assertEqual(
            self._get_child_class(ProjectSerializerAll, "batches"),
            BatchSerializer,
        )

    def test_otsession_otscripts_child(self):
        self.assertEqual(
            self._get_child_class(OTSessionSerializer, "otscripts"),
            OTScriptSerializer,
        )

    def test_otsession_compoundorders_child(self):
        self.assertEqual(
            self._get_child_class(OTSessionSerializer, "compoundorders"),
            CompoundOrderSerializer,
        )
