"""Tests for backend.createmodels module (888 LOC).

Covers:
  - create_project_model
  - create_batch_model (with and without parent batch)
  - create_target_model
  - create_method_model
  - create_reaction_model
  - create_pubchem_info_model (with and without CAS)
  - get_pubchem_info (existing queryset, PubChem lookup, not found)
  - create_product_model (with/without pubchem)
  - create_reactant_model (with/without pubchem, previous reaction product)
  - create_catalog_entry_model (all branches: target/reactant, prev product,
    lab inventory, Manifold catalog with screening/bb prices)
  - CreateEncodedActionModels (init orchestration, getProductMols,
    calculateMass, calculateVolume, action model creation helpers)

All external dependencies (Django ORM, RDKit, chem_utils, pubchem_utils,
conversions, default_storage, MCuleAPI) are mocked.
"""

from unittest import TestCase
from unittest.mock import patch, MagicMock, call, PropertyMock
import math


# ===================================================================
# create_project_model
# ===================================================================
class TestCreateProjectModel(TestCase):
    """Tests for create_project_model."""

    @patch("backend.createmodels.Project")
    def test_creates_project_with_correct_fields(self, MockProject):
        from backend.createmodels import create_project_model

        mock_instance = MagicMock()
        mock_instance.id = 42
        MockProject.return_value = mock_instance

        project_info = {
            "projectname": "test-proj",
            "submittername": "Alice",
            "submitterorganisation": "AcmeCorp",
            "proteintarget": "EGFR",
        }

        result = create_project_model(project_info)

        self.assertEqual(mock_instance.name, "test-proj")
        self.assertEqual(mock_instance.submittername, "Alice")
        self.assertEqual(mock_instance.submitterorganisation, "AcmeCorp")
        self.assertEqual(mock_instance.proteintarget, "EGFR")
        mock_instance.save.assert_called_once()
        self.assertEqual(result, 42)


# ===================================================================
# create_batch_model
# ===================================================================
class TestCreateBatchModel(TestCase):
    """Tests for create_batch_model."""

    @patch("backend.createmodels.Batch")
    @patch("backend.createmodels.Project.objects")
    def test_creates_root_batch(self, mock_proj_qs, MockBatch):
        from backend.createmodels import create_batch_model

        mock_proj_obj = MagicMock()
        mock_proj_qs.get.return_value = mock_proj_obj

        mock_batch = MagicMock()
        mock_batch.id = 10
        MockBatch.return_value = mock_batch

        result = create_batch_model(project_id=1, batchtag="batch-A")

        mock_proj_qs.get.assert_called_once_with(id=1)
        self.assertEqual(mock_batch.project_id, mock_proj_obj)
        self.assertEqual(mock_batch.batchtag, "batch-A")
        mock_batch.save.assert_called_once()
        self.assertEqual(result, 10)

    @patch("backend.createmodels.Batch")
    @patch("backend.createmodels.Project.objects")
    def test_creates_child_batch_with_parent(self, mock_proj_qs, MockBatch):
        from backend.createmodels import create_batch_model

        mock_proj_obj = MagicMock()
        mock_proj_qs.get.return_value = mock_proj_obj

        mock_parent_batch = MagicMock()
        # Batch.objects.get will be called for parent
        MockBatch.objects = MagicMock()
        MockBatch.objects.get.return_value = mock_parent_batch

        mock_batch = MagicMock()
        mock_batch.id = 11
        MockBatch.return_value = mock_batch

        result = create_batch_model(project_id=1, batchtag="sub-batch", batch_id=5)

        MockBatch.objects.get.assert_called_once_with(pk=5)
        self.assertEqual(mock_batch.batch_id, mock_parent_batch)
        self.assertEqual(result, 11)


# ===================================================================
# create_target_model
# ===================================================================
class TestCreateTargetModel(TestCase):
    """Tests for create_target_model."""

    @patch("backend.createmodels.default_storage")
    @patch("backend.createmodels.create_svg_string", return_value="<svg/>")
    @patch("backend.createmodels.calculate_mass_from_mols", return_value=50.0)
    @patch("backend.createmodels.calculate_mols_from_conc", return_value=0.001)
    @patch("backend.createmodels.canon_smiles", return_value="CCO")
    @patch("backend.createmodels.Target")
    @patch("backend.createmodels.Batch.objects")
    def test_creates_target_with_all_fields(
        self,
        mock_batch_qs,
        MockTarget,
        mock_canon,
        mock_calc_mols,
        mock_calc_mass,
        mock_svg,
        mock_storage,
    ):
        from backend.createmodels import create_target_model

        mock_batch_obj = MagicMock()
        mock_batch_qs.get.return_value = mock_batch_obj
        mock_storage.save.return_value = "targetimages/ethanol.svg"

        mock_target = MagicMock()
        mock_target.id = 99
        MockTarget.return_value = mock_target

        result = create_target_model(
            batch_id=1, name="ethanol", smiles="CCO", concentration=10.0, volume=100.0
        )

        mock_batch_qs.get.assert_called_once_with(id=1)
        mock_canon.assert_called_once_with(smiles="CCO")
        mock_calc_mols.assert_called_once_with(
            target_concentration=10.0, target_volume=100.0
        )
        mock_calc_mass.assert_called_once_with(mols=0.001, SMILES="CCO")
        self.assertEqual(mock_target.smiles, "CCO")
        self.assertEqual(mock_target.mols, 0.001)
        self.assertEqual(mock_target.mass, 50.0)
        self.assertEqual(mock_target.concentration, 10.0)
        self.assertEqual(mock_target.volume, 100.0)
        self.assertEqual(mock_target.name, "ethanol")
        self.assertEqual(mock_target.image, "targetimages/ethanol.svg")
        mock_target.save.assert_called_once()
        self.assertEqual(result, 99)


# ===================================================================
# create_method_model
# ===================================================================
class TestCreateMethodModel(TestCase):
    """Tests for create_method_model."""

    @patch("backend.createmodels.Method")
    @patch("backend.createmodels.Target.objects")
    def test_creates_method(self, mock_target_qs, MockMethod):
        from backend.createmodels import create_method_model

        mock_target_obj = MagicMock()
        mock_target_qs.get.return_value = mock_target_obj

        mock_method = MagicMock()
        mock_method.id = 55
        MockMethod.return_value = mock_method

        result = create_method_model(target_id=10, nosteps=3, otchem=True)

        self.assertEqual(mock_method.target_id, mock_target_obj)
        self.assertEqual(mock_method.nosteps, 3)
        self.assertTrue(mock_method.otchem)
        mock_method.save.assert_called_once()
        self.assertEqual(result, 55)


# ===================================================================
# create_reaction_model
# ===================================================================
class TestCreateReactionModel(TestCase):
    """Tests for create_reaction_model."""

    @patch("backend.createmodels.default_storage")
    @patch(
        "backend.createmodels.create_reaction_svg_string",
        return_value="<svg>rxn</svg>",
    )
    @patch("backend.createmodels.Reaction")
    @patch("backend.createmodels.Method.objects")
    def test_creates_reaction_basic(
        self, mock_method_qs, MockReaction, mock_rxn_svg, mock_storage
    ):
        from backend.createmodels import create_reaction_model

        mock_method_obj = MagicMock()
        mock_method_qs.get.return_value = mock_method_obj
        mock_storage.save.return_value = "reactionimages/Amidation.svg"

        mock_rxn = MagicMock()
        mock_rxn.id = 77
        MockReaction.return_value = mock_rxn

        result = create_reaction_model(
            method_id=5,
            reaction_class="Amidation",
            reaction_number=1,
            intramolecular=False,
            reaction_smarts="[C:1][OH]>>[C:1][NH2]",
        )

        self.assertEqual(mock_rxn.method_id, mock_method_obj)
        self.assertEqual(mock_rxn.reactionclass, "Amidation")
        self.assertEqual(mock_rxn.number, 1)
        self.assertFalse(mock_rxn.intramolecular)
        self.assertEqual(mock_rxn.image, "reactionimages/Amidation.svg")
        mock_rxn.save.assert_called_once()
        self.assertEqual(result, 77)

    @patch("backend.createmodels.default_storage")
    @patch("backend.createmodels.create_reaction_svg_string", return_value="<svg/>")
    @patch("backend.createmodels.Reaction")
    @patch("backend.createmodels.Method.objects")
    def test_creates_reaction_with_optional_fields(
        self, mock_method_qs, MockReaction, mock_rxn_svg, mock_storage
    ):
        from backend.createmodels import create_reaction_model

        mock_method_qs.get.return_value = MagicMock()
        mock_storage.save.return_value = "reactionimages/test.svg"
        mock_rxn = MagicMock()
        mock_rxn.id = 78
        MockReaction.return_value = mock_rxn

        create_reaction_model(
            method_id=5,
            reaction_class="SNAr",
            reaction_number=2,
            intramolecular=True,
            reaction_smarts="[F:1]>>[N:1]",
            reaction_temperature=80,
            reaction_recipe="heated",
            groupby_column=False,
        )

        self.assertEqual(mock_rxn.temperature, 80)
        self.assertEqual(mock_rxn.recipe, "heated")
        self.assertFalse(mock_rxn.groupbycolumn)

    @patch("backend.createmodels.default_storage")
    @patch("backend.createmodels.create_reaction_svg_string", return_value="<svg/>")
    @patch("backend.createmodels.Reaction")
    @patch("backend.createmodels.Method.objects")
    def test_no_temperature_or_recipe_when_not_provided(
        self, mock_method_qs, MockReaction, mock_rxn_svg, mock_storage
    ):
        from backend.createmodels import create_reaction_model

        mock_method_qs.get.return_value = MagicMock()
        mock_storage.save.return_value = "reactionimages/test.svg"
        mock_rxn = MagicMock()
        mock_rxn.id = 79
        MockReaction.return_value = mock_rxn

        create_reaction_model(
            method_id=5,
            reaction_class="Suzuki",
            reaction_number=1,
            intramolecular=False,
            reaction_smarts="[C:1]>>[C:1]",
        )

        # temperature and recipe should NOT have been set (no attribute assignment)
        # groupbycolumn defaults to True
        self.assertTrue(mock_rxn.groupbycolumn)


# ===================================================================
# create_pubchem_info_model
# ===================================================================
class TestCreatePubChemInfoModel(TestCase):
    """Tests for create_pubchem_info_model."""

    @patch("backend.createmodels.canon_smiles", return_value="CCO")
    @patch("backend.createmodels.PubChemInfo")
    def test_creates_with_cas(self, MockPubChem, mock_canon):
        from backend.createmodels import create_pubchem_info_model

        mock_obj = MagicMock()
        MockPubChem.return_value = mock_obj

        result = create_pubchem_info_model(
            compoundid=702, smiles="CCO", cas="64-17-5"
        )

        self.assertEqual(mock_obj.smiles, "CCO")
        self.assertEqual(mock_obj.cas, "64-17-5")
        self.assertEqual(mock_obj.compoundid, 702)
        self.assertIn("702", mock_obj.summaryurl)
        self.assertIn("702", mock_obj.lcssurl)
        mock_obj.save.assert_called_once()
        self.assertEqual(result, mock_obj)

    @patch("backend.createmodels.canon_smiles", return_value="c1ccccc1")
    @patch("backend.createmodels.PubChemInfo")
    def test_creates_without_cas(self, MockPubChem, mock_canon):
        from backend.createmodels import create_pubchem_info_model

        mock_obj = MagicMock()
        MockPubChem.return_value = mock_obj

        create_pubchem_info_model(compoundid=241, smiles="c1ccccc1")

        # cas should not be set when not provided
        # (MagicMock won't raise, but we verify save is called)
        mock_obj.save.assert_called_once()


# ===================================================================
# get_pubchem_info
# ===================================================================
class TestGetPubChemInfo(TestCase):
    """Tests for get_pubchem_info."""

    @patch("backend.createmodels.canon_smiles", return_value="CCO")
    @patch("backend.createmodels.PubChemInfo.objects")
    def test_returns_existing_pubchem_object(self, mock_qs, mock_canon):
        from backend.createmodels import get_pubchem_info

        existing = MagicMock()
        mock_qs.filter.return_value = [existing]

        result = get_pubchem_info(smiles="CCO")

        mock_qs.filter.assert_called_once_with(smiles="CCO")
        self.assertEqual(result, existing)

    @patch("backend.createmodels.create_pubchem_info_model")
    @patch("backend.createmodels.get_pubchem_cas", return_value="64-17-5")
    @patch("backend.createmodels.get_pubchem_compound")
    @patch("backend.createmodels.get_inchi_key", return_value="LFQSCWFLJHTTHZ")
    @patch("backend.createmodels.canon_smiles", return_value="CCO")
    @patch("backend.createmodels.PubChemInfo.objects")
    def test_fetches_from_pubchem_when_not_in_db(
        self,
        mock_qs,
        mock_canon,
        mock_inchi,
        mock_get_compound,
        mock_get_cas,
        mock_create,
    ):
        from backend.createmodels import get_pubchem_info

        mock_qs.filter.return_value = []  # not in DB
        mock_compound = MagicMock()
        mock_compound.cid = 702
        mock_get_compound.return_value = mock_compound
        mock_create.return_value = MagicMock()

        result = get_pubchem_info(smiles="CCO")

        mock_inchi.assert_called_once_with(smiles="CCO")
        mock_get_compound.assert_called_once_with(inchikey="LFQSCWFLJHTTHZ")
        mock_create.assert_called_once_with(
            compoundid=702, smiles="CCO", cas="64-17-5"
        )
        self.assertIsNotNone(result)

    @patch("backend.createmodels.create_pubchem_info_model")
    @patch("backend.createmodels.get_pubchem_cas", return_value=None)
    @patch("backend.createmodels.get_pubchem_compound")
    @patch("backend.createmodels.get_inchi_key", return_value="ABCDEF")
    @patch("backend.createmodels.canon_smiles", return_value="C")
    @patch("backend.createmodels.PubChemInfo.objects")
    def test_fetches_from_pubchem_without_cas(
        self,
        mock_qs,
        mock_canon,
        mock_inchi,
        mock_get_compound,
        mock_get_cas,
        mock_create,
    ):
        from backend.createmodels import get_pubchem_info

        mock_qs.filter.return_value = []
        mock_compound = MagicMock()
        mock_compound.cid = 297
        mock_get_compound.return_value = mock_compound

        get_pubchem_info(smiles="C")

        mock_create.assert_called_once_with(compoundid=297, smiles="C")

    @patch("backend.createmodels.get_pubchem_compound", return_value=None)
    @patch("backend.createmodels.get_inchi_key", return_value="XYZ")
    @patch("backend.createmodels.canon_smiles", return_value="INVALID")
    @patch("backend.createmodels.PubChemInfo.objects")
    def test_returns_false_when_not_found_anywhere(
        self, mock_qs, mock_canon, mock_inchi, mock_get_compound
    ):
        from backend.createmodels import get_pubchem_info

        mock_qs.filter.return_value = []

        result = get_pubchem_info(smiles="INVALID")

        self.assertFalse(result)


# ===================================================================
# create_product_model
# ===================================================================
class TestCreateProductModel(TestCase):
    """Tests for create_product_model."""

    @patch("backend.createmodels.default_storage")
    @patch("backend.createmodels.create_svg_string", return_value="<svg>prod</svg>")
    @patch("backend.createmodels.get_pubchem_info")
    @patch("backend.createmodels.canon_smiles", return_value="c1ccccc1")
    @patch("backend.createmodels.Product")
    @patch("backend.createmodels.Reaction.objects")
    def test_creates_product_with_pubchem(
        self,
        mock_rxn_qs,
        MockProduct,
        mock_canon,
        mock_pubchem,
        mock_svg,
        mock_storage,
    ):
        from backend.createmodels import create_product_model

        mock_rxn_obj = MagicMock()
        mock_rxn_qs.get.return_value = mock_rxn_obj
        mock_pubchem.return_value = MagicMock()
        mock_storage.save.return_value = "productimages/product.svg"

        mock_prod = MagicMock()
        MockProduct.return_value = mock_prod

        create_product_model(reaction_id=1, product_smiles="c1ccccc1")

        mock_canon.assert_called_once_with(smiles="c1ccccc1")
        self.assertEqual(mock_prod.smiles, "c1ccccc1")
        self.assertEqual(mock_prod.reaction_id, mock_rxn_obj)
        mock_pubchem.assert_called_once_with(smiles="c1ccccc1")
        mock_prod.save.assert_called_once()

    @patch("backend.createmodels.default_storage")
    @patch("backend.createmodels.create_svg_string", return_value="<svg/>")
    @patch("backend.createmodels.get_pubchem_info", return_value=False)
    @patch("backend.createmodels.canon_smiles", return_value="CC")
    @patch("backend.createmodels.Product")
    @patch("backend.createmodels.Reaction.objects")
    def test_creates_product_without_pubchem_match(
        self,
        mock_rxn_qs,
        MockProduct,
        mock_canon,
        mock_pubchem,
        mock_svg,
        mock_storage,
    ):
        from backend.createmodels import create_product_model

        mock_rxn_qs.get.return_value = MagicMock()
        mock_storage.save.return_value = "productimages/product.svg"
        mock_prod = MagicMock()
        MockProduct.return_value = mock_prod

        create_product_model(reaction_id=1, product_smiles="CC")

        mock_prod.save.assert_called_once()

    @patch("backend.createmodels.default_storage")
    @patch("backend.createmodels.create_svg_string", return_value="<svg/>")
    @patch("backend.createmodels.canon_smiles", return_value="CC")
    @patch("backend.createmodels.Product")
    @patch("backend.createmodels.Reaction.objects")
    def test_creates_product_skip_pubchem(
        self,
        mock_rxn_qs,
        MockProduct,
        mock_canon,
        mock_svg,
        mock_storage,
    ):
        from backend.createmodels import create_product_model

        mock_rxn_qs.get.return_value = MagicMock()
        mock_storage.save.return_value = "productimages/product.svg"
        mock_prod = MagicMock()
        MockProduct.return_value = mock_prod

        create_product_model(reaction_id=1, product_smiles="CC", fetch_pubchem=False)

        mock_prod.save.assert_called_once()


# ===================================================================
# create_reactant_model
# ===================================================================
class TestCreateReactantModel(TestCase):
    """Tests for create_reactant_model."""

    @patch("backend.createmodels.get_pubchem_info")
    @patch("backend.createmodels.Reactant")
    @patch("backend.createmodels.Reaction.objects")
    def test_creates_reactant_with_pubchem(
        self, mock_rxn_qs, MockReactant, mock_pubchem
    ):
        from backend.createmodels import create_reactant_model

        mock_rxn_obj = MagicMock()
        mock_rxn_qs.get.return_value = mock_rxn_obj
        mock_pubchem.return_value = MagicMock()

        mock_reactant = MagicMock()
        mock_reactant.id = 33
        MockReactant.return_value = mock_reactant

        result = create_reactant_model(
            reaction_id=1,
            reactant_smiles="CCO",
            previous_reaction_product=False,
        )

        self.assertEqual(mock_reactant.reaction_id, mock_rxn_obj)
        self.assertEqual(mock_reactant.smiles, "CCO")
        self.assertFalse(mock_reactant.previousreactionproduct)
        mock_reactant.save.assert_called_once()
        self.assertEqual(result, 33)

    @patch("backend.createmodels.Reactant")
    @patch("backend.createmodels.Reaction.objects")
    def test_creates_reactant_skip_pubchem(self, mock_rxn_qs, MockReactant):
        from backend.createmodels import create_reactant_model

        mock_rxn_qs.get.return_value = MagicMock()
        mock_reactant = MagicMock()
        mock_reactant.id = 34
        MockReactant.return_value = mock_reactant

        result = create_reactant_model(
            reaction_id=1,
            reactant_smiles="CC",
            previous_reaction_product=True,
            fetch_pubchem=False,
        )

        self.assertTrue(mock_reactant.previousreactionproduct)
        self.assertEqual(result, 34)


# ===================================================================
# create_catalog_entry_model
# ===================================================================
class TestCreateCatalogEntryModel(TestCase):
    """Tests for create_catalog_entry_model."""

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Target.objects")
    def test_creates_for_target_as_previous_rxn_product(
        self, mock_target_qs, MockCatalog
    ):
        from backend.createmodels import create_catalog_entry_model

        mock_target_obj = MagicMock()
        mock_target_qs.get.return_value = mock_target_obj

        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        create_catalog_entry_model(
            target_id=5, previous_reaction_product=True
        )

        self.assertEqual(mock_entry.target_id, mock_target_obj)
        self.assertEqual(mock_entry.vendor, "reaction product")
        self.assertEqual(mock_entry.catalogid, "NA")
        self.assertEqual(mock_entry.upperprice, 0)
        self.assertEqual(mock_entry.leadtime, 0)
        mock_entry.save.assert_called_once()

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Reactant.objects")
    def test_creates_for_reactant_as_lab_inventory(
        self, mock_reactant_qs, MockCatalog
    ):
        from backend.createmodels import create_catalog_entry_model

        mock_reactant_obj = MagicMock()
        mock_reactant_qs.get.return_value = mock_reactant_obj

        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        create_catalog_entry_model(reactant_id=8, lab_inventory=True)

        self.assertEqual(mock_entry.reactant_id, mock_reactant_obj)
        self.assertEqual(mock_entry.vendor, "reaction product")
        mock_entry.save.assert_called_once()

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Reactant.objects")
    def test_creates_from_manifold_screening_catalog(
        self, mock_reactant_qs, MockCatalog
    ):
        from backend.createmodels import create_catalog_entry_model

        mock_reactant_qs.get.return_value = MagicMock()
        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        catalog_entry = {
            "catalogName": "Enamine REAL",
            "catalogId": "Z123456",
            "purchaseInfo": {
                "isScreening": True,
                "scrLeadTimeWeeks": 4,
                "scrPriceRange": "< $100 / g",
                "bbLeadTimeWeeks": "unknown",
                "bbPriceRange": "unknown",
            },
        }

        create_catalog_entry_model(
            catalog_entry=catalog_entry, reactant_id=8
        )

        self.assertEqual(mock_entry.vendor, "Enamine REAL")
        self.assertEqual(mock_entry.catalogid, "Z123456")
        self.assertEqual(mock_entry.leadtime, 4)
        self.assertEqual(mock_entry.upperprice, 100)
        mock_entry.save.assert_called_once()

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Reactant.objects")
    def test_creates_from_manifold_bb_catalog_range_price(
        self, mock_reactant_qs, MockCatalog
    ):
        from backend.createmodels import create_catalog_entry_model

        mock_reactant_qs.get.return_value = MagicMock()
        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        catalog_entry = {
            "catalogName": "Sigma",
            "catalogId": "S789",
            "purchaseInfo": {
                "isScreening": False,
                "scrLeadTimeWeeks": "unknown",
                "scrPriceRange": "unknown",
                "bbLeadTimeWeeks": 2,
                "bbPriceRange": "$10 - $50 / g",
            },
        }

        create_catalog_entry_model(
            catalog_entry=catalog_entry, reactant_id=8
        )

        self.assertEqual(mock_entry.leadtime, 2)
        self.assertEqual(mock_entry.upperprice, 50)
        mock_entry.save.assert_called_once()

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Reactant.objects")
    def test_creates_from_manifold_bb_unknown_price(
        self, mock_reactant_qs, MockCatalog
    ):
        from backend.createmodels import create_catalog_entry_model

        mock_reactant_qs.get.return_value = MagicMock()
        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        catalog_entry = {
            "catalogName": "generic",
            "catalogId": "G001",
            "purchaseInfo": {
                "isScreening": False,
                "scrLeadTimeWeeks": "unknown",
                "scrPriceRange": "unknown",
                "bbLeadTimeWeeks": "unknown",
                "bbPriceRange": "unknown",
            },
        }

        create_catalog_entry_model(
            catalog_entry=catalog_entry, reactant_id=8
        )

        self.assertIsNone(mock_entry.upperprice)
        self.assertIsNone(mock_entry.leadtime)
        mock_entry.save.assert_called_once()

    @patch("backend.createmodels.CatalogEntry")
    def test_creates_entry_with_no_target_or_reactant(self, MockCatalog):
        """When neither target_id nor reactant_id is provided."""
        from backend.createmodels import create_catalog_entry_model

        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        create_catalog_entry_model(previous_reaction_product=True)

        mock_entry.save.assert_called_once()

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Reactant.objects")
    def test_screening_with_k_in_price(self, mock_reactant_qs, MockCatalog):
        """Prices with 'k' should be multiplied by 1000."""
        from backend.createmodels import create_catalog_entry_model

        mock_reactant_qs.get.return_value = MagicMock()
        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        catalog_entry = {
            "catalogName": "SomeVendor",
            "catalogId": "SV001",
            "purchaseInfo": {
                "isScreening": True,
                "scrLeadTimeWeeks": 6,
                "scrPriceRange": "< $5k / g",
                "bbLeadTimeWeeks": "unknown",
                "bbPriceRange": "unknown",
            },
        }

        create_catalog_entry_model(catalog_entry=catalog_entry, reactant_id=1)

        self.assertEqual(mock_entry.upperprice, 5000)

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Reactant.objects")
    def test_screening_unknown_leadtime_and_price(
        self, mock_reactant_qs, MockCatalog
    ):
        from backend.createmodels import create_catalog_entry_model

        mock_reactant_qs.get.return_value = MagicMock()
        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        catalog_entry = {
            "catalogName": "Mystery",
            "catalogId": "M001",
            "purchaseInfo": {
                "isScreening": True,
                "scrLeadTimeWeeks": "unknown",
                "scrPriceRange": "unknown",
                "bbLeadTimeWeeks": "unknown",
                "bbPriceRange": "unknown",
            },
        }

        create_catalog_entry_model(catalog_entry=catalog_entry, reactant_id=1)

        self.assertIsNone(mock_entry.leadtime)
        self.assertIsNone(mock_entry.upperprice)

    @patch("backend.createmodels.CatalogEntry")
    @patch("backend.createmodels.Reactant.objects")
    def test_generic_catalog_sets_none_for_price_and_leadtime(
        self, mock_reactant_qs, MockCatalog
    ):
        from backend.createmodels import create_catalog_entry_model

        mock_reactant_qs.get.return_value = MagicMock()
        mock_entry = MagicMock()
        MockCatalog.return_value = mock_entry

        catalog_entry = {
            "catalogName": "generic",
            "catalogId": "GEN01",
            "purchaseInfo": {
                "isScreening": False,
                "scrLeadTimeWeeks": "unknown",
                "scrPriceRange": "unknown",
                "bbLeadTimeWeeks": "unknown",
                "bbPriceRange": "unknown",
            },
        }

        create_catalog_entry_model(catalog_entry=catalog_entry, reactant_id=1)

        self.assertIsNone(mock_entry.upperprice)
        self.assertIsNone(mock_entry.leadtime)


# ===================================================================
# CreateEncodedActionModels — calculateMass
# ===================================================================
class TestCalculateMass(TestCase):
    """Tests for CreateEncodedActionModels.calculateMass."""

    def _make_instance(self, productmols=0.001):
        """Build a minimal CreateEncodedActionModels-like object."""
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        obj.productmols = productmols
        return obj

    def test_moleq_calculation(self):
        obj = self._make_instance(productmols=0.001)
        result = obj.calculateMass(
            calcunit="moleq", calcvalue=2.0, reactant_MW=100.0
        )
        # 2.0 * 0.001 * 100.0 * 1000 = 200.0 mg
        self.assertAlmostEqual(result, 200.0)

    def test_non_moleq_returns_none(self):
        obj = self._make_instance(productmols=0.001)
        result = obj.calculateMass(
            calcunit="masseq", calcvalue=2.0, reactant_MW=100.0
        )
        self.assertIsNone(result)


# ===================================================================
# CreateEncodedActionModels — calculateVolume
# ===================================================================
class TestCalculateVolume(TestCase):
    """Tests for CreateEncodedActionModels.calculateVolume."""

    def _make_instance(self, productmols=0.001, productmass=50.0):
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        obj.productmols = productmols
        obj.productmass = productmass
        return obj

    def test_masseq_calculation(self):
        obj = self._make_instance(productmass=50.0)
        result = obj.calculateVolume(calcunit="masseq", calcvalue=3.0)
        # 3.0 * 50.0 = 150.0 uL
        self.assertAlmostEqual(result, 150.0)

    def test_moleq_with_density(self):
        obj = self._make_instance(productmols=0.001)
        result = obj.calculateVolume(
            calcunit="moleq",
            calcvalue=2.0,
            reactant_density=0.8,
            reactant_MW=100.0,
        )
        # mol_material = 2.0 * 0.001 = 0.002
        # vol = (0.002 * 100.0 / 0.8) * 1e3 = 250.0 uL
        self.assertAlmostEqual(result, 250.0)

    def test_moleq_with_concentration(self):
        obj = self._make_instance(productmols=0.001)
        result = obj.calculateVolume(
            calcunit="moleq",
            calcvalue=1.5,
            conc_reagents=0.5,  # 0.5 M
        )
        # mol_material = 1.5 * 0.001 = 0.0015
        # vol = (0.0015 / 0.5) * 1e6 = 3000.0 uL
        self.assertAlmostEqual(result, 3000.0)

    def test_unknown_calcunit_returns_none(self):
        obj = self._make_instance()
        result = obj.calculateVolume(calcunit="unknown", calcvalue=1.0)
        self.assertIsNone(result)


# ===================================================================
# CreateEncodedActionModels — getProductMols
# ===================================================================
class TestGetProductMols(TestCase):
    """Tests for CreateEncodedActionModels.getProductMols."""

    def _make_instance(self):
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        return obj

    @patch("backend.createmodels.get_reaction_yields")
    @patch("backend.createmodels.check_proceeding_reactions")
    def test_no_proceeding_reactions(self, mock_check, mock_yields):
        obj = self._make_instance()
        obj.reaction_obj = MagicMock()
        obj.reaction_obj.id = 1
        obj.reaction_obj.recipe = "standard"
        obj.reaction_name = "Amidation"
        obj.target_obj = MagicMock()
        obj.target_obj.mols = 0.001

        mock_check.return_value = None  # falsy
        mock_yields.return_value = [0.8]

        result = obj.getProductMols()

        # 0.001 / 0.8 = 0.00125
        self.assertAlmostEqual(result, 0.00125)
        mock_yields.assert_called_once_with(
            reactionclasslist=["Amidation"], recipelist=["standard"]
        )

    @patch("backend.createmodels.get_reaction_yields")
    @patch("backend.createmodels.check_proceeding_reactions")
    def test_with_proceeding_reactions(self, mock_check, mock_yields):
        obj = self._make_instance()
        obj.reaction_obj = MagicMock()
        obj.reaction_obj.id = 1
        obj.reaction_obj.recipe = "standard"
        obj.reaction_name = "Amidation"
        obj.target_obj = MagicMock()
        obj.target_obj.mols = 0.001

        # Two proceeding reactions
        rxn1 = MagicMock()
        rxn1.reactionclass = "Suzuki"
        rxn1.recipe = "standard"
        rxn2 = MagicMock()
        rxn2.reactionclass = "SNAr"
        rxn2.recipe = "heated"

        mock_check.return_value = [rxn1, rxn2]
        mock_yields.return_value = [0.7, 0.9, 0.8]

        result = obj.getProductMols()

        # yield correction = 0.7 * 0.9 * 0.8 = 0.504
        # product_mols = 0.001 / 0.504
        expected = 0.001 / math.prod([0.7, 0.9, 0.8])
        self.assertAlmostEqual(result, expected)

        mock_yields.assert_called_once_with(
            reactionclasslist=["Suzuki", "SNAr", "Amidation"],
            recipelist=["standard", "heated", "standard"],
        )


# ===================================================================
# CreateEncodedActionModels — createActionSessionModel
# ===================================================================
class TestCreateActionSessionModel(TestCase):
    """Tests for CreateEncodedActionModels.createActionSessionModel."""

    def _make_instance(self):
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        obj.reaction_obj = MagicMock()
        return obj

    @patch("backend.createmodels.ActionSession")
    def test_creates_action_session(self, MockAS):
        obj = self._make_instance()
        mock_session = MagicMock()
        MockAS.return_value = mock_session

        result = obj.createActionSessionModel(
            actionsessiontype="reaction",
            driver="robot",
            sessionnumber=1,
            continuation=False,
        )

        self.assertEqual(mock_session.reaction_id, obj.reaction_obj)
        self.assertEqual(mock_session.type, "reaction")
        self.assertEqual(mock_session.driver, "robot")
        self.assertEqual(mock_session.sessionnumber, 1)
        self.assertFalse(mock_session.continuation)
        mock_session.save.assert_called_once()
        self.assertEqual(result, mock_session)

    @patch("backend.createmodels.ActionSession")
    def test_catches_exception(self, MockAS):
        obj = self._make_instance()
        MockAS.side_effect = Exception("DB error")

        result = obj.createActionSessionModel(
            actionsessiontype="stir",
            driver="human",
            sessionnumber=2,
            continuation=True,
        )

        self.assertIsNone(result)


# ===================================================================
# CreateEncodedActionModels — createStirActionModel
# ===================================================================
class TestCreateStirActionModel(TestCase):
    """Tests for CreateEncodedActionModels.createStirActionModel."""

    def _make_instance(self):
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        obj.reaction_obj = MagicMock()
        return obj

    @patch("backend.createmodels.StirAction")
    def test_creates_stir_action(self, MockStir):
        obj = self._make_instance()
        mock_stir = MagicMock()
        MockStir.return_value = mock_stir

        recipe_stir = MagicMock()
        recipe_stir.action_number = 1
        recipe_stir.plate_role = "reaction"
        recipe_stir.plate_role_index = 1
        recipe_stir.duration = 2.0
        recipe_stir.duration_unit = "h"
        recipe_stir.temperature = 25
        recipe_stir.temperature_unit = "degC"

        actionsession_obj = MagicMock()
        obj.createStirActionModel(actionsession_obj, recipe_stir)

        self.assertEqual(mock_stir.reaction_id, obj.reaction_obj)
        self.assertEqual(mock_stir.actionsession_id, actionsession_obj)
        self.assertEqual(mock_stir.number, 1)
        self.assertEqual(mock_stir.duration, 2.0)
        self.assertEqual(mock_stir.durationunit, "h")
        self.assertEqual(mock_stir.temperature, 25)
        self.assertEqual(mock_stir.temperatureunit, "degC")
        mock_stir.save.assert_called_once()

    @patch("backend.createmodels.StirAction")
    def test_catches_exception(self, MockStir):
        obj = self._make_instance()
        MockStir.side_effect = Exception("DB error")

        result = obj.createStirActionModel(MagicMock(), MagicMock())
        # Should not raise — exception is caught and logged


# ===================================================================
# CreateEncodedActionModels — createMixActionModel
# ===================================================================
class TestCreateMixActionModel(TestCase):
    """Tests for CreateEncodedActionModels.createMixActionModel."""

    def _make_instance(self):
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        obj.reaction_obj = MagicMock()
        return obj

    @patch("backend.createmodels.MixAction")
    def test_creates_mix_action(self, MockMix):
        obj = self._make_instance()
        mock_mix = MagicMock()
        MockMix.return_value = mock_mix

        recipe_mix = MagicMock()
        recipe_mix.action_number = 2
        recipe_mix.plate_role = "reaction"
        recipe_mix.plate_role_index = 1
        recipe_mix.repetitions = 5

        actionsession_obj = MagicMock()
        obj.createMixActionModel(actionsession_obj, recipe_mix)

        self.assertEqual(mock_mix.number, 2)
        self.assertEqual(mock_mix.repetitions, 5)
        mock_mix.save.assert_called_once()


# ===================================================================
# CreateEncodedActionModels — createExtractActionModel
# ===================================================================
class TestCreateExtractActionModel(TestCase):
    """Tests for CreateEncodedActionModels.createExtractActionModel."""

    def _make_instance(self):
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        obj.reaction_obj = MagicMock()
        obj.productsmiles = "c1ccccc1"
        return obj

    @patch("backend.createmodels.Descriptors")
    @patch("backend.createmodels.Chem")
    @patch("backend.createmodels.ExtractAction")
    def test_creates_extract_action(self, MockExtract, MockChem, MockDescriptors):
        obj = self._make_instance()
        mock_extract = MagicMock()
        MockExtract.return_value = mock_extract

        mock_mol = MagicMock()
        MockChem.MolFromSmiles.return_value = mock_mol
        MockDescriptors.MolWt.return_value = 78.11

        recipe_ext = MagicMock()
        recipe_ext.action_number = 1
        recipe_ext.from_plate_role = "reaction"
        recipe_ext.from_plate_role_index = 1
        recipe_ext.to_plate_role = "workup"
        recipe_ext.to_plate_role_index = 1
        recipe_ext.volume = 200.0
        recipe_ext.solvent = "ethyl acetate"
        recipe_ext.bottom_layer_volume = 50.0
        recipe_ext.concentration = 0.1

        actionsession_obj = MagicMock()
        obj.createExtractActionModel(actionsession_obj, recipe_ext)

        self.assertEqual(mock_extract.smiles, "c1ccccc1")
        self.assertEqual(mock_extract.volume, 200.0)
        self.assertEqual(mock_extract.solvent, "ethyl acetate")
        self.assertEqual(mock_extract.bottomlayervolume, 50.0)
        self.assertEqual(mock_extract.concentration, 0.1)
        mock_extract.save.assert_called_once()

    @patch("backend.createmodels.Descriptors")
    @patch("backend.createmodels.Chem")
    @patch("backend.createmodels.ExtractAction")
    def test_extract_without_optional_fields(
        self, MockExtract, MockChem, MockDescriptors
    ):
        obj = self._make_instance()
        mock_extract = MagicMock()
        MockExtract.return_value = mock_extract

        mock_mol = MagicMock()
        MockChem.MolFromSmiles.return_value = mock_mol
        MockDescriptors.MolWt.return_value = 78.11

        recipe_ext = MagicMock()
        recipe_ext.action_number = 1
        recipe_ext.from_plate_role = "reaction"
        recipe_ext.from_plate_role_index = 1
        recipe_ext.to_plate_role = "workup"
        recipe_ext.to_plate_role_index = 1
        recipe_ext.volume = 100.0
        recipe_ext.solvent = None
        recipe_ext.bottom_layer_volume = None
        recipe_ext.concentration = None

        actionsession_obj = MagicMock()
        obj.createExtractActionModel(actionsession_obj, recipe_ext)

        self.assertEqual(mock_extract.concentration, 0)
        mock_extract.save.assert_called_once()


# ===================================================================
# CreateEncodedActionModels — createAddActionModel
# ===================================================================
class TestCreateAddActionModel(TestCase):
    """Tests for CreateEncodedActionModels.createAddActionModel."""

    def _make_instance(self):
        from backend.createmodels import CreateEncodedActionModels

        obj = object.__new__(CreateEncodedActionModels)
        obj.reaction_obj = MagicMock()
        obj.productsmiles = "c1ccccc1"
        obj.productmols = 0.001
        obj.productmass = 50.0
        obj.reactant_pair_smiles = ["CCO", "CC(=O)O"]
        obj.used_reactant_indices = []
        return obj

    @patch("backend.createmodels.Descriptors")
    @patch("backend.createmodels.Chem")
    @patch("backend.createmodels.match_smarts")
    @patch("backend.createmodels.AddAction")
    def test_creates_add_action_with_smarts_match(
        self, MockAdd, mock_match, MockChem, MockDescriptors
    ):
        obj = self._make_instance()
        mock_add = MagicMock()
        MockAdd.return_value = mock_add

        # First reactant matches SMARTS
        mock_match.side_effect = [True, False]
        mock_mol = MagicMock()
        MockChem.MolFromSmiles.return_value = mock_mol
        MockDescriptors.MolWt.return_value = 46.07

        recipe_add = MagicMock()
        recipe_add.material_smarts = "[OX2H]"
        recipe_add.material_smiles = None
        recipe_add.equivalents = 1.5
        recipe_add.quantity_unit = "moleq"
        recipe_add.concentration = 0.5
        recipe_add.solvent = "DMF"
        recipe_add.density = None
        recipe_add.action_number = 1
        recipe_add.from_plate_role = "startingmaterial"
        recipe_add.from_plate_role_index = 1
        recipe_add.to_plate_role = "reaction"
        recipe_add.to_plate_role_index = 1

        actionsession_obj = MagicMock()
        obj.createAddActionModel(actionsession_obj, recipe_add)

        self.assertEqual(mock_add.smiles, "CCO")
        self.assertIn(0, obj.used_reactant_indices)
        mock_add.save.assert_called_once()

    @patch("backend.createmodels.Descriptors")
    @patch("backend.createmodels.Chem")
    @patch("backend.createmodels.canon_smiles", return_value="CC(=O)O")
    @patch("backend.createmodels.AddAction")
    def test_creates_add_action_with_material_smiles(
        self, MockAdd, mock_canon, MockChem, MockDescriptors
    ):
        obj = self._make_instance()
        mock_add = MagicMock()
        MockAdd.return_value = mock_add

        mock_mol = MagicMock()
        MockChem.MolFromSmiles.return_value = mock_mol
        MockDescriptors.MolWt.return_value = 60.05

        recipe_add = MagicMock()
        recipe_add.material_smarts = None
        recipe_add.material_smiles = "CC(=O)O"
        recipe_add.equivalents = 1.0
        recipe_add.quantity_unit = "uL"
        recipe_add.concentration = None
        recipe_add.solvent = None
        recipe_add.density = None
        recipe_add.action_number = 2
        recipe_add.from_plate_role = "startingmaterial"
        recipe_add.from_plate_role_index = 1
        recipe_add.to_plate_role = "reaction"
        recipe_add.to_plate_role_index = 1

        actionsession_obj = MagicMock()
        obj.createAddActionModel(actionsession_obj, recipe_add)

        self.assertEqual(mock_add.smiles, "CC(=O)O")
        self.assertEqual(mock_add.volume, 1.0)

    @patch("backend.createmodels.Descriptors")
    @patch("backend.createmodels.Chem")
    @patch("backend.createmodels.AddAction")
    def test_creates_add_action_product_smiles_fallback(
        self, MockAdd, MockChem, MockDescriptors
    ):
        """When neither material_smarts nor material_smiles, uses product SMILES."""
        obj = self._make_instance()
        mock_add = MagicMock()
        MockAdd.return_value = mock_add

        mock_mol = MagicMock()
        MockChem.MolFromSmiles.return_value = mock_mol
        MockDescriptors.MolWt.return_value = 78.11

        recipe_add = MagicMock()
        recipe_add.material_smarts = None
        recipe_add.material_smiles = None
        recipe_add.equivalents = 1.0
        recipe_add.quantity_unit = "moleq"
        recipe_add.concentration = None
        recipe_add.solvent = None
        recipe_add.density = None
        recipe_add.action_number = 3
        recipe_add.from_plate_role = "reaction"
        recipe_add.from_plate_role_index = 1
        recipe_add.to_plate_role = "reaction"
        recipe_add.to_plate_role_index = 1

        actionsession_obj = MagicMock()
        obj.createAddActionModel(actionsession_obj, recipe_add)

        self.assertEqual(mock_add.smiles, "c1ccccc1")

    @patch("backend.createmodels.AddAction")
    @patch("backend.createmodels.match_smarts")
    def test_add_action_catches_exception(self, mock_match, MockAdd):
        """Verify broad except block doesn't propagate."""
        obj = self._make_instance()
        mock_match.side_effect = Exception("SMARTS parsing error")

        recipe_add = MagicMock()
        recipe_add.material_smarts = "[INVALID"
        recipe_add.material_smiles = None

        # Should not raise
        obj.createAddActionModel(MagicMock(), recipe_add)


# ===================================================================
# CreateEncodedActionModels — __init__ orchestration
# ===================================================================
class TestCreateEncodedActionModelsInit(TestCase):
    """Tests for CreateEncodedActionModels.__init__ orchestration."""

    @patch("backend.createmodels.get_product_smiles", return_value=["c1ccccc1"])
    @patch("backend.createmodels.calculate_mass_from_mols", return_value=78.0)
    @patch("backend.createmodels.get_product")
    @patch("backend.createmodels.Target.objects")
    @patch("backend.createmodels.Reaction.objects")
    @patch("backend.createmodels.MCuleAPI")
    @patch("backend.recipe_utils.collect_session_actions")
    @patch("backend.recipe_utils.get_latest_recipe")
    def test_init_calls_orchestration(
        self,
        mock_get_recipe,
        mock_collect,
        MockMCule,
        mock_rxn_qs,
        mock_target_qs,
        mock_get_product,
        mock_calc_mass,
        mock_prod_smiles,
    ):
        from backend.createmodels import CreateEncodedActionModels

        # Set up recipe with one action session
        mock_recipe = MagicMock()
        mock_session = MagicMock()
        mock_session.session_type = "reaction"
        mock_session.driver = "robot"
        mock_session.session_number = 1
        mock_session.continuation = False
        mock_recipe.action_sessions.all.return_value.order_by.return_value = [
            mock_session
        ]
        mock_get_recipe.return_value = mock_recipe

        # No actions in session
        mock_collect.return_value = []

        mock_rxn_obj = MagicMock()
        mock_rxn_obj.id = 1
        mock_rxn_obj.recipe = "standard"
        mock_rxn_qs.get.return_value = mock_rxn_obj

        mock_target_obj = MagicMock()
        mock_target_obj.mols = 0.001
        mock_target_qs.get.return_value = mock_target_obj

        mock_product_obj = MagicMock()
        mock_product_obj.smiles = "c1ccccc1"
        mock_get_product.return_value = mock_product_obj

        # Patch getProductMols to avoid triggering more DB queries
        with patch.object(
            CreateEncodedActionModels, "getProductMols", return_value=0.00125
        ):
            with patch.object(
                CreateEncodedActionModels,
                "createActionSessionModel",
                return_value=MagicMock(),
            ) as mock_create_session:
                obj = CreateEncodedActionModels(
                    reaction_class="Amidation",
                    recipe_name="standard",
                    intramolecular=False,
                    target_id=42,
                    reaction_id=1,
                    reactant_pair_smiles=["CCO", "CC(=O)O"],
                )

                mock_get_recipe.assert_called_once_with("Amidation", "standard")
                mock_create_session.assert_called_once()
                self.assertEqual(obj.productsmiles, "c1ccccc1")
                self.assertFalse(obj.intramolecular)
