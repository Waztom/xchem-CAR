"""Tests for CAR's celery upload tasks.

Tests for:
  - upload_manifold_reaction
  - upload_custom_reaction
  - upload_combi_custom_reaction
  - canonicalize_smiles

All external dependencies (model creation, API calls, file I/O) are mocked.
These tests verify the orchestration logic: correct wiring of data through
groupby loops, conditional branching (otchem True/False), and proper
cleanup on both success and failure paths.
"""

from unittest import TestCase
from unittest.mock import patch, call, MagicMock, ANY

from backend.tasks import (
    upload_manifold_reaction,
    upload_custom_reaction,
    upload_combi_custom_reaction,
    canonicalize_smiles,
    create_ot_script,
    _process_reaction_step_sessions,
    _get_custom_sm_csv_path,
    ZipOTBatchProtocol,
)

# ---------------------------------------------------------------------------
# Shared constants & helpers
# ---------------------------------------------------------------------------

PROJECT_INFO = {
    "submittername": "Test User",
    "submitterorganisation": "Test Org",
    "proteintarget": "SARS-CoV-2 Mpro",
    "projectname": "Test Project",
}

# Mock SMILES for building blocks and products used across tests.
# These are orchestration tests, so the SMILES don't need to represent real chemistry.
SMILES_ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
SMILES_IBUPROFEN = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"
SMILES_REACTANT_A = "CC(=O)Cl"
SMILES_REACTANT_B = "Oc1ccccc1C(=O)O"
SMILES_REACTANT_C = "NCC1CCOCC1"
SMILES_REACTANT_D = "c1ccc(Br)cc1"
SMILES_PRODUCT_A = "CC(=O)Oc1ccccc1C(=O)O"
SMILES_PRODUCT_B = "CCO"


def _make_retro_validate_output(
    validated=True,
    smiles=None,
    names=None,
    concentrations=None,
    volumes=None,
    batch_tags=None,
    csv_fp="/tmp/test.csv",
):
    """Build a validate_output tuple mimicking retro-API validation result."""
    if smiles is None:
        smiles = [SMILES_ASPIRIN, SMILES_IBUPROFEN]
    if names is None:
        names = ["Target-{}".format(i + 1) for i in range(len(smiles))]
    if concentrations is None:
        concentrations = [10.0] * len(smiles)
    if volumes is None:
        volumes = [100.0] * len(smiles)
    if batch_tags is None:
        batch_tags = ["batch-1"] * len(smiles)

    uploaded_dict = {
        "target-SMILES": smiles,
        "target-names": names,
        "concentration-required-mM": concentrations,
        "amount-required-uL": volumes,
        "batch-tag": batch_tags,
    }
    validate_dict = {"field": [], "warning_string": []}
    return (
        "retro-API",
        validate_dict,
        validated,
        dict(PROJECT_INFO),
        csv_fp,
        uploaded_dict,
    )


def _make_manifold_route(
    reactions,
    molecules=None,
):
    """Build a single Manifold route dict.

    Parameters
    ----------
    reactions : list[dict]
        Each dict has keys: name, reactantSmiles, productSmiles
    molecules : list[dict] or None
        If None, auto-generates from reactant SMILES with a dummy catalog entry.
    """
    if molecules is None:
        molecules = []
        seen = set()
        for rxn in reactions:
            for smi in rxn["reactantSmiles"]:
                if smi not in seen:
                    seen.add(smi)
                    molecules.append(
                        {
                            "smiles": smi,
                            "isBuildingBlock": True,
                            "catalogEntries": [
                                {"vendorName": "TestVendor", "smiles": smi}
                            ],
                        }
                    )
    return {"reactions": reactions, "molecules": molecules}


def _make_manifold_response(routes_per_target):
    """Build a full Manifold retrosynthesis batch response.

    Parameters
    ----------
    routes_per_target : list[list[dict]]
        Outer list = per target, inner list = routes for that target.
        Each route is as returned by _make_manifold_route.
    """
    results = []
    for target_routes in routes_per_target:
        results.append({"routes": target_routes})
    return {"results": results}


# ---------------------------------------------------------------------------
# Helpers to build custom-chem and combi-custom-chem validate outputs
# ---------------------------------------------------------------------------


def _make_custom_chem_validate_output(
    validated=True,
    csv_fp="/tmp/test_custom.csv",
    n_targets=1,
):
    """Build validate_output for a single-step custom-chem upload."""
    smiles = [SMILES_PRODUCT_A] * n_targets
    names = ["Exp-{}".format(i + 1) for i in range(n_targets)]
    concentrations = [10.0] * n_targets
    volumes = [100.0] * n_targets
    batch_tags = ["batch-1"] * n_targets
    no_steps = [1] * n_targets
    # Each value is a tuple-of-tuples with one step
    reactant_pair_smiles = [((SMILES_REACTANT_A, SMILES_REACTANT_B),)] * n_targets
    reaction_groupby_columns = [("A",)] * n_targets
    reaction_names = [("Acylation",)] * n_targets
    reaction_recipes = [("standard",)] * n_targets
    product_smiles = [(SMILES_PRODUCT_A,)] * n_targets

    uploaded_dict = {
        "target-SMILES": smiles,
        "target-names": names,
        "concentration-required-mM": concentrations,
        "amount-required-uL": volumes,
        "batch-tag": batch_tags,
        "no-steps": no_steps,
        "reactant-pair-smiles": reactant_pair_smiles,
        "reaction-groupby-column": reaction_groupby_columns,
        "reaction-name": reaction_names,
        "reaction-recipe": reaction_recipes,
        "product-smiles": product_smiles,
    }
    validate_dict = {"field": [], "warning_string": []}
    return (
        "custom-chem",
        validate_dict,
        validated,
        dict(PROJECT_INFO),
        csv_fp,
        uploaded_dict,
    )


def _make_combi_validate_output(
    validated=True,
    csv_fp="/tmp/test_combi.csv",
    n_targets=2,
):
    """Build validate_output for a single-step combi-custom-chem upload."""
    smiles = [SMILES_PRODUCT_A] * n_targets
    names = ["Combi-{}".format(i + 1) for i in range(n_targets)]
    concentrations = [10.0] * n_targets
    volumes = [100.0] * n_targets
    batch_tags = ["batch-1"] * n_targets
    no_steps = [1] * n_targets
    reactant_pair_smiles = [((SMILES_REACTANT_A, SMILES_REACTANT_C),)] * n_targets
    reaction_names = [("SNAr",)] * n_targets
    reaction_recipes = [("standard",)] * n_targets
    product_smiles = [(SMILES_PRODUCT_A,)] * n_targets

    uploaded_dict = {
        "target-SMILES": smiles,
        "target-names": names,
        "concentration-required-mM": concentrations,
        "amount-required-uL": volumes,
        "batch-tag": batch_tags,
        "no-steps": no_steps,
        "reactant-pair-smiles": reactant_pair_smiles,
        "reaction-name": reaction_names,
        "reaction-recipe": reaction_recipes,
        "product-smiles": product_smiles,
    }
    validate_dict = {"field": [], "warning_string": []}
    return (
        "combi-custom-chem",
        validate_dict,
        validated,
        dict(PROJECT_INFO),
        csv_fp,
        uploaded_dict,
    )


# ===========================================================================
# upload_manifold_reaction tests
# ===========================================================================

# All patches target backend.tasks where the names are imported
MANIFOLD_PATCHES = [
    "backend.tasks.delete_tmp_file",
    "backend.tasks.create_catalog_entry_model",
    "backend.tasks.create_reactant_model",
    "backend.tasks.create_product_model",
    "backend.tasks.create_reaction_model",
    "backend.tasks.create_method_model",
    "backend.tasks.create_target_model",
    "backend.tasks.create_batch_model",
    "backend.tasks.create_project_model",
    "backend.tasks.check_previous_reaction_products",
    "backend.tasks.recipe_exists",
    "backend.tasks.get_recipe_intramolecular",
    "backend.tasks.get_manifold_retrosynthesis_batch",
]


class TestUploadManifoldNotValidated(TestCase):
    """When validated=False the task should short-circuit."""

    @patch("backend.tasks.delete_tmp_file")
    def test_returns_early_and_deletes_file(self, mock_del):
        validate_output = _make_retro_validate_output(validated=False)
        result = upload_manifold_reaction(validate_output)

        validate_dict, validated, project_info = result
        self.assertFalse(validated)
        mock_del.assert_called_once_with("/tmp/test.csv")

    @patch("backend.tasks.delete_tmp_file")
    @patch("backend.tasks.create_project_model")
    def test_no_models_created(self, mock_project, mock_del):
        validate_output = _make_retro_validate_output(validated=False)
        upload_manifold_reaction(validate_output)
        mock_project.assert_not_called()


class TestUploadManifoldHappyPath(TestCase):
    """Happy path: one target, one route, all reactions encoded → otchem=True."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_recipe = patch("backend.tasks.recipe_exists", return_value=True)
        self.patcher_intra = patch(
            "backend.tasks.get_recipe_intramolecular", return_value=False
        )
        self.patcher_manifold = patch(
            "backend.tasks.get_manifold_retrosynthesis_batch"
        )

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_recipe = self.patcher_recipe.start()
        self.mock_intra = self.patcher_intra.start()
        self.mock_manifold = self.patcher_manifold.start()

    def tearDown(self):
        patch.stopall()

    def _single_target_single_route(self):
        """Set up a single target with one building-block route + one encoded route."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": True,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [
                        {"vendorName": "Sigma", "smiles": SMILES_ASPIRIN}
                    ],
                }
            ],
            "reactions": [],
        }
        encoded_route = _make_manifold_route(
            reactions=[
                {
                    "name": "Acylation",
                    "reactantSmiles": [SMILES_REACTANT_A, SMILES_REACTANT_B],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
        )

        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, encoded_route]]
        )
        return _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN],
            names=["Aspirin"],
        )

    def test_creates_project(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        self.mock_project.assert_called_once()
        call_args = self.mock_project.call_args[0][0]
        self.assertEqual(call_args["projectname"], "Test Project")

    def test_creates_batch(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        self.mock_batch.assert_called_once_with(
            project_id=1,
            batchtag="batch-1",
        )

    def test_creates_target(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        self.mock_target.assert_called_once_with(
            batch_id=10,
            name="Aspirin",
            smiles=SMILES_ASPIRIN,
            concentration=10.0,
            volume=100.0,
        )

    def test_creates_building_block_catalog_entry(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        # Building block catalog entry from routes[0]
        self.mock_catalog.assert_any_call(
            catalog_entry={"vendorName": "Sigma", "smiles": SMILES_ASPIRIN},
            target_id=100,
        )

    def test_encoded_route_creates_method_with_otchem_true(self):
        """All reactions found encoded → otchem=True."""
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        self.mock_method.assert_called_with(
            target_id=100,
            nosteps=1,
            otchem=True,
        )

    def test_creates_reaction_with_standard_recipe(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        self.mock_reaction.assert_called_once()
        call_kwargs = self.mock_reaction.call_args[1]
        self.assertEqual(call_kwargs["reaction_class"], "Acylation")
        self.assertEqual(call_kwargs["reaction_number"], 1)
        self.assertFalse(call_kwargs["intramolecular"])
        self.assertEqual(call_kwargs["reaction_recipe"], "standard")

    def test_creates_product(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        self.mock_product.assert_called_once_with(
            reaction_id=5000,
            product_smiles=SMILES_PRODUCT_A,
            fetch_pubchem=True,
        )

    def test_creates_reactants_with_catalog_entries(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)

        # 2 reactants, both are new (no previous reaction product)
        self.assertEqual(self.mock_reactant.call_count, 2)
        # Each new reactant gets a catalog entry from the Manifold molecule data
        for c in self.mock_reactant.call_args_list:
            self.assertFalse(c[1]["previous_reaction_product"])

    def test_deletes_tmp_file(self):
        validate_output = self._single_target_single_route()
        upload_manifold_reaction(validate_output)
        self.mock_del.assert_called_once_with("/tmp/test.csv")

    def test_return_tuple(self):
        validate_output = self._single_target_single_route()
        result = upload_manifold_reaction(validate_output)

        validate_dict, validated, project_info = result
        self.assertTrue(validated)
        self.assertIn("project_id", project_info)
        self.assertEqual(project_info["project_id"], 1)


class TestUploadManifoldNotEncoded(TestCase):
    """When NOT all reactions are encoded → otchem=False path."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        # recipe_exists returns False → not encoded
        self.patcher_recipe = patch("backend.tasks.recipe_exists", return_value=False)
        self.patcher_intra = patch(
            "backend.tasks.get_recipe_intramolecular", return_value=False
        )
        self.patcher_manifold = patch(
            "backend.tasks.get_manifold_retrosynthesis_batch"
        )

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_recipe = self.patcher_recipe.start()
        self.mock_intra = self.patcher_intra.start()
        self.mock_manifold = self.patcher_manifold.start()

    def tearDown(self):
        patch.stopall()

    def test_not_encoded_creates_method_with_otchem_false(self):
        """When recipe_exists returns False, otchem should be False."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        unencoded_route = _make_manifold_route(
            reactions=[
                {
                    "name": "NovelReaction",
                    "reactantSmiles": [SMILES_REACTANT_A, SMILES_REACTANT_B],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
        )
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, unencoded_route]]
        )

        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN],
            names=["Aspirin"],
        )
        upload_manifold_reaction(validate_output)

        self.mock_method.assert_called_with(
            target_id=100,
            nosteps=1,
            otchem=False,
        )

    def test_not_encoded_no_recipe_kwarg_on_reaction(self):
        """Non-encoded reactions should NOT pass reaction_recipe."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        unencoded_route = _make_manifold_route(
            reactions=[
                {
                    "name": "NovelReaction",
                    "reactantSmiles": [SMILES_REACTANT_A, SMILES_REACTANT_B],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
        )
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, unencoded_route]]
        )

        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN],
            names=["Aspirin"],
        )
        upload_manifold_reaction(validate_output)

        call_kwargs = self.mock_reaction.call_args[1]
        self.assertNotIn("reaction_recipe", call_kwargs)


class TestUploadManifoldEdgeCases(TestCase):
    """Edge cases for Manifold upload."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_recipe = patch("backend.tasks.recipe_exists", return_value=True)
        self.patcher_intra = patch(
            "backend.tasks.get_recipe_intramolecular", return_value=False
        )
        self.patcher_manifold = patch(
            "backend.tasks.get_manifold_retrosynthesis_batch"
        )

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_recipe = self.patcher_recipe.start()
        self.mock_intra = self.patcher_intra.start()
        self.mock_manifold = self.patcher_manifold.start()

    def tearDown(self):
        patch.stopall()

    def test_no_routes_found(self):
        """When Manifold returns no routes, no target/method/reaction created."""
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[]]  # empty routes list
        )
        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"]
        )
        upload_manifold_reaction(validate_output)

        self.mock_target.assert_not_called()
        self.mock_method.assert_not_called()

    def test_no_results_key_in_response(self):
        """When Manifold response has no 'results' key, nothing happens."""
        self.mock_manifold.return_value = {"error": "API failure"}
        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"]
        )
        upload_manifold_reaction(validate_output)

        self.mock_target.assert_not_called()

    def test_multiple_batch_tags_create_multiple_batches(self):
        """Two different batch-tags should create two batch models."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": True,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [{"vendorName": "Sigma", "smiles": SMILES_ASPIRIN}],
                }
            ],
            "reactions": [],
        }
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[
                [building_block_route],
                [building_block_route],
            ]
        )

        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN, SMILES_IBUPROFEN],
            names=["Aspirin", "Ibuprofen"],
            batch_tags=["batch-A", "batch-B"],
        )
        upload_manifold_reaction(validate_output)

        self.assertEqual(self.mock_batch.call_count, 2)
        batch_tags_used = [c[1]["batchtag"] for c in self.mock_batch.call_args_list]
        self.assertIn("batch-A", batch_tags_used)
        self.assertIn("batch-B", batch_tags_used)

    def test_previous_reaction_product_detected(self):
        """When check_previous_reaction_products returns a truthy value,
        the reactant should be marked as previous_reaction_product=True."""
        self.mock_prev.return_value = True  # truthy queryset

        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        encoded_route = _make_manifold_route(
            reactions=[
                {
                    "name": "Acylation",
                    "reactantSmiles": [SMILES_REACTANT_A],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
        )
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, encoded_route]]
        )
        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"]
        )
        upload_manifold_reaction(validate_output)

        # Reactant created with previous_reaction_product=True
        self.mock_reactant.assert_called_once()
        call_kwargs = self.mock_reactant.call_args[1]
        self.assertTrue(call_kwargs["previous_reaction_product"])

        # Catalog entry marked as previous reaction product
        self.mock_catalog.assert_any_call(
            reactant_id=9000,
            previous_reaction_product=True,
        )

    def test_intramolecular_single_reactant(self):
        """Single reactant in non-encoded route → intramolecular=True."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        route = _make_manifold_route(
            reactions=[
                {
                    "name": "Cyclisation",
                    "reactantSmiles": [SMILES_REACTANT_A],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
        )
        self.mock_recipe.return_value = False  # not encoded
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, route]]
        )
        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"]
        )
        upload_manifold_reaction(validate_output)

        call_kwargs = self.mock_reaction.call_args[1]
        self.assertTrue(call_kwargs["intramolecular"])

    def test_intramolecular_encoded_route_checks_recipe(self):
        """Single reactant in encoded route, intramolecular depends on get_recipe_intramolecular."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        route = _make_manifold_route(
            reactions=[
                {
                    "name": "Cyclisation",
                    "reactantSmiles": [SMILES_REACTANT_A],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
        )
        self.mock_recipe.return_value = True  # encoded
        self.mock_intra.return_value = True  # recipe supports intramolecular
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, route]]
        )
        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"]
        )
        upload_manifold_reaction(validate_output)

        call_kwargs = self.mock_reaction.call_args[1]
        self.assertTrue(call_kwargs["intramolecular"])

    def test_two_step_route_creates_two_reactions(self):
        """A 2-step route should create 2 reactions."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        two_step_route = _make_manifold_route(
            reactions=[
                {
                    "name": "Acylation",
                    "reactantSmiles": [SMILES_REACTANT_A, SMILES_REACTANT_B],
                    "productSmiles": SMILES_PRODUCT_A,
                },
                {
                    "name": "Reduction",
                    "reactantSmiles": [SMILES_PRODUCT_A, SMILES_REACTANT_C],
                    "productSmiles": "CCO",
                },
            ],
        )
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, two_step_route]]
        )
        validate_output = _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"]
        )
        upload_manifold_reaction(validate_output)

        self.assertEqual(self.mock_reaction.call_count, 2)
        # Check reaction numbers are 1 and 2
        reaction_numbers = [
            c[1]["reaction_number"] for c in self.mock_reaction.call_args_list
        ]
        self.assertEqual(reaction_numbers, [1, 2])


# ===========================================================================
# upload_custom_reaction tests
# ===========================================================================


class TestUploadCustomNotValidated(TestCase):
    """When validated=False the custom upload task should short-circuit."""

    @patch("backend.tasks.delete_tmp_file")
    def test_returns_early_and_deletes_file(self, mock_del):
        validate_output = _make_custom_chem_validate_output(validated=False)
        result = upload_custom_reaction(validate_output)

        validate_dict, validated, project_info = result
        self.assertFalse(validated)
        mock_del.assert_called_once_with("/tmp/test_custom.csv")

    @patch("backend.tasks.delete_tmp_file")
    @patch("backend.tasks.create_project_model")
    def test_no_models_created(self, mock_project, mock_del):
        validate_output = _make_custom_chem_validate_output(validated=False)
        upload_custom_reaction(validate_output)
        mock_project.assert_not_called()


class TestUploadCustomHappyPath(TestCase):
    """Happy path for custom-chem upload: single target, single step."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_temp = patch(
            "backend.tasks.get_recipe_stir_temperature", return_value=25.0
        )
        self.patcher_exact = patch("backend.tasks.get_exact_search")

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_temp = self.patcher_temp.start()
        self.mock_exact = self.patcher_exact.start()

    def tearDown(self):
        patch.stopall()

    def test_creates_project(self):
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        self.mock_project.assert_called_once()

    def test_creates_batch_per_tag(self):
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        self.mock_batch.assert_called_once_with(
            project_id=1,
            batchtag="batch-1",
        )

    def test_creates_target(self):
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        self.mock_target.assert_called_once_with(
            batch_id=10,
            name="Exp-1",
            smiles=SMILES_PRODUCT_A,
            concentration=10.0,
            volume=100.0,
        )

    def test_method_always_otchem_true(self):
        """Custom-chem uploads always set otchem=True."""
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        self.mock_method.assert_called_once_with(
            target_id=100,
            nosteps=1,
            otchem=True,
        )

    def test_creates_reaction_with_temperature_and_groupby(self):
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        self.mock_reaction.assert_called_once()
        call_kwargs = self.mock_reaction.call_args[1]
        self.assertEqual(call_kwargs["reaction_class"], "Acylation")
        self.assertEqual(call_kwargs["reaction_recipe"], "standard")
        self.assertEqual(call_kwargs["reaction_temperature"], 25.0)
        self.assertEqual(call_kwargs["groupby_column"], "A")
        self.assertEqual(call_kwargs["reaction_number"], 1)
        self.assertFalse(call_kwargs["intramolecular"])

    def test_creates_product(self):
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        self.mock_product.assert_called_once_with(
            reaction_id=5000,
            product_smiles=SMILES_PRODUCT_A,
            fetch_pubchem=True,
        )

    def test_reactant_not_previous_creates_lab_inventory_entry(self):
        """New reactants get a lab_inventory=True catalog entry."""
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        # 2 reactants
        self.assertEqual(self.mock_reactant.call_count, 2)
        # Each gets a lab inventory catalog entry
        self.mock_catalog.assert_any_call(
            reactant_id=9000,
            previous_reaction_product=False,
            lab_inventory=True,
        )

    def test_fetch_catalogue_calls_exact_search(self):
        """When fetch_catalogue=True, get_exact_search is called per reactant."""
        self.mock_exact.return_value = {
            "results": [{"vendorName": "Enamine", "smiles": SMILES_REACTANT_A}]
        }
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output, fetch_catalogue=True)

        self.assertEqual(self.mock_exact.call_count, 2)

    def test_fetch_catalogue_false_skips_exact_search(self):
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output, fetch_catalogue=False)

        self.mock_exact.assert_not_called()

    def test_deletes_tmp_file(self):
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)
        self.mock_del.assert_called_once_with("/tmp/test_custom.csv")

    def test_return_tuple(self):
        validate_output = _make_custom_chem_validate_output()
        result = upload_custom_reaction(validate_output)

        validate_dict, validated, project_info = result
        self.assertTrue(validated)
        self.assertIn("project_id", project_info)

    def test_multiple_targets(self):
        """Two targets should create two target models and two methods."""
        validate_output = _make_custom_chem_validate_output(n_targets=2)
        upload_custom_reaction(validate_output)

        self.assertEqual(self.mock_target.call_count, 2)
        self.assertEqual(self.mock_method.call_count, 2)

    def test_previous_reaction_product(self):
        """When a reactant matches a previous reaction product."""
        self.mock_prev.return_value = True  # truthy
        validate_output = _make_custom_chem_validate_output()
        upload_custom_reaction(validate_output)

        # All reactant calls should have previous_reaction_product=True
        for c in self.mock_reactant.call_args_list:
            self.assertTrue(c[1]["previous_reaction_product"])
        # And catalog entries marked as previous product
        self.mock_catalog.assert_any_call(
            reactant_id=9000,
            previous_reaction_product=True,
        )


# ===========================================================================
# upload_combi_custom_reaction tests
# ===========================================================================


class TestUploadCombiNotValidated(TestCase):
    """When validated=False the combi upload task should short-circuit."""

    @patch("backend.tasks.delete_tmp_file")
    def test_returns_early_and_deletes_file(self, mock_del):
        validate_output = _make_combi_validate_output(validated=False)
        result = upload_combi_custom_reaction(validate_output)

        validate_dict, validated, project_info = result
        self.assertFalse(validated)
        mock_del.assert_called_once_with("/tmp/test_combi.csv")

    @patch("backend.tasks.delete_tmp_file")
    @patch("backend.tasks.create_project_model")
    def test_no_models_created(self, mock_project, mock_del):
        validate_output = _make_combi_validate_output(validated=False)
        upload_combi_custom_reaction(validate_output)
        mock_project.assert_not_called()


class TestUploadCombiHappyPath(TestCase):
    """Happy path for combi-custom-chem upload."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_temp = patch(
            "backend.tasks.get_recipe_stir_temperature", return_value=80.0
        )
        self.patcher_exact = patch("backend.tasks.get_exact_search")

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_temp = self.patcher_temp.start()
        self.mock_exact = self.patcher_exact.start()

    def tearDown(self):
        patch.stopall()

    def test_creates_project_and_batch(self):
        validate_output = _make_combi_validate_output()
        upload_combi_custom_reaction(validate_output)

        self.mock_project.assert_called_once()
        self.mock_batch.assert_called_once_with(
            project_id=1,
            batchtag="batch-1",
        )

    def test_creates_targets_per_combi_row(self):
        validate_output = _make_combi_validate_output(n_targets=2)
        upload_combi_custom_reaction(validate_output)

        self.assertEqual(self.mock_target.call_count, 2)

    def test_method_always_otchem_true(self):
        validate_output = _make_combi_validate_output(n_targets=1)
        upload_combi_custom_reaction(validate_output)

        self.mock_method.assert_called_once_with(
            target_id=100,
            nosteps=1,
            otchem=True,
        )

    def test_reaction_has_temperature_but_no_groupby(self):
        """Combi upload creates reactions with temperature but without groupby_column."""
        validate_output = _make_combi_validate_output(n_targets=1)
        upload_combi_custom_reaction(validate_output)

        self.mock_reaction.assert_called_once()
        call_kwargs = self.mock_reaction.call_args[1]
        self.assertEqual(call_kwargs["reaction_class"], "SNAr")
        self.assertEqual(call_kwargs["reaction_temperature"], 80.0)
        # No groupby_column for combi uploads
        self.assertNotIn("groupby_column", call_kwargs)

    def test_creates_product(self):
        validate_output = _make_combi_validate_output(n_targets=1)
        upload_combi_custom_reaction(validate_output)

        self.mock_product.assert_called_once_with(
            reaction_id=5000,
            product_smiles=SMILES_PRODUCT_A,
            fetch_pubchem=True,
        )

    def test_creates_reactants_with_lab_inventory(self):
        """New reactants get lab_inventory=True catalog entries."""
        validate_output = _make_combi_validate_output(n_targets=1)
        upload_combi_custom_reaction(validate_output)

        self.assertEqual(self.mock_reactant.call_count, 2)
        self.mock_catalog.assert_any_call(
            reactant_id=9000,
            previous_reaction_product=False,
            lab_inventory=True,
        )

    def test_fetch_catalogue_calls_exact_search(self):
        self.mock_exact.return_value = {
            "results": [{"vendorName": "Enamine", "smiles": SMILES_REACTANT_C}]
        }
        validate_output = _make_combi_validate_output(n_targets=1)
        upload_combi_custom_reaction(validate_output, fetch_catalogue=True)

        self.assertEqual(self.mock_exact.call_count, 2)

    def test_fetch_catalogue_false_skips_search(self):
        validate_output = _make_combi_validate_output(n_targets=1)
        upload_combi_custom_reaction(validate_output, fetch_catalogue=False)

        self.mock_exact.assert_not_called()

    def test_deletes_tmp_file(self):
        validate_output = _make_combi_validate_output()
        upload_combi_custom_reaction(validate_output)
        self.mock_del.assert_called_once_with("/tmp/test_combi.csv")

    def test_return_tuple(self):
        validate_output = _make_combi_validate_output()
        result = upload_combi_custom_reaction(validate_output)

        validate_dict, validated, project_info = result
        self.assertTrue(validated)
        self.assertIn("project_id", project_info)

    def test_previous_reaction_product(self):
        self.mock_prev.return_value = True
        validate_output = _make_combi_validate_output(n_targets=1)
        upload_combi_custom_reaction(validate_output)

        for c in self.mock_reactant.call_args_list:
            self.assertTrue(c[1]["previous_reaction_product"])
        self.mock_catalog.assert_any_call(
            reactant_id=9000,
            previous_reaction_product=True,
        )


# ===========================================================================
# Multi-step helper variants
# ===========================================================================


def _make_custom_chem_multistep_output(
    validated=True,
    csv_fp="/tmp/test_custom_multi.csv",
):
    """Build validate_output for a 2-step custom-chem upload (1 target)."""
    uploaded_dict = {
        "target-SMILES": [SMILES_PRODUCT_B],
        "target-names": ["MultiStep-1"],
        "concentration-required-mM": [10.0],
        "amount-required-uL": [100.0],
        "batch-tag": ["batch-1"],
        "no-steps": [2],
        "reactant-pair-smiles": [
            (
                (SMILES_REACTANT_A, SMILES_REACTANT_B),  # step 1
                (SMILES_PRODUCT_A, SMILES_REACTANT_C),   # step 2
            )
        ],
        "reaction-groupby-column": [("A", "B")],
        "reaction-name": [("Acylation", "Reduction")],
        "reaction-recipe": [("standard", "standard")],
        "product-smiles": [(SMILES_PRODUCT_A, SMILES_PRODUCT_B)],
    }
    validate_dict = {"field": [], "warning_string": []}
    return (
        "custom-chem",
        validate_dict,
        validated,
        dict(PROJECT_INFO),
        csv_fp,
        uploaded_dict,
    )


def _make_combi_multistep_output(
    validated=True,
    csv_fp="/tmp/test_combi_multi.csv",
):
    """Build validate_output for a 2-step combi-custom-chem upload (1 target)."""
    uploaded_dict = {
        "target-SMILES": [SMILES_PRODUCT_B],
        "target-names": ["CombiMulti-1"],
        "concentration-required-mM": [10.0],
        "amount-required-uL": [100.0],
        "batch-tag": ["batch-1"],
        "no-steps": [2],
        "reactant-pair-smiles": [
            (
                (SMILES_REACTANT_A, SMILES_REACTANT_B),
                (SMILES_PRODUCT_A, SMILES_REACTANT_C),
            )
        ],
        "reaction-name": [("SNAr", "Amidation")],
        "reaction-recipe": [("standard", "standard")],
        "product-smiles": [(SMILES_PRODUCT_A, SMILES_PRODUCT_B)],
    }
    validate_dict = {"field": [], "warning_string": []}
    return (
        "combi-custom-chem",
        validate_dict,
        validated,
        dict(PROJECT_INFO),
        csv_fp,
        uploaded_dict,
    )


# ===========================================================================
# Multi-step custom-chem tests
# ===========================================================================


class TestUploadCustomMultiStep(TestCase):
    """Verify 2-step custom-chem upload creates correct reactions per step."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", side_effect=[5001, 5002]
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_temp = patch(
            "backend.tasks.get_recipe_stir_temperature", return_value=25.0
        )
        self.patcher_exact = patch("backend.tasks.get_exact_search")

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_temp = self.patcher_temp.start()
        self.mock_exact = self.patcher_exact.start()

    def tearDown(self):
        patch.stopall()

    def test_creates_two_reactions(self):
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)
        self.assertEqual(self.mock_reaction.call_count, 2)

    def test_reaction_numbers_sequential(self):
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)

        reaction_numbers = [
            c[1]["reaction_number"] for c in self.mock_reaction.call_args_list
        ]
        self.assertEqual(reaction_numbers, [1, 2])

    def test_reaction_classes_match_steps(self):
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)

        classes = [
            c[1]["reaction_class"] for c in self.mock_reaction.call_args_list
        ]
        self.assertEqual(classes, ["Acylation", "Reduction"])

    def test_groupby_columns_per_step(self):
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)

        groupby_cols = [
            c[1]["groupby_column"] for c in self.mock_reaction.call_args_list
        ]
        self.assertEqual(groupby_cols, ["A", "B"])

    def test_products_per_step(self):
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)

        product_calls = self.mock_product.call_args_list
        self.assertEqual(len(product_calls), 2)
        self.assertEqual(
            product_calls[0][1]["product_smiles"], SMILES_PRODUCT_A,
        )
        self.assertEqual(
            product_calls[1][1]["product_smiles"], SMILES_PRODUCT_B,
        )

    def test_method_nosteps_is_two(self):
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)

        self.mock_method.assert_called_once_with(
            target_id=100, nosteps=2, otchem=True,
        )

    def test_reactants_across_both_steps(self):
        """2 reactants per step × 2 steps = 4 reactant models."""
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)
        self.assertEqual(self.mock_reactant.call_count, 4)

    def test_recipes_per_step(self):
        validate_output = _make_custom_chem_multistep_output()
        upload_custom_reaction(validate_output)

        recipes = [
            c[1]["reaction_recipe"] for c in self.mock_reaction.call_args_list
        ]
        self.assertEqual(recipes, ["standard", "standard"])


# ===========================================================================
# Multi-step combi-custom-chem tests
# ===========================================================================


class TestUploadCombiMultiStep(TestCase):
    """Verify 2-step combi-custom-chem upload creates correct reactions per step."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", side_effect=[5001, 5002]
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_temp = patch(
            "backend.tasks.get_recipe_stir_temperature", return_value=80.0
        )
        self.patcher_exact = patch("backend.tasks.get_exact_search")

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_temp = self.patcher_temp.start()
        self.mock_exact = self.patcher_exact.start()

    def tearDown(self):
        patch.stopall()

    def test_creates_two_reactions(self):
        validate_output = _make_combi_multistep_output()
        upload_combi_custom_reaction(validate_output)
        self.assertEqual(self.mock_reaction.call_count, 2)

    def test_reaction_numbers_sequential(self):
        validate_output = _make_combi_multistep_output()
        upload_combi_custom_reaction(validate_output)

        numbers = [
            c[1]["reaction_number"] for c in self.mock_reaction.call_args_list
        ]
        self.assertEqual(numbers, [1, 2])

    def test_reaction_classes_match_steps(self):
        validate_output = _make_combi_multistep_output()
        upload_combi_custom_reaction(validate_output)

        classes = [
            c[1]["reaction_class"] for c in self.mock_reaction.call_args_list
        ]
        self.assertEqual(classes, ["SNAr", "Amidation"])

    def test_no_groupby_column_on_combi(self):
        """Combi reactions should not have groupby_column, even multi-step."""
        validate_output = _make_combi_multistep_output()
        upload_combi_custom_reaction(validate_output)

        for c in self.mock_reaction.call_args_list:
            self.assertNotIn("groupby_column", c[1])

    def test_products_per_step(self):
        validate_output = _make_combi_multistep_output()
        upload_combi_custom_reaction(validate_output)

        product_calls = self.mock_product.call_args_list
        self.assertEqual(len(product_calls), 2)
        self.assertEqual(
            product_calls[0][1]["product_smiles"], SMILES_PRODUCT_A,
        )
        self.assertEqual(
            product_calls[1][1]["product_smiles"], SMILES_PRODUCT_B,
        )

    def test_method_nosteps_is_two(self):
        validate_output = _make_combi_multistep_output()
        upload_combi_custom_reaction(validate_output)

        self.mock_method.assert_called_once_with(
            target_id=100, nosteps=2, otchem=True,
        )

    def test_reactants_across_both_steps(self):
        validate_output = _make_combi_multistep_output()
        upload_combi_custom_reaction(validate_output)
        self.assertEqual(self.mock_reactant.call_count, 4)


# ===========================================================================
# Multi-batch custom/combi tests
# ===========================================================================


def _make_custom_chem_multibatch_output(
    csv_fp="/tmp/test_custom_mb.csv",
):
    """Build validate_output for custom-chem with 2 targets in different batches."""
    uploaded_dict = {
        "target-SMILES": [SMILES_PRODUCT_A, SMILES_PRODUCT_A],
        "target-names": ["Exp-1", "Exp-2"],
        "concentration-required-mM": [10.0, 20.0],
        "amount-required-uL": [100.0, 200.0],
        "batch-tag": ["batch-X", "batch-Y"],
        "no-steps": [1, 1],
        "reactant-pair-smiles": [
            ((SMILES_REACTANT_A, SMILES_REACTANT_B),),
            ((SMILES_REACTANT_A, SMILES_REACTANT_C),),
        ],
        "reaction-groupby-column": [("A",), ("A",)],
        "reaction-name": [("Acylation",), ("SNAr",)],
        "reaction-recipe": [("standard",), ("standard",)],
        "product-smiles": [(SMILES_PRODUCT_A,), (SMILES_PRODUCT_A,)],
    }
    validate_dict = {"field": [], "warning_string": []}
    return (
        "custom-chem",
        validate_dict,
        True,
        dict(PROJECT_INFO),
        csv_fp,
        uploaded_dict,
    )


def _make_combi_multibatch_output(
    csv_fp="/tmp/test_combi_mb.csv",
):
    """Build validate_output for combi-custom-chem with 2 targets in different batches."""
    uploaded_dict = {
        "target-SMILES": [SMILES_PRODUCT_A, SMILES_PRODUCT_A],
        "target-names": ["Combi-1", "Combi-2"],
        "concentration-required-mM": [10.0, 20.0],
        "amount-required-uL": [100.0, 200.0],
        "batch-tag": ["batch-X", "batch-Y"],
        "no-steps": [1, 1],
        "reactant-pair-smiles": [
            ((SMILES_REACTANT_A, SMILES_REACTANT_B),),
            ((SMILES_REACTANT_A, SMILES_REACTANT_C),),
        ],
        "reaction-name": [("Acylation",), ("SNAr",)],
        "reaction-recipe": [("standard",), ("standard",)],
        "product-smiles": [(SMILES_PRODUCT_A,), (SMILES_PRODUCT_A,)],
    }
    validate_dict = {"field": [], "warning_string": []}
    return (
        "combi-custom-chem",
        validate_dict,
        True,
        dict(PROJECT_INFO),
        csv_fp,
        uploaded_dict,
    )


class TestUploadCustomMultiBatch(TestCase):
    """Custom-chem with targets across different batch tags."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch(
            "backend.tasks.create_batch_model", side_effect=[10, 20]
        )
        self.patcher_target = patch(
            "backend.tasks.create_target_model", side_effect=[100, 200]
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_temp = patch(
            "backend.tasks.get_recipe_stir_temperature", return_value=25.0
        )
        self.patcher_exact = patch("backend.tasks.get_exact_search")

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_temp = self.patcher_temp.start()
        self.mock_exact = self.patcher_exact.start()

    def tearDown(self):
        patch.stopall()

    def test_creates_two_batches(self):
        validate_output = _make_custom_chem_multibatch_output()
        upload_custom_reaction(validate_output)

        self.assertEqual(self.mock_batch.call_count, 2)
        tags = sorted(c[1]["batchtag"] for c in self.mock_batch.call_args_list)
        self.assertEqual(tags, ["batch-X", "batch-Y"])

    def test_each_target_in_correct_batch(self):
        validate_output = _make_custom_chem_multibatch_output()
        upload_custom_reaction(validate_output)

        self.assertEqual(self.mock_target.call_count, 2)
        # First call should use batch_id=10, second batch_id=20
        batch_ids = [c[1]["batch_id"] for c in self.mock_target.call_args_list]
        self.assertEqual(sorted(batch_ids), [10, 20])

    def test_single_project_created(self):
        validate_output = _make_custom_chem_multibatch_output()
        upload_custom_reaction(validate_output)
        self.mock_project.assert_called_once()


class TestUploadCombiMultiBatch(TestCase):
    """Combi-custom-chem with targets across different batch tags."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch(
            "backend.tasks.create_batch_model", side_effect=[10, 20]
        )
        self.patcher_target = patch(
            "backend.tasks.create_target_model", side_effect=[100, 200]
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", return_value=9000
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_temp = patch(
            "backend.tasks.get_recipe_stir_temperature", return_value=80.0
        )
        self.patcher_exact = patch("backend.tasks.get_exact_search")

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_temp = self.patcher_temp.start()
        self.mock_exact = self.patcher_exact.start()

    def tearDown(self):
        patch.stopall()

    def test_creates_two_batches(self):
        validate_output = _make_combi_multibatch_output()
        upload_combi_custom_reaction(validate_output)

        self.assertEqual(self.mock_batch.call_count, 2)
        tags = sorted(c[1]["batchtag"] for c in self.mock_batch.call_args_list)
        self.assertEqual(tags, ["batch-X", "batch-Y"])

    def test_each_target_in_correct_batch(self):
        validate_output = _make_combi_multibatch_output()
        upload_combi_custom_reaction(validate_output)

        self.assertEqual(self.mock_target.call_count, 2)
        batch_ids = [c[1]["batch_id"] for c in self.mock_target.call_args_list]
        self.assertEqual(sorted(batch_ids), [10, 20])

    def test_single_project_created(self):
        validate_output = _make_combi_multibatch_output()
        upload_combi_custom_reaction(validate_output)
        self.mock_project.assert_called_once()


# ===========================================================================
# Manifold catalog entry wiring – per-reactant correctness
# ===========================================================================


class TestManifoldCatalogWiring(TestCase):
    """Verify catalog entries from route['molecules'] are wired to the correct reactant."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        # Return distinct IDs per reactant so we can tell them apart
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", side_effect=[9001, 9002]
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        self.patcher_recipe = patch("backend.tasks.recipe_exists", return_value=True)
        self.patcher_intra = patch(
            "backend.tasks.get_recipe_intramolecular", return_value=False
        )
        self.patcher_manifold = patch(
            "backend.tasks.get_manifold_retrosynthesis_batch"
        )

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_recipe = self.patcher_recipe.start()
        self.mock_intra = self.patcher_intra.start()
        self.mock_manifold = self.patcher_manifold.start()

    def tearDown(self):
        patch.stopall()

    def _two_reactant_route(self):
        """Route with 2 reactants, each having a different vendor catalog entry."""
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        encoded_route = {
            "reactions": [
                {
                    "name": "Acylation",
                    "reactantSmiles": [SMILES_REACTANT_A, SMILES_REACTANT_B],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
            "molecules": [
                {
                    "smiles": SMILES_REACTANT_A,
                    "isBuildingBlock": True,
                    "catalogEntries": [
                        {"vendorName": "Sigma", "smiles": SMILES_REACTANT_A}
                    ],
                },
                {
                    "smiles": SMILES_REACTANT_B,
                    "isBuildingBlock": True,
                    "catalogEntries": [
                        {"vendorName": "Enamine", "smiles": SMILES_REACTANT_B}
                    ],
                },
            ],
        }
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, encoded_route]]
        )
        return _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"],
        )

    def test_reactant_a_gets_sigma_catalog(self):
        """Reactant A (id=9001) should get Sigma catalog entry."""
        validate_output = self._two_reactant_route()
        upload_manifold_reaction(validate_output)

        self.mock_catalog.assert_any_call(
            catalog_entry={"vendorName": "Sigma", "smiles": SMILES_REACTANT_A},
            reactant_id=9001,
        )

    def test_reactant_b_gets_enamine_catalog(self):
        """Reactant B (id=9002) should get Enamine catalog entry."""
        validate_output = self._two_reactant_route()
        upload_manifold_reaction(validate_output)

        self.mock_catalog.assert_any_call(
            catalog_entry={"vendorName": "Enamine", "smiles": SMILES_REACTANT_B},
            reactant_id=9002,
        )

    def test_no_cross_wiring(self):
        """Sigma entry should NOT be sent to reactant 9002, and vice versa."""
        validate_output = self._two_reactant_route()
        upload_manifold_reaction(validate_output)

        # Collect all catalog calls that have a reactant_id
        reactant_catalog_calls = [
            c for c in self.mock_catalog.call_args_list
            if "reactant_id" in c[1]
        ]
        for c in reactant_catalog_calls:
            rid = c[1]["reactant_id"]
            entry = c[1].get("catalog_entry", {})
            vendor = entry.get("vendorName", "")
            if rid == 9001 and vendor:
                self.assertEqual(vendor, "Sigma")
            if rid == 9002 and vendor:
                self.assertEqual(vendor, "Enamine")


class TestManifoldCatalogWiringNotEncoded(TestCase):
    """Same catalog wiring test for the not-encoded branch."""

    def setUp(self):
        self.patcher_del = patch("backend.tasks.delete_tmp_file")
        self.patcher_project = patch(
            "backend.tasks.create_project_model", return_value=1
        )
        self.patcher_batch = patch("backend.tasks.create_batch_model", return_value=10)
        self.patcher_target = patch(
            "backend.tasks.create_target_model", return_value=100
        )
        self.patcher_method = patch(
            "backend.tasks.create_method_model", return_value=1000
        )
        self.patcher_reaction = patch(
            "backend.tasks.create_reaction_model", return_value=5000
        )
        self.patcher_product = patch("backend.tasks.create_product_model")
        self.patcher_reactant = patch(
            "backend.tasks.create_reactant_model", side_effect=[9001, 9002]
        )
        self.patcher_catalog = patch("backend.tasks.create_catalog_entry_model")
        self.patcher_prev = patch(
            "backend.tasks.check_previous_reaction_products", return_value=None
        )
        # Not encoded
        self.patcher_recipe = patch("backend.tasks.recipe_exists", return_value=False)
        self.patcher_intra = patch(
            "backend.tasks.get_recipe_intramolecular", return_value=False
        )
        self.patcher_manifold = patch(
            "backend.tasks.get_manifold_retrosynthesis_batch"
        )

        self.mock_del = self.patcher_del.start()
        self.mock_project = self.patcher_project.start()
        self.mock_batch = self.patcher_batch.start()
        self.mock_target = self.patcher_target.start()
        self.mock_method = self.patcher_method.start()
        self.mock_reaction = self.patcher_reaction.start()
        self.mock_product = self.patcher_product.start()
        self.mock_reactant = self.patcher_reactant.start()
        self.mock_catalog = self.patcher_catalog.start()
        self.mock_prev = self.patcher_prev.start()
        self.mock_recipe = self.patcher_recipe.start()
        self.mock_intra = self.patcher_intra.start()
        self.mock_manifold = self.patcher_manifold.start()

    def tearDown(self):
        patch.stopall()

    def _two_reactant_not_encoded_route(self):
        building_block_route = {
            "molecules": [
                {
                    "isBuildingBlock": False,
                    "smiles": SMILES_ASPIRIN,
                    "catalogEntries": [],
                }
            ],
            "reactions": [],
        }
        not_encoded_route = {
            "reactions": [
                {
                    "name": "NovelReaction",
                    "reactantSmiles": [SMILES_REACTANT_A, SMILES_REACTANT_B],
                    "productSmiles": SMILES_PRODUCT_A,
                }
            ],
            "molecules": [
                {
                    "smiles": SMILES_REACTANT_A,
                    "isBuildingBlock": True,
                    "catalogEntries": [
                        {"vendorName": "Alfa", "smiles": SMILES_REACTANT_A}
                    ],
                },
                {
                    "smiles": SMILES_REACTANT_B,
                    "isBuildingBlock": True,
                    "catalogEntries": [
                        {"vendorName": "TCI", "smiles": SMILES_REACTANT_B}
                    ],
                },
            ],
        }
        self.mock_manifold.return_value = _make_manifold_response(
            routes_per_target=[[building_block_route, not_encoded_route]]
        )
        return _make_retro_validate_output(
            smiles=[SMILES_ASPIRIN], names=["Aspirin"],
        )

    def test_reactant_a_gets_alfa_catalog(self):
        validate_output = self._two_reactant_not_encoded_route()
        upload_manifold_reaction(validate_output)

        self.mock_catalog.assert_any_call(
            catalog_entry={"vendorName": "Alfa", "smiles": SMILES_REACTANT_A},
            reactant_id=9001,
        )

    def test_reactant_b_gets_tci_catalog(self):
        validate_output = self._two_reactant_not_encoded_route()
        upload_manifold_reaction(validate_output)

        self.mock_catalog.assert_any_call(
            catalog_entry={"vendorName": "TCI", "smiles": SMILES_REACTANT_B},
            reactant_id=9002,
        )

    def test_no_cross_wiring(self):
        validate_output = self._two_reactant_not_encoded_route()
        upload_manifold_reaction(validate_output)

        reactant_catalog_calls = [
            c for c in self.mock_catalog.call_args_list
            if "reactant_id" in c[1]
        ]
        for c in reactant_catalog_calls:
            rid = c[1]["reactant_id"]
            entry = c[1].get("catalog_entry", {})
            vendor = entry.get("vendorName", "")
            if rid == 9001 and vendor:
                self.assertEqual(vendor, "Alfa")
            if rid == 9002 and vendor:
                self.assertEqual(vendor, "TCI")


# ===========================================================================
# canonicalize_smiles tests
# ===========================================================================


class TestCanonicalizeSmilesList(TestCase):
    """Test canonicalize_smiles with a list of SMILES (no CSV)."""

    @patch("backend.tasks.canon_smiles", side_effect=lambda s: s)
    def test_valid_smiles_returns_validated_true(self, mock_canon):
        validated, result = canonicalize_smiles(smiles=["CCO", "CC"])
        self.assertTrue(validated)

    @patch("backend.tasks.canon_smiles", side_effect=lambda s: s)
    def test_valid_smiles_returns_canonicalized_list(self, mock_canon):
        validated, result = canonicalize_smiles(smiles=["CCO", "CC"])
        self.assertIsInstance(result, list)
        self.assertEqual(len(result), 2)

    def test_invalid_smiles_returns_validated_false(self):
        validated, result = canonicalize_smiles(smiles=["CCO", "INVALID_SMILES_XYZ"])
        self.assertFalse(validated)

    def test_invalid_smiles_returns_error_with_index(self):
        validated, result = canonicalize_smiles(smiles=["CCO", "INVALID_SMILES_XYZ"])
        self.assertIn("1", str(result))  # index 1 is the bad one

    @patch("backend.tasks.canon_smiles", side_effect=lambda s: s)
    def test_single_smiles_valid(self, mock_canon):
        validated, result = canonicalize_smiles(smiles=["CCO"])
        self.assertTrue(validated)
        self.assertEqual(len(result), 1)

    def test_multiple_invalid_smiles_reports_all_indices(self):
        validated, result = canonicalize_smiles(
            smiles=["INVALID_1", "CCO", "INVALID_2"]
        )
        self.assertFalse(validated)
        self.assertIn("0", str(result))
        self.assertIn("2", str(result))


class TestCanonicalizeSmilesCSV(TestCase):
    """Test canonicalize_smiles reading from a CSV file."""

    @patch("backend.tasks.canon_smiles", side_effect=lambda s: s)
    @patch("backend.tasks.delete_tmp_file")
    @patch("backend.tasks.pd.read_csv")
    def test_reads_csv_and_deletes_file(self, mock_csv, mock_del, mock_canon):
        import pandas as pd

        mock_csv.return_value = pd.DataFrame({"SMILES": ["CCO", "CC"]})
        validated, result = canonicalize_smiles(csvfile="/tmp/test_smiles.csv")

        self.assertTrue(validated)
        mock_csv.assert_called_once_with("/tmp/test_smiles.csv", encoding="utf8")
        mock_del.assert_called_once_with("/tmp/test_smiles.csv")

    @patch("backend.tasks.canon_smiles", side_effect=lambda s: s)
    @patch("backend.tasks.delete_tmp_file")
    @patch("backend.tasks.pd.read_csv")
    def test_csv_with_invalid_smiles(self, mock_csv, mock_del, mock_canon):
        import pandas as pd

        mock_csv.return_value = pd.DataFrame({"SMILES": ["CCO", "NOT_A_MOLECULE"]})
        validated, result = canonicalize_smiles(csvfile="/tmp/bad.csv")

        self.assertFalse(validated)
        mock_del.assert_called_once_with("/tmp/bad.csv")


# ===========================================================================
# _get_custom_sm_csv_path tests
# ===========================================================================


class TestGetCustomSmCsvPath(TestCase):
    """Unit tests for the _get_custom_sm_csv_path helper."""

    def test_returns_path_when_key_present(self):
        result = _get_custom_sm_csv_path({"42": "/tmp/sm.csv"}, 42)
        self.assertEqual(result, "/tmp/sm.csv")

    def test_returns_none_when_key_missing(self):
        result = _get_custom_sm_csv_path({"99": "/tmp/sm.csv"}, 42)
        self.assertIsNone(result)

    def test_returns_none_when_dict_is_none(self):
        result = _get_custom_sm_csv_path(None, 42)
        self.assertIsNone(result)

    def test_returns_none_when_dict_is_empty(self):
        result = _get_custom_sm_csv_path({}, 42)
        self.assertIsNone(result)


# ===========================================================================
# _process_reaction_step_sessions tests
# ===========================================================================

OT_STEP_PATCHES = [
    "backend.tasks.create_multiple_ot_sessions",
    "backend.tasks.get_grouped_action_session_types",
    "backend.tasks.get_action_session_types",
    "backend.tasks.get_grouped_action_session_sequences",
    "backend.tasks.get_action_session_sequence_numbers",
    "backend.tasks.get_action_session_query_set",
]


class TestProcessReactionStepSessions(TestCase):
    """Tests for the _process_reaction_step_sessions helper."""

    @patch(*OT_STEP_PATCHES[:1])
    @patch(*OT_STEP_PATCHES[1:2])
    @patch(*OT_STEP_PATCHES[2:3])
    @patch(*OT_STEP_PATCHES[3:4])
    @patch(*OT_STEP_PATCHES[4:5])
    @patch(*OT_STEP_PATCHES[5:6])
    def _run(
        self,
        reaction_ids,
        mock_get_as_qs,
        mock_seq_nums,
        mock_grouped_seqs,
        mock_as_types,
        mock_grouped_types,
        mock_create_sessions,
        robot_count=1,
        human_count=0,
        custom_sm_csv_path=None,
    ):
        """Helper to set up a single-sequence, single-type scenario."""
        # Build a fake queryset-like object with a .filter method
        from unittest.mock import MagicMock

        mock_qs = MagicMock()
        mock_get_as_qs.return_value = mock_qs
        mock_seq_nums.return_value = [1]

        # Build the type-level queryset
        type_qs = MagicMock()
        robot_qs = MagicMock()
        robot_qs.__bool__ = lambda self: robot_count > 0
        human_qs = MagicMock()
        human_qs.__bool__ = lambda self: human_count > 0
        type_qs.filter.side_effect = lambda driver=None: (
            robot_qs if driver == "robot" else human_qs
        )

        # One sequence group containing one type group
        mock_grouped_seqs.return_value = [mock_qs]
        mock_as_types.return_value = ["robot"]
        mock_grouped_types.return_value = [type_qs]

        protocol_obj = MagicMock()
        _process_reaction_step_sessions(
            reaction_ids=reaction_ids,
            reactionstep=1,
            otbatchprotocolobj=protocol_obj,
            custom_sm_csv_path=custom_sm_csv_path,
            batchtag="batch-1",
        )
        return mock_create_sessions, robot_qs, protocol_obj

    def test_calls_create_multiple_ot_sessions_for_robot(self):
        mock_create, robot_qs, protocol_obj = self._run([10, 20])
        mock_create.assert_called_once_with(
            reactionstep=1,
            otbatchprotocolobj=protocol_obj,
            actionsessionqueryset=robot_qs,
            customSMcsvpath=None,
            batchtag="batch-1",
        )

    def test_passes_custom_sm_csv_path(self):
        mock_create, _, _ = self._run(
            [10], custom_sm_csv_path="/tmp/custom.csv"
        )
        self.assertEqual(
            mock_create.call_args.kwargs["customSMcsvpath"], "/tmp/custom.csv"
        )

    def test_skips_when_no_robot_sessions(self):
        mock_create, _, _ = self._run([10], robot_count=0)
        mock_create.assert_not_called()


# ===========================================================================
# create_ot_script tests
# ===========================================================================

# Common patches for create_ot_script tests
OT_SCRIPT_PATCHES = [
    "backend.tasks.ZipOTBatchProtocol",
    "backend.tasks._process_reaction_step_sessions",
    "backend.tasks.group_reactions",
    "backend.tasks.get_max_reaction_number",
    "backend.tasks._get_custom_sm_csv_path",
    "backend.tasks.get_ot_batch_protocol_query_set",
    "backend.tasks.get_batch_tag",
    "backend.tasks.CreateEncodedActionModels",
    "backend.tasks.get_recipe_intramolecular",
    "backend.tasks.get_action_session_query_set",
    "backend.tasks.get_batch_reactions",
    "backend.tasks.get_reactions_to_do",
    "backend.tasks.current_task",
    "backend.tasks.Batch",
    "backend.tasks.OTProject",
    "backend.tasks.OTBatchProtocol",
]


def _make_reaction_obj(
    reaction_id=1,
    reactionclass="Amidation",
    recipe="standard",
    target_id=100,
    reactant_smiles=None,
):
    """Build a mock reaction ORM object."""
    from unittest.mock import MagicMock

    if reactant_smiles is None:
        reactant_smiles = [SMILES_REACTANT_A, SMILES_REACTANT_B]
    obj = MagicMock()
    obj.id = reaction_id
    obj.reactionclass = reactionclass
    obj.recipe = recipe
    obj.method_id.target_id.id = target_id
    obj.reactants.all.return_value.values_list.return_value = reactant_smiles
    return obj


@patch("backend.tasks.OTBatchProtocol")
@patch("backend.tasks.OTProject")
@patch("backend.tasks.Batch")
@patch("backend.tasks.current_task")
@patch("backend.tasks.get_reactions_to_do")
@patch("backend.tasks.get_batch_reactions")
@patch("backend.tasks.get_action_session_query_set")
@patch("backend.tasks.get_recipe_intramolecular")
@patch("backend.tasks.CreateEncodedActionModels")
@patch("backend.tasks.get_batch_tag")
@patch("backend.tasks.get_ot_batch_protocol_query_set")
@patch("backend.tasks._get_custom_sm_csv_path")
@patch("backend.tasks.get_max_reaction_number")
@patch("backend.tasks.group_reactions")
@patch("backend.tasks._process_reaction_step_sessions")
@patch("backend.tasks.ZipOTBatchProtocol")
class TestCreateOTScriptHappyPath(TestCase):
    """Single-batch, single-step, action sessions already exist."""

    def _setup_mocks(
        self,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        # OTProject mock
        mock_otproject = mock_OTProject.return_value
        mock_otproject.id = 500

        # Batch
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch

        # current_task
        mock_current_task.request.id = "celery-task-id-123"

        # Reactions already have action sessions
        rxn = _make_reaction_obj(reaction_id=10)
        mock_get_batch_rxns.return_value = [rxn]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: True)

        # Protocol doesn't exist yet
        mock_get_protocol.return_value = []

        # Batch tag
        mock_get_tag.return_value = "batch-1"

        # Grouping: single step with one reaction
        mock_max_rxn.return_value = 1
        mock_group_rxns.return_value = [[rxn]]

        # Custom SM path
        mock_custom_sm.return_value = None

        # Zip succeeds
        mock_zip.return_value.errors = []

        return mock_otproject

    def test_returns_task_summary_and_project_id(self, *mocks):
        mock_otproject = self._setup_mocks(*mocks)
        result = create_ot_script(batchids=[1], protocol_name="P1")
        task_summary, project_id = result
        self.assertEqual(project_id, 500)
        self.assertIn(1, task_summary)

    def test_creates_ot_project(self, *mocks):
        self._setup_mocks(*mocks)
        mock_OTProject = mocks[-2]
        create_ot_script(batchids=[1], protocol_name="P1")
        mock_OTProject.return_value.save.assert_called()
        self.assertEqual(mock_OTProject.return_value.name, "P1")

    def test_creates_batch_protocol(self, *mocks):
        self._setup_mocks(*mocks)
        mock_OTBatchProtocol = mocks[-1]
        create_ot_script(batchids=[1], protocol_name="P1")
        mock_OTBatchProtocol.return_value.save.assert_called()

    def test_calls_process_reaction_step(self, *mocks):
        self._setup_mocks(*mocks)
        mock_process_step = mocks[1]
        create_ot_script(batchids=[1], protocol_name="P1")
        mock_process_step.assert_called_once()
        kwargs = mock_process_step.call_args.kwargs
        self.assertEqual(kwargs["reactionstep"], 1)
        self.assertEqual(kwargs["reaction_ids"], [10])
        self.assertEqual(kwargs["batchtag"], "batch-1")

    def test_creates_zip_protocol(self, *mocks):
        self._setup_mocks(*mocks)
        mock_zip = mocks[0]
        create_ot_script(batchids=[1], protocol_name="P1")
        mock_zip.assert_called_once()

    def test_does_not_call_create_encoded_when_actions_exist(self, *mocks):
        self._setup_mocks(*mocks)
        mock_create_encoded = mocks[8]
        create_ot_script(batchids=[1], protocol_name="P1")
        mock_create_encoded.assert_not_called()

    def test_task_summary_true_when_zip_no_errors(self, *mocks):
        self._setup_mocks(*mocks)
        mock_zip = mocks[0]
        mock_zip.return_value.errors = []
        result = create_ot_script(batchids=[1], protocol_name="P1")
        self.assertTrue(result[0][1])


@patch("backend.tasks.OTBatchProtocol")
@patch("backend.tasks.OTProject")
@patch("backend.tasks.Batch")
@patch("backend.tasks.current_task")
@patch("backend.tasks.get_reactions_to_do")
@patch("backend.tasks.get_batch_reactions")
@patch("backend.tasks.get_action_session_query_set")
@patch("backend.tasks.get_recipe_intramolecular")
@patch("backend.tasks.CreateEncodedActionModels")
@patch("backend.tasks.get_batch_tag")
@patch("backend.tasks.get_ot_batch_protocol_query_set")
@patch("backend.tasks._get_custom_sm_csv_path")
@patch("backend.tasks.get_max_reaction_number")
@patch("backend.tasks.group_reactions")
@patch("backend.tasks._process_reaction_step_sessions")
@patch("backend.tasks.ZipOTBatchProtocol")
class TestCreateOTScriptNoActionSessions(TestCase):
    """When action sessions don't exist, CreateEncodedActionModels is called."""

    def test_creates_encoded_action_models(
        self,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn = _make_reaction_obj(reaction_id=10, reactionclass="Amidation")
        mock_get_batch_rxns.return_value = [rxn]
        # No action sessions exist yet
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: False)
        mock_get_protocol.return_value = []
        mock_get_tag.return_value = "batch-1"
        mock_max_rxn.return_value = 1
        mock_group_rxns.return_value = [[rxn]]
        mock_custom_sm.return_value = None
        mock_zip.return_value.errors = []
        mock_intramolecular.return_value = False

        create_ot_script(batchids=[1], protocol_name="P1")

        mock_create_encoded.assert_called_once_with(
            reaction_class="Amidation",
            recipe_name="standard",
            intramolecular=False,
            target_id=100,
            reaction_id=10,
            reactant_pair_smiles=[SMILES_REACTANT_A, SMILES_REACTANT_B],
        )

    def test_intramolecular_true_for_single_reactant(
        self,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn = _make_reaction_obj(
            reaction_id=10,
            reactionclass="Cyclization",
            reactant_smiles=[SMILES_REACTANT_A],  # single reactant
        )
        mock_get_batch_rxns.return_value = [rxn]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: False)
        mock_get_protocol.return_value = []
        mock_get_tag.return_value = "batch-1"
        mock_max_rxn.return_value = 1
        mock_group_rxns.return_value = [[rxn]]
        mock_custom_sm.return_value = None
        mock_zip.return_value.errors = []
        mock_intramolecular.return_value = True

        create_ot_script(batchids=[1], protocol_name="P1")

        # Single reactant + intramolecular recipe → intramolecular=True
        enc_call = mock_create_encoded.call_args
        self.assertTrue(enc_call.kwargs["intramolecular"])


@patch("backend.tasks.OTBatchProtocol")
@patch("backend.tasks.OTProject")
@patch("backend.tasks.Batch")
@patch("backend.tasks.current_task")
@patch("backend.tasks.get_reactions_to_do")
@patch("backend.tasks.get_batch_reactions")
@patch("backend.tasks.get_action_session_query_set")
@patch("backend.tasks.get_recipe_intramolecular")
@patch("backend.tasks.CreateEncodedActionModels")
@patch("backend.tasks.get_batch_tag")
@patch("backend.tasks.get_ot_batch_protocol_query_set")
@patch("backend.tasks._get_custom_sm_csv_path")
@patch("backend.tasks.get_max_reaction_number")
@patch("backend.tasks.group_reactions")
@patch("backend.tasks._process_reaction_step_sessions")
@patch("backend.tasks.ZipOTBatchProtocol")
class TestCreateOTScriptExistingProtocol(TestCase):
    """When a completed protocol already exists, it should be re-used."""

    def test_reuses_existing_protocol(
        self,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn = _make_reaction_obj(reaction_id=10)
        mock_get_batch_rxns.return_value = [rxn]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: True)
        mock_get_tag.return_value = "batch-1"

        # Protocol already exists with a zipfile
        existing_protocol = MagicMock()
        existing_protocol.zipfile = "/existing/zip.zip"
        mock_get_protocol.return_value = [existing_protocol]

        result = create_ot_script(batchids=[1], protocol_name="P1")

        # Should NOT create a new protocol or call step processing
        mock_process_step.assert_not_called()
        mock_zip.assert_not_called()
        mock_group_rxns.assert_not_called()

        # Should mark as True
        self.assertTrue(result[0][1])

        # Existing protocol should be updated and saved
        existing_protocol.save.assert_called()

    def test_no_new_otbatchprotocol_created(
        self,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn = _make_reaction_obj()
        mock_get_batch_rxns.return_value = [rxn]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: True)
        mock_get_tag.return_value = "batch-1"

        existing_protocol = MagicMock()
        existing_protocol.zipfile = "/existing/zip.zip"
        mock_get_protocol.return_value = [existing_protocol]

        create_ot_script(batchids=[1], protocol_name="P1")

        # The OTBatchProtocol constructor should only be called for the
        # class-level mock setup, not for creating a new instance
        # (the function uses the existing one, not OTBatchProtocol())
        mock_OTBatchProtocol.assert_not_called()


@patch("backend.tasks.OTBatchProtocol")
@patch("backend.tasks.OTProject")
@patch("backend.tasks.Batch")
@patch("backend.tasks.current_task")
@patch("backend.tasks.get_reactions_to_do")
@patch("backend.tasks.get_batch_reactions")
@patch("backend.tasks.get_action_session_query_set")
@patch("backend.tasks.get_recipe_intramolecular")
@patch("backend.tasks.CreateEncodedActionModels")
@patch("backend.tasks.get_batch_tag")
@patch("backend.tasks.get_ot_batch_protocol_query_set")
@patch("backend.tasks._get_custom_sm_csv_path")
@patch("backend.tasks.get_max_reaction_number")
@patch("backend.tasks.group_reactions")
@patch("backend.tasks._process_reaction_step_sessions")
@patch("backend.tasks.ZipOTBatchProtocol")
class TestCreateOTScriptMultiStep(TestCase):
    """Two-step reaction: step 1 always runs, step 2 uses QC filter."""

    def _setup(
        self,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
        step2_reactions_remain=True,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn1 = _make_reaction_obj(reaction_id=10)
        rxn2 = _make_reaction_obj(reaction_id=20)
        mock_get_batch_rxns.return_value = [rxn1, rxn2]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: True)
        mock_get_protocol.return_value = []
        mock_get_tag.return_value = "batch-1"
        mock_max_rxn.return_value = 2
        mock_custom_sm.return_value = None
        mock_zip.return_value.errors = []

        # Two steps: step 0 has rxn1, step 1 has rxn2
        step1_group = [rxn1]
        step2_group = [rxn2]
        mock_group_rxns.return_value = [step1_group, step2_group]

        if step2_reactions_remain:
            mock_rxns_to_do.return_value = [rxn2]
        else:
            mock_rxns_to_do.return_value = []

        return mock_process_step, mock_rxns_to_do

    def test_both_steps_processed(self, *mocks):
        mock_process_step, _ = self._setup(*mocks)
        create_ot_script(batchids=[1], protocol_name="P1")
        self.assertEqual(mock_process_step.call_count, 2)

        # Step 1: reactionstep=1, reaction_ids=[10]
        call1 = mock_process_step.call_args_list[0]
        self.assertEqual(call1.kwargs["reactionstep"], 1)
        self.assertEqual(call1.kwargs["reaction_ids"], [10])

        # Step 2: reactionstep=2, reaction_ids=[20]
        call2 = mock_process_step.call_args_list[1]
        self.assertEqual(call2.kwargs["reactionstep"], 2)
        self.assertEqual(call2.kwargs["reaction_ids"], [20])

    def test_qc_filter_applied_only_to_step_2(self, *mocks):
        """get_reactions_to_do is called for step 2 but not step 1."""
        _, mock_rxns_to_do = self._setup(*mocks)
        create_ot_script(batchids=[1], protocol_name="P1")
        mock_rxns_to_do.assert_called_once()

    def test_step2_skipped_when_no_reactions_pass_qc(self, *mocks):
        mock_process_step, mock_rxns_to_do = self._setup(
            *mocks, step2_reactions_remain=False
        )
        create_ot_script(batchids=[1], protocol_name="P1")

        # Only step 1 processed
        self.assertEqual(mock_process_step.call_count, 1)
        call1 = mock_process_step.call_args_list[0]
        self.assertEqual(call1.kwargs["reactionstep"], 1)

    def test_three_steps_breaks_at_empty_step2(self, *mocks):
        """If step 2 has no reactions to do, step 3 never runs."""
        mock_process_step = mocks[1]
        mock_group_rxns = mocks[2]
        mock_rxns_to_do = mocks[11]  # get_reactions_to_do

        # Build 3-step scenario
        mock_OTProject = mocks[-2]
        mock_Batch = mocks[-3]
        mock_current_task = mocks[-4]

        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn1 = _make_reaction_obj(reaction_id=10)
        rxn2 = _make_reaction_obj(reaction_id=20)
        rxn3 = _make_reaction_obj(reaction_id=30)
        mocks[10].return_value = [rxn1, rxn2, rxn3]  # get_batch_reactions
        mocks[7].return_value = MagicMock(__bool__=lambda self: True)  # get_as_qs
        mocks[5].return_value = []  # get_ot_batch_protocol_query_set
        mocks[6].return_value = "batch-1"  # get_batch_tag
        mocks[3].return_value = 3  # get_max_reaction_number
        mocks[4].return_value = None  # custom_sm_csv_path
        mocks[0].return_value.errors = []  # ZipOTBatchProtocol

        mock_group_rxns.return_value = [[rxn1], [rxn2], [rxn3]]
        mock_rxns_to_do.return_value = []  # step 2 fails QC → empty

        create_ot_script(batchids=[1], protocol_name="P1")

        # Only step 1 runs; step 2 QC check returns empty → break
        self.assertEqual(mock_process_step.call_count, 1)


@patch("backend.tasks.OTBatchProtocol")
@patch("backend.tasks.OTProject")
@patch("backend.tasks.Batch")
@patch("backend.tasks.current_task")
@patch("backend.tasks.get_reactions_to_do")
@patch("backend.tasks.get_batch_reactions")
@patch("backend.tasks.get_action_session_query_set")
@patch("backend.tasks.get_recipe_intramolecular")
@patch("backend.tasks.CreateEncodedActionModels")
@patch("backend.tasks.get_batch_tag")
@patch("backend.tasks.get_ot_batch_protocol_query_set")
@patch("backend.tasks._get_custom_sm_csv_path")
@patch("backend.tasks.get_max_reaction_number")
@patch("backend.tasks.group_reactions")
@patch("backend.tasks._process_reaction_step_sessions")
@patch("backend.tasks.ZipOTBatchProtocol")
class TestCreateOTScriptTempFileCleanup(TestCase):
    """Custom SM files should be deleted after the transaction."""

    @patch("backend.tasks.os.path.exists", return_value=True)
    @patch("backend.tasks.os.remove")
    def test_custom_sm_files_cleaned_up(
        self,
        mock_remove,
        mock_exists,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn = _make_reaction_obj()
        mock_get_batch_rxns.return_value = [rxn]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: True)
        mock_get_protocol.return_value = []
        mock_get_tag.return_value = "batch-1"
        mock_max_rxn.return_value = 1
        mock_group_rxns.return_value = [[rxn]]
        mock_custom_sm.return_value = "/tmp/custom.csv"
        mock_zip.return_value.errors = []

        create_ot_script(
            batchids=[1],
            protocol_name="P1",
            custom_SM_files={"1": "/tmp/custom.csv"},
        )

        mock_remove.assert_called_once_with("/tmp/custom.csv")

    @patch("backend.tasks.os.path.exists", return_value=False)
    @patch("backend.tasks.os.remove")
    def test_no_remove_when_file_missing(
        self,
        mock_remove,
        mock_exists,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn = _make_reaction_obj()
        mock_get_batch_rxns.return_value = [rxn]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: True)
        mock_get_protocol.return_value = []
        mock_get_tag.return_value = "batch-1"
        mock_max_rxn.return_value = 1
        mock_group_rxns.return_value = [[rxn]]
        mock_custom_sm.return_value = None
        mock_zip.return_value.errors = []

        create_ot_script(
            batchids=[1],
            protocol_name="P1",
            custom_SM_files={"1": "/tmp/gone.csv"},
        )

        mock_remove.assert_not_called()


@patch("backend.tasks.OTBatchProtocol")
@patch("backend.tasks.OTProject")
@patch("backend.tasks.Batch")
@patch("backend.tasks.current_task")
@patch("backend.tasks.get_reactions_to_do")
@patch("backend.tasks.get_batch_reactions")
@patch("backend.tasks.get_action_session_query_set")
@patch("backend.tasks.get_recipe_intramolecular")
@patch("backend.tasks.CreateEncodedActionModels")
@patch("backend.tasks.get_batch_tag")
@patch("backend.tasks.get_ot_batch_protocol_query_set")
@patch("backend.tasks._get_custom_sm_csv_path")
@patch("backend.tasks.get_max_reaction_number")
@patch("backend.tasks.group_reactions")
@patch("backend.tasks._process_reaction_step_sessions")
@patch("backend.tasks.ZipOTBatchProtocol")
class TestCreateOTScriptZipErrors(TestCase):
    """When ZipOTBatchProtocol reports errors, task_summary should be False."""

    def test_task_summary_false_on_zip_errors(
        self,
        mock_zip,
        mock_process_step,
        mock_group_rxns,
        mock_max_rxn,
        mock_custom_sm,
        mock_get_protocol,
        mock_get_tag,
        mock_create_encoded,
        mock_intramolecular,
        mock_get_as_qs,
        mock_get_batch_rxns,
        mock_rxns_to_do,
        mock_current_task,
        mock_Batch,
        mock_OTProject,
        mock_OTBatchProtocol,
    ):
        mock_OTProject.return_value.id = 500
        mock_batch = MagicMock()
        mock_batch.project_id = MagicMock()
        mock_Batch.objects.get.return_value = mock_batch
        mock_current_task.request.id = "task-1"

        rxn = _make_reaction_obj()
        mock_get_batch_rxns.return_value = [rxn]
        mock_get_as_qs.return_value = MagicMock(__bool__=lambda self: True)
        mock_get_protocol.return_value = []
        mock_get_tag.return_value = "batch-1"
        mock_max_rxn.return_value = 1
        mock_group_rxns.return_value = [[rxn]]
        mock_custom_sm.return_value = None

        # Zip reports errors
        mock_zip.return_value.errors = [
            {"function": "someFunc", "errorwarning": "something went wrong"},
        ]

        result = create_ot_script(batchids=[1], protocol_name="P1")
        # Errors list is truthy → not errors = False
        self.assertFalse(result[0][1])


# =========================================================================
# Tests for ZipOTBatchProtocol (direct class tests)
# =========================================================================


def _make_file_obj(name):
    """Return a mock with a .name attribute (mimics Django FileField)."""
    obj = MagicMock()
    obj.name = name
    return obj


def _make_artifact_obj(field_name, file_path):
    """Return a mock ORM object whose *field_name* attr is a FileField-like."""
    obj = MagicMock()
    file_field = MagicMock()
    file_field.name = file_path
    setattr(obj, field_name, file_field)
    return obj


@patch("backend.tasks.os.remove")
@patch("builtins.open", MagicMock())
@patch("backend.tasks.ZipFile")
@patch("backend.tasks.OTScript.objects.filter")
@patch("backend.tasks.CompoundOrder.objects.filter")
@patch("backend.tasks.SolventPrep.objects.filter")
@patch("backend.tasks.OTSession.objects.filter")
@patch("backend.tasks.settings")
class TestZipOTBatchProtocolHappyPath(TestCase):
    """Full construction with one session that has one of each artifact type."""

    def _setup(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"

        session = MagicMock()
        mock_ot_filter.return_value = [session]

        sp = _make_artifact_obj("solventprepcsv", "solventprep/sp1.csv")
        co = _make_artifact_obj("ordercsv", "compoundorders/co1.csv")
        ot = _make_artifact_obj("otscript", "otscripts/run1.py")

        def related_side_effect(otsession_id=None):
            return None  # overridden per-model below

        mock_sp_filter.return_value = [sp]
        mock_co_filter.return_value = [co]
        mock_os_filter.return_value = [ot]

        return mock_zipfile, mock_remove

    def test_no_errors_on_success(self, *args):
        self._setup(*args)
        proto_obj = MagicMock()
        z = ZipOTBatchProtocol(otbatchprotocolobj=proto_obj, batchtag="B1")
        self.assertEqual(z.errors, [])

    def test_zipfn_uses_batchtag(self, *args):
        self._setup(*args)
        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="X42")
        self.assertEqual(z.zipfn, "batch-X42-protocol.zip")

    def test_zip_write_called_three_times(self, *args):
        mock_zipfile, _ = self._setup(*args)
        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        # One artifact per type → 3 write calls
        zip_instance = mock_zipfile.return_value.__enter__.return_value
        self.assertEqual(zip_instance.write.call_count, 3)

    def test_zip_arcnames_use_correct_subdirs(self, *args):
        mock_zipfile, _ = self._setup(*args)
        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        zip_instance = mock_zipfile.return_value.__enter__.return_value
        arcnames = [c.kwargs["arcname"] for c in zip_instance.write.call_args_list]
        self.assertIn("solventprep/sp1.csv", arcnames)
        self.assertIn("compoundorders/co1.csv", arcnames)
        self.assertIn("otscripts/run1.py", arcnames)

    def test_zip_filenames_use_mediaroot(self, *args):
        mock_zipfile, _ = self._setup(*args)
        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        zip_instance = mock_zipfile.return_value.__enter__.return_value
        filenames = [c.kwargs["filename"] for c in zip_instance.write.call_args_list]
        for fn in filenames:
            self.assertTrue(fn.startswith("/media/"))

    def test_saves_zip_to_model(self, *args):
        self._setup(*args)
        proto_obj = MagicMock()
        z = ZipOTBatchProtocol(otbatchprotocolobj=proto_obj, batchtag="B1")
        proto_obj.zipfile.save.assert_called_once()
        proto_obj.save.assert_called_once()

    def test_deletes_tmp_zip(self, *args):
        _, mock_remove = self._setup(*args)
        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        mock_remove.assert_called_once_with("/media/tmp/batchprotocoltmp.zip")


@patch("backend.tasks.os.remove")
@patch("builtins.open", MagicMock())
@patch("backend.tasks.ZipFile")
@patch("backend.tasks.OTScript.objects.filter")
@patch("backend.tasks.CompoundOrder.objects.filter")
@patch("backend.tasks.SolventPrep.objects.filter")
@patch("backend.tasks.OTSession.objects.filter")
@patch("backend.tasks.settings")
class TestZipOTBatchProtocolNoSessions(TestCase):
    """When no OT sessions exist, zip is not created and a warning is logged."""

    def test_records_warning_when_no_sessions(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"
        mock_ot_filter.return_value = []  # empty queryset

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        self.assertEqual(len(z.errors), 1)
        self.assertIn("_get_ot_sessions", z.errors[0]["function"])

    def test_no_zip_created_when_no_sessions(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"
        mock_ot_filter.return_value = []

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        mock_zipfile.assert_not_called()
        mock_remove.assert_not_called()


@patch("backend.tasks.os.remove")
@patch("builtins.open", MagicMock())
@patch("backend.tasks.ZipFile")
@patch("backend.tasks.OTScript.objects.filter")
@patch("backend.tasks.CompoundOrder.objects.filter")
@patch("backend.tasks.SolventPrep.objects.filter")
@patch("backend.tasks.OTSession.objects.filter")
@patch("backend.tasks.settings")
class TestZipOTBatchProtocolMissingArtifacts(TestCase):
    """When some artifact types are missing, warnings are logged but
    other artifacts still get zipped."""

    def test_missing_solventprep_logs_warning(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"
        session = MagicMock()
        mock_ot_filter.return_value = [session]
        mock_sp_filter.return_value = []  # no solvent preps
        mock_co_filter.return_value = [
            _make_artifact_obj("ordercsv", "compoundorders/co1.csv")
        ]
        mock_os_filter.return_value = [
            _make_artifact_obj("otscript", "otscripts/run1.py")
        ]

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        warning_funcs = [e["function"] for e in z.errors]
        self.assertTrue(
            any("SolventPrep" in f for f in warning_funcs),
            f"Expected SolventPrep warning, got: {warning_funcs}",
        )
        # Other two artifacts still written
        zip_instance = mock_zipfile.return_value.__enter__.return_value
        self.assertEqual(zip_instance.write.call_count, 2)

    def test_all_artifacts_missing_logs_three_warnings(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"
        session = MagicMock()
        mock_ot_filter.return_value = [session]
        mock_sp_filter.return_value = []
        mock_co_filter.return_value = []
        mock_os_filter.return_value = []

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        self.assertEqual(len(z.errors), 3)
        zip_instance = mock_zipfile.return_value.__enter__.return_value
        zip_instance.write.assert_not_called()


@patch("backend.tasks.os.remove")
@patch("builtins.open", MagicMock())
@patch("backend.tasks.ZipFile")
@patch("backend.tasks.OTScript.objects.filter")
@patch("backend.tasks.CompoundOrder.objects.filter")
@patch("backend.tasks.SolventPrep.objects.filter")
@patch("backend.tasks.OTSession.objects.filter")
@patch("backend.tasks.settings")
class TestZipOTBatchProtocolMultipleSessions(TestCase):
    """Multiple OT sessions each contribute their own artifacts."""

    def test_two_sessions_both_processed(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"
        s1, s2 = MagicMock(), MagicMock()
        mock_ot_filter.return_value = [s1, s2]

        # Only OTScript for each session; others empty
        mock_sp_filter.return_value = []
        mock_co_filter.return_value = []
        ot1 = _make_artifact_obj("otscript", "otscripts/run1.py")
        ot2 = _make_artifact_obj("otscript", "otscripts/run2.py")
        mock_os_filter.side_effect = [[ot1], [ot2]]

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        zip_instance = mock_zipfile.return_value.__enter__.return_value
        # 2 otscripts written (SolventPrep + CompoundOrder empty → warnings only)
        self.assertEqual(zip_instance.write.call_count, 2)

    def test_multiple_artifacts_per_type(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"
        session = MagicMock()
        mock_ot_filter.return_value = [session]

        sp1 = _make_artifact_obj("solventprepcsv", "solventprep/sp1.csv")
        sp2 = _make_artifact_obj("solventprepcsv", "solventprep/sp2.csv")
        mock_sp_filter.return_value = [sp1, sp2]
        mock_co_filter.return_value = []
        mock_os_filter.return_value = []

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        zip_instance = mock_zipfile.return_value.__enter__.return_value
        # 2 solvent preps written
        self.assertEqual(zip_instance.write.call_count, 2)


@patch("backend.tasks.os.remove")
@patch("builtins.open", MagicMock())
@patch("backend.tasks.ZipFile")
@patch("backend.tasks.OTScript.objects.filter")
@patch("backend.tasks.CompoundOrder.objects.filter")
@patch("backend.tasks.SolventPrep.objects.filter")
@patch("backend.tasks.OTSession.objects.filter")
@patch("backend.tasks.settings")
class TestZipOTBatchProtocolAddWarning(TestCase):
    """Test the add_warning method directly."""

    def test_add_warning_appends_to_errors(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        mock_settings.MEDIA_ROOT = "/media"
        mock_ot_filter.return_value = []  # triggers a warning via _get_ot_sessions

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        # Already has 1 warning from no sessions; add another manually
        z.add_warning(function="test_func", errorwarning="test warning")
        self.assertEqual(len(z.errors), 2)
        self.assertEqual(z.errors[1]["function"], "test_func")
        self.assertEqual(z.errors[1]["errorwarning"], "test warning")

    def test_errors_is_falsy_when_no_warnings(
        self,
        mock_settings,
        mock_ot_filter,
        mock_sp_filter,
        mock_co_filter,
        mock_os_filter,
        mock_zipfile,
        mock_remove,
    ):
        """Validates the bug fix: empty errors list is falsy (was a truthy
        dict before refactoring)."""
        mock_settings.MEDIA_ROOT = "/media"
        session = MagicMock()
        mock_ot_filter.return_value = [session]
        mock_sp_filter.return_value = [
            _make_artifact_obj("solventprepcsv", "solventprep/sp1.csv")
        ]
        mock_co_filter.return_value = [
            _make_artifact_obj("ordercsv", "compoundorders/co1.csv")
        ]
        mock_os_filter.return_value = [
            _make_artifact_obj("otscript", "otscripts/run1.py")
        ]

        z = ZipOTBatchProtocol(otbatchprotocolobj=MagicMock(), batchtag="B1")
        self.assertFalse(z.errors)  # [] is falsy
        self.assertTrue(not z.errors)  # mirrors create_ot_script check
