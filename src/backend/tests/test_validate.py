"""Tests for the ValidateFile class in backend/validate.py

Uses unittest.TestCase for all tests — no Django DB required because
chemistry DB calls (get_recipe_smarts) are mocked.

Three upload types are tested:
  - retro-API          (postera-ver1.5.csv format)
  - custom-chem        (custom-chem-ver1.5.csv format)
  - combi-custom-chem  (custom-combi-batch-ver1.5.csv format)
"""

import io
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import pandas as pd
from rdkit import Chem

from backend.validate import ValidateFile


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"

# ---------------------------------------------------------------------------
# Real SMILES data taken from the actual test fixture CSVs
# ---------------------------------------------------------------------------

# postera-ver1.5.csv
POSTERA_SMILES = [
    "Fc1ccc(Oc2ccc(F)c(NCC3CCOCC3)n2)cc1",
    "Nc1nc(Oc2ccc(F)cc2)ccc1F",
    "Fc1ccc(Oc2ccc(F)cn2)cc1",
]
POSTERA_NAMES = ["Frag-1", "Frag-2", "Frag-3"]
POSTERA_CONCENTRATIONS = [100.0, 100.0, 100.0]
POSTERA_AMOUNTS = [100.0, 100.0, 100.0]
POSTERA_BATCHTAGS = ["to-check", "to-check", "to-check"]

# custom-chem-ver1.5.csv  (3-step synthesis, single row)
CUSTOM_CHEM_REACTANT_1_STEP1 = "CC(C)(C)OC(=O)N1CCCNCC1"
CUSTOM_CHEM_REACTANT_2_STEP1 = "Cc1cc(Cl)nc(C)n1"
CUSTOM_CHEM_PRODUCT_STEP1 = "Cc1cc(N2CCCN(C(=O)OC(C)(C)C)CC2)nc(C)n1"
CUSTOM_CHEM_REACTION_NAME_STEP1 = "N-nucleophilic aromatic substitution"
CUSTOM_CHEM_RECIPE_STEP1 = "standard-NMP"

CUSTOM_CHEM_REACTANT_1_STEP2 = "Cc1cc(N2CCCN(C(=O)OC(C)(C)C)CC2)nc(C)n1"
CUSTOM_CHEM_PRODUCT_STEP2 = "Cc1cc(N2CCCNCC2)nc(C)n1"
CUSTOM_CHEM_REACTION_NAME_STEP2 = "Boc deprotection"
CUSTOM_CHEM_RECIPE_STEP2 = "standard"

CUSTOM_CHEM_REACTANT_1_STEP3 = "Cc1cc(N2CCCNCC2)nc(C)n1"
CUSTOM_CHEM_REACTANT_2_STEP3 = "COc1ccccc1C(C)=O"
CUSTOM_CHEM_PRODUCT_STEP3 = "COc1ccccc1C(C)N1CCCN(c2cc(C)nc(C)n2)CC1"
CUSTOM_CHEM_REACTION_NAME_STEP3 = "Reductive amination"
CUSTOM_CHEM_RECIPE_STEP3 = "standard"


# ---------------------------------------------------------------------------
# Helpers to create in-memory CSV buffers
# ---------------------------------------------------------------------------


def _make_csv(data: dict) -> io.BytesIO:
    """Create a BytesIO CSV from a dict of columns."""
    df = pd.DataFrame(data)
    buf = io.BytesIO()
    df.to_csv(buf, index=False)
    buf.seek(0)
    return buf


def _make_retro_csv(
    smiles=None,
    names=None,
    concentrations=None,
    amounts=None,
    batch_tags=None,
):
    """Build a retro-API CSV with defaults from the postera fixture."""
    return _make_csv(
        {
            "target-SMILES": smiles or POSTERA_SMILES,
            "target-names": names or POSTERA_NAMES,
            "concentration-required-mM": concentrations or POSTERA_CONCENTRATIONS,
            "amount-required-uL": amounts or POSTERA_AMOUNTS,
            "batch-tag": batch_tags or POSTERA_BATCHTAGS,
        }
    )


def _make_custom_chem_csv_1step(
    target_names=None,
    no_steps=None,
    concentrations=None,
    amounts=None,
    batch_tags=None,
    reactant_1=None,
    reactant_2=None,
    product_smiles=None,
    reaction_names=None,
    reaction_recipes=None,
    groupby_columns=None,
    n_rows=2,
):
    """Build a single-step custom-chem CSV (includes groupby column)."""
    target_names = target_names or [f"target-{i}" for i in range(n_rows)]
    no_steps = no_steps or [1] * n_rows
    concentrations = concentrations or [10.0] * n_rows
    amounts = amounts or [100.0] * n_rows
    batch_tags = batch_tags or ["batch-1"] * n_rows
    reactant_1 = reactant_1 or ["CC(=O)Cl"] * n_rows
    reactant_2 = reactant_2 or ["NCC1CCOCC1"] * n_rows
    product_smiles = product_smiles or ["CC(=O)NCC1CCOCC1"] * n_rows
    reaction_names = reaction_names or ["Amidation"] * n_rows
    reaction_recipes = reaction_recipes or ["standard"] * n_rows
    groupby_columns = groupby_columns or [True] * n_rows

    return _make_csv(
        {
            "target-names": target_names,
            "no-steps": no_steps,
            "concentration-required-mM": concentrations,
            "amount-required-uL": amounts,
            "batch-tag": batch_tags,
            "reaction-groupby-column-1": groupby_columns,
            "reactant-1-1": reactant_1,
            "reactant-2-1": reactant_2,
            "reaction-product-smiles-1": product_smiles,
            "reaction-name-1": reaction_names,
            "reaction-recipe-1": reaction_recipes,
        }
    )


def _make_combi_csv_1step(
    combi_groups=None,
    no_steps=None,
    concentrations=None,
    amounts=None,
    batch_tags=None,
    reactant_1=None,
    reactant_2=None,
    reaction_names=None,
    reaction_recipes=None,
):
    """Build a single-step combi-custom-chem CSV."""
    combi_groups = combi_groups or [1, 1]
    no_steps = no_steps or [1, 1]
    concentrations = concentrations or [10.0, 10.0]
    amounts = amounts or [100.0, 100.0]
    batch_tags = batch_tags or ["batch-1", "batch-1"]
    reactant_1 = reactant_1 or ["CC(=O)Cl", ""]
    reactant_2 = reactant_2 or ["NCC1CCOCC1", "NCc1ccccc1"]
    reaction_names = reaction_names or ["Amidation", "Amidation"]
    reaction_recipes = reaction_recipes or ["standard", "standard"]

    return _make_csv(
        {
            "combi-group": combi_groups,
            "no-steps": no_steps,
            "concentration-required-mM": concentrations,
            "amount-required-uL": amounts,
            "batch-tag": batch_tags,
            "reactant-1-1": reactant_1,
            "reactant-2-1": reactant_2,
            "reaction-name-1": reaction_names,
            "reaction-recipe-1": reaction_recipes,
        }
    )


def _open_fixture(filename: str) -> io.BytesIO:
    """Open a fixture CSV as BytesIO for ValidateFile consumption."""
    path = FIXTURES_DIR / filename
    with open(path, "rb") as f:
        data = f.read()
    return io.BytesIO(data)


# =========================================================================
# Tests for add_warning
# =========================================================================


class TestAddWarning(TestCase):
    """Tests for the add_warning helper method."""

    def _make_vf(self):
        with patch.object(ValidateFile, "__init__", lambda self, *a, **kw: None):
            vf = ValidateFile.__new__(ValidateFile)
            vf.validate_dict = {"field": [], "warning_string": []}
            vf.validated = True
            return vf

    def test_add_warning_appends_field_and_string(self):
        vf = self._make_vf()
        vf.add_warning("test_field", "something went wrong")
        self.assertEqual(vf.validate_dict["field"], ["test_field"])
        self.assertEqual(vf.validate_dict["warning_string"], ["something went wrong"])

    def test_add_warning_accumulates(self):
        vf = self._make_vf()
        vf.add_warning("f1", "w1")
        vf.add_warning("f2", "w2")
        self.assertEqual(len(vf.validate_dict["field"]), 2)
        self.assertEqual(len(vf.validate_dict["warning_string"]), 2)


# =========================================================================
# Tests for checkColumnNames
# =========================================================================


class TestCheckColumnNames(TestCase):
    """Tests for the checkColumnNames method."""

    def _make_vf(self, columns):
        with patch.object(ValidateFile, "__init__", lambda self, *a, **kw: None):
            vf = ValidateFile.__new__(ValidateFile)
            vf.df_columns = pd.Index(columns)
            vf.validate_dict = {"field": [], "warning_string": []}
            vf.validated = True
            return vf

    def test_matching_columns_passes(self):
        cols = ["a", "b", "c"]
        vf = self._make_vf(cols)
        vf.checkColumnNames(expected_column_names=cols)
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.validate_dict["field"]), 0)

    def test_matching_columns_different_order_passes(self):
        vf = self._make_vf(["c", "a", "b"])
        vf.checkColumnNames(expected_column_names=["a", "b", "c"])
        self.assertTrue(vf.validated)

    def test_mismatched_columns_fails(self):
        vf = self._make_vf(["a", "b", "x"])
        vf.checkColumnNames(expected_column_names=["a", "b", "c"])
        self.assertFalse(vf.validated)
        self.assertIn("name_columns", vf.validate_dict["field"])

    def test_extra_column_fails(self):
        vf = self._make_vf(["a", "b", "c", "d"])
        vf.checkColumnNames(expected_column_names=["a", "b", "c"])
        self.assertFalse(vf.validated)

    def test_missing_column_fails(self):
        vf = self._make_vf(["a", "b"])
        vf.checkColumnNames(expected_column_names=["a", "b", "c"])
        self.assertFalse(vf.validated)

    def test_retro_api_expected_columns(self):
        expected = [
            "target-SMILES",
            "target-names",
            "concentration-required-mM",
            "amount-required-uL",
            "batch-tag",
        ]
        vf = self._make_vf(expected)
        vf.checkColumnNames(expected_column_names=expected)
        self.assertTrue(vf.validated)


# =========================================================================
# Tests for checkNumberColumns
# =========================================================================


class TestCheckNumberColumns(TestCase):
    """Tests for the checkNumberColumns method."""

    def _make_vf(self, n_cols):
        with patch.object(ValidateFile, "__init__", lambda self, *a, **kw: None):
            vf = ValidateFile.__new__(ValidateFile)
            vf.no_df_columns = n_cols
            vf.validate_dict = {"field": [], "warning_string": []}
            vf.validated = True
            return vf

    def test_correct_number_passes(self):
        vf = self._make_vf(5)
        vf.checkNumberColumns(expected_no_columns=5, expected_column_names=["a"] * 5)
        self.assertTrue(vf.validated)

    def test_too_few_columns_fails(self):
        vf = self._make_vf(3)
        vf.checkNumberColumns(expected_no_columns=5, expected_column_names=["a"] * 5)
        self.assertFalse(vf.validated)
        self.assertIn("number_columns", vf.validate_dict["field"])
        self.assertIn("Found 3 columns", vf.validate_dict["warning_string"][0])

    def test_too_many_columns_fails(self):
        vf = self._make_vf(7)
        vf.checkNumberColumns(expected_no_columns=5, expected_column_names=["a"] * 5)
        self.assertFalse(vf.validated)

    def test_retro_api_column_count(self):
        vf = self._make_vf(5)
        expected_column_names = [
            "target-SMILES",
            "target-names",
            "concentration-required-mM",
            "amount-required-uL",
            "batch-tag",
        ]
        vf.checkNumberColumns(
            expected_no_columns=5, expected_column_names=expected_column_names
        )
        self.assertTrue(vf.validated)


# =========================================================================
# Tests for checkSMILES
# =========================================================================


class TestCheckSMILES(TestCase):
    """Tests for the checkSMILES method."""

    def _make_vf(self):
        with patch.object(ValidateFile, "__init__", lambda self, *a, **kw: None):
            vf = ValidateFile.__new__(ValidateFile)
            vf.validate_dict = {"field": [], "warning_string": []}
            vf.validated = True
            return vf

    def test_valid_postera_smiles(self):
        vf = self._make_vf()
        vf.checkSMILES(
            df_rows_index=list(range(len(POSTERA_SMILES))),
            smiles=POSTERA_SMILES,
            smiles_type="target",
        )
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.validate_dict["field"]), 0)

    def test_valid_single_smiles(self):
        vf = self._make_vf()
        vf.checkSMILES(
            df_rows_index=[0],
            smiles=["CC(=O)Oc1ccccc1C(=O)O"],  # Aspirin
            smiles_type="target",
        )
        self.assertTrue(vf.validated)

    def test_invalid_smiles_fails(self):
        vf = self._make_vf()
        vf.checkSMILES(
            df_rows_index=[0, 1],
            smiles=[POSTERA_SMILES[0], "INVALID_SMILES_XYZ"],
            smiles_type="target",
        )
        self.assertFalse(vf.validated)
        self.assertIn("check_smiles", vf.validate_dict["field"])

    def test_empty_smiles_is_valid_in_rdkit(self):
        """RDKit treats empty string as a valid (no-atom) molecule."""
        vf = self._make_vf()
        vf.checkSMILES(
            df_rows_index=[0],
            smiles=[""],
            smiles_type="target",
        )
        self.assertTrue(vf.validated)

    def test_multiple_invalid_smiles_all_reported(self):
        vf = self._make_vf()
        vf.checkSMILES(
            df_rows_index=[0, 1, 2],
            smiles=["BAD1", POSTERA_SMILES[0], "BAD2"],
            smiles_type="target",
        )
        self.assertFalse(vf.validated)
        smiles_warnings = [f for f in vf.validate_dict["field"] if f == "check_smiles"]
        self.assertEqual(len(smiles_warnings), 2)

    def test_smiles_with_stereochemistry(self):
        vf = self._make_vf()
        stereo_smiles = [
            "C[C@H](O)c1ccccc1",
            "C[C@@H](O)c1ccccc1",
        ]
        vf.checkSMILES(df_rows_index=[0, 1], smiles=stereo_smiles, smiles_type="target")
        self.assertTrue(vf.validated)

    def test_reactant_smiles_from_custom_chem(self):
        """Reactant SMILES from custom-chem-ver1.5.csv are valid."""
        vf = self._make_vf()
        vf.checkSMILES(
            df_rows_index=[0, 1],
            smiles=[CUSTOM_CHEM_REACTANT_1_STEP1, CUSTOM_CHEM_REACTANT_2_STEP1],
            smiles_type="reactant",
        )
        self.assertTrue(vf.validated)


# =========================================================================
# Tests for checkIsNumber
# =========================================================================


class TestCheckIsNumber(TestCase):
    """Tests for the checkIsNumber method."""

    def _make_vf(self, n_rows=2):
        with patch.object(ValidateFile, "__init__", lambda self, *a, **kw: None):
            vf = ValidateFile.__new__(ValidateFile)
            vf.validate_dict = {"field": [], "warning_string": []}
            vf.validated = True
            vf.index_df_rows = range(n_rows)
            return vf

    def test_valid_concentrations(self):
        vf = self._make_vf(3)
        vf.checkIsNumber(values=POSTERA_CONCENTRATIONS)
        self.assertTrue(vf.validated)

    def test_valid_amounts(self):
        vf = self._make_vf(3)
        vf.checkIsNumber(values=POSTERA_AMOUNTS)
        self.assertTrue(vf.validated)

    def test_valid_integers(self):
        vf = self._make_vf(2)
        vf.checkIsNumber(values=[10, 20])
        self.assertTrue(vf.validated)

    def test_valid_floats(self):
        vf = self._make_vf(2)
        vf.checkIsNumber(values=[10.5, 20.3])
        self.assertTrue(vf.validated)

    def test_mixed_int_float(self):
        vf = self._make_vf(2)
        vf.checkIsNumber(values=[10, 20.5])
        self.assertTrue(vf.validated)

    def test_string_value_fails(self):
        vf = self._make_vf(2)
        vf.checkIsNumber(values=[10.0, "not_a_number"])
        self.assertFalse(vf.validated)
        self.assertIn("check_number", vf.validate_dict["field"])

    def test_none_value_fails(self):
        vf = self._make_vf(1)
        vf.checkIsNumber(values=[None])
        self.assertFalse(vf.validated)

    def test_empty_list_passes(self):
        vf = self._make_vf(0)
        vf.checkIsNumber(values=[])
        self.assertTrue(vf.validated)

    def test_negative_numbers_pass(self):
        vf = self._make_vf(2)
        vf.checkIsNumber(values=[-5.0, -10])
        self.assertTrue(vf.validated)

    def test_zero_passes(self):
        vf = self._make_vf(1)
        vf.checkIsNumber(values=[0])
        self.assertTrue(vf.validated)


# =========================================================================
# Tests for checkIsString
# =========================================================================


class TestCheckIsString(TestCase):
    """Tests for the checkIsString method."""

    def _make_vf(self, batch_tags):
        with patch.object(ValidateFile, "__init__", lambda self, *a, **kw: None):
            vf = ValidateFile.__new__(ValidateFile)
            vf.validate_dict = {"field": [], "warning_string": []}
            vf.validated = True
            n = len(batch_tags)
            vf.index_df_rows = range(n)
            vf.df = pd.DataFrame({"batch-tag": batch_tags})
            return vf

    def test_valid_strings_pass(self):
        vf = self._make_vf(POSTERA_BATCHTAGS)
        vf.checkIsString()
        self.assertTrue(vf.validated)

    def test_strips_whitespace(self):
        vf = self._make_vf(["  batch-1  ", " batch-2"])
        vf.checkIsString()
        self.assertTrue(vf.validated)
        self.assertEqual(vf.batchtags, ["batch-1", "batch-2"])


# =========================================================================
# Tests for checkReaction
# =========================================================================


class TestCheckReaction(TestCase):
    """Tests for the checkReaction method."""

    def _make_vf(self):
        with patch.object(ValidateFile, "__init__", lambda self, *a, **kw: None):
            vf = ValidateFile.__new__(ValidateFile)
            vf.validate_dict = {"field": [], "warning_string": []}
            vf.validated = True
            return vf

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_valid_reaction_returns_product(self, mock_check, mock_smarts):
        """Mock a successful reaction and verify product is returned."""
        vf = self._make_vf()
        mock_smarts.return_value = ["[C:1](=O)Cl.[N:2]>>[C:1](=O)[N:2]"]
        product_mol = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mock_check.return_value = [product_mol]

        result = vf.checkReaction(
            reactant_pair_smiles=[("CC(=O)Cl", "NCC1CCOCC1")],
            reaction_names=["Amidation"],
            reaction_recipes=["standard"],
        )
        self.assertTrue(vf.validated)
        self.assertEqual(len(result), 1)
        self.assertIsNotNone(Chem.MolFromSmiles(result[0]))

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_invalid_reaction_sets_validated_false(self, mock_check, mock_smarts):
        vf = self._make_vf()
        mock_smarts.return_value = ["[C:1](=O)Cl.[N:2]>>[C:1](=O)[N:2]"]
        mock_check.return_value = []  # no products → failure

        vf.checkReaction(
            reactant_pair_smiles=[("BAD", "SMILES")],
            reaction_names=["Amidation"],
            reaction_recipes=["standard"],
        )
        self.assertFalse(vf.validated)
        self.assertIn("check_reaction", vf.validate_dict["field"])

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_reaction_with_explicit_product_smiles(self, mock_check, mock_smarts):
        """When product_smiles are provided, the explicit product is used."""
        vf = self._make_vf()
        mock_smarts.return_value = ["SMARTS"]
        product_mol = Chem.MolFromSmiles(CUSTOM_CHEM_PRODUCT_STEP1)
        mock_check.return_value = [product_mol]

        result = vf.checkReaction(
            reactant_pair_smiles=[
                (CUSTOM_CHEM_REACTANT_1_STEP1, CUSTOM_CHEM_REACTANT_2_STEP1)
            ],
            reaction_names=[CUSTOM_CHEM_REACTION_NAME_STEP1],
            reaction_recipes=[CUSTOM_CHEM_RECIPE_STEP1],
            product_smiles=[CUSTOM_CHEM_PRODUCT_STEP1],
        )
        self.assertTrue(vf.validated)
        expected = Chem.MolToSmiles(Chem.MolFromSmiles(CUSTOM_CHEM_PRODUCT_STEP1))
        self.assertEqual(result[0], expected)

    @patch("backend.validate.get_recipe_smarts")
    def test_no_smarts_falls_back_to_product_smiles(self, mock_smarts):
        """When get_recipe_smarts returns falsy, provided product_smiles are used."""
        vf = self._make_vf()
        mock_smarts.return_value = (None, None)

        result = vf.checkReaction(
            reactant_pair_smiles=[
                (CUSTOM_CHEM_REACTANT_1_STEP1, CUSTOM_CHEM_REACTANT_2_STEP1)
            ],
            reaction_names=["Unknown-reaction"],
            reaction_recipes=["unknown-recipe"],
            product_smiles=[CUSTOM_CHEM_PRODUCT_STEP1],
        )
        self.assertTrue(vf.validated)
        self.assertEqual(result, [CUSTOM_CHEM_PRODUCT_STEP1])

    @patch("backend.validate.get_recipe_smarts")
    def test_reaction_exception_sets_validated_false(self, mock_smarts):
        vf = self._make_vf()
        mock_smarts.side_effect = Exception("SMARTS lookup failed")

        vf.checkReaction(
            reactant_pair_smiles=[("CC(=O)Cl", "NCC1CCOCC1")],
            reaction_names=["Amidation"],
            reaction_recipes=["standard"],
        )
        self.assertFalse(vf.validated)
        self.assertIn("check_reaction", vf.validate_dict["field"])

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_multiple_reactions_validated(self, mock_check, mock_smarts):
        """Multiple reaction pairs all succeed."""
        vf = self._make_vf()
        mock_smarts.return_value = ["SMARTS"]
        mol1 = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mol2 = Chem.MolFromSmiles("CC(=O)NCc1ccccc1")
        mock_check.side_effect = [[mol1], [mol2]]

        result = vf.checkReaction(
            reactant_pair_smiles=[
                ("CC(=O)Cl", "NCC1CCOCC1"),
                ("CC(=O)Cl", "NCc1ccccc1"),
            ],
            reaction_names=["Amidation", "Amidation"],
            reaction_recipes=["standard", "standard"],
        )
        self.assertTrue(vf.validated)
        self.assertEqual(len(result), 2)

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_second_reaction_invalid_first_valid(self, mock_check, mock_smarts):
        vf = self._make_vf()
        mock_smarts.return_value = ["SMARTS"]
        mol1 = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mock_check.side_effect = [[mol1], []]  # second fails

        vf.checkReaction(
            reactant_pair_smiles=[("CC(=O)Cl", "NCC1CCOCC1"), ("BAD", "SMILES")],
            reaction_names=["Amidation", "Amidation"],
            reaction_recipes=["standard", "standard"],
        )
        self.assertFalse(vf.validated)


# =========================================================================
# Tests for retro-API upload type (integration via constructor)
# =========================================================================


class TestRetroAPIValidation(TestCase):
    """Tests for the retro-API validation path through __init__."""

    def test_valid_retro_csv(self):
        csv = _make_retro_csv()
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.validate_dict["field"]), 0)

    def test_retro_csv_from_fixture_file(self):
        """Use the actual postera-ver1.5.csv fixture file."""
        csv = _open_fixture("postera-ver1.5.csv")
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.target_smiles), 3)

    def test_retro_csv_wrong_number_columns(self):
        csv = _make_csv(
            {
                "target-SMILES": [POSTERA_SMILES[0]],
                "target-names": [POSTERA_NAMES[0]],
                "concentration-required-mM": [100.0],
                # missing amount-required-uL and batch-tag
            }
        )
        vf = ValidateFile(csv, "retro-API")
        self.assertFalse(vf.validated)
        self.assertIn("number_columns", vf.validate_dict["field"])

    def test_retro_csv_wrong_column_names(self):
        csv = _make_csv(
            {
                "WRONG-SMILES": [POSTERA_SMILES[0]],
                "target-names": [POSTERA_NAMES[0]],
                "concentration-required-mM": [100.0],
                "amount-required-uL": [100.0],
                "batch-tag": ["to-check"],
            }
        )
        vf = ValidateFile(csv, "retro-API")
        self.assertFalse(vf.validated)
        self.assertIn("name_columns", vf.validate_dict["field"])

    def test_retro_csv_invalid_smiles(self):
        csv = _make_retro_csv(
            smiles=[POSTERA_SMILES[0], "NOT_VALID!!!", POSTERA_SMILES[2]]
        )
        vf = ValidateFile(csv, "retro-API")
        self.assertFalse(vf.validated)
        self.assertIn("check_smiles", vf.validate_dict["field"])

    def test_retro_csv_strips_whitespace_from_smiles(self):
        padded = [f"  {s}  " for s in POSTERA_SMILES]
        csv = _make_retro_csv(smiles=padded)
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)
        self.assertEqual(vf.target_smiles, POSTERA_SMILES)

    def test_retro_csv_single_row(self):
        csv = _make_retro_csv(
            smiles=[POSTERA_SMILES[0]],
            names=[POSTERA_NAMES[0]],
            concentrations=[100.0],
            amounts=[100.0],
            batch_tags=["to-check"],
        )
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)

    def test_retro_csv_extra_column_fails(self):
        csv = _make_csv(
            {
                "target-SMILES": [POSTERA_SMILES[0]],
                "target-names": [POSTERA_NAMES[0]],
                "concentration-required-mM": [100.0],
                "amount-required-uL": [100.0],
                "batch-tag": ["to-check"],
                "extra-col": ["oops"],
            }
        )
        vf = ValidateFile(csv, "retro-API")
        self.assertFalse(vf.validated)

    def test_retro_csv_all_invalid_smiles(self):
        csv = _make_retro_csv(smiles=["INVALID1", "INVALID2", "INVALID3"])
        vf = ValidateFile(csv, "retro-API")
        self.assertFalse(vf.validated)
        smiles_warnings = [
            f for f in vf.validate_dict["field"] if f == "check_smiles"
        ]
        self.assertEqual(len(smiles_warnings), 3)

    def test_retro_csv_empty_rows(self):
        """CSV with headers but no data rows should pass."""
        csv = _make_csv(
            {
                "target-SMILES": [],
                "target-names": [],
                "concentration-required-mM": [],
                "amount-required-uL": [],
                "batch-tag": [],
            }
        )
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)

    def test_retro_csv_unicode_names(self):
        csv = _make_retro_csv(names=["target-α", "target-β", "target-γ"])
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)

    def test_retro_csv_large_numbers(self):
        csv = _make_retro_csv(
            concentrations=[1e6, 1e-6, 1e3], amounts=[1e8, 1e-8, 1e4]
        )
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)

    def test_retro_csv_complex_drug_smiles(self):
        complex_smiles = [
            "C[C@H](O)c1ccccc1",
            "C[C@@H](O)c1ccccc1",
            "c1ccc2[nH]ccc2c1",  # indole
        ]
        csv = _make_retro_csv(smiles=complex_smiles)
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)

    def test_retro_csv_windows_line_endings(self):
        content = (
            "target-SMILES,target-names,concentration-required-mM,"
            "amount-required-uL,batch-tag\r\n"
            f"{POSTERA_SMILES[0]},{POSTERA_NAMES[0]},100,100,to-check\r\n"
        )
        buf = io.BytesIO(content.encode("utf-8"))
        vf = ValidateFile(buf, "retro-API")
        self.assertTrue(vf.validated)

    def test_multiple_instances_independent(self):
        csv1 = _make_retro_csv()
        csv2 = _make_retro_csv(smiles=["INVALID", "INVALID", "INVALID"])
        vf1 = ValidateFile(csv1, "retro-API")
        vf2 = ValidateFile(csv2, "retro-API")
        self.assertTrue(vf1.validated)
        self.assertFalse(vf2.validated)
        self.assertEqual(len(vf1.validate_dict["field"]), 0)
        self.assertGreater(len(vf2.validate_dict["field"]), 0)


# =========================================================================
# Tests for custom-chem upload type
# =========================================================================


class TestCustomChemValidation(TestCase):
    """Tests for the custom-chem validation path."""

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_valid_custom_chem_1step(self, mock_check, mock_smarts):
        """A valid single-step custom-chem CSV passes validation."""
        mock_smarts.return_value = ["SMARTS"]
        mock_check.return_value = [Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")]

        csv = _make_custom_chem_csv_1step(n_rows=1)
        vf = ValidateFile(csv, "custom-chem")
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.validate_dict["field"]), 0)

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_custom_chem_builds_dataframe(self, mock_check, mock_smarts):
        """Output DataFrame has the expected columns after validation."""
        mock_smarts.return_value = ["SMARTS"]
        mol = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mock_check.return_value = [mol]

        csv = _make_custom_chem_csv_1step(
            target_names=["t1"],
            batch_tags=["b1"],
            concentrations=[5.0],
            amounts=[50.0],
            n_rows=1,
        )
        vf = ValidateFile(csv, "custom-chem")
        self.assertTrue(vf.validated)
        for col in [
            "target-names", "batch-tag", "target-SMILES",
            "reactant-pair-smiles", "reaction-name", "reaction-recipe",
            "product-smiles",
        ]:
            self.assertIn(col, vf.df.columns)
        self.assertEqual(list(vf.df["target-names"]), ["t1"])
        self.assertEqual(list(vf.df["batch-tag"]), ["b1"])
        self.assertEqual(list(vf.df["concentration-required-mM"]), [5.0])
        self.assertEqual(list(vf.df["amount-required-uL"]), [50.0])

    def test_custom_chem_wrong_columns_raises(self):
        """Completely wrong columns cause a KeyError on 'no-steps' access."""
        csv = _make_csv({"wrong-col-1": ["a"], "wrong-col-2": ["b"]})
        with self.assertRaises(KeyError):
            ValidateFile(csv, "custom-chem")

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_custom_chem_reaction_failure_invalidates(self, mock_check, mock_smarts):
        mock_smarts.return_value = ["SMARTS"]
        mock_check.return_value = []  # no products

        csv = _make_custom_chem_csv_1step(n_rows=1)
        vf = ValidateFile(csv, "custom-chem")
        self.assertFalse(vf.validated)
        self.assertIn("check_reaction", vf.validate_dict["field"])

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_custom_chem_string_concentration_fails(self, mock_check, mock_smarts):
        mock_smarts.return_value = ["SMARTS"]
        mock_check.return_value = [Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")]

        csv = _make_custom_chem_csv_1step(concentrations=["not_num"], n_rows=1)
        vf = ValidateFile(csv, "custom-chem")
        self.assertFalse(vf.validated)
        self.assertIn("check_number", vf.validate_dict["field"])

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_custom_chem_stores_reaction_info(self, mock_check, mock_smarts):
        mock_smarts.return_value = ["SMARTS"]
        mol = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mock_check.return_value = [mol]

        csv = _make_custom_chem_csv_1step(n_rows=1)
        vf = ValidateFile(csv, "custom-chem")
        self.assertTrue(vf.validated)
        for col in ["reactant-pair-smiles", "reaction-name", "reaction-recipe", "product-smiles"]:
            self.assertIn(col, vf.df.columns)

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    def test_custom_chem_two_rows(self, mock_check, mock_smarts):
        """Two-row custom-chem CSV processes both rows."""
        mock_smarts.return_value = ["SMARTS"]
        mol1 = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mol2 = Chem.MolFromSmiles("CC(=O)NCc1ccccc1")
        mock_check.side_effect = [[mol1], [mol2]]

        csv = _make_custom_chem_csv_1step(
            target_names=["t1", "t2"],
            reactant_2=["NCC1CCOCC1", "NCc1ccccc1"],
            product_smiles=["CC(=O)NCC1CCOCC1", "CC(=O)NCc1ccccc1"],
            n_rows=2,
        )
        vf = ValidateFile(csv, "custom-chem")
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.df), 2)

    def test_custom_chem_fixture_file_column_mismatch(self):
        """The actual fixture CSV lacks reaction-groupby-column-* columns,
        so it will fail column validation (expected behaviour)."""
        csv = _open_fixture("custom-chem-ver1.5.csv")
        vf = ValidateFile(csv, "custom-chem")
        self.assertFalse(vf.validated)


# =========================================================================
# Tests for combi-custom-chem upload type
# =========================================================================


class TestCombiCustomChemValidation(TestCase):
    """Tests for the combi-custom-chem validation path."""

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    @patch("backend.validate.combi_chem")
    def test_valid_combi_csv(self, mock_combi, mock_check, mock_smarts):
        mock_smarts.return_value = ["SMARTS"]
        mock_combi.return_value = [
            ("CC(=O)Cl", "NCC1CCOCC1"),
            ("CC(=O)Cl", "NCc1ccccc1"),
        ]
        mol1 = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mol2 = Chem.MolFromSmiles("CC(=O)NCc1ccccc1")
        mock_check.side_effect = [[mol1], [mol2]]

        csv = _make_combi_csv_1step()
        vf = ValidateFile(csv, "combi-custom-chem")
        self.assertTrue(vf.validated)

    def test_combi_csv_wrong_columns_raises(self):
        """Completely wrong columns cause a KeyError on 'no-steps' access."""
        csv = _make_csv({"wrong": ["a"], "columns": ["b"]})
        with self.assertRaises(KeyError):
            ValidateFile(csv, "combi-custom-chem")

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    @patch("backend.validate.combi_chem")
    def test_combi_csv_builds_output_dataframe(self, mock_combi, mock_check, mock_smarts):
        mock_smarts.return_value = ["SMARTS"]
        mock_combi.return_value = [
            ("CC(=O)Cl", "NCC1CCOCC1"),
            ("CC(=O)Cl", "NCc1ccccc1"),
        ]
        mol1 = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mol2 = Chem.MolFromSmiles("CC(=O)NCc1ccccc1")
        mock_check.side_effect = [[mol1], [mol2]]

        csv = _make_combi_csv_1step()
        vf = ValidateFile(csv, "combi-custom-chem")
        self.assertTrue(vf.validated)
        for col in [
            "batch-tag", "target-names", "target-SMILES",
            "concentration-required-mM", "amount-required-uL",
            "no-steps", "reactant-pair-smiles", "reaction-name",
            "reaction-recipe", "product-smiles",
        ]:
            self.assertIn(col, vf.df.columns, f"Missing column: {col}")

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    @patch("backend.validate.combi_chem")
    def test_combi_csv_produces_correct_number_of_targets(
        self, mock_combi, mock_check, mock_smarts
    ):
        mock_smarts.return_value = ["SMARTS"]
        mock_combi.return_value = [
            ("CC(=O)Cl", "NCC1CCOCC1"),
            ("CC(=O)Cl", "NCc1ccccc1"),
        ]
        mol1 = Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")
        mol2 = Chem.MolFromSmiles("CC(=O)NCc1ccccc1")
        mock_check.side_effect = [[mol1], [mol2]]

        csv = _make_combi_csv_1step()
        vf = ValidateFile(csv, "combi-custom-chem")
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.target_smiles), 2)

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    @patch("backend.validate.combi_chem")
    def test_combi_csv_reaction_failure(self, mock_combi, mock_check, mock_smarts):
        mock_smarts.return_value = ["SMARTS"]
        mock_combi.return_value = [("CC(=O)Cl", "NCC1CCOCC1")]
        mock_check.return_value = []  # fail

        csv = _make_combi_csv_1step()
        vf = ValidateFile(csv, "combi-custom-chem")
        self.assertFalse(vf.validated)

    @patch("backend.validate.get_recipe_smarts")
    @patch("backend.validate.check_reactant_smarts")
    @patch("backend.validate.combi_chem")
    def test_combi_csv_non_numeric_concentration_raises(self, mock_combi, mock_check, mock_smarts):
        """Non-numeric concentration causes a ValueError during float() conversion."""
        mock_smarts.return_value = ["SMARTS"]
        mock_combi.return_value = [("CC(=O)Cl", "NCC1CCOCC1")]
        mock_check.return_value = [Chem.MolFromSmiles("CC(=O)NCC1CCOCC1")]

        csv = _make_combi_csv_1step(concentrations=["NaN_str", "NaN_str"])
        with self.assertRaises(ValueError):
            ValidateFile(csv, "combi-custom-chem")


# =========================================================================
# Tests for unknown upload type
# =========================================================================


class TestUnknownUploadType(TestCase):
    """An unrecognised upload_type passes through without validation."""

    def test_unknown_type_no_validation(self):
        csv = _make_retro_csv()
        vf = ValidateFile(csv, "unknown-type")
        self.assertTrue(vf.validated)
        self.assertEqual(vf.upload_type, "unknown-type")
        self.assertEqual(len(vf.validate_dict["field"]), 0)


# =========================================================================
# Tests for validate_dict structure
# =========================================================================


class TestValidateDict(TestCase):
    """validate_dict should always be well-formed."""

    def test_keys_always_present(self):
        csv = _make_retro_csv()
        vf = ValidateFile(csv, "retro-API")
        self.assertIn("field", vf.validate_dict)
        self.assertIn("warning_string", vf.validate_dict)
        self.assertEqual(
            len(vf.validate_dict["field"]), len(vf.validate_dict["warning_string"])
        )

    def test_validate_dict_on_failure(self):
        csv = _make_retro_csv(smiles=["INVALID", "INVALID", "INVALID"])
        vf = ValidateFile(csv, "retro-API")
        self.assertFalse(vf.validated)
        self.assertGreater(len(vf.validate_dict["field"]), 0)
        self.assertEqual(
            len(vf.validate_dict["field"]), len(vf.validate_dict["warning_string"])
        )


# =========================================================================
# Edge cases
# =========================================================================


class TestEdgeCases(TestCase):
    """Edge-case tests for ValidateFile."""

    def test_retro_csv_with_rings(self):
        ring_smiles = ["c1ccccc1", "C1CCCCC1", "c1ccc2[nH]ccc2c1"]
        csv = _make_retro_csv(smiles=ring_smiles)
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)

    def test_retro_csv_with_charged_atoms(self):
        charged = ["[NH4+]", "CC(=O)[O-]", "c1cc[nH+]cc1"]
        csv = _make_retro_csv(smiles=charged)
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)


# =========================================================================
# Integration tests using actual fixture files
# =========================================================================


class TestFixtureFileIntegration(TestCase):
    """Integration tests loading real fixture CSV files."""

    def test_postera_fixture_retro_api(self):
        """postera-ver1.5.csv should pass retro-API validation."""
        csv = _open_fixture("postera-ver1.5.csv")
        vf = ValidateFile(csv, "retro-API")
        self.assertTrue(vf.validated)
        self.assertEqual(len(vf.target_smiles), 3)
        # Verify SMILES match what we expect
        self.assertEqual(vf.target_smiles, POSTERA_SMILES)

    def test_custom_chem_fixture_lacks_groupby_columns(self):
        """custom-chem-ver1.5.csv does not include reaction-groupby-column-*
        so it should fail the column-count / column-name checks."""
        csv = _open_fixture("custom-chem-ver1.5.csv")
        vf = ValidateFile(csv, "custom-chem")
        self.assertFalse(vf.validated)

    def test_combi_fixture_column_structure(self):
        """custom-combi-batch-ver1.5.csv has the correct columns for
        combi-custom-chem validation (column check should pass)."""
        csv = _open_fixture("custom-combi-batch-ver1.5.csv")
        # Column and number checks should pass; any failure would be
        # in checkReaction which hits the DB — we just check columns here
        df = pd.read_csv(csv)
        expected_cols = {
            "combi-group", "no-steps", "concentration-required-mM",
            "amount-required-uL", "batch-tag", "reactant-1-1",
            "reactant-2-1", "reaction-name-1", "reaction-recipe-1",
            "reactant-1-2", "reaction-name-2", "reaction-recipe-2",
            "reactant-1-3", "reaction-name-3", "reaction-recipe-3",
        }
        self.assertEqual(set(df.columns), expected_cols)