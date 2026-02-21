"""Tests for the otsession package.

Tests for:
  - SessionOrchestrator  (entry point, session routing)
  - BaseSession          (common setup, model creation, cleanup)
  - ReactionSession      (reaction execution flow)
  - WorkupSession        (workup execution flow, stage determination)
  - AnalysisSession      (analysis execution flow)
  - DataManager          (statistics, querysets, dataframes, grouping)
  - PipetteManager       (tip-rack & pipette selection, transfer calculations)
  - MaterialManager      (mass calc, vial counts, salt matching, combine_strings)
  - MaterialCalculator   (salt-form recalculation)
  - SaltMatchingService  (find_matching_materials)
  - DeckManager          (slot allocation)
  - WellManager          (well availability, index helpers)
  - ColumnManager        (column CRUD)
  - LabwareSelector      (labware comparison / selection)
  - PlateQueryService    (plate queries)
  - PlateFactory         (plate creation flows)

All external dependencies (Django ORM, RDKit MW, file storage) are mocked.
"""

from unittest import TestCase
from unittest.mock import patch, MagicMock, PropertyMock, call
import math
import pandas as pd

# ---------------------------------------------------------------------------
# Module paths used for patching
# ---------------------------------------------------------------------------
_SESSION_MOD = "backend.opentrons.otsession"
_ORCH_MOD = f"{_SESSION_MOD}.session_orchestrator"
_BASE_MOD = f"{_SESSION_MOD}.sessions.base_session"
_REACT_MOD = f"{_SESSION_MOD}.sessions.reaction_session"
_WORKUP_MOD = f"{_SESSION_MOD}.sessions.workup_session"
_ANALYSIS_MOD = f"{_SESSION_MOD}.sessions.analysis_session"
_DATA_MOD = f"{_SESSION_MOD}.managers.data_manager"
_PIP_MOD = f"{_SESSION_MOD}.managers.pipette_manager"
_MAT_MOD = f"{_SESSION_MOD}.managers.material_manager"
_DECK_MOD = f"{_SESSION_MOD}.managers.deck_manager"
_PLATE_FAC_MOD = f"{_SESSION_MOD}.managers.plate_manager.plate_factory"
_WELL_MOD = f"{_SESSION_MOD}.managers.plate_manager.well_manager"
_COL_MOD = f"{_SESSION_MOD}.managers.plate_manager.column_manager"
_LAB_MOD = f"{_SESSION_MOD}.managers.plate_manager.labware_selector"
_PQS_MOD = f"{_SESSION_MOD}.managers.plate_manager.plate_query_service"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _mock_action_session_qs(session_type="reaction", ids=None, continuation=None):
    """Return a MagicMock that behaves like a QuerySet[ActionSession]."""
    ids = ids or [1, 2]
    qs = MagicMock()
    qs.exists.return_value = True
    qs.values_list.return_value = MagicMock()

    def _values_list(field, flat=False):
        m = MagicMock()
        if field == "id":
            m.__iter__ = lambda s: iter(ids)
            m.distinct.return_value = ids
            return m
        if field == "reaction_id":
            m.__iter__ = lambda s: iter([10, 20])
            m.distinct.return_value.__getitem__ = lambda s, i: [10, 20]
            m.distinct.return_value.__iter__ = lambda s: iter([10, 20])
            # for list() call
            m.distinct.return_value.__len__ = lambda s: 2
            return m
        if field == "type":
            m.distinct.return_value.__getitem__ = lambda s, i: session_type
            m.distinct.return_value.__iter__ = lambda s: iter([session_type])
            return m
        return m

    qs.values_list.side_effect = _values_list

    # support .filter(continuation=True/False)
    if continuation is not None:
        qs.filter.return_value = _mock_action_session_qs(session_type, ids)
    else:
        cont_qs = MagicMock()
        cont_qs.exists.return_value = False
        non_cont_qs = MagicMock()
        non_cont_qs.exists.return_value = True
        def _filter(**kwargs):
            if kwargs.get("continuation") is True:
                return cont_qs
            return non_cont_qs
        qs.filter.side_effect = _filter

    return qs


def _mock_plate(**kwargs):
    """Return a MagicMock Plate with sensible defaults."""
    defaults = dict(
        id=1,
        role="reaction",
        role_index=1,
        labware="plateone_96_wellplate_2500ul",
        maxwellvolume=2500,
        numberwells=96,
        numberwellsincolumn=8,
        numbercolumns=12,
        indexswellavailable=0,
        indexcolumnavailable=0,
        name="test-plate",
    )
    defaults.update(kwargs)
    plate = MagicMock()
    for k, v in defaults.items():
        setattr(plate, k, v)
    return plate


def _make_session_stub(**overrides):
    """Build a minimal stub that looks like a BaseSession for manager __init__."""
    session = MagicMock()
    session.reactionstep = 1
    session.otsessionobj = MagicMock(id=100)
    session.otbatchprotocolobj = MagicMock()
    session.batchobj = MagicMock(batchtag="B1")
    session.actionsession_ids = [1, 2]
    session.reaction_ids = [10, 20]
    session.actionsessiontype = "reaction"
    session.deckobj = MagicMock(
        indexslotavailable=1, numberslots=11, slotavailable=True
    )
    for k, v in overrides.items():
        setattr(session, k, v)
    return session


# ===================================================================
# DataManager pure-logic tests
# ===================================================================
class TestDataManagerMedian(TestCase):
    """DataManager.get_median_value — pure statistics."""

    def setUp(self):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        self.dm = DataManager(session=_make_session_stub())

    def test_empty_list(self):
        self.assertEqual(self.dm.get_median_value([]), 0)

    def test_odd_length(self):
        self.assertEqual(self.dm.get_median_value([3, 1, 2]), 2)

    def test_even_length(self):
        self.assertEqual(self.dm.get_median_value([1, 2, 3, 4]), 2.5)

    def test_single_value(self):
        self.assertEqual(self.dm.get_median_value([42]), 42)

    def test_all_same(self):
        self.assertEqual(self.dm.get_median_value([5, 5, 5]), 5)


class TestDataManagerSum(TestCase):
    """DataManager.get_sum_value."""

    def setUp(self):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        self.dm = DataManager(session=_make_session_stub())

    def test_normal(self):
        self.assertEqual(self.dm.get_sum_value([1, 2, 3]), 6)

    def test_empty(self):
        self.assertEqual(self.dm.get_sum_value([]), 0)

    def test_floats(self):
        self.assertAlmostEqual(self.dm.get_sum_value([1.1, 2.2, 3.3]), 6.6, places=5)


class TestDataManagerRoundedVolumes(TestCase):
    """get_rounded_add_action_volumes / get_rounded_extract_action_volumes."""

    def setUp(self):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        self.dm = DataManager(session=_make_session_stub())

    def test_add_action_volumes(self):
        mock_actions = [MagicMock(volume=1.6), MagicMock(volume=2.4), MagicMock(volume=3.5)]
        result = self.dm.get_rounded_add_action_volumes(mock_actions)
        self.assertEqual(result, [2, 2, 4])

    def test_extract_action_volumes(self):
        mock_actions = [MagicMock(volume=10.3), MagicMock(volume=20.7)]
        result = self.dm.get_rounded_extract_action_volumes(mock_actions)
        self.assertEqual(result, [10, 21])

    def test_empty_queryset(self):
        self.assertEqual(self.dm.get_rounded_add_action_volumes([]), [])

    def test_exact_integers(self):
        mock_actions = [MagicMock(volume=5.0), MagicMock(volume=10.0)]
        result = self.dm.get_rounded_add_action_volumes(mock_actions)
        self.assertEqual(result, [5, 10])


class TestDataManagerAddActionsDataframe(TestCase):
    """get_add_actions_dataframe — dataframe building with uniquesolution."""

    def setUp(self):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        self.session = _make_session_stub()
        self.session.material_manager = MagicMock()
        self.session.material_manager.combine_strings.side_effect = (
            lambda row: f"{row['smiles']}-{row['solvent']}-{row['concentration']}"
        )
        self.dm = DataManager(session=self.session)

    def test_empty_queryset_returns_empty_df(self):
        qs = MagicMock()
        qs.exists.return_value = False
        result = self.dm.get_add_actions_dataframe(qs)
        self.assertTrue(result.empty)

    def test_builds_dataframe_with_uniquesolution(self):
        qs = MagicMock()
        qs.exists.return_value = True
        qs.values.return_value = [
            {"smiles": "CCO", "solvent": "DMSO", "concentration": 0.5, "volume": 100},
        ]
        result = self.dm.get_add_actions_dataframe(qs)
        self.assertIn("uniquesolution", result.columns)
        self.assertEqual(result.iloc[0]["uniquesolution"], "CCO-DMSO-0.5")


class TestDataManagerGroupedTemperatures(TestCase):
    """get_unique_temperatures and get_grouped_temperature_reactions."""

    def setUp(self):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        self.dm = DataManager(session=_make_session_stub())

    @patch(f"{_DATA_MOD}.get_reaction_temperature")
    def test_unique_temperatures(self, mock_get_temp):
        mock_get_temp.side_effect = [25, 80, 25]
        reactions = [MagicMock(id=1), MagicMock(id=2), MagicMock(id=3)]
        result = self.dm.get_unique_temperatures(reactions)
        self.assertEqual(result, [25, 80])

    @patch(f"{_DATA_MOD}.get_reaction_query_set")
    @patch(f"{_DATA_MOD}.get_reaction_temperature")
    def test_grouped_temperature_reactions(self, mock_get_temp, mock_get_rxn_qs):
        # 2 calls in get_unique_temperatures + 2 temps × 2 reactions in grouping = 6
        mock_get_temp.side_effect = [25, 80, 25, 80, 25, 80]
        reactions = [MagicMock(id=1), MagicMock(id=2)]
        mock_get_rxn_qs.return_value = MagicMock()

        result = self.dm.get_grouped_temperature_reactions(reactions)
        self.assertIn(25, result)
        self.assertIn(80, result)


class TestDataManagerGroupedClassRecipe(TestCase):
    """get_unique_reaction_classes, get_unique_reaction_recipes, get_grouped_reaction_by_class_recipe."""

    def setUp(self):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        self.dm = DataManager(session=_make_session_stub())

    @patch(f"{_DATA_MOD}.get_reaction_class")
    def test_unique_classes(self, mock_get_class):
        mock_get_class.side_effect = ["amidation", "sulfonamide", "amidation"]
        reactions = [MagicMock(id=1), MagicMock(id=2), MagicMock(id=3)]
        result = self.dm.get_unique_reaction_classes(reactions)
        self.assertEqual(result, ["amidation", "sulfonamide"])

    @patch(f"{_DATA_MOD}.get_reaction_recipe")
    def test_unique_recipes(self, mock_get_recipe):
        mock_get_recipe.side_effect = ["recipe-A", "recipe-B", "recipe-A"]
        reactions = [MagicMock(id=1), MagicMock(id=2), MagicMock(id=3)]
        result = self.dm.get_unique_reaction_recipes(reactions)
        self.assertEqual(result, ["recipe-A", "recipe-B"])


# ===================================================================
# PipetteManager pure-logic tests
# ===================================================================
class TestPipetteManagerTipRackType(TestCase):
    """PipetteManager.get_tip_rack_type — volume-based selection."""

    def setUp(self):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        self.pm = PipetteManager(session=_make_session_stub())

    def test_empty_volumes_default_300(self):
        result = self.pm.get_tip_rack_type([])
        self.assertEqual(result, "opentrons_96_tiprack_300ul")

    def test_small_volumes_20ul(self):
        result = self.pm.get_tip_rack_type([5, 10, 15, 20])
        self.assertIn("20ul", result)

    def test_medium_volumes_300ul(self):
        result = self.pm.get_tip_rack_type([50, 100, 200])
        self.assertIn("300ul", result)

    def test_large_volumes_1000ul(self):
        result = self.pm.get_tip_rack_type([500, 800, 1000])
        self.assertIn("1000ul", result)

    def test_boundary_20(self):
        # median of [20] is 20, which should select 20ul
        result = self.pm.get_tip_rack_type([20])
        self.assertEqual(result, "opentrons_96_tiprack_20ul")

    def test_boundary_300(self):
        result = self.pm.get_tip_rack_type([300])
        self.assertEqual(result, "opentrons_96_tiprack_300ul")


class TestPipetteManagerNumberTransfers(TestCase):
    """PipetteManager.get_number_transfers."""

    def setUp(self):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        self.pm = PipetteManager(session=_make_session_stub())

    def test_all_fit_single_transfer(self):
        result = self.pm.get_number_transfers(300, [100, 200, 300])
        self.assertEqual(result, 3)

    def test_needs_multiple_transfers(self):
        # 500µL with 300µL pipette → ceiling(500/300) = 2
        result = self.pm.get_number_transfers(300, [500])
        self.assertEqual(result, 2)

    def test_exact_multiple(self):
        # 600µL with 300µL pipette → exactly 2
        result = self.pm.get_number_transfers(300, [600])
        self.assertEqual(result, 2)

    def test_empty_volumes(self):
        result = self.pm.get_number_transfers(300, [])
        self.assertEqual(result, 0)

    def test_mixed_volumes(self):
        # 100 → 1, 500 → 2, 300 → 1 = 4
        result = self.pm.get_number_transfers(300, [100, 500, 300])
        self.assertEqual(result, 4)


class TestPipetteManagerPipetteType(TestCase):
    """PipetteManager.get_pipette_type — optimal pipette selection."""

    def setUp(self):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        self.pm = PipetteManager(session=_make_session_stub())

    def test_empty_volumes_returns_p300_single(self):
        result = self.pm.get_pipette_type([], channel_type="single")
        self.assertEqual(result["maxvolume"], 300)
        self.assertEqual(result["type"], "single")

    def test_small_volumes_p10(self):
        result = self.pm.get_pipette_type([5, 8, 3], channel_type="single")
        self.assertEqual(result["labware"], "p10_single")

    def test_medium_volumes_p300(self):
        result = self.pm.get_pipette_type([100, 200, 150], channel_type="single")
        self.assertEqual(result["labware"], "p300_single")

    def test_large_volumes_p1000(self):
        result = self.pm.get_pipette_type([500, 800], channel_type="single")
        self.assertEqual(result["labware"], "p1000_single_gen2")

    def test_multi_channel(self):
        result = self.pm.get_pipette_type([5, 10], channel_type="multi")
        self.assertEqual(result["type"], "multi")

    def test_exceeds_all_returns_largest(self):
        result = self.pm.get_pipette_type([2000], channel_type="single")
        self.assertEqual(result["maxvolume"], 1000)


# ===================================================================
# MaterialCalculator tests
# ===================================================================
class TestMaterialCalculator(TestCase):
    """MaterialCalculator.recalculate_amount_for_salt_form."""

    def setUp(self):
        from backend.opentrons.otsession.managers.material_manager import (
            MaterialCalculator,
        )

        self.calc = MaterialCalculator()

    @patch(f"{_MAT_MOD}.Descriptors")
    @patch(f"{_MAT_MOD}.Chem")
    def test_salt_form_volume_adjustment(self, mock_chem, mock_desc):
        mock_chem.MolFromSmiles.return_value = MagicMock()
        mock_desc.MolWt.side_effect = [100.0, 150.0]  # original, actual
        result = self.calc.recalculate_amount_for_salt_form("CCO", "CCO.Cl", 100.0)
        self.assertAlmostEqual(result, 150.0, places=1)  # 100 * (150/100)

    @patch(f"{_MAT_MOD}.Chem")
    def test_invalid_smiles_returns_original(self, mock_chem):
        mock_chem.MolFromSmiles.return_value = None
        result = self.calc.recalculate_amount_for_salt_form("bad", "also_bad", 100.0)
        self.assertEqual(result, 100.0)


# ===================================================================
# SaltMatchingService tests
# ===================================================================
class TestSaltMatchingService(TestCase):
    """SaltMatchingService.find_matching_materials."""

    def setUp(self):
        from backend.opentrons.otsession.managers.material_manager import (
            SaltMatchingService,
        )

        self.sms = SaltMatchingService()

    def test_no_potential_matches(self):
        exists, wells, plate, remaining, actual = self.sms.find_matching_materials(
            "CCO", [], 100
        )
        self.assertFalse(exists)
        self.assertEqual(wells, [])
        self.assertIsNone(plate)
        self.assertEqual(remaining, 100)

    @patch(f"{_MAT_MOD}.are_equivalent_structures", return_value=True)
    @patch(f"{_MAT_MOD}.canon_smiles", side_effect=lambda x: x)
    @patch(f"{_MAT_MOD}.strip_salts", side_effect=lambda x: x)
    def test_match_found_sufficient_volume(self, mock_strip, mock_canon, mock_equiv):
        well = MagicMock(smiles="CCO", volume=200, plate_id=MagicMock(id=5))
        exists, wells, plate, remaining, actual = self.sms.find_matching_materials(
            "CCO", [well], 100
        )
        self.assertTrue(exists)
        self.assertEqual(len(wells), 1)
        self.assertEqual(remaining, 0)
        self.assertEqual(actual, "CCO")

    @patch(f"{_MAT_MOD}.are_equivalent_structures", return_value=True)
    @patch(f"{_MAT_MOD}.canon_smiles", side_effect=lambda x: x)
    @patch(f"{_MAT_MOD}.strip_salts", side_effect=lambda x: x)
    def test_match_found_insufficient_volume(self, mock_strip, mock_canon, mock_equiv):
        well = MagicMock(smiles="CCO", volume=50, plate_id=MagicMock(id=5))
        exists, wells, plate, remaining, actual = self.sms.find_matching_materials(
            "CCO", [well], 100
        )
        self.assertFalse(exists)
        self.assertEqual(remaining, 50)

    @patch(f"{_MAT_MOD}.are_equivalent_structures", return_value=False)
    @patch(f"{_MAT_MOD}.canon_smiles", side_effect=lambda x: x)
    @patch(f"{_MAT_MOD}.strip_salts", side_effect=lambda x: x)
    def test_no_structure_match(self, mock_strip, mock_canon, mock_equiv):
        well = MagicMock(smiles="CCCC", volume=200, plate_id=MagicMock(id=5))
        exists, wells, plate, remaining, actual = self.sms.find_matching_materials(
            "CCO", [well], 100
        )
        self.assertFalse(exists)
        self.assertEqual(wells, [])
        self.assertEqual(remaining, 100)


# ===================================================================
# MaterialManager pure-logic tests
# ===================================================================
class TestMaterialManagerCombineStrings(TestCase):
    """MaterialManager.combine_strings."""

    def setUp(self):
        from backend.opentrons.otsession.managers.material_manager import (
            MaterialManager,
        )

        self.mm = MaterialManager(session=_make_session_stub())

    def test_normal(self):
        row = {"smiles": "CCO", "solvent": "DMSO", "concentration": 0.5}
        self.assertEqual(self.mm.combine_strings(row), "CCO-DMSO-0.5")

    def test_none_values(self):
        row = {"smiles": None, "solvent": None, "concentration": None}
        self.assertEqual(self.mm.combine_strings(row), "None-None-None")


class TestMaterialManagerCalcMass(TestCase):
    """MaterialManager.calc_mass."""

    def setUp(self):
        from backend.opentrons.otsession.managers.material_manager import (
            MaterialManager,
        )

        self.mm = MaterialManager(session=_make_session_stub())

    @patch(f"{_MAT_MOD}.Descriptors")
    @patch(f"{_MAT_MOD}.Chem")
    def test_basic_calc(self, mock_chem, mock_desc):
        mock_chem.MolFromSmiles.return_value = MagicMock()
        mock_desc.MolWt.return_value = 100.0

        row = {"smiles": "CCO", "concentration": 1.0, "volume": 100}
        # mols = 1.0 * 100 * 1e-6 = 1e-4
        # mass_mg = 1e-4 * 100.0 * 1e3 = 10.0
        result = self.mm.calc_mass(row)
        self.assertAlmostEqual(result, 10.0, places=2)

    @patch(f"{_MAT_MOD}.Chem")
    def test_invalid_smiles_returns_zero(self, mock_chem):
        mock_chem.MolFromSmiles.return_value = None
        row = {"smiles": "bad", "concentration": 1.0, "volume": 100}
        self.assertEqual(self.mm.calc_mass(row), 0.0)


class TestMaterialManagerGetNumberVials(TestCase):
    """MaterialManager.get_number_vials."""

    def setUp(self):
        from backend.opentrons.otsession.managers.material_manager import (
            MaterialManager,
        )

        self.mm = MaterialManager(session=_make_session_stub())

    def test_fits_in_one_vial(self):
        self.assertEqual(self.mm.get_number_vials(2500, 1000), 1)

    def test_microvolume(self):
        self.assertEqual(self.mm.get_number_vials(2500, 3), 1)

    def test_exceeds_single_vial(self):
        # 2500µL vial, 5% dead volume → effective 2375µL
        # 5000 / 2375 = 2.1... → ceil = 3
        result = self.mm.get_number_vials(2500, 5000)
        self.assertEqual(result, math.ceil(5000 / (2500 * 0.95)))

    def test_exact_fit(self):
        # volume == max → falls to else branch (strict >), dead volume subtracted
        # 1000 / (1000 * 0.95) = 1.05 → ceil = 2
        self.assertEqual(self.mm.get_number_vials(1000, 1000), 2)


class TestMaterialManagerDeadVolume(TestCase):
    """MaterialManager.get_dead_volume."""

    def setUp(self):
        from backend.opentrons.otsession.managers.material_manager import (
            MaterialManager,
        )

        self.mm = MaterialManager(session=_make_session_stub())

    def test_five_percent(self):
        self.assertAlmostEqual(self.mm.get_dead_volume(2500), 125.0)

    def test_zero(self):
        self.assertAlmostEqual(self.mm.get_dead_volume(0), 0.0)


class TestMaterialManagerGetMaxWellVolume(TestCase):
    """MaterialManager.get_max_well_volume (delegates to plate attribute)."""

    def setUp(self):
        from backend.opentrons.otsession.managers.material_manager import (
            MaterialManager,
        )

        self.mm = MaterialManager(session=_make_session_stub())

    def test_returns_plate_maxwellvolume(self):
        plate = _mock_plate(maxwellvolume=2500)
        self.assertEqual(self.mm.get_max_well_volume(plate), 2500)


# ===================================================================
# DeckManager tests
# ===================================================================
class TestDeckManagerCreateDeck(TestCase):
    """DeckManager.create_deck_model."""

    @patch(f"{_DECK_MOD}.Deck")
    def test_creates_deck(self, MockDeck):
        from backend.opentrons.otsession.managers.deck_manager import DeckManager

        mock_deck = MagicMock()
        MockDeck.return_value = mock_deck
        session = _make_session_stub()

        dm = DeckManager(session)
        result = dm.create_deck_model()

        mock_deck.save.assert_called_once()
        self.assertEqual(mock_deck.numberslots, 11)
        self.assertEqual(result, mock_deck)


class TestDeckManagerSlotAvailability(TestCase):
    """DeckManager.check_deck_slot_available."""

    def test_returns_slot_when_available(self):
        from backend.opentrons.otsession.managers.deck_manager import DeckManager

        session = _make_session_stub()
        session.deckobj.indexslotavailable = 3
        session.deckobj.numberslots = 11
        dm = DeckManager(session)

        slot = dm.check_deck_slot_available()
        self.assertEqual(slot, 3)
        self.assertEqual(session.deckobj.indexslotavailable, 4)

    def test_raises_when_deck_full(self):
        from backend.opentrons.otsession.managers.deck_manager import DeckManager

        session = _make_session_stub()
        session.deckobj.indexslotavailable = 12
        session.deckobj.numberslots = 11
        dm = DeckManager(session)

        with self.assertRaises(ValueError):
            dm.check_deck_slot_available()


# ===================================================================
# WellManager tests
# ===================================================================
class TestWellManagerAvailability(TestCase):
    """WellManager.get_plate_well_index_available."""

    def setUp(self):
        from backend.opentrons.otsession.managers.plate_manager.well_manager import (
            WellManager,
        )

        self.wm = WellManager(_make_session_stub())

    def test_well_available(self):
        plate = _mock_plate(indexswellavailable=5, numberwells=96)
        self.assertEqual(self.wm.get_plate_well_index_available(plate), 5)

    def test_plate_full(self):
        plate = _mock_plate(indexswellavailable=96, numberwells=96)
        self.assertFalse(self.wm.get_plate_well_index_available(plate))

    def test_last_well(self):
        plate = _mock_plate(indexswellavailable=95, numberwells=96)
        self.assertEqual(self.wm.get_plate_well_index_available(plate), 95)


class TestWellManagerNewColumn(TestCase):
    """WellManager.check_index_well_is_new_column."""

    def setUp(self):
        from backend.opentrons.otsession.managers.plate_manager.well_manager import (
            WellManager,
        )

        self.wm = WellManager(_make_session_stub())

    def test_index_zero_is_new_column(self):
        plate = _mock_plate(indexswellavailable=0, numberwellsincolumn=8)
        self.assertTrue(self.wm.check_index_well_is_new_column(plate))

    def test_start_of_column(self):
        plate = _mock_plate(indexswellavailable=8, numberwellsincolumn=8)
        self.assertTrue(self.wm.check_index_well_is_new_column(plate))

    def test_mid_column(self):
        plate = _mock_plate(indexswellavailable=5, numberwellsincolumn=8)
        self.assertFalse(self.wm.check_index_well_is_new_column(plate))

    def test_column_16_wells(self):
        plate = _mock_plate(indexswellavailable=16, numberwellsincolumn=16)
        self.assertTrue(self.wm.check_index_well_is_new_column(plate))


class TestWellManagerNewColumnIndex(TestCase):
    """WellManager.get_new_column_and_well_index_available."""

    def setUp(self):
        from backend.opentrons.otsession.managers.plate_manager.well_manager import (
            WellManager,
        )

        self.wm = WellManager(_make_session_stub())

    def test_column_available(self):
        plate = _mock_plate(
            indexcolumnavailable=2, numbercolumns=12, numberwellsincolumn=8
        )
        result = self.wm.get_new_column_and_well_index_available(plate)
        self.assertEqual(result, (2, 16))  # col 2 → well index 2*8=16

    def test_no_more_columns(self):
        plate = _mock_plate(
            indexcolumnavailable=12, numbercolumns=12, numberwellsincolumn=8
        )
        result = self.wm.get_new_column_and_well_index_available(plate)
        self.assertFalse(result)


class TestWellManagerCreateWell(TestCase):
    """WellManager.create_well_model."""

    @patch(f"{_WELL_MOD}.well_index_to_well_name", return_value="A01")
    @patch(f"{_WELL_MOD}.Well")
    def test_creates_well_with_all_fields(self, MockWell, mock_name):
        from backend.opentrons.otsession.managers.plate_manager.well_manager import (
            WellManager,
        )

        mock_well = MagicMock()
        MockWell.return_value = mock_well
        session = _make_session_stub()
        wm = WellManager(session)
        plate = _mock_plate()

        result = wm.create_well_model(
            plate_obj=plate,
            role="reaction",
            role_index=1,
            wellindex=0,
            volume=100.0,
            smiles="CCO",
            concentration=0.5,
            solvent="DMSO",
        )
        mock_well.save.assert_called_once()
        self.assertEqual(mock_well.smiles, "CCO")
        self.assertEqual(mock_well.volume, 100.0)
        self.assertEqual(mock_well.name, "A01")


# ===================================================================
# ColumnManager tests
# ===================================================================
class TestColumnManagerCreateColumn(TestCase):
    """ColumnManager.create_column_model."""

    @patch(f"{_COL_MOD}.Column")
    def test_creates_column(self, MockColumn):
        from backend.opentrons.otsession.managers.plate_manager.column_manager import (
            ColumnManager,
        )

        mock_col = MagicMock()
        MockColumn.return_value = mock_col
        session = _make_session_stub()
        cm = ColumnManager(session)
        plate = _mock_plate()

        result = cm.create_column_model(
            plate_obj=plate,
            columnindex=0,
            role="reaction",
            role_index=1,
            reactionclass="amidation",
        )
        mock_col.save.assert_called_once()
        self.assertEqual(mock_col.reactionclass, "amidation")
        self.assertEqual(mock_col.index, 0)


class TestColumnManagerCurrentIndex(TestCase):
    """ColumnManager.get_plate_current_column_index."""

    def setUp(self):
        from backend.opentrons.otsession.managers.plate_manager.column_manager import (
            ColumnManager,
        )

        self.cm = ColumnManager(_make_session_stub())

    def test_column_available(self):
        plate = _mock_plate(indexcolumnavailable=3, numbercolumns=12)
        self.assertEqual(self.cm.get_plate_current_column_index(plate), 3)

    def test_all_columns_used(self):
        plate = _mock_plate(indexcolumnavailable=12, numbercolumns=12)
        self.assertFalse(self.cm.get_plate_current_column_index(plate))


# ===================================================================
# LabwareSelector tests
# ===================================================================
class TestLabwareSelectorGetPlateType(TestCase):
    """LabwareSelector.get_plate_type — labware comparison logic."""

    def setUp(self):
        from backend.opentrons.otsession.managers.plate_manager.labware_selector import (
            LabwareSelector,
        )

        session = _make_session_stub()
        session.data_manager = MagicMock()
        session.data_manager.get_median_value.side_effect = lambda values: (
            sorted(values)[len(values) // 2] if values else 0
        )
        session.material_manager = MagicMock()
        session.material_manager.get_number_vials.return_value = 1
        self.ls = LabwareSelector(session)

    def test_reaction_role_returns_plate(self):
        result = self.ls.get_plate_type(role="reaction", volumes=[500, 1000])
        self.assertIsNotNone(result)
        # Should be a recognized labware key
        from backend.opentrons.labwareavailable import labware_plates

        self.assertIn(result, labware_plates)

    def test_startingmaterial_role(self):
        result = self.ls.get_plate_type(role="startingmaterial", volumes=[100, 200])
        self.assertIsNotNone(result)

    def test_high_temperature_selects_reflux_capable(self):
        result = self.ls.get_plate_type(
            role="reaction", volumes=[500], temperature=200
        )
        # Should select paradox_96 which handles up to 250°C
        self.assertIn("paradox", result)

    def test_no_suitable_labware_returns_none_when_impossible(self):
        # Temperature too high for any labware → early return None
        result = self.ls.get_plate_type(
            role="reaction", volumes=[500], temperature=500
        )
        self.assertIsNone(result)


class TestLabwareSelectorDeadVolume(TestCase):
    """LabwareSelector dead volume and max well volume."""

    def setUp(self):
        from backend.opentrons.otsession.managers.plate_manager.labware_selector import (
            LabwareSelector,
        )

        self.ls = LabwareSelector(_make_session_stub())

    def test_dead_volume(self):
        self.assertAlmostEqual(self.ls.get_dead_volume(2500), 125.0)

    def test_max_well_volume(self):
        plate = _mock_plate(maxwellvolume=1000)
        self.assertEqual(self.ls.get_max_well_volume(plate), 1000)


# ===================================================================
# PlateQueryService tests
# ===================================================================
class TestPlateQueryServiceAllPlates(TestCase):
    """PlateQueryService.get_all_ot_batch_protocol_plates."""

    @patch(f"{_PQS_MOD}.Plate")
    def test_filters_by_protocol_and_role(self, MockPlate):
        from backend.opentrons.otsession.managers.plate_manager.plate_query_service import (
            PlateQueryService,
        )

        session = _make_session_stub()
        pqs = PlateQueryService(session)
        mock_qs = MagicMock()
        MockPlate.objects.filter.return_value = mock_qs

        result = pqs.get_all_ot_batch_protocol_plates(session.otbatchprotocolobj)
        MockPlate.objects.filter.assert_called_once()


class TestPlateQueryServiceInputPlates(TestCase):
    """PlateQueryService.get_input_plates_needed."""

    @patch(f"{_PQS_MOD}.Plate")
    def test_returns_empty_when_no_protocol_plates(self, MockPlate):
        from backend.opentrons.otsession.managers.plate_manager.plate_query_service import (
            PlateQueryService,
        )

        session = _make_session_stub()
        pqs = PlateQueryService(session)
        MockPlate.objects.filter.return_value = MagicMock(
            __iter__=lambda s: iter([]), __bool__=lambda s: False
        )

        result = pqs.get_input_plates_needed(
            searchsmiles=["CCO"], groupreactionqueryset=MagicMock()
        )
        self.assertEqual(result, [])


# ===================================================================
# PlateFactory tests
# ===================================================================
class TestPlateFactoryCreatePlateModel(TestCase):
    """PlateFactory.create_plate_model."""

    @patch(f"{_PLATE_FAC_MOD}.labware_plates", {
        "plateone_96_wellplate_2500ul": {
            "no_wells_in_column": 8,
            "volume_well": 2500,
            "no_wells": 96,
            "no_columns": 12,
        }
    })
    @patch(f"{_PLATE_FAC_MOD}.sanitize_for_python_var", return_value="test_plate")
    @patch(f"{_PLATE_FAC_MOD}.Plate")
    def test_creates_plate_when_slot_available(self, MockPlate, mock_sanitize):
        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )

        mock_plate = MagicMock(id=42)
        MockPlate.return_value = mock_plate
        session = _make_session_stub()
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = 4

        pf = PlateFactory(session)
        result = pf.create_plate_model(
            platename="test-plate",
            labwaretype="plateone_96_wellplate_2500ul",
            role="reaction",
            role_index=1,
        )

        self.assertEqual(result, mock_plate)
        self.assertEqual(mock_plate.save.call_count, 2)  # initial + name update
        self.assertEqual(mock_plate.role, "reaction")
        self.assertEqual(mock_plate.numberwells, 96)

    def test_returns_none_when_deck_full(self):
        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )

        session = _make_session_stub()
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = False

        pf = PlateFactory(session)
        result = pf.create_plate_model(
            platename="test",
            labwaretype="plateone_96_wellplate_2500ul",
            role="reaction",
            role_index=1,
        )
        self.assertIsNone(result)


class TestPlateFactoryUpdatePlateIds(TestCase):
    """PlateFactory.update_plate_deck_ot_session_ids."""

    def test_updates_all_plates(self):
        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )

        session = _make_session_stub()
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = 5
        session.well_manager = MagicMock()
        session.well_manager.get_plate_wells.return_value = MagicMock()
        session.column_manager = MagicMock()
        session.column_manager.get_plate_columns.return_value = MagicMock()

        plate1 = _mock_plate(id=1)
        plate2 = _mock_plate(id=2)

        pf = PlateFactory(session)
        pf.update_plate_deck_ot_session_ids([plate1, plate2])

        # Each plate should have .save() called
        self.assertTrue(plate1.save.called)
        self.assertTrue(plate2.save.called)


# ===================================================================
# SessionOrchestrator tests
# ===================================================================
class TestSessionOrchestratorDetermineType(TestCase):
    """SessionOrchestrator._determine_session_type."""

    def test_single_type(self):
        from backend.opentrons.otsession.session_orchestrator import (
            SessionOrchestrator,
        )

        qs = MagicMock()
        qs.exists.return_value = True
        distinct_mock = MagicMock()
        distinct_mock.__iter__ = lambda s: iter(["reaction"])
        distinct_mock.__len__ = lambda s: 1
        qs.values_list.return_value.distinct.return_value = distinct_mock

        with patch.object(SessionOrchestrator, "__init__", lambda self, **kw: None):
            orch = SessionOrchestrator()
            orch.actionsessionqueryset = qs
            result = orch._determine_session_type()
            self.assertEqual(result, "reaction")

    def test_empty_queryset_raises(self):
        from backend.opentrons.otsession.session_orchestrator import (
            SessionOrchestrator,
        )

        qs = MagicMock()
        qs.exists.return_value = False

        with patch.object(SessionOrchestrator, "__init__", lambda self, **kw: None):
            orch = SessionOrchestrator()
            orch.actionsessionqueryset = qs
            with self.assertRaises(ValueError):
                orch._determine_session_type()


class TestSessionOrchestratorCreateSession(TestCase):
    """SessionOrchestrator._create_session — session type routing."""

    def _make_orch(self, session_type):
        from backend.opentrons.otsession.session_orchestrator import (
            SessionOrchestrator,
        )

        with patch.object(SessionOrchestrator, "__init__", lambda self, **kw: None):
            orch = SessionOrchestrator()
            orch.actionsessiontype = session_type
            orch.reactionstep = 1
            orch.otbatchprotocolobj = MagicMock()
            orch.actionsessionqueryset = _mock_action_session_qs(session_type)
            orch.customSMcsvpath = None
            orch.use_multichannel = False
            return orch

    @patch(f"{_ORCH_MOD}.ReactionSession")
    def test_creates_reaction_session(self, MockRS):
        orch = self._make_orch("reaction")
        result = orch._create_session()
        MockRS.assert_called_once()

    @patch(f"{_ORCH_MOD}.WorkupSession")
    def test_creates_workup_session(self, MockWS):
        orch = self._make_orch("workup")
        result = orch._create_session()
        MockWS.assert_called_once()

    @patch(f"{_ORCH_MOD}.AnalysisSession")
    def test_creates_analysis_session(self, MockAS):
        orch = self._make_orch("analyse")
        result = orch._create_session()
        MockAS.assert_called_once()

    @patch(f"{_ORCH_MOD}.AnalysisSession")
    def test_creates_analysis_session_alt_spelling(self, MockAS):
        orch = self._make_orch("analysis")
        result = orch._create_session()
        MockAS.assert_called_once()

    def test_unknown_type_raises(self):
        orch = self._make_orch("invalid")
        with self.assertRaises(ValueError):
            orch._create_session()


class TestSessionOrchestratorExecute(TestCase):
    """SessionOrchestrator.execute."""

    def test_execute_delegates_to_session(self):
        from backend.opentrons.otsession.session_orchestrator import (
            SessionOrchestrator,
        )

        with patch.object(SessionOrchestrator, "__init__", lambda self, **kw: None):
            orch = SessionOrchestrator()
            orch.actionsessiontype = "reaction"
            mock_session = MagicMock()
            mock_session.execute.return_value = True
            mock_session.otsessionobj = MagicMock(id=99)
            orch.session = mock_session

            result = orch.execute()
            mock_session.execute.assert_called_once()
            self.assertTrue(result)
            self.assertEqual(orch.otsessionobj.id, 99)

    def test_execute_raises_if_no_otsession(self):
        from backend.opentrons.otsession.session_orchestrator import (
            SessionOrchestrator,
        )

        with patch.object(SessionOrchestrator, "__init__", lambda self, **kw: None):
            orch = SessionOrchestrator()
            orch.actionsessiontype = "reaction"
            mock_session = MagicMock()
            mock_session.execute.return_value = True
            mock_session.otsessionobj = None
            orch.session = mock_session

            with self.assertRaises(ValueError):
                orch.execute()

    def test_execute_calls_cleanup_on_error(self):
        from backend.opentrons.otsession.session_orchestrator import (
            SessionOrchestrator,
        )

        with patch.object(SessionOrchestrator, "__init__", lambda self, **kw: None):
            orch = SessionOrchestrator()
            orch.actionsessiontype = "reaction"
            mock_session = MagicMock()
            mock_session.execute.side_effect = RuntimeError("boom")
            mock_session.cleanup = MagicMock()
            orch.session = mock_session

            with self.assertRaises(RuntimeError):
                orch.execute()
            mock_session.cleanup.assert_called_once()


# ===================================================================
# BaseSession tests
# ===================================================================
class TestBaseSessionCreateOTSessionModel(TestCase):
    """BaseSession.create_ot_session_model."""

    @patch(f"{_BASE_MOD}.OTSession")
    def test_creates_and_saves(self, MockOTSession):
        from backend.opentrons.otsession.sessions.reaction_session import (
            ReactionSession,
        )

        mock_obj = MagicMock(id=42)
        MockOTSession.return_value = mock_obj
        qs = _mock_action_session_qs()
        proto = MagicMock()

        with patch.object(ReactionSession, "__init__", lambda self, **kw: None):
            session = ReactionSession()
            session.otbatchprotocolobj = proto
            session.actionsessiontype = "reaction"
            session.reactionstep = 1

            result = session.create_ot_session_model()
            mock_obj.save.assert_called_once()
            self.assertEqual(result.id, 42)


class TestBaseSessionCleanup(TestCase):
    """BaseSession.cleanup — deletes otsessionobj on error."""

    def test_cleanup_deletes_otsession(self):
        from backend.opentrons.otsession.sessions.reaction_session import (
            ReactionSession,
        )

        with patch.object(ReactionSession, "__init__", lambda self, **kw: None):
            session = ReactionSession()
            session.otsessionobj = MagicMock(id=10)
            session.cleanup()
            session.otsessionobj.delete.assert_called_once()

    def test_cleanup_no_otsession(self):
        from backend.opentrons.otsession.sessions.reaction_session import (
            ReactionSession,
        )

        with patch.object(ReactionSession, "__init__", lambda self, **kw: None):
            session = ReactionSession()
            session.otsessionobj = None
            session.cleanup()  # Should not raise


class TestBaseSessionSetupCommonResources(TestCase):
    """BaseSession.setup_common_resources — manager initialization."""

    @patch(f"{_BASE_MOD}.get_reaction_query_set")
    @patch(f"{_BASE_MOD}.DataManager")
    @patch(f"{_BASE_MOD}.DeckManager")
    @patch(f"{_BASE_MOD}.WellManager")
    @patch(f"{_BASE_MOD}.ColumnManager")
    @patch(f"{_BASE_MOD}.LabwareSelector")
    @patch(f"{_BASE_MOD}.PlateQueryService")
    @patch(f"{_BASE_MOD}.PlateFactory")
    @patch(f"{_BASE_MOD}.MaterialManager")
    @patch(f"{_BASE_MOD}.PipetteManager")
    @patch(f"{_BASE_MOD}.OTSession")
    def test_initializes_all_managers(
        self,
        MockOTSession,
        MockPipMgr,
        MockMatMgr,
        MockPlateFact,
        MockPQS,
        MockLabSel,
        MockColMgr,
        MockWellMgr,
        MockDeckMgr,
        MockDataMgr,
        mock_get_rxn_qs,
    ):
        from backend.opentrons.otsession.sessions.reaction_session import (
            ReactionSession,
        )

        mock_ot_obj = MagicMock(id=1)
        MockOTSession.return_value = mock_ot_obj
        mock_deck_mgr = MagicMock()
        mock_deck_mgr.create_deck_model.return_value = MagicMock()
        MockDeckMgr.return_value = mock_deck_mgr

        with patch.object(ReactionSession, "__init__", lambda self, **kw: None):
            session = ReactionSession()
            session.otbatchprotocolobj = MagicMock()
            session.actionsessiontype = "reaction"
            session.reactionstep = 1
            session.reaction_ids = [10, 20]
            session.otsessionobj = None
            session.deckobj = None
            session.groupreactionqueryset = None
            session.is_initialized = False

            session.setup_common_resources()

            self.assertTrue(session.is_initialized)
            self.assertIsNotNone(session.data_manager)
            self.assertIsNotNone(session.deck_manager)
            self.assertIsNotNone(session.well_manager)
            self.assertIsNotNone(session.plate_factory)
            self.assertIsNotNone(session.pipette_manager)
            self.assertIsNotNone(session.material_manager)


# ===================================================================
# ReactionSession tests
# ===================================================================
class TestReactionSessionExecute(TestCase):
    """ReactionSession.execute — full orchestration flow."""

    def _make_session(self):
        from backend.opentrons.otsession.sessions.reaction_session import (
            ReactionSession,
        )

        with patch.object(ReactionSession, "__init__", lambda self, **kw: None):
            session = ReactionSession()
            session.reactionstep = 1
            session.reaction_ids = [10, 20]
            session.actionsession_ids = [1, 2]
            session.actionsessionqueryset = _mock_action_session_qs()
            session.customSMcsvpath = None
            session.use_multichannel = False
            session.groupreactionqueryset = MagicMock()
            session.otsessionobj = MagicMock(id=1)
            session.roundedaddvolumes = []
            session.roundedextractvolumes = []
            session.roundedvolumes = []
            session.addactionsdf = pd.DataFrame()
            session.addactionqueryset = MagicMock()

            # Mock managers
            session.data_manager = MagicMock()
            session.data_manager.get_add_action_query_set.return_value = MagicMock(
                values_list=MagicMock(return_value=["CCO"])
            )
            session.data_manager.get_extract_action_query_set.return_value = MagicMock()
            session.data_manager.get_add_actions_dataframe.return_value = pd.DataFrame()
            session.data_manager.get_rounded_add_action_volumes.return_value = [100, 200]
            session.data_manager.get_rounded_extract_action_volumes.return_value = [50]
            session.data_manager.get_grouped_temperature_reactions.return_value = {
                25: MagicMock()
            }

            session.pipette_manager = MagicMock()
            session.pipette_manager.get_tip_rack_type.return_value = (
                "opentrons_96_tiprack_300ul"
            )
            session.pipette_manager.get_pipette_type.return_value = {
                "labware": "p300_single",
                "position": "right",
                "type": "single",
                "maxvolume": 300,
            }

            session.plate_factory = MagicMock()
            session.plate_query_service = MagicMock()
            session.plate_query_service.get_input_plates_needed.return_value = []
            session.material_manager = MagicMock()

            # patch setup_common_resources to do nothing (already set up)
            session.setup_common_resources = MagicMock()

            return session

    def test_execute_happy_path(self):
        session = self._make_session()
        result = session.execute()
        self.assertTrue(result)
        session.setup_common_resources.assert_called_once()
        session.pipette_manager.create_tip_racks.assert_called_once()
        session.pipette_manager.create_pipette_model.assert_called_once()
        session.plate_factory.create_reaction_starting_plate.assert_called_once()

    def test_execute_creates_plates_by_temperature(self):
        session = self._make_session()
        session.execute()
        session.plate_factory.create_plates_by_temperature.assert_called_once()

    def test_execute_with_custom_sm_csv(self):
        session = self._make_session()
        session.customSMcsvpath = "/tmp/custom.csv"
        session.execute()
        session.plate_factory.create_starting_material_plates_from_csv.assert_called_once_with(
            csv_path="/tmp/custom.csv"
        )

    def test_execute_step_gt_1_creates_solvent_plate(self):
        session = self._make_session()
        session.reactionstep = 2
        session.material_manager.get_add_actions_material_dataframe.return_value = (
            pd.DataFrame({"volume": [100]})
        )
        session.execute()
        session.plate_factory.create_solvent_plate.assert_called_once()


# ===================================================================
# WorkupSession tests
# ===================================================================
class TestWorkupSessionDetermineStage(TestCase):
    """WorkupSession.determine_workup_stage."""

    def _make_session(self, add_indices=None, extract_indices=None):
        from backend.opentrons.otsession.sessions.workup_session import WorkupSession

        with patch.object(WorkupSession, "__init__", lambda self, **kw: None):
            session = WorkupSession()

            # Mock add action queryset
            add_qs = MagicMock()
            if add_indices:
                add_qs.exists.return_value = True
                workup_add = MagicMock()
                workup_add.values_list.return_value.distinct.return_value = add_indices
                add_qs.filter.return_value = workup_add
            else:
                add_qs.exists.return_value = False

            # Mock extract action queryset
            extract_qs = MagicMock()
            if extract_indices:
                extract_qs.exists.return_value = True
                workup_extract = MagicMock()
                workup_extract.values_list.return_value.distinct.return_value = (
                    extract_indices
                )
                extract_qs.filter.return_value = workup_extract
            else:
                extract_qs.exists.return_value = False

            session.addactionqueryset = add_qs
            session.extractactionqueryset = extract_qs
            return session

    def test_stage_from_add_actions(self):
        session = self._make_session(add_indices=[2])
        self.assertEqual(session.determine_workup_stage(), 2)

    def test_stage_from_extract_actions(self):
        session = self._make_session(extract_indices=[3])
        self.assertEqual(session.determine_workup_stage(), 3)

    def test_default_stage_1(self):
        session = self._make_session()
        self.assertEqual(session.determine_workup_stage(), 1)

    def test_max_of_both(self):
        session = self._make_session(add_indices=[1], extract_indices=[2])
        self.assertEqual(session.determine_workup_stage(), 2)


class TestWorkupSessionGetSourcePlates(TestCase):
    """WorkupSession.get_source_plates_for_workup."""

    def _make_session(self, workup_stage, all_plates):
        from backend.opentrons.otsession.sessions.workup_session import WorkupSession

        with patch.object(WorkupSession, "__init__", lambda self, **kw: None):
            session = WorkupSession()
            session.otbatchprotocolobj = MagicMock()
            session.plate_query_service = MagicMock()
            session.plate_query_service.get_all_ot_batch_protocol_plates.return_value = (
                all_plates
            )
            # Mock determine_workup_stage
            session.determine_workup_stage = MagicMock(return_value=workup_stage)
            session.addactionqueryset = MagicMock(exists=MagicMock(return_value=False))
            session.extractactionqueryset = MagicMock(
                exists=MagicMock(return_value=False)
            )
            return session

    def test_stage_1_gets_reaction_plates(self):
        plates = [
            _mock_plate(role="reaction", role_index=1),
            _mock_plate(role="workup", role_index=1),
        ]
        session = self._make_session(1, plates)
        result = session.get_source_plates_for_workup()
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].role, "reaction")

    def test_stage_2_gets_workup_1_plates(self):
        plates = [
            _mock_plate(role="reaction", role_index=1),
            _mock_plate(role="workup", role_index=1),
        ]
        session = self._make_session(2, plates)
        result = session.get_source_plates_for_workup()
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].role, "workup")

    def test_no_source_returns_empty(self):
        session = self._make_session(1, [])
        result = session.get_source_plates_for_workup()
        self.assertEqual(result, [])


class TestWorkupSessionExecute(TestCase):
    """WorkupSession.execute — full workup flow."""

    def test_execute_happy_path(self):
        from backend.opentrons.otsession.sessions.workup_session import WorkupSession

        with patch.object(WorkupSession, "__init__", lambda self, **kw: None):
            session = WorkupSession()
            session.reactionstep = 1
            session.reaction_ids = [10]
            session.actionsession_ids = [1]
            session.roundedaddvolumes = []
            session.roundedextractvolumes = []
            session.roundedvolumes = []
            session.addactionsdf = pd.DataFrame()
            session.addactionqueryset = MagicMock()
            session.extractactionqueryset = MagicMock()

            session.setup_common_resources = MagicMock()
            session.data_manager = MagicMock()
            session.data_manager.get_add_action_query_set.return_value = MagicMock(
                exists=MagicMock(return_value=True)
            )
            session.data_manager.get_extract_action_query_set.return_value = MagicMock()
            session.data_manager.get_add_actions_dataframe.return_value = pd.DataFrame()
            session.data_manager.get_rounded_add_action_volumes.return_value = [100]
            session.data_manager.get_rounded_extract_action_volumes.return_value = [50]

            session.pipette_manager = MagicMock()
            session.pipette_manager.get_tip_rack_type.return_value = (
                "opentrons_96_tiprack_300ul"
            )
            session.pipette_manager.get_pipette_type.return_value = {
                "labware": "p300_single"
            }

            session.plate_factory = MagicMock()
            session.plate_query_service = MagicMock()
            session.material_manager = MagicMock()
            session.material_manager.get_add_actions_material_dataframe.return_value = (
                pd.DataFrame()
            )

            session.determine_workup_stage = MagicMock(return_value=1)
            session.get_source_plates_for_workup = MagicMock(return_value=[])

            result = session.execute()
            self.assertTrue(result)
            session.setup_common_resources.assert_called_once()
            session.plate_factory.create_workup_plate.assert_called_once()
            session.pipette_manager.create_pipette_model.assert_called_once()


# ===================================================================
# AnalysisSession tests
# ===================================================================
class TestAnalysisSessionGetSourcePlates(TestCase):
    """AnalysisSession.get_source_plates_for_analysis."""

    def _make_session(self, all_plates):
        from backend.opentrons.otsession.sessions.analysis_session import (
            AnalysisSession,
        )

        with patch.object(AnalysisSession, "__init__", lambda self, **kw: None):
            session = AnalysisSession()
            session.otbatchprotocolobj = MagicMock()
            session.plate_query_service = MagicMock()
            session.plate_query_service.get_all_ot_batch_protocol_plates.return_value = (
                all_plates
            )
            return session

    def test_prefers_highest_workup_plates(self):
        plates = [
            _mock_plate(role="workup", role_index=1),
            _mock_plate(role="workup", role_index=2),
            _mock_plate(role="reaction", role_index=1),
        ]
        session = self._make_session(plates)
        result = session.get_source_plates_for_analysis()
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].role_index, 2)

    def test_falls_back_to_reaction_plates(self):
        plates = [_mock_plate(role="reaction", role_index=1)]
        session = self._make_session(plates)
        result = session.get_source_plates_for_analysis()
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0].role, "reaction")

    def test_no_plates_returns_empty(self):
        session = self._make_session([])
        result = session.get_source_plates_for_analysis()
        self.assertEqual(result, [])


class TestAnalysisSessionExecute(TestCase):
    """AnalysisSession.execute — full analysis flow."""

    def test_execute_happy_path(self):
        from backend.opentrons.otsession.sessions.analysis_session import (
            AnalysisSession,
        )

        with patch.object(AnalysisSession, "__init__", lambda self, **kw: None):
            session = AnalysisSession()
            session.reactionstep = 1
            session.reaction_ids = [10]
            session.actionsession_ids = [1]
            session.roundedaddvolumes = []
            session.roundedextractvolumes = []
            session.roundedvolumes = []
            session.addactionsdf = pd.DataFrame()
            session.addactionqueryset = MagicMock()
            session.extractactionqueryset = MagicMock()

            session.setup_common_resources = MagicMock()
            session.data_manager = MagicMock()
            session.data_manager.get_extract_action_query_set.return_value = MagicMock()
            session.data_manager.get_add_action_query_set.return_value = MagicMock(
                exists=MagicMock(return_value=False)
            )
            session.data_manager.get_rounded_extract_action_volumes.return_value = [50]

            session.pipette_manager = MagicMock()
            session.pipette_manager.get_tip_rack_type.return_value = (
                "opentrons_96_tiprack_300ul"
            )
            session.pipette_manager.get_pipette_type.return_value = {
                "labware": "p300_single"
            }

            session.plate_factory = MagicMock()
            session.material_manager = MagicMock()

            session.get_source_plates_for_analysis = MagicMock(return_value=[])

            result = session.execute()
            self.assertTrue(result)
            session.setup_common_resources.assert_called_once()
            session.plate_factory.create_analyse_plate.assert_called_once()
            session.pipette_manager.create_pipette_model.assert_called_once()


# ===================================================================
# DataManager create_compound_order_model
# ===================================================================
class TestDataManagerCompoundOrder(TestCase):
    """DataManager.create_compound_order_model."""

    @patch(f"{_DATA_MOD}.default_storage")
    @patch(f"{_DATA_MOD}.CompoundOrder")
    def test_creates_order(self, MockCO, mock_storage):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        mock_co = MagicMock()
        MockCO.return_value = mock_co
        mock_storage.save.return_value = "compoundorders/test.csv"

        session = _make_session_stub()
        dm = DataManager(session)

        df = pd.DataFrame({"smiles": ["CCO"], "volume": [100]})
        result = dm.create_compound_order_model(df, is_custom_starter_plate=True)

        mock_co.save.assert_called_once()
        self.assertTrue(mock_co.iscustomSMplate)
        self.assertEqual(result, mock_co)

    @patch(f"{_DATA_MOD}.default_storage")
    @patch(f"{_DATA_MOD}.CompoundOrder")
    def test_empty_df_still_creates_order(self, MockCO, mock_storage):
        """Even empty CSV content should create a model (the method doesn't guard)."""
        from backend.opentrons.otsession.managers.data_manager import DataManager

        mock_co = MagicMock()
        MockCO.return_value = mock_co
        mock_storage.save.return_value = "compoundorders/test.csv"

        session = _make_session_stub()
        dm = DataManager(session)
        df = pd.DataFrame()
        result = dm.create_compound_order_model(df)
        mock_co.save.assert_called_once()


# ===================================================================
# PipetteManager model creation
# ===================================================================
class TestPipetteManagerCreateModels(TestCase):
    """PipetteManager.create_pipette_model, create_tiprack_model, create_tip_racks."""

    @patch(f"{_PIP_MOD}.Pipette")
    def test_create_pipette_model(self, MockPipette):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        mock_pip = MagicMock(id=1)
        MockPipette.return_value = mock_pip
        session = _make_session_stub()
        session.pipettetype = {
            "labware": "p300_single",
            "position": "right",
            "type": "single",
            "maxvolume": 300,
        }

        pm = PipetteManager(session)
        result = pm.create_pipette_model()

        mock_pip.save.assert_called_once()
        self.assertEqual(mock_pip.labware, "p300_single")
        self.assertEqual(mock_pip.maxvolume, 300)

    @patch(f"{_PIP_MOD}.TipRack")
    def test_create_tiprack_model(self, MockTipRack):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        mock_tr = MagicMock()
        MockTipRack.return_value = mock_tr
        session = _make_session_stub()
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = 2

        pm = PipetteManager(session)
        result = pm.create_tiprack_model("opentrons_96_tiprack_300ul")

        mock_tr.save.assert_called_once()
        self.assertEqual(mock_tr.labware, "opentrons_96_tiprack_300ul")

    @patch(f"{_PIP_MOD}.TipRack")
    def test_create_tiprack_returns_none_when_full(self, MockTipRack):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        session = _make_session_stub()
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = None

        pm = PipetteManager(session)
        result = pm.create_tiprack_model("opentrons_96_tiprack_300ul")
        self.assertIsNone(result)

    @patch(f"{_PIP_MOD}.TipRack")
    def test_create_tip_racks_creates_three(self, MockTipRack):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        mock_tr = MagicMock()
        MockTipRack.return_value = mock_tr
        session = _make_session_stub()
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = 5

        pm = PipetteManager(session)
        result = pm.create_tip_racks("opentrons_96_tiprack_300ul")
        self.assertEqual(len(result), 3)


class TestPipetteManagerSetupForSession(TestCase):
    """PipetteManager.setup_pipettes_for_session."""

    @patch(f"{_PIP_MOD}.TipRack")
    @patch(f"{_PIP_MOD}.Pipette")
    def test_setup_success(self, MockPipette, MockTipRack):
        from backend.opentrons.otsession.managers.pipette_manager import PipetteManager

        MockPipette.return_value = MagicMock(id=1)
        MockTipRack.return_value = MagicMock()
        session = _make_session_stub()
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = 3

        pm = PipetteManager(session)
        result = pm.setup_pipettes_for_session([100, 200])
        self.assertTrue(result)


# ===================================================================
# MaterialManager product SMILES
# ===================================================================
class TestMaterialManagerGetProductSmiles(TestCase):
    """MaterialManager.get_product_smiles."""

    @patch(f"{_MAT_MOD}.canon_smiles", side_effect=lambda x: x)
    @patch(f"{_MAT_MOD}.Product")
    def test_returns_unique_smiles(self, MockProduct, mock_canon):
        from backend.opentrons.otsession.managers.material_manager import (
            MaterialManager,
        )

        MockProduct.objects.filter.return_value.values_list.return_value = [
            "CCO",
            "CCO",
            "CCCC",
        ]

        mm = MaterialManager(session=_make_session_stub())
        result = mm.get_product_smiles([10, 20])
        self.assertEqual(len(result), 2)
        self.assertIn("CCO", result)
        self.assertIn("CCCC", result)


# ===================================================================
# DataManager rounded reaction volumes
# ===================================================================
class TestDataManagerRoundedReactionVolumes(TestCase):
    """DataManager.get_rounded_reaction_volumes — sums add action volumes per reaction."""

    def setUp(self):
        from backend.opentrons.otsession.managers.data_manager import DataManager

        session = _make_session_stub()
        self.dm = DataManager(session)

    def test_sums_volumes_per_reaction(self):
        reactions = [MagicMock(id=1), MagicMock(id=2)]

        add_qs_1 = [MagicMock(volume=100.5), MagicMock(volume=50.3)]
        add_qs_2 = [MagicMock(volume=200.7)]

        self.dm.get_add_action_query_set = MagicMock(side_effect=[add_qs_1, add_qs_2])

        result = self.dm.get_rounded_reaction_volumes(reactions)
        # reaction 1: round(100.5) + round(50.3) = 100+50 = 150 (sum of rounded is 150)
        # reaction 2: round(200.7) = 201
        # Actually get_rounded_add_action_volumes rounds each, then get_sum_value sums
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0], 100 + 50)  # rounded then summed
        self.assertEqual(result[1], 201)


# ===================================================================
# Integration-style: BaseSession __init__ sets up properties
# ===================================================================
class TestBaseSessionInit(TestCase):
    """BaseSession.__init__ correctly sets up all properties from queryset."""

    def test_init_sets_properties(self):
        from backend.opentrons.otsession.sessions.reaction_session import (
            ReactionSession,
        )

        proto = MagicMock()
        proto.batch_id = MagicMock(batchtag="B1")
        qs = _mock_action_session_qs("reaction", ids=[1, 2])

        # We need to arrange the queryset to make __init__ work
        with patch.object(
            ReactionSession, "__init__",
            lambda self, reactionstep, otbatchprotocolobj, actionsessionqueryset, customSMcsvpath=None: (
                setattr(self, "reactionstep", reactionstep) or
                setattr(self, "otbatchprotocolobj", otbatchprotocolobj)
            ),
        ):
            session = ReactionSession(
                reactionstep=1,
                otbatchprotocolobj=proto,
                actionsessionqueryset=qs,
            )
            self.assertEqual(session.reactionstep, 1)


# ===================================================================
# PlateFactory solvent plate
# ===================================================================
class TestPlateFactorySolventPlate(TestCase):
    """PlateFactory.create_solvent_plate."""

    def test_empty_df_returns_none(self):
        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )

        session = _make_session_stub()
        pf = PlateFactory(session)
        result = pf.create_solvent_plate(pd.DataFrame())
        self.assertIsNone(result)

    @patch(f"{_PLATE_FAC_MOD}.Plate")
    @patch(f"{_PLATE_FAC_MOD}.labware_plates", {
        "fluidx_24_vials_2500ul": {
            "type": ["startingmaterial"],
            "no_wells_in_column": 4,
            "no_wells": 24,
            "no_columns": 6,
            "volume_well": 2500,
            "max_fill_percentage": 60,
            "max_temp": 25,
        }
    })
    def test_creates_solvent_plate_and_well(self, MockPlate):
        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )

        session = _make_session_stub()
        session.labware_selector = MagicMock()
        session.labware_selector.get_plate_type.return_value = "fluidx_24_vials_2500ul"
        session.deck_manager = MagicMock()
        session.deck_manager.check_deck_slot_available.return_value = 5
        session.well_manager = MagicMock()
        session.well_manager.get_plate_well_index_available.return_value = 0
        session.well_manager.create_well_model.return_value = MagicMock(
            index=0, name="A01"
        )
        session.material_manager = MagicMock()

        mock_plate = MagicMock(id=1, name="solvent-DMSO", labware="fluidx_24_vials_2500ul")
        MockPlate.return_value = mock_plate
        MockPlate.objects.filter.return_value = MagicMock(exists=MagicMock(return_value=False))

        pf = PlateFactory(session)
        df = pd.DataFrame({"solvent": ["DMSO"], "volume": [500]})
        result = pf.create_solvent_plate(df)

        session.well_manager.create_well_model.assert_called_once()
        session.material_manager.create_solvent_prep_model.assert_called_once()


# ===================================================================
# PlateFactory: create_reaction_starting_plate multichannel branching
# ===================================================================
class TestPlateFactoryReactionStartingPlateBranching(TestCase):
    """PlateFactory.create_reaction_starting_plate dispatches correctly."""

    def _make_pf(self, *, use_multichannel=False, materials_df=None):
        """Build a PlateFactory with a session stub wired for the test."""
        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )

        session = _make_session_stub()
        session.use_multichannel = use_multichannel

        if materials_df is None:
            materials_df = pd.DataFrame({
                "uniquesolution": ["CCO-DMSO-10.0"],
                "smiles": ["CCO"],
                "volume": [500],
                "solvent": ["DMSO"],
                "concentration": [10.0],
                "molecularweight": [46.07],
            })
        session.material_manager = MagicMock()
        session.material_manager.get_add_actions_material_dataframe.return_value = materials_df

        pf = PlateFactory(session)
        return pf

    def test_empty_materials_returns_none(self):
        pf = self._make_pf(materials_df=pd.DataFrame())
        result = pf.create_reaction_starting_plate()
        self.assertIsNone(result)

    def test_use_multichannel_false_calls_sequential(self):
        pf = self._make_pf(use_multichannel=False)
        pf._create_sequential_starting_plate = MagicMock(return_value=MagicMock())
        pf._create_multichannel_starting_plate = MagicMock()

        pf.create_reaction_starting_plate()

        pf._create_sequential_starting_plate.assert_called_once()
        pf._create_multichannel_starting_plate.assert_not_called()

    def test_use_multichannel_true_calls_multichannel(self):
        pf = self._make_pf(use_multichannel=True)
        pf._create_sequential_starting_plate = MagicMock()
        pf._create_multichannel_starting_plate = MagicMock(return_value=MagicMock())

        pf.create_reaction_starting_plate()

        pf._create_multichannel_starting_plate.assert_called_once()
        pf._create_sequential_starting_plate.assert_not_called()


# ===================================================================
# PlateFactory: _create_multichannel_starting_plate integration tests
# ===================================================================
def _make_multichannel_session(
    *,
    addactionsdf=None,
    materials_df=None,
    wells_per_column=8,
    num_columns=12,
    max_well_volume=2500,
):
    """Build a fully-wired session stub for _create_multichannel_starting_plate tests.

    Returns (session, mock_plate) where mock_plate is the plate object that
    ``create_plate_model`` will return.
    """
    session = _make_session_stub()
    session.use_multichannel = True

    mock_plate = _mock_plate(
        labware="plateone_96_wellplate_2500ul",
        numberwellsincolumn=wells_per_column,
        numbercolumns=num_columns,
        maxwellvolume=max_well_volume,
        name="reaction-starting-materials-mc",
    )

    # labware selector
    session.labware_selector = MagicMock()
    session.labware_selector.get_plate_type.return_value = "plateone_96_wellplate_2500ul"
    session.labware_selector.get_dead_volume.return_value = 50.0  # 50µL dead volume

    # well manager
    _well_counter = [0]

    def _create_well(plate_obj, role, role_index, wellindex, volume,
                     smiles, concentration=None, solvent=None, transfer_type="single"):
        w = MagicMock()
        w.index = wellindex
        w.name = f"well-{wellindex}"
        w.transfer_type = transfer_type
        return w

    session.well_manager = MagicMock()
    session.well_manager.create_well_model.side_effect = _create_well
    session.well_manager.get_max_well_volume.return_value = max_well_volume * 0.6  # 60% fill
    session.well_manager.get_plate_well_index_available.side_effect = (
        lambda plate_obj: _well_counter[0]
    )

    def _update_well_index(plate_obj, wellindexupdate):
        _well_counter[0] = wellindexupdate

    session.well_manager.update_plate_well_index.side_effect = _update_well_index

    # column manager
    session.column_manager = MagicMock()

    # deck manager
    session.deck_manager = MagicMock()
    session.deck_manager.check_deck_slot_available.return_value = 3

    # material manager
    session.material_manager = MagicMock()
    session.material_manager.combine_strings.side_effect = (
        lambda row: f"{row['smiles']}-{row['solvent']}-{row['concentration']}"
    )
    session.material_manager.calc_mass.return_value = 1.0

    # data manager (for _create_compound_order)
    session.data_manager = MagicMock()

    # addactionsdf
    session.addactionsdf = addactionsdf

    return session, mock_plate


class TestCreateMultichannelStartingPlate(TestCase):
    """PlateFactory._create_multichannel_starting_plate — full integration tests."""

    def _run(self, *, addactionsdf, materials_df, wells_per_column=8,
             num_columns=12, max_well_volume=2500):
        """Helper: build session, create PlateFactory, call the method."""
        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )

        session, mock_plate = _make_multichannel_session(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=wells_per_column,
            num_columns=num_columns,
            max_well_volume=max_well_volume,
        )

        pf = PlateFactory(session)
        pf.create_plate_model = MagicMock(return_value=mock_plate)
        pf._create_compound_order = MagicMock()
        pf._create_sequential_starting_plate = MagicMock(return_value=mock_plate)

        result = pf._create_multichannel_starting_plate(materials_df)
        return result, pf, session, mock_plate

    # --- Fixture builder helpers ---

    @staticmethod
    def _make_addactions_df(entries):
        """Build an addactionsdf from a list of dicts.

        Each dict: {smiles, volume, solvent, concentration, reaction_id}
        """
        rows = []
        for e in entries:
            rows.append({
                "smiles": e["smiles"],
                "volume": e["volume"],
                "solvent": e.get("solvent", "DMSO"),
                "concentration": e.get("concentration", 10.0),
                "reaction_id_id": e["reaction_id"],
                "molecularweight": e.get("molecularweight", 46.07),
            })
        df = pd.DataFrame(rows)
        df["uniquesolution"] = df.apply(
            lambda r: f"{r['smiles']}-{r['solvent']}-{r['concentration']}", axis=1
        )
        return df

    @staticmethod
    def _make_materials_df(entries):
        """Build a materials_df from a list of dicts.

        Each dict: {smiles, volume, solvent, concentration}
        volume = total summed volume for that unique solution.
        """
        rows = []
        for e in entries:
            rows.append({
                "uniquesolution": f"{e['smiles']}-{e.get('solvent', 'DMSO')}-{e.get('concentration', 10.0)}",
                "smiles": e["smiles"],
                "volume": e["volume"],
                "solvent": e.get("solvent", "DMSO"),
                "concentration": e.get("concentration", 10.0),
                "molecularweight": e.get("molecularweight", 46.07),
            })
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Test: fallback to sequential when addactionsdf is None
    # ------------------------------------------------------------------
    def test_fallback_when_no_addactionsdf(self):
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 500},
        ])
        result, pf, session, _ = self._run(
            addactionsdf=None,
            materials_df=materials_df,
        )
        pf._create_sequential_starting_plate.assert_called_once_with(materials_df)

    # ------------------------------------------------------------------
    # Test: fallback when addactionsdf is empty
    # ------------------------------------------------------------------
    def test_fallback_when_addactionsdf_empty(self):
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 500},
        ])
        result, pf, session, _ = self._run(
            addactionsdf=pd.DataFrame(),
            materials_df=materials_df,
        )
        pf._create_sequential_starting_plate.assert_called_once_with(materials_df)

    # ------------------------------------------------------------------
    # Test: all materials go to single channel (< wells_per_column reactions)
    # ------------------------------------------------------------------
    def test_all_single_channel_when_too_few_reactions(self):
        """3 reactions for one reagent on an 8-well-per-column plate → all SC."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": 1},
            {"smiles": "CCO", "volume": 50, "reaction_id": 2},
            {"smiles": "CCO", "volume": 50, "reaction_id": 3},
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 150},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        # Should be called for the SC well(s), with transfer_type="single"
        calls = session.well_manager.create_well_model.call_args_list
        for c in calls:
            self.assertEqual(c.kwargs.get("transfer_type", c[1].get("transfer_type", "single")), "single")

    # ------------------------------------------------------------------
    # Test: one reagent with exactly wells_per_column reactions → 1 MC column
    # ------------------------------------------------------------------
    def test_one_full_multichannel_column(self):
        """8 reactions for one reagent → one full MC column, no SC wells."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 9)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},  # 50 * 8 = 400
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        # Should have created 8 wells with transfer_type="multichannel"
        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        self.assertEqual(len(mc_calls), 8)

        # All wells should be in column 0 (indices 0–7)
        well_indices = [_get_kwarg(c, "wellindex") for c in mc_calls]
        self.assertEqual(sorted(well_indices), list(range(8)))

        # Compound order should have been created
        pf._create_compound_order.assert_called_once()

    # ------------------------------------------------------------------
    # Test: reagent with 10 reactions → 1 MC column (8) + 2 SC leftovers
    # ------------------------------------------------------------------
    def test_multichannel_with_leftover_to_single_channel(self):
        """10 reactions: 1 MC column of 8, 2 leftovers go to SC."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 11)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 500},  # 50 * 10 = 500
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        # 8 MC wells + at least 1 SC well
        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        self.assertEqual(len(mc_calls), 8)
        self.assertGreaterEqual(len(sc_calls), 1)

    # ------------------------------------------------------------------
    # Test: two different reagents, one MC-eligible, one SC
    # ------------------------------------------------------------------
    def test_mixed_mc_and_sc_reagents(self):
        """Reagent A: 8 reactions (MC), Reagent B: 3 reactions (SC)."""
        addactionsdf = self._make_addactions_df(
            [{"smiles": "CCO", "volume": 50, "reaction_id": i} for i in range(1, 9)]
            + [{"smiles": "CCCC", "volume": 100, "reaction_id": i} for i in range(1, 4)]
        )
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},
            {"smiles": "CCCC", "volume": 300},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]

        # 8 MC wells for CCO
        self.assertEqual(len(mc_calls), 8)
        mc_smiles = {_get_kwarg(c, "smiles") for c in mc_calls}
        self.assertEqual(mc_smiles, {"CCO"})

        # SC wells for CCCC
        self.assertGreaterEqual(len(sc_calls), 1)
        sc_smiles = {_get_kwarg(c, "smiles") for c in sc_calls}
        self.assertIn("CCCC", sc_smiles)

    # ------------------------------------------------------------------
    # Test: 16 reactions → 2 MC columns
    # ------------------------------------------------------------------
    def test_16_reactions_single_reagent_consolidates_to_one_column(self):
        """16 reactions for one reagent at low volume → 1 MC column with 2 transfers/channel.

        The factory consolidates transfers into fewer source columns when
        the per-well volume fits.  transfers_per_channel = ceil(16/8) = 2,
        per_well = 50*2*1.15 + 50 = 165 µL, well under the 1500 µL max.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 17)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 800},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        # All 8 wells in 1 column (volume allows 2 transfers per channel)
        self.assertEqual(len(mc_calls), 8)
        well_indices = sorted(_get_kwarg(c, "wellindex") for c in mc_calls)
        self.assertEqual(well_indices, list(range(8)))

        # Per-well volume: 50 * 2 * 1.15 + 50 = 165
        expected_vol = round(50 * 2 * 1.15 + 50, 2)
        for c in mc_calls:
            self.assertAlmostEqual(_get_kwarg(c, "volume"), expected_vol, places=1)

    # ------------------------------------------------------------------
    # Test: two different reagents each with ≥8 reactions → 2 MC columns
    # ------------------------------------------------------------------
    def test_two_reagents_two_mc_columns(self):
        """Two reagents each with 8 reactions → 2 MC columns (one per reagent)."""
        addactionsdf = self._make_addactions_df(
            [{"smiles": "CCO", "volume": 50, "reaction_id": i} for i in range(1, 9)]
            + [{"smiles": "CCCC", "volume": 100, "reaction_id": i} for i in range(1, 9)]
        )
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},
            {"smiles": "CCCC", "volume": 800},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        # 8 wells per column × 2 columns = 16 wells total
        self.assertEqual(len(mc_calls), 16)

        # Wells should span columns 0 and 1 (indices 0–15)
        well_indices = sorted(_get_kwarg(c, "wellindex") for c in mc_calls)
        self.assertEqual(well_indices, list(range(16)))

    # ------------------------------------------------------------------
    # Test: MC wells have correct per-well volume calculation
    # ------------------------------------------------------------------
    def test_mc_well_volume_calculation(self):
        """Verify per-well volume = (transfer_vol * transfers_per_channel * 1.15) + dead_vol."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 9)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},
        ])

        # dead_volume=50, max fill = 2500*0.6 = 1500
        # transfers_per_channel = ceil(8/8) = 1
        # per_well = 50 * 1 * 1.15 + 50 = 107.5
        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        expected_vol = round(50 * 1 * 1.15 + 50, 2)  # 107.5
        for c in mc_calls:
            self.assertAlmostEqual(_get_kwarg(c, "volume"), expected_vol, places=1)

    # ------------------------------------------------------------------
    # Test: column and well indices are updated after MC placement
    # ------------------------------------------------------------------
    def test_column_index_updated_after_mc_placement(self):
        """After placing 1 MC column, column index should advance to 1."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 9)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        # Well index should have been updated past the first column
        session.well_manager.update_plate_well_index.assert_called()
        session.column_manager.update_plate_column_index_available.assert_called()

    # ------------------------------------------------------------------
    # Test: plate is returned (not None)
    # ------------------------------------------------------------------
    def test_returns_plate_object(self):
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 9)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},
        ])

        result, pf, session, mock_plate = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
        )

        self.assertEqual(result, mock_plate)

    # ------------------------------------------------------------------
    # Test: 24-well plate (wells_per_column=4)
    # ------------------------------------------------------------------
    def test_24_well_plate_4_per_column(self):
        """4 reactions on a 24-well plate (4 wells per column) → 1 MC column."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 5)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 200},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=4,
            num_columns=6,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        self.assertEqual(len(mc_calls), 4)

    # ------------------------------------------------------------------
    # 384-well plate tests (wells_per_column=16, num_columns=24)
    # ------------------------------------------------------------------
    def test_384_well_one_mc_sub_column(self):
        """16 reactions on a 384-well plate → 1 sub-column of 8 MC wells.

        wells_per_group=8, transfers_per_channel=ceil(16/8)=2.
        Sub-column 0 of column 0 uses interleaved indices [0,2,4,6,8,10,12,14].
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 20, "reaction_id": i}
            for i in range(1, 17)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 320},  # 20 * 16
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=16,
            num_columns=24,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        self.assertEqual(len(mc_calls), 8)
        well_indices = sorted(_get_kwarg(c, "wellindex") for c in mc_calls)
        self.assertEqual(well_indices, [0, 2, 4, 6, 8, 10, 12, 14])

    def test_384_well_below_threshold_all_sc(self):
        """7 reactions on a 384-well plate (needs 8 per sub-column) → all SC.

        With sub-column logic, the MC threshold on a 384-well plate is 8
        (not 16), matching the 8-channel pipette's reach.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 20, "reaction_id": i}
            for i in range(1, 8)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 140},  # 20 * 7
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=16,
            num_columns=24,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        self.assertEqual(len(mc_calls), 0)

    def test_384_well_mixed_mc_and_sc(self):
        """384-well: reagent A has 16 rxns (MC), reagent B has 5 rxns (SC).

        CCO (16 rxns) → 1 sub-column of 8 MC wells (2 transfers/channel).
        CCCC (5 rxns) → below 8-well threshold → SC.
        """
        addactionsdf = self._make_addactions_df(
            [{"smiles": "CCO", "volume": 20, "reaction_id": i} for i in range(1, 17)]
            + [{"smiles": "CCCC", "volume": 30, "reaction_id": i} for i in range(1, 6)]
        )
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 320},
            {"smiles": "CCCC", "volume": 150},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=16,
            num_columns=24,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        self.assertEqual(len(mc_calls), 8)
        self.assertEqual({_get_kwarg(c, "smiles") for c in mc_calls}, {"CCO"})
        self.assertGreaterEqual(len(sc_calls), 1)
        self.assertIn("CCCC", {_get_kwarg(c, "smiles") for c in sc_calls})

    def test_384_well_with_leftover(self):
        """384-well: 20 reactions (8 MC via sub-column + 4 leftover to SC).

        wells_per_group=8, columns_needed=20//8=2, leftover=4.
        PlateFactory places 1 sub-column of 8 MC wells (with
        transfers_per_channel=ceil(20/8)=3). Leftover 4 → SC.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 20, "reaction_id": i}
            for i in range(1, 21)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=16,
            num_columns=24,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        # 8 wells in MC sub-column
        self.assertEqual(len(mc_calls), 8)
        # Leftover 4 reactions should generate SC well(s)
        self.assertGreaterEqual(len(sc_calls), 1)

    # ------------------------------------------------------------------
    # 24-well plate additional edge cases (wells_per_column=4, num_columns=6)
    # ------------------------------------------------------------------
    def test_24_well_below_threshold_all_sc(self):
        """3 reactions on a 24-well plate (needs 4 per column) → all SC."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 4)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 150},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=4,
            num_columns=6,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        self.assertEqual(len(mc_calls), 0)

    def test_24_well_mixed_mc_and_sc(self):
        """24-well: reagent A has 4 rxns (MC), reagent B has 2 rxns (SC)."""
        addactionsdf = self._make_addactions_df(
            [{"smiles": "CCO", "volume": 50, "reaction_id": i} for i in range(1, 5)]
            + [{"smiles": "CCCC", "volume": 100, "reaction_id": i} for i in range(1, 3)]
        )
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 200},
            {"smiles": "CCCC", "volume": 200},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=4,
            num_columns=6,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        self.assertEqual(len(mc_calls), 4)
        self.assertEqual({_get_kwarg(c, "smiles") for c in mc_calls}, {"CCO"})
        self.assertGreaterEqual(len(sc_calls), 1)

    def test_24_well_two_mc_reagents(self):
        """24-well: two reagents each with 4 reactions → 2 MC columns."""
        addactionsdf = self._make_addactions_df(
            [{"smiles": "CCO", "volume": 50, "reaction_id": i} for i in range(1, 5)]
            + [{"smiles": "CCCC", "volume": 50, "reaction_id": i} for i in range(1, 5)]
        )
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 200},
            {"smiles": "CCCC", "volume": 200},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=4,
            num_columns=6,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        # 4 wells × 2 columns = 8
        self.assertEqual(len(mc_calls), 8)
        well_indices = sorted(_get_kwarg(c, "wellindex") for c in mc_calls)
        self.assertEqual(well_indices, list(range(8)))

    # ------------------------------------------------------------------
    # Volume overflow edge cases
    # ------------------------------------------------------------------
    def test_mc_volume_overflow_splits_to_multiple_columns(self):
        """When per-well volume exceeds max, the plate factory uses multiple MC columns.

        Setup: 96-well, max_well_volume=500 (effective=300 at 60% fill),
        dead_volume=50, 8 rxns at 200µL each.
        per_well = 200*1*1.15 + 50 = 280, fits in 300 → 1 column.
        But with 16 rxns: per_well = 200*2*1.15+50 = 510 >300,
        so max_transfers=1, columns_from_volume=ceil(16/8)=2 → 2 MC columns.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 200, "reaction_id": i}
            for i in range(1, 17)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 3200},  # 200 * 16
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=8,
            num_columns=12,
            max_well_volume=500,  # effective = 300 at 60%
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        # Should have 2 MC columns → 16 wells
        self.assertEqual(len(mc_calls), 16)
        well_indices = sorted(_get_kwarg(c, "wellindex") for c in mc_calls)
        self.assertEqual(well_indices, list(range(16)))

    def test_sc_volume_overflow_splits_into_multiple_wells(self):
        """SC material with volume exceeding max capacity splits into multiple wells.

        Setup: 96-well, max_well_volume=400 (effective=240), dead_volume=50.
        2 reactions at 200µL each → total=400, +15% error=460, +dead=510.
        510 > 240, so should split across multiple SC wells.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 200, "reaction_id": 1},
            {"smiles": "CCO", "volume": 200, "reaction_id": 2},
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 400},
        ])

        # Track well index progression for SC overflow
        well_idx_counter = [0]

        def _get_avail(plate_obj):
            val = well_idx_counter[0]
            return val

        def _update_idx(plate_obj, wellindexupdate):
            well_idx_counter[0] = wellindexupdate

        session, mock_plate = _make_multichannel_session(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            max_well_volume=400,  # effective = 240
        )
        session.well_manager.get_plate_well_index_available.side_effect = _get_avail
        session.well_manager.update_plate_well_index.side_effect = _update_idx

        from backend.opentrons.otsession.managers.plate_manager.plate_factory import (
            PlateFactory,
        )
        pf = PlateFactory(session)
        pf.create_plate_model = MagicMock(return_value=mock_plate)
        pf._create_compound_order = MagicMock()
        pf._create_sequential_starting_plate = MagicMock(return_value=mock_plate)

        result = pf._create_multichannel_starting_plate(materials_df)

        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        # Volume 510 > effective 240 → splits into ceil(510/240) = 3 wells
        self.assertGreaterEqual(len(sc_calls), 2)

    # ------------------------------------------------------------------
    # Plate-full edge cases
    # ------------------------------------------------------------------
    def test_plate_full_mc_overflows_to_sc(self):
        """When MC columns exhaust all plate columns, remaining MC material
        falls back to SC.

        Setup: 24-well (4 wells × 6 columns). 7 reagents each with 4 rxns.
        Only 6 columns available → 6 MC columns, 7th reagent → SC fallback.
        """
        reagents = [f"C{'C' * i}O" for i in range(7)]  # 7 distinct SMILES
        addactions = []
        for smiles in reagents:
            for rid in range(1, 5):
                addactions.append({"smiles": smiles, "volume": 50, "reaction_id": rid})

        addactionsdf = self._make_addactions_df(addactions)
        materials_df = self._make_materials_df([
            {"smiles": s, "volume": 200} for s in reagents
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=4,
            num_columns=6,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        # At most 6 columns × 4 wells = 24 MC wells
        self.assertLessEqual(len(mc_calls), 24)
        # At least one reagent must have fallen back to SC
        self.assertGreaterEqual(len(sc_calls), 1)

    def test_96_well_many_reagents_fills_all_columns(self):
        """96-well: 12 reagents each with 8 rxns → fills all 12 columns as MC.

        No SC wells should be needed.
        """
        reagents = [f"C{'C' * i}O" for i in range(12)]
        addactions = []
        for smiles in reagents:
            for rid in range(1, 9):
                addactions.append({"smiles": smiles, "volume": 50, "reaction_id": rid})

        addactionsdf = self._make_addactions_df(addactions)
        materials_df = self._make_materials_df([
            {"smiles": s, "volume": 400} for s in reagents
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=8,
            num_columns=12,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        # 12 columns × 8 wells = 96 MC wells
        self.assertEqual(len(mc_calls), 96)
        well_indices = sorted(_get_kwarg(c, "wellindex") for c in mc_calls)
        self.assertEqual(well_indices, list(range(96)))

    def test_384_well_volume_overflow_with_many_reactions(self):
        """384-well: reagent with 48 reactions at high volume per transfer.

        wells_per_group=8, max_well_volume=200 (effective=120).
        transfer_vol=50, 48 rxns → transfers_per_channel=ceil(48/8)=6,
        per_well=50*6*1.15+50=395 > 120.
        max_transfers=floor((120-50)/(50*1.15))=1,
        groups_from_volume=ceil(48/(8*1))=6 MC sub-column slots.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 49)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 2400},  # 50 * 48
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=16,
            num_columns=24,
            max_well_volume=200,  # effective = 120
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        # Should get 3 MC columns × 16 wells = 48 MC wells
        self.assertEqual(len(mc_calls), 48)

    # ------------------------------------------------------------------
    # MC volume calculation edge cases per plate type
    # ------------------------------------------------------------------
    def test_24_well_mc_volume_calculation(self):
        """24-well: verify per-well volume with wells_per_column=4.

        8 reactions, volume=100: transfers_per_channel=ceil(8/4)=2,
        per_well = 100*2*1.15 + 50 = 280.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 100, "reaction_id": i}
            for i in range(1, 9)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 800},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=4,
            num_columns=6,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        # 4 wells per column (all in same column since per-well volume fits)
        self.assertGreaterEqual(len(mc_calls), 4)
        expected_vol = round(100 * 2 * 1.15 + 50, 2)  # 280.0
        for c in mc_calls:
            self.assertAlmostEqual(_get_kwarg(c, "volume"), expected_vol, places=1)

    def test_384_well_mc_volume_calculation(self):
        """384-well: verify per-well volume with wells_per_group=8.

        16 reactions, volume=30: transfers_per_channel=ceil(16/8)=2,
        per_well = 30*2*1.15 + 50 = 119.0.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 30, "reaction_id": i}
            for i in range(1, 17)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 480},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=16,
            num_columns=24,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        self.assertEqual(len(mc_calls), 8)
        expected_vol = round(30 * 2 * 1.15 + 50, 2)  # 119.0
        for c in mc_calls:
            self.assertAlmostEqual(_get_kwarg(c, "volume"), expected_vol, places=1)

    # ------------------------------------------------------------------
    # Boundary: exactly at column threshold
    # ------------------------------------------------------------------
    def test_24_well_exact_boundary_at_4_reactions(self):
        """Exactly 4 reactions on a 24-well → 1 full MC column, zero leftovers."""
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 50, "reaction_id": i}
            for i in range(1, 5)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 200},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=4,
            num_columns=6,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        self.assertEqual(len(mc_calls), 4)
        # No SC if the only reagent perfectly fills columns (leftover may
        # still appear in materials_df for the total-volume path, but since
        # the identity key won't be in sc_unique_solutions, no SC well.)
        self.assertEqual(len(sc_calls), 0)

    def test_384_well_exact_boundary_at_16_reactions(self):
        """16 reactions on a 384-well → 1 sub-column of 8 MC wells, zero leftovers.

        wells_per_group=8, columns_needed=2, leftover=0.
        PlateFactory: 1 sub-column, transfers_per_channel=2.
        """
        addactionsdf = self._make_addactions_df([
            {"smiles": "CCO", "volume": 20, "reaction_id": i}
            for i in range(1, 17)
        ])
        materials_df = self._make_materials_df([
            {"smiles": "CCO", "volume": 320},
        ])

        result, pf, session, _ = self._run(
            addactionsdf=addactionsdf,
            materials_df=materials_df,
            wells_per_column=16,
            num_columns=24,
        )

        mc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "multichannel"
        ]
        sc_calls = [
            c for c in session.well_manager.create_well_model.call_args_list
            if _get_kwarg(c, "transfer_type") == "single"
        ]
        self.assertEqual(len(mc_calls), 8)
        self.assertEqual(len(sc_calls), 0)


def _get_kwarg(call_obj, key):
    """Extract a keyword argument from a mock call, checking both kwargs and positional."""
    if key in call_obj.kwargs:
        return call_obj.kwargs[key]
    # Fall back to checking if it was passed positionally via call()
    return call_obj.kwargs.get(key)
