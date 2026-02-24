"""Tests for the MultichannelAnalyzer module.

Tests for:
  - WellInfo / MultichannelGroup / CherryPickWell / MaterialGroup dataclasses
  - MultichannelAnalyzer.make_identity_key
  - MultichannelAnalyzer.get_column_index_for_well
  - MultichannelAnalyzer.get_well_indices_for_column
  - MultichannelAnalyzer.analyze_plate (96-, 384-, 24-well)
  - MultichannelAnalyzer._analyze_column (eligibility rules)
  - MultichannelAnalyzer.prioritize_materials
  - MultichannelAnalyzer.plan_plate_layout
  - MultichannelAnalyzer.analyze_custom_plate_csv
  - Edge cases: empty plates, boundary counts, overflow, mixed solvents

No external dependencies are required — MultichannelAnalyzer is pure Python.
"""

from unittest import TestCase
import math

from backend.opentrons.otsession.managers.multichannel_analyzer import (
    MultichannelAnalyzer,
    WellInfo,
    MultichannelGroup,
    CherryPickWell,
    ColumnAnalysis,
    MaterialGroup,
    PlateAnalysisResult,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_uniform_column(
    col_index,
    smiles="CCO",
    volume=50.0,
    wells_per_column=8,
    solvent=None,
    concentration=None,
):
    """Create a list of WellInfo filling one full column with identical reagent."""
    start = col_index * wells_per_column
    return [
        WellInfo(
            index=start + i,
            smiles=smiles,
            volume=volume,
            solvent=solvent,
            concentration=concentration,
        )
        for i in range(wells_per_column)
    ]


def _make_partial_column(col_index, count, smiles="CCO", volume=50.0, wells_per_column=8):
    """Create a partial column with *count* wells (< wells_per_column)."""
    start = col_index * wells_per_column
    return [
        WellInfo(index=start + i, smiles=smiles, volume=volume)
        for i in range(count)
    ]


def _make_uniform_sub_column(
    col_index,
    sub_col_index,
    smiles="CCO",
    volume=50.0,
    wells_per_column=16,
    pipette_channels=8,
):
    """Create WellInfo entries for one sub-column of a 384-well plate.

    sub_col_index 0 → even row indices (0, 2, 4, …)
    sub_col_index 1 → odd  row indices (1, 3, 5, …)
    """
    sub_cols = max(1, wells_per_column // pipette_channels)
    group_size = min(wells_per_column, pipette_channels)
    col_start = col_index * wells_per_column
    return [
        WellInfo(
            index=col_start + sub_col_index + i * sub_cols,
            smiles=smiles,
            volume=volume,
        )
        for i in range(group_size)
    ]


# ===================================================================
# Identity key tests
# ===================================================================
class TestMakeIdentityKey(TestCase):
    """MultichannelAnalyzer.make_identity_key static method."""

    def test_basic_key(self):
        key = MultichannelAnalyzer.make_identity_key("CCO", "DMSO", 10.0)
        self.assertEqual(key, "CCO-DMSO-10.0")

    def test_none_solvent_and_concentration(self):
        key = MultichannelAnalyzer.make_identity_key("CCO", None, None)
        self.assertEqual(key, "CCO--")

    def test_same_smiles_different_solvent(self):
        k1 = MultichannelAnalyzer.make_identity_key("CCO", "DMSO", 10.0)
        k2 = MultichannelAnalyzer.make_identity_key("CCO", "water", 10.0)
        self.assertNotEqual(k1, k2)

    def test_same_smiles_different_concentration(self):
        k1 = MultichannelAnalyzer.make_identity_key("CCO", "DMSO", 10.0)
        k2 = MultichannelAnalyzer.make_identity_key("CCO", "DMSO", 20.0)
        self.assertNotEqual(k1, k2)


# ===================================================================
# Column / well index helpers
# ===================================================================
class TestColumnHelpers(TestCase):
    """get_column_index_for_well and get_well_indices_for_column."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=8, num_columns=12)

    # --- get_column_index_for_well ---
    def test_first_well_column_0(self):
        self.assertEqual(self.analyzer.get_column_index_for_well(0), 0)

    def test_last_well_column_0(self):
        self.assertEqual(self.analyzer.get_column_index_for_well(7), 0)

    def test_first_well_column_1(self):
        self.assertEqual(self.analyzer.get_column_index_for_well(8), 1)

    def test_last_well_on_plate(self):
        self.assertEqual(self.analyzer.get_column_index_for_well(95), 11)

    def test_384_well_plate(self):
        a384 = MultichannelAnalyzer(wells_per_column=16, num_columns=24)
        self.assertEqual(a384.get_column_index_for_well(0), 0)
        self.assertEqual(a384.get_column_index_for_well(15), 0)
        self.assertEqual(a384.get_column_index_for_well(16), 1)
        self.assertEqual(a384.get_column_index_for_well(383), 23)

    def test_384_sub_column_indices(self):
        """384-well sub-columns: even/odd row interleaving."""
        a384 = MultichannelAnalyzer(wells_per_column=16, num_columns=24)
        # Sub-column 0 of column 0 → even rows: 0,2,4,6,8,10,12,14
        self.assertEqual(
            a384.get_sub_column_well_indices(0, 0),
            [0, 2, 4, 6, 8, 10, 12, 14],
        )
        # Sub-column 1 of column 0 → odd rows: 1,3,5,7,9,11,13,15
        self.assertEqual(
            a384.get_sub_column_well_indices(0, 1),
            [1, 3, 5, 7, 9, 11, 13, 15],
        )
        # Sub-column 0 of column 1 → 16,18,20,22,24,26,28,30
        self.assertEqual(
            a384.get_sub_column_well_indices(1, 0),
            [16, 18, 20, 22, 24, 26, 28, 30],
        )

    def test_96_well_sub_column_same_as_full_column(self):
        """96-well: sub-column 0 == full column (sub_columns_per_column=1)."""
        a96 = MultichannelAnalyzer(wells_per_column=8, num_columns=12)
        self.assertEqual(
            a96.get_sub_column_well_indices(0, 0),
            list(range(0, 8)),
        )
        self.assertEqual(
            a96.get_sub_column_well_indices(2, 0),
            list(range(16, 24)),
        )

    # --- get_well_indices_for_column ---
    def test_column_0_indices(self):
        self.assertEqual(
            self.analyzer.get_well_indices_for_column(0), list(range(0, 8))
        )

    def test_column_11_indices(self):
        self.assertEqual(
            self.analyzer.get_well_indices_for_column(11), list(range(88, 96))
        )

    def test_24_well_plate(self):
        a24 = MultichannelAnalyzer(wells_per_column=4, num_columns=6)
        self.assertEqual(a24.get_well_indices_for_column(0), [0, 1, 2, 3])
        self.assertEqual(a24.get_well_indices_for_column(5), [20, 21, 22, 23])


# ===================================================================
# analyze_plate — 96-well standard cases
# ===================================================================
class TestAnalyzePlate96Well(TestCase):
    """analyze_plate on a standard 96-well plate (8 rows × 12 columns)."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=8, num_columns=12)

    # --- Full uniform column → multichannel eligible ---
    def test_single_full_uniform_column(self):
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(result.multichannel_groups[0].column_index, 0)
        self.assertEqual(result.multichannel_groups[0].smiles, "CCO")
        self.assertEqual(result.multichannel_groups[0].volume, 50.0)
        self.assertEqual(result.multichannel_well_count, 8)
        self.assertEqual(result.cherry_pick_well_count, 0)
        self.assertAlmostEqual(result.efficiency, 1.0)

    def test_two_uniform_columns_different_reagent(self):
        wells = (
            _make_uniform_column(0, smiles="CCO", volume=50.0)
            + _make_uniform_column(1, smiles="CCC", volume=100.0)
        )
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 2)
        self.assertEqual(result.multichannel_well_count, 16)
        self.assertEqual(result.total_wells, 16)
        self.assertAlmostEqual(result.efficiency, 1.0)

    def test_uniform_column_with_solvent_and_concentration(self):
        wells = _make_uniform_column(
            0, smiles="CCO", volume=50.0, solvent="DMSO", concentration=10.0
        )
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(result.multichannel_groups[0].solvent, "DMSO")
        self.assertEqual(result.multichannel_groups[0].concentration, 10.0)

    # --- Partial column → cherry picks ---
    def test_partial_column(self):
        wells = _make_partial_column(0, count=5)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 5)
        for cp in result.cherry_pick_wells:
            self.assertEqual(cp.reason, "incomplete_column")

    def test_single_well(self):
        wells = [WellInfo(index=0, smiles="CCO", volume=50.0)]
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 1)

    # --- Mixed reagents in column → cherry picks ---
    def test_mixed_smiles_in_column(self):
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0)
        wells[3] = WellInfo(index=3, smiles="DIFFERENT", volume=50.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 8)
        self.assertEqual(result.column_analyses[0].reason, "mixed_reagents")

    def test_mixed_solvent_in_column(self):
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0, solvent="DMSO")
        wells[7] = WellInfo(index=7, smiles="CCO", volume=50.0, solvent="water")
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.column_analyses[0].reason, "mixed_reagents")

    def test_mixed_concentration_in_column(self):
        wells = _make_uniform_column(
            0, smiles="CCO", volume=50.0, concentration=10.0
        )
        wells[0] = WellInfo(
            index=0, smiles="CCO", volume=50.0, concentration=20.0
        )
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.column_analyses[0].reason, "mixed_reagents")

    # --- Mixed volumes in column → cherry picks ---
    def test_mixed_volumes_in_column(self):
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0)
        wells[4] = WellInfo(index=4, smiles="CCO", volume=75.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.column_analyses[0].reason, "mixed_volumes")

    def test_volume_within_tolerance(self):
        """Very small volume differences should still be considered uniform."""
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0)
        wells[0] = WellInfo(index=0, smiles="CCO", volume=50.005)  # 0.01% diff
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 1)

    def test_volume_outside_tolerance(self):
        """Volumes differing by >0.1% should fail."""
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0)
        wells[0] = WellInfo(index=0, smiles="CCO", volume=50.5)  # 1% diff
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.column_analyses[0].reason, "mixed_volumes")

    # --- None / empty SMILES ---
    def test_none_smiles_in_full_column(self):
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0)
        wells[7] = WellInfo(index=7, smiles=None, volume=50.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.column_analyses[0].reason, "empty_well_contents")

    def test_empty_string_smiles_in_full_column(self):
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0)
        wells[0] = WellInfo(index=0, smiles="", volume=50.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.column_analyses[0].reason, "empty_well_contents")

    # --- Empty plate ---
    def test_empty_plate(self):
        result = self.analyzer.analyze_plate([])
        self.assertEqual(result.total_wells, 0)
        self.assertAlmostEqual(result.efficiency, 0.0)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(len(result.cherry_pick_wells), 0)

    # --- Mixed: some columns multichannel, some cherry-pick ---
    def test_mixed_plate(self):
        """Two MC columns + one partial column."""
        wells = (
            _make_uniform_column(0, smiles="CCO", volume=50.0)
            + _make_uniform_column(1, smiles="CCC", volume=100.0)
            + _make_partial_column(2, count=4, smiles="CCCC", volume=30.0)
        )
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 2)
        self.assertEqual(result.multichannel_well_count, 16)
        self.assertEqual(result.cherry_pick_well_count, 4)
        self.assertEqual(result.total_wells, 20)
        self.assertAlmostEqual(result.efficiency, 16 / 20)

    def test_full_plate_all_multichannel(self):
        """All 12 columns filled identically → 100% multichannel."""
        wells = []
        for col in range(12):
            wells += _make_uniform_column(col, smiles=f"SMILES_{col}", volume=50.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 12)
        self.assertEqual(result.multichannel_well_count, 96)
        self.assertAlmostEqual(result.efficiency, 1.0)

    def test_full_plate_no_multichannel(self):
        """96 wells, every one with a different SMILES → 0% multichannel."""
        wells = [
            WellInfo(index=i, smiles=f"S{i}", volume=50.0) for i in range(96)
        ]
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 96)
        self.assertAlmostEqual(result.efficiency, 0.0)


# ===================================================================
# analyze_plate — 384- and 24-well plates
# ===================================================================
class TestAnalyzePlate384Well(TestCase):
    """analyze_plate on a 384-well plate (16 rows × 24 columns).

    The 8-channel pipette reaches every other row, so each physical
    column contains two sub-columns of 8 wells.
    """

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=16, num_columns=24)

    def test_full_uniform_column_produces_two_mc_groups(self):
        """All 16 wells identical → 2 MC groups (one per sub-column)."""
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0, wells_per_column=16)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 2)
        self.assertEqual(len(result.multichannel_groups[0].well_indices), 8)
        self.assertEqual(len(result.multichannel_groups[1].well_indices), 8)
        self.assertEqual(result.multichannel_well_count, 16)
        self.assertAlmostEqual(result.efficiency, 1.0)

    def test_sub_column_indices_are_interleaved(self):
        """MC groups from a 384-well column use interleaved indices."""
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0, wells_per_column=16)
        result = self.analyzer.analyze_plate(wells)
        g0 = result.multichannel_groups[0]
        g1 = result.multichannel_groups[1]
        self.assertEqual(g0.well_indices, [0, 2, 4, 6, 8, 10, 12, 14])
        self.assertEqual(g1.well_indices, [1, 3, 5, 7, 9, 11, 13, 15])
        self.assertEqual(g0.sub_column_index, 0)
        self.assertEqual(g1.sub_column_index, 1)

    def test_only_even_sub_column_filled(self):
        """8 wells at even rows → 1 MC group; odd rows empty → nothing."""
        wells = _make_uniform_sub_column(0, sub_col_index=0, smiles="CCO", volume=50.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(result.multichannel_groups[0].sub_column_index, 0)
        self.assertEqual(len(result.multichannel_groups[0].well_indices), 8)

    def test_only_odd_sub_column_filled(self):
        """8 wells at odd rows → 1 MC group."""
        wells = _make_uniform_sub_column(0, sub_col_index=1, smiles="CCO", volume=50.0)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(result.multichannel_groups[0].sub_column_index, 1)

    def test_even_uniform_odd_mixed_reagents(self):
        """Even sub-column uniform → MC; odd sub-column mixed → cherry-picks."""
        even_wells = _make_uniform_sub_column(0, sub_col_index=0, smiles="CCO", volume=50.0)
        odd_wells = _make_uniform_sub_column(0, sub_col_index=1, smiles="CCO", volume=50.0)
        # Break odd sub-column uniformity
        odd_wells[3] = WellInfo(index=odd_wells[3].index, smiles="DIFFERENT", volume=50.0)
        result = self.analyzer.analyze_plate(even_wells + odd_wells)
        self.assertEqual(len(result.multichannel_groups), 1)  # even only
        self.assertEqual(result.multichannel_groups[0].sub_column_index, 0)
        self.assertEqual(result.cherry_pick_well_count, 8)  # odd sub-column

    def test_partial_sub_column(self):
        """Fewer than 8 wells in a sub-column → incomplete, cherry-picks."""
        # Only 5 wells at even positions in column 0
        col_start = 0
        wells = [
            WellInfo(index=col_start + i * 2, smiles="CCO", volume=50.0)
            for i in range(5)
        ]
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 5)

    def test_partial_column_12_wells(self):
        """12 of 16 wells → potentially one sub-col full, one partial."""
        # Wells 0-11 of column 0
        wells = _make_partial_column(0, count=12, wells_per_column=16)
        result = self.analyzer.analyze_plate(wells)
        # Even sub-col: indices 0,2,4,6,8,10 = 6 wells (incomplete)
        # Odd sub-col: indices 1,3,5,7,9,11 = 6 wells (incomplete)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 12)


class TestAnalyzePlate24Well(TestCase):
    """analyze_plate on a 24-well plate (4 rows × 6 columns)."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=4, num_columns=6)

    def test_full_column_24(self):
        wells = _make_uniform_column(0, smiles="CCO", volume=50.0, wells_per_column=4)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(len(result.multichannel_groups[0].well_indices), 4)

    def test_partial_column_24_three_wells(self):
        wells = _make_partial_column(0, count=3, smiles="CCO", volume=50.0, wells_per_column=4)
        result = self.analyzer.analyze_plate(wells)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 3)


# ===================================================================
# prioritize_materials
# ===================================================================
class TestPrioritizeMaterials(TestCase):
    """MultichannelAnalyzer.prioritize_materials."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=8, num_columns=12)

    def test_single_material_fills_one_column(self):
        materials = [{"smiles": "CCO", "volume": 50.0, "reaction_count": 8}]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc), 1)
        self.assertEqual(mc[0].columns_needed, 1)
        self.assertEqual(mc[0].leftover_wells, 0)
        self.assertEqual(len(sc), 0)

    def test_single_material_too_few_reactions(self):
        materials = [{"smiles": "CCO", "volume": 50.0, "reaction_count": 7}]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc), 0)
        self.assertEqual(len(sc), 1)

    def test_material_with_leftover(self):
        materials = [{"smiles": "CCO", "volume": 50.0, "reaction_count": 20}]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc), 1)
        self.assertEqual(mc[0].columns_needed, 2)  # 16 in columns
        self.assertEqual(mc[0].leftover_wells, 4)  # 4 remaining

    def test_multiple_materials_mixed(self):
        """Some materials fill columns, others don't."""
        materials = [
            {"smiles": "A", "volume": 50.0, "reaction_count": 20},   # 2 cols + 4
            {"smiles": "B", "volume": 50.0, "reaction_count": 16},   # 2 cols
            {"smiles": "C", "volume": 50.0, "reaction_count": 3},    # SC
            {"smiles": "D", "volume": 100.0, "reaction_count": 5},   # SC
        ]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc), 2)  # A and B
        self.assertEqual(len(sc), 2)  # C and D

    def test_priority_ordering(self):
        """Materials with more full columns come first."""
        materials = [
            {"smiles": "A", "volume": 50.0, "reaction_count": 16},   # 2 cols
            {"smiles": "B", "volume": 50.0, "reaction_count": 24},   # 3 cols
            {"smiles": "C", "volume": 50.0, "reaction_count": 8},    # 1 col
        ]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(mc[0].smiles, "B")  # 3 cols first
        self.assertEqual(mc[1].smiles, "A")  # 2 cols second
        self.assertEqual(mc[2].smiles, "C")  # 1 col third

    def test_same_smiles_different_volume_separate_groups(self):
        """Same SMILES but different volumes → separate groups."""
        materials = [
            {"smiles": "CCO", "volume": 50.0, "reaction_count": 8},
            {"smiles": "CCO", "volume": 100.0, "reaction_count": 8},
        ]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc), 2)
        volumes = {m.volume for m in mc}
        self.assertEqual(volumes, {50.0, 100.0})

    def test_same_smiles_same_volume_aggregated(self):
        """Same SMILES + same volume from multiple entries → merged."""
        materials = [
            {"smiles": "CCO", "volume": 50.0, "reaction_count": 5},
            {"smiles": "CCO", "volume": 50.0, "reaction_count": 5},
        ]
        mc, sc = self.analyzer.prioritize_materials(materials)
        # 5 + 5 = 10 reactions → 1 full column (8) + 2 leftover
        self.assertEqual(len(mc), 1)
        self.assertEqual(mc[0].reaction_count, 10)
        self.assertEqual(mc[0].columns_needed, 1)
        self.assertEqual(mc[0].leftover_wells, 2)

    def test_solvent_differentiates_groups(self):
        """Same SMILES+volume but different solvents → separate groups."""
        materials = [
            {"smiles": "CCO", "volume": 50.0, "reaction_count": 8, "solvent": "DMSO"},
            {"smiles": "CCO", "volume": 50.0, "reaction_count": 8, "solvent": "water"},
        ]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc), 2)

    def test_all_below_threshold(self):
        """No materials fill a column → all single-channel."""
        materials = [
            {"smiles": "A", "volume": 50.0, "reaction_count": 3},
            {"smiles": "B", "volume": 50.0, "reaction_count": 5},
        ]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc), 0)
        self.assertEqual(len(sc), 2)
        # Single-channel sorted by reaction_count descending
        self.assertEqual(sc[0].smiles, "B")
        self.assertEqual(sc[1].smiles, "A")

    def test_empty_materials_list(self):
        mc, sc = self.analyzer.prioritize_materials([])
        self.assertEqual(len(mc), 0)
        self.assertEqual(len(sc), 0)

    def test_exact_column_boundary(self):
        """Exactly wells_per_column reactions: 1 column, 0 leftovers."""
        materials = [{"smiles": "CCO", "volume": 50.0, "reaction_count": 8}]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(mc[0].columns_needed, 1)
        self.assertEqual(mc[0].leftover_wells, 0)

    def test_one_more_than_column(self):
        """wells_per_column + 1 reactions: 1 column + 1 leftover."""
        materials = [{"smiles": "CCO", "volume": 50.0, "reaction_count": 9}]
        mc, sc = self.analyzer.prioritize_materials(materials)
        self.assertEqual(mc[0].columns_needed, 1)
        self.assertEqual(mc[0].leftover_wells, 1)


# ===================================================================
# plan_plate_layout
# ===================================================================
class TestPlanPlateLayout(TestCase):
    """MultichannelAnalyzer.plan_plate_layout."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=8, num_columns=12)

    def test_only_multichannel(self):
        mc_materials = [
            MaterialGroup(
                identity_key="CCO--", smiles="CCO", volume=50.0,
                columns_needed=2, leftover_wells=0, reaction_count=16,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(mc_materials, [])
        self.assertEqual(len(mc_groups), 2)
        self.assertEqual(mc_groups[0].column_index, 0)
        self.assertEqual(mc_groups[1].column_index, 1)
        self.assertEqual(len(cp_wells), 0)

    def test_only_single_channel(self):
        sc_materials = [
            MaterialGroup(
                identity_key="CCO--", smiles="CCO", volume=50.0,
                columns_needed=0, leftover_wells=0, reaction_count=3,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout([], sc_materials)
        self.assertEqual(len(mc_groups), 0)
        self.assertEqual(len(cp_wells), 3)
        # Sequential filling starts at index 0 when no MC columns
        self.assertEqual(cp_wells[0].well_index, 0)
        self.assertEqual(cp_wells[1].well_index, 1)
        self.assertEqual(cp_wells[2].well_index, 2)

    def test_mixed_mc_and_sc(self):
        """MC columns first, then SC wells fill remaining space."""
        mc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=2, leftover_wells=3, reaction_count=19,
            ),
        ]
        sc_materials = [
            MaterialGroup(
                identity_key="B--", smiles="B", volume=100.0,
                columns_needed=0, leftover_wells=0, reaction_count=5,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(
            mc_materials, sc_materials
        )
        self.assertEqual(len(mc_groups), 2)
        # 3 leftover from A + 5 from B = 8 cherry picks
        self.assertEqual(len(cp_wells), 8)
        # Cherry picks start after MC columns (2 columns × 8 = index 16)
        self.assertEqual(cp_wells[0].well_index, 16)
        # First 3 cherry picks should be A (leftover)
        for cp in cp_wells[:3]:
            self.assertEqual(cp.smiles, "A")
        # Next 5 should be B
        for cp in cp_wells[3:]:
            self.assertEqual(cp.smiles, "B")

    def test_column_overflow_handled(self):
        """When MC materials require more columns than available."""
        mc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=15, leftover_wells=0, reaction_count=120,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(mc_materials, [])
        # Only 12 columns fit on a 96-well plate
        self.assertEqual(len(mc_groups), 12)

    def test_plate_capacity_limits_cherry_picks(self):
        """SC wells stop when plate is full."""
        sc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=0, leftover_wells=0, reaction_count=200,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout([], sc_materials)
        # 96-well plate → max 96 wells
        self.assertEqual(len(cp_wells), 96)

    def test_empty_inputs(self):
        mc_groups, cp_wells = self.analyzer.plan_plate_layout([], [])
        self.assertEqual(len(mc_groups), 0)
        self.assertEqual(len(cp_wells), 0)

    def test_multi_material_column_assignment(self):
        """Multiple MC materials each get their columns in priority order."""
        mc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=3, leftover_wells=0, reaction_count=24,
            ),
            MaterialGroup(
                identity_key="B--", smiles="B", volume=100.0,
                columns_needed=2, leftover_wells=0, reaction_count=16,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(mc_materials, [])
        self.assertEqual(len(mc_groups), 5)
        # A gets columns 0-2
        for i in range(3):
            self.assertEqual(mc_groups[i].smiles, "A")
            self.assertEqual(mc_groups[i].column_index, i)
        # B gets columns 3-4
        for i in range(3, 5):
            self.assertEqual(mc_groups[i].smiles, "B")
            self.assertEqual(mc_groups[i].column_index, i)


class TestPlanPlateLayout384Well(TestCase):
    """plan_plate_layout on 384-well plates (16 rows, 2 sub-columns)."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(
            wells_per_column=16, num_columns=24
        )

    def test_even_mc_groups_384(self):
        """Two MC groups fill both sub-columns of the first physical column."""
        mc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=2, leftover_wells=0, reaction_count=16,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(mc_materials, [])
        self.assertEqual(len(mc_groups), 2)
        # Both groups in physical column 0, sub-columns 0 and 1
        self.assertEqual(mc_groups[0].column_index, 0)
        self.assertEqual(mc_groups[0].sub_column_index, 0)
        self.assertEqual(mc_groups[1].column_index, 0)
        self.assertEqual(mc_groups[1].sub_column_index, 1)
        self.assertEqual(len(cp_wells), 0)

    def test_odd_mc_groups_no_sc_overlap(self):
        """Three MC groups: 1 physical column full + 1 partial.

        SC wells must start AFTER the partial column, never in it.
        This is the key regression test for the Phase B overlap bug.
        """
        mc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=3, leftover_wells=2, reaction_count=26,
            ),
        ]
        sc_materials = [
            MaterialGroup(
                identity_key="B--", smiles="B", volume=100.0,
                columns_needed=0, leftover_wells=0, reaction_count=4,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(
            mc_materials, sc_materials
        )
        self.assertEqual(len(mc_groups), 3)
        # Columns 0-sub0, 0-sub1, 1-sub0
        self.assertEqual(mc_groups[2].column_index, 1)
        self.assertEqual(mc_groups[2].sub_column_index, 0)

        # SC wells must start at column 2 (index 32), NOT column 1
        mc_well_indices = set()
        for g in mc_groups:
            mc_well_indices.update(g.well_indices)

        for cp in cp_wells:
            self.assertNotIn(
                cp.well_index, mc_well_indices,
                f"SC well {cp.well_index} overlaps with MC wells {sorted(mc_well_indices)}",
            )
        # First SC well should be at column 2 (index 32)
        self.assertEqual(cp_wells[0].well_index, 32)

    def test_single_mc_group_384(self):
        """A single MC group occupies sub-column 0 of column 0.

        SC must start at column 1 (index 16), not column 0.
        """
        mc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=1, leftover_wells=3, reaction_count=11,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(mc_materials, [])
        self.assertEqual(len(mc_groups), 1)
        self.assertEqual(mc_groups[0].column_index, 0)
        # 3 leftover → first SC well at column 1 (index 16)
        self.assertEqual(len(cp_wells), 3)
        self.assertEqual(cp_wells[0].well_index, 16)

    def test_no_mc_groups_384(self):
        """Pure SC on 384-well — starts at well 0."""
        sc_materials = [
            MaterialGroup(
                identity_key="X--", smiles="X", volume=50.0,
                columns_needed=0, leftover_wells=0, reaction_count=5,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout([], sc_materials)
        self.assertEqual(len(mc_groups), 0)
        self.assertEqual(cp_wells[0].well_index, 0)

    def test_mc_sc_well_indices_never_overlap(self):
        """No MC well index should ever appear in the SC well list."""
        mc_materials = [
            MaterialGroup(
                identity_key="A--", smiles="A", volume=50.0,
                columns_needed=5, leftover_wells=4, reaction_count=44,
            ),
            MaterialGroup(
                identity_key="B--", smiles="B", volume=75.0,
                columns_needed=1, leftover_wells=6, reaction_count=14,
            ),
        ]
        sc_materials = [
            MaterialGroup(
                identity_key="C--", smiles="C", volume=25.0,
                columns_needed=0, leftover_wells=0, reaction_count=10,
            ),
        ]
        mc_groups, cp_wells = self.analyzer.plan_plate_layout(
            mc_materials, sc_materials
        )
        mc_indices = set()
        for g in mc_groups:
            mc_indices.update(g.well_indices)
        sc_indices = {cp.well_index for cp in cp_wells}
        overlap = mc_indices & sc_indices
        self.assertEqual(
            overlap, set(),
            f"MC and SC well indices overlap: {overlap}",
        )


# ===================================================================
# analyze_custom_plate_csv
# ===================================================================
class TestAnalyzeCustomPlateCsv(TestCase):
    """MultichannelAnalyzer.analyze_custom_plate_csv."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=8, num_columns=12)

    def _make_csv_rows(self, col_index, smiles, volume, count=8, solvent=None):
        """Simulate parsed CSV rows for one column."""
        start = col_index * 8
        rows = []
        for i in range(count):
            row = {
                "well-index": str(start + i),
                "SMILES": smiles,
                "amount-uL": str(volume),
            }
            if solvent:
                row["solvent"] = solvent
            rows.append(row)
        return rows

    def test_full_multichannel_column(self):
        rows = self._make_csv_rows(0, "CCO", 50.0)
        result = self.analyzer.analyze_custom_plate_csv(rows)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(result.multichannel_well_count, 8)

    def test_partial_column_csv(self):
        rows = self._make_csv_rows(0, "CCO", 50.0, count=5)
        result = self.analyzer.analyze_custom_plate_csv(rows)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.cherry_pick_well_count, 5)

    def test_mixed_csv(self):
        """One MC column + one partial column."""
        rows = (
            self._make_csv_rows(0, "CCO", 50.0, count=8)
            + self._make_csv_rows(1, "CCC", 100.0, count=5)
        )
        result = self.analyzer.analyze_custom_plate_csv(rows)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(result.cherry_pick_well_count, 5)
        self.assertEqual(result.total_wells, 13)

    def test_csv_with_concentration(self):
        rows = self._make_csv_rows(0, "CCO", 50.0)
        for r in rows:
            r["concentration"] = "10.0"
        result = self.analyzer.analyze_custom_plate_csv(rows)
        self.assertEqual(len(result.multichannel_groups), 1)
        self.assertEqual(result.multichannel_groups[0].concentration, 10.0)

    def test_csv_with_mixed_concentration(self):
        rows = self._make_csv_rows(0, "CCO", 50.0)
        rows[0]["concentration"] = "10.0"
        rows[1]["concentration"] = "20.0"
        result = self.analyzer.analyze_custom_plate_csv(rows)
        self.assertEqual(len(result.multichannel_groups), 0)
        self.assertEqual(result.column_analyses[0].reason, "mixed_reagents")

    def test_csv_override_plate_dimensions(self):
        """Override to 384-well dimensions.

        All 16 wells in column 0 → 2 MC groups (one per sub-column of 8).
        """
        rows = [
            {"well-index": str(i), "SMILES": "CCO", "amount-uL": "50.0"}
            for i in range(16)
        ]
        result = self.analyzer.analyze_custom_plate_csv(
            rows, wells_per_column=16, num_columns=24
        )
        self.assertEqual(len(result.multichannel_groups), 2)
        self.assertEqual(len(result.multichannel_groups[0].well_indices), 8)
        self.assertEqual(len(result.multichannel_groups[1].well_indices), 8)

    def test_csv_restores_original_dimensions(self):
        """After calling with overrides, original dimensions must be restored."""
        self.analyzer.analyze_custom_plate_csv(
            [{"well-index": "0", "SMILES": "CCO", "amount-uL": "50.0"}],
            wells_per_column=16,
            num_columns=24,
        )
        self.assertEqual(self.analyzer.wells_per_column, 8)
        self.assertEqual(self.analyzer.num_columns, 12)
        # Derived sub-column properties must also be restored
        self.assertEqual(self.analyzer.wells_per_group, 8)
        self.assertEqual(self.analyzer.sub_columns_per_column, 1)

    def test_csv_identical_to_wellinfo_analysis(self):
        """CSV analysis must produce the same result as direct WellInfo analysis."""
        csv_rows = self._make_csv_rows(0, "CCO", 50.0) + self._make_csv_rows(
            1, "CCC", 100.0, count=4
        )
        csv_result = self.analyzer.analyze_custom_plate_csv(csv_rows)

        # Build equivalent WellInfo list
        well_infos = (
            _make_uniform_column(0, smiles="CCO", volume=50.0)
            + _make_partial_column(1, count=4, smiles="CCC", volume=100.0)
        )
        direct_result = self.analyzer.analyze_plate(well_infos)

        self.assertEqual(
            csv_result.multichannel_well_count, direct_result.multichannel_well_count
        )
        self.assertEqual(
            csv_result.cherry_pick_well_count, direct_result.cherry_pick_well_count
        )
        self.assertAlmostEqual(csv_result.efficiency, direct_result.efficiency)


# ===================================================================
# End-to-end: prioritize → plan → analyze round-trip
# ===================================================================
class TestEndToEndPipeline(TestCase):
    """Verifies the full pipeline: prioritize → plan → analyze round-trip."""

    def setUp(self):
        self.analyzer = MultichannelAnalyzer(wells_per_column=8, num_columns=12)

    def test_round_trip(self):
        """Materials → prioritise → plan layout → analyze the result."""
        materials = [
            {"smiles": "A", "volume": 50.0, "reaction_count": 20},   # 2 cols + 4
            {"smiles": "B", "volume": 50.0, "reaction_count": 16},   # 2 cols
            {"smiles": "C", "volume": 100.0, "reaction_count": 3},   # SC
        ]
        mc_materials, sc_materials = self.analyzer.prioritize_materials(materials)

        mc_groups, cp_wells = self.analyzer.plan_plate_layout(
            mc_materials, sc_materials
        )

        # Build WellInfo from the planned layout
        well_infos = []
        for g in mc_groups:
            for idx in g.well_indices:
                well_infos.append(
                    WellInfo(index=idx, smiles=g.smiles, volume=g.volume)
                )
        for cp in cp_wells:
            well_infos.append(
                WellInfo(index=cp.well_index, smiles=cp.smiles, volume=cp.volume)
            )

        # Analyze the resulting plate
        result = self.analyzer.analyze_plate(well_infos)

        # The MC columns should round-trip correctly
        self.assertEqual(result.multichannel_well_count, 4 * 8)  # 2+2 cols = 32
        # A leftover (4) + C (3) = 7 cherry picks
        self.assertEqual(result.cherry_pick_well_count, 7)
        self.assertEqual(result.total_wells, 39)

    def test_round_trip_all_single_channel(self):
        """When nothing fills a column, the round-trip still works."""
        materials = [
            {"smiles": "X", "volume": 50.0, "reaction_count": 3},
            {"smiles": "Y", "volume": 75.0, "reaction_count": 5},
        ]
        mc_materials, sc_materials = self.analyzer.prioritize_materials(materials)
        self.assertEqual(len(mc_materials), 0)

        mc_groups, cp_wells = self.analyzer.plan_plate_layout(
            mc_materials, sc_materials
        )
        self.assertEqual(len(mc_groups), 0)
        self.assertEqual(len(cp_wells), 8)

        well_infos = [
            WellInfo(index=cp.well_index, smiles=cp.smiles, volume=cp.volume)
            for cp in cp_wells
        ]
        result = self.analyzer.analyze_plate(well_infos)
        # All wells are cherry-picks (8 wells spread across a column
        # won't form a uniform column because they have different SMILES)
        self.assertEqual(result.total_wells, 8)

    def test_realistic_amide_coupling_scenario(self):
        """Simulates a real amide coupling where one amine is used in many
        reactions with the same coupling conditions.

        Scenario: 24 reactions using the same amine + coupling reagent,
        3 reactions with a different amine, 2 one-off reagents.
        """
        materials = [
            # Common amine — used in 24 reactions → 3 full columns
            {"smiles": "NCC1=CC=CC=C1", "volume": 50.0, "reaction_count": 24},
            # Secondary amine — only 3 reactions
            {"smiles": "N(C)CC1=CC=CC=C1", "volume": 50.0, "reaction_count": 3},
            # HATU coupling reagent — used in all 27 reactions
            {"smiles": "O=C(ON1N=NC2=CC=CC=N21)N(C)C", "volume": 75.0, "reaction_count": 27},
            # DIPEA base — used in all 27 reactions
            {"smiles": "CCN(CC)CC", "volume": 100.0, "reaction_count": 27},
        ]
        mc_materials, sc_materials = self.analyzer.prioritize_materials(materials)

        # HATU (3 cols + 3 leftover), DIPEA (3 cols + 3 leftover),
        # amine-1 (3 cols + 0 leftover) = 9 MC cols
        # amine-2 (3 rxns) = SC
        total_mc_cols = sum(m.columns_needed for m in mc_materials)
        self.assertEqual(total_mc_cols, 9)
        self.assertGreaterEqual(len(sc_materials), 1)

        mc_groups, cp_wells = self.analyzer.plan_plate_layout(
            mc_materials, sc_materials
        )
        self.assertEqual(len(mc_groups), 9)
        # Leftover: HATU(3) + DIPEA(3) + amine-2(3) = 9 cherry picks
        self.assertEqual(len(cp_wells), 9)

        # Columns 0-8 are MC, SC starts at index 72 (9 * 8)
        self.assertEqual(cp_wells[0].well_index, 72)


# ===================================================================
# Dataclass defaults / construction
# ===================================================================
class TestDataclasses(TestCase):
    """Sanity checks on dataclass defaults."""

    def test_well_info_defaults(self):
        w = WellInfo(index=0)
        self.assertIsNone(w.smiles)
        self.assertEqual(w.volume, 0.0)
        self.assertIsNone(w.concentration)
        self.assertIsNone(w.solvent)

    def test_multichannel_group_defaults(self):
        g = MultichannelGroup(column_index=0)
        self.assertEqual(g.well_indices, [])
        self.assertEqual(g.smiles, "")
        self.assertEqual(g.volume, 0.0)

    def test_cherry_pick_well_defaults(self):
        cp = CherryPickWell(well_index=5)
        self.assertEqual(cp.smiles, "")
        self.assertEqual(cp.reason, "")

    def test_plate_analysis_result_defaults(self):
        r = PlateAnalysisResult()
        self.assertEqual(r.total_wells, 0)
        self.assertAlmostEqual(r.efficiency, 0.0)

    def test_material_group_defaults(self):
        mg = MaterialGroup(identity_key="test", smiles="CCO", volume=50.0)
        self.assertEqual(mg.reaction_count, 0)
        self.assertEqual(mg.columns_needed, 0)
        self.assertEqual(mg.leftover_wells, 0)
