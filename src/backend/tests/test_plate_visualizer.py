"""Tests for the session visualizer and transfer record modules.

Tests for:
  - TransferRecord     (transfer_mode field, dataclass defaults)
  - TransferLedger     (recording with transfer_mode passthrough)
  - SessionVisualizer  (multichannel annotations on wells, legend, transfer
                        table, transfer summaries, detail panel)

All external dependencies (Django ORM) are mocked.
"""

import re
from unittest import TestCase
from unittest.mock import patch, MagicMock


# ===================================================================
# Transfer Record – transfer_mode field
# ===================================================================

class TestTransferRecordTransferMode(TestCase):
    """Tests for the ``transfer_mode`` field on TransferRecord."""

    def test_default_transfer_mode_is_single(self):
        from backend.opentrons.otwriter.transfer_record import TransferRecord

        rec = TransferRecord(
            action_type="add",
            source_plate_name="sm-plate-1",
            source_plate_role="startingmaterial",
            source_well_index=0,
            source_well_name="A1",
            dest_plate_name="rxn-plate-1",
            dest_plate_role="reaction",
            dest_well_index=0,
            dest_well_name="A1",
            volume=10.0,
        )
        self.assertEqual(rec.transfer_mode, "single")

    def test_multichannel_transfer_mode(self):
        from backend.opentrons.otwriter.transfer_record import TransferRecord

        rec = TransferRecord(
            action_type="add",
            source_plate_name="sm-plate-1",
            source_plate_role="startingmaterial",
            source_well_index=0,
            source_well_name="A1",
            dest_plate_name="rxn-plate-1",
            dest_plate_role="reaction",
            dest_well_index=0,
            dest_well_name="A1",
            volume=10.0,
            transfer_mode="multichannel",
        )
        self.assertEqual(rec.transfer_mode, "multichannel")

    def test_ledger_record_passes_transfer_mode(self):
        from backend.opentrons.otwriter.transfer_record import TransferLedger

        ledger = TransferLedger()
        rec = ledger.record(
            action_type="add",
            source_plate_name="sm-plate-1",
            source_plate_role="startingmaterial",
            source_well_index=0,
            source_well_name="A1",
            dest_plate_name="rxn-plate-1",
            dest_plate_role="reaction",
            dest_well_index=0,
            dest_well_name="A1",
            volume=10.0,
            transfer_mode="multichannel",
        )
        self.assertEqual(rec.transfer_mode, "multichannel")
        self.assertEqual(len(ledger), 1)


# ===================================================================
# Helpers
# ===================================================================

def _make_visualizer(wells, ledger_records=None):
    """Build a SessionVisualizer with mocked dependencies."""
    from backend.opentrons.otwriter.session_visualizer import SessionVisualizer
    from backend.opentrons.otwriter.transfer_record import TransferLedger

    plate = MagicMock()
    plate.id = 10
    plate.name = "sm-plate-1"
    plate.role = "startingmaterial"
    plate.role_index = 1
    plate.labware = "plateone_96_wellplate_2500ul"
    plate.numberwells = 96
    plate.numberwellsincolumn = 8
    plate.numbercolumns = 12
    plate.index = 1

    sg = MagicMock()
    sg.platequeryset = [plate]
    sg.otsessiontype = "reaction"
    sg.protocolname = "test-protocol"

    ledger = TransferLedger()
    if ledger_records:
        for kwargs in ledger_records:
            ledger.record(**kwargs)
    sg.transfer_ledger = ledger

    # Patch Well query inside __init__
    with patch("backend.models.Well") as MockWell:
        MockWell.objects.filter.return_value.order_by.return_value = wells
        viz = SessionVisualizer(sg)

    return viz


# ===================================================================
# SessionVisualizer – multichannel annotations
# ===================================================================

class TestSessionVisualizerMultichannel(TestCase):
    """Tests for MC annotations in SessionVisualizer."""

    def test_mc_well_gets_dashed_border_class(self):
        """Starting material wells with transfer_type='multichannel'
        should have the 'mc-well' CSS class."""
        well = MagicMock()
        well.id = 1
        well.index = 0
        well.smiles = "CCO"
        well.volume = 100.0
        well.solvent = "DMSO"
        well.concentration = 0.5
        well.name = "A1"
        well.transfer_type = "multichannel"

        viz = _make_visualizer(wells=[well])
        html = viz.generate_html()
        self.assertIn("mc-well", html)
        self.assertIn("mc-badge", html)

    def test_single_well_no_mc_class(self):
        """Wells with transfer_type='single' should NOT have mc-well on
        any well div (the CSS definition will still contain the class name)."""
        well = MagicMock()
        well.id = 1
        well.index = 0
        well.smiles = "CCO"
        well.volume = 100.0
        well.solvent = "DMSO"
        well.concentration = 0.5
        well.name = "A1"
        well.transfer_type = "single"

        viz = _make_visualizer(wells=[well])
        html = viz.generate_html()
        # No well div should carry the mc-well class
        mc_well_divs = re.findall(r'<div class="well[^"]*mc-well', html)
        self.assertEqual(mc_well_divs, [])

    def test_mc_legend_entry_present(self):
        """Legend should contain the multichannel entry."""
        viz = _make_visualizer(wells=[])
        html = viz.generate_html()
        self.assertIn("mc-legend-swatch", html)
        self.assertIn("Multichannel transfer", html)

    def test_mc_transfer_mode_in_transfer_table(self):
        """Transfer table should show MC tag for multichannel transfers."""
        records = [
            dict(
                action_type="add",
                source_plate_name="sm-plate-1",
                source_plate_role="startingmaterial",
                source_well_index=0,
                source_well_name="A1",
                dest_plate_name="rxn-plate-1",
                dest_plate_role="reaction",
                dest_well_index=0,
                dest_well_name="A1",
                volume=10.0,
                smiles="CCO",
                transfer_mode="multichannel",
            ),
        ]
        viz = _make_visualizer(wells=[], ledger_records=records)
        html = viz.generate_html()
        self.assertIn("mc-tag", html)
        self.assertIn("Mode", html)

    def test_single_transfer_shows_sc(self):
        """Transfer table should show 'SC' for single-channel transfers."""
        records = [
            dict(
                action_type="add",
                source_plate_name="sm-plate-1",
                source_plate_role="startingmaterial",
                source_well_index=0,
                source_well_name="A1",
                dest_plate_name="rxn-plate-1",
                dest_plate_role="reaction",
                dest_well_index=0,
                dest_well_name="A1",
                volume=10.0,
                smiles="CCO",
            ),
        ]
        viz = _make_visualizer(wells=[], ledger_records=records)
        html = viz.generate_html()
        # Should contain SC (for single channel) but not the MC tag
        self.assertIn(">SC<", html)

    def test_reaction_well_annotated_via_incoming_mc_transfer(self):
        """Reaction plate wells receiving MC transfers should get annotated
        even if their model transfer_type is 'single'."""
        # Reaction plate well — transfer_type is 'single'
        well = MagicMock()
        well.id = 2
        well.index = 0
        well.smiles = "CCO"
        well.volume = 100.0
        well.solvent = "DMSO"
        well.concentration = 0.5
        well.name = "A1"
        well.transfer_type = "single"

        # MC transfer arriving at this well
        records = [
            dict(
                action_type="add",
                source_plate_name="sm-plate-1",
                source_plate_role="startingmaterial",
                source_well_index=0,
                source_well_name="A1",
                dest_plate_name="sm-plate-1",  # same plate as viz
                dest_plate_role="reaction",
                dest_well_index=0,
                dest_well_name="A1",
                volume=10.0,
                smiles="CCO",
                transfer_mode="multichannel",
            ),
        ]
        viz = _make_visualizer(wells=[well], ledger_records=records)
        html = viz.generate_html()
        # Well should be annotated as MC via incoming transfer
        self.assertIn("mc-well", html)
        self.assertIn("mc-badge", html)

    def test_detail_panel_shows_multichannel_mode(self):
        """The JS detail panel should display 'Multichannel' label for
        MC wells via data-transfer-type attribute."""
        well = MagicMock()
        well.id = 1
        well.index = 0
        well.smiles = "CCO"
        well.volume = 100.0
        well.solvent = "DMSO"
        well.concentration = 0.5
        well.name = "A1"
        well.transfer_type = "multichannel"

        viz = _make_visualizer(wells=[well])
        html = viz.generate_html()
        self.assertIn('data-transfer-type="multichannel"', html)


# ===================================================================
# Transfer summaries – MC tags
# ===================================================================

class TestTransferSummariesMultichannel(TestCase):
    """Tests for MC tags in transfer summary list items."""

    def test_mc_tag_in_summary(self):
        from backend.opentrons.otwriter.session_visualizer import SessionVisualizer
        from backend.opentrons.otwriter.transfer_record import TransferRecord

        rec = TransferRecord(
            action_type="add",
            source_plate_name="sm-plate-1",
            source_plate_role="startingmaterial",
            source_well_index=0,
            source_well_name="A1",
            dest_plate_name="rxn-plate-1",
            dest_plate_role="reaction",
            dest_well_index=0,
            dest_well_name="A1",
            volume=10.0,
            smiles="CCO",
            transfer_mode="multichannel",
        )
        result = SessionVisualizer._transfer_summaries([rec], "in")
        self.assertIn("mc-tag", result)
        self.assertIn("MC", result)

    def test_no_mc_tag_for_single(self):
        from backend.opentrons.otwriter.session_visualizer import SessionVisualizer
        from backend.opentrons.otwriter.transfer_record import TransferRecord

        rec = TransferRecord(
            action_type="add",
            source_plate_name="sm-plate-1",
            source_plate_role="startingmaterial",
            source_well_index=0,
            source_well_name="A1",
            dest_plate_name="rxn-plate-1",
            dest_plate_role="reaction",
            dest_well_index=0,
            dest_well_name="A1",
            volume=10.0,
            smiles="CCO",
        )
        result = SessionVisualizer._transfer_summaries([rec], "in")
        self.assertNotIn("mc-tag", result)
