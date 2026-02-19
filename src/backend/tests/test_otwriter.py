"""Tests for the otwriter package.

Tests for:
  - CommandGenerator   (protocol command string generation)
  - FileManager        (file path creation, content writing, OTScript model)
  - VolumeManager      (volume checks, dead volume, aspirate height, reactant status)
  - WellFinder         (reaction well, solvent wells, starting material wells)
  - SessionHandler     (base handler: add_command, get_session_number)
  - QueryService       (action/mix/extract querysets, plates, pipettes, columns)
  - ScriptGenerator    (orchestration, setup, session routing)
  - ReactionSessionHandler  (process session, add/mix actions, dilution)
  - WorkupSessionHandler    (process session, extract/mix actions)
  - AnalysisSessionHandler  (process session, add actions)

All external dependencies (Django ORM, file I/O, storage) are mocked.
"""

import math
import os
from unittest import TestCase
from unittest.mock import patch, MagicMock, PropertyMock, call, mock_open

# ---------------------------------------------------------------------------
# Module paths used for patching
# ---------------------------------------------------------------------------
_WRITER_MOD = "backend.opentrons.otwriter"
_SG_MOD = f"{_WRITER_MOD}.script_generator"
_CG_MOD = f"{_WRITER_MOD}.command_generator"
_FM_MOD = f"{_WRITER_MOD}.file_manager"
_BH_MOD = f"{_WRITER_MOD}.session_handlers.base_handler"
_RH_MOD = f"{_WRITER_MOD}.session_handlers.reaction_handler"
_WH_MOD = f"{_WRITER_MOD}.session_handlers.workup_handler"
_AH_MOD = f"{_WRITER_MOD}.session_handlers.analysis_handler"
_VM_MOD = f"{_WRITER_MOD}.utils.volume_manager"
_WF_MOD = f"{_WRITER_MOD}.utils.well_finder"
_QS_MOD = f"{_WRITER_MOD}.utils.query_service"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_script_generator_stub(**overrides):
    """Return a MagicMock with the interface expected by helper components."""
    sg = MagicMock()
    sg.otsession_id = overrides.get("otsession_id", 1)
    sg.actionsession_ids = overrides.get("actionsession_ids", [1, 2])
    sg.otsessiontype = overrides.get("otsessiontype", "reaction")
    sg.reactionstep = overrides.get("reactionstep", 1)
    sg.batchtag = overrides.get("batchtag", "ABC123")
    sg.apiLevel = overrides.get("apiLevel", "2.9")
    sg.protocolname = overrides.get(
        "protocolname", "reaction-bABC123-r1-s1"
    )
    sg.pipettename = overrides.get("pipettename", "p300_single")
    sg.content = []
    sg.otsessionobj = MagicMock()
    sg.otsessionobj.id = sg.otsession_id

    # Sub-components
    sg.command_generator = MagicMock()
    sg.query_service = MagicMock()
    sg.volume_manager = MagicMock()
    sg.well_finder = MagicMock()
    sg.file_manager = MagicMock()

    def _add_command(cmds):
        if isinstance(cmds, str):
            sg.content.append(cmds)
        elif isinstance(cmds, list):
            sg.content.extend(cmds)

    sg.add_command.side_effect = _add_command
    return sg


def _make_well_stub(**overrides):
    """Return a MagicMock with the interface expected for a Well object."""
    w = MagicMock()
    w.id = overrides.get("id", 1)
    w.index = overrides.get("index", 0)
    w.volume = overrides.get("volume", 100.0)
    w.available = overrides.get("available", True)
    w.smiles = overrides.get("smiles", "CCO")
    w.solvent = overrides.get("solvent", "DMSO")
    w.concentration = overrides.get("concentration", 0.5)
    w.role = overrides.get("role", "startingmaterial")
    w.reactantfornextstep = overrides.get("reactantfornextstep", False)
    plate = MagicMock()
    plate.id = overrides.get("plate_id", 10)
    plate.name = overrides.get("plate_name", "plate_1")
    plate.maxwellvolume = overrides.get("maxwellvolume", 2500.0)
    w.plate_id = plate
    return w


def _make_plate_stub(**overrides):
    """Return a MagicMock with the interface expected for a Plate object."""
    p = MagicMock()
    p.id = overrides.get("id", 10)
    p.name = overrides.get("name", "plate_1")
    p.role = overrides.get("role", "reaction")
    p.role_index = overrides.get("role_index", 1)
    p.labware = overrides.get("labware", "plateone_96_wellplate_2500ul")
    p.maxwellvolume = overrides.get("maxwellvolume", 2500.0)
    return p


def _make_action_session_qs(session_number=1, action_sessions=None):
    """Return a MagicMock queryset for action sessions in handler tests."""
    if action_sessions is None:
        as1 = MagicMock()
        as1.reaction_id.id = 1
        action_sessions = [as1]
    qs = MagicMock()
    qs.exists.return_value = True
    qs.count.return_value = len(action_sessions)
    qs.__iter__ = lambda s: iter(action_sessions)
    distinct_result = MagicMock()
    distinct_result.__len__ = lambda s: 1
    distinct_result.__getitem__ = lambda s, i: session_number
    vl = MagicMock()
    vl.distinct.return_value = distinct_result
    qs.values_list.return_value = vl
    return qs


# ===================================================================
# CommandGenerator tests
# ===================================================================
class TestCommandGeneratorComment(TestCase):
    """Tests for CommandGenerator.comment()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.sg = _make_script_generator_stub()
        self.cg = CommandGenerator(self.sg)

    def test_default_newline_and_level(self):
        result = self.cg.comment("hello world")
        self.assertEqual(result, "\n\t# hello world")

    def test_no_newlines(self):
        result = self.cg.comment("test", num_newlines=0)
        self.assertEqual(result, "\t# test")

    def test_double_newline(self):
        result = self.cg.comment("test", num_newlines=2)
        self.assertEqual(result, "\n\n\t# test")

    def test_deeper_indent(self):
        result = self.cg.comment("test", level=3)
        self.assertEqual(result, "\n\t\t\t# test")


class TestCommandGeneratorIndent(TestCase):
    """Tests for CommandGenerator.indent()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_default_level(self):
        self.assertEqual(self.cg.indent(), "\t")

    def test_level_zero(self):
        self.assertEqual(self.cg.indent(0), "")

    def test_level_four(self):
        self.assertEqual(self.cg.indent(4), "\t\t\t\t")


class TestCommandGeneratorScriptSetup(TestCase):
    """Tests for CommandGenerator.get_script_setup()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_returns_list(self):
        result = self.cg.get_script_setup("myprotocol", "2.9")
        self.assertIsInstance(result, list)

    def test_contains_import(self):
        result = self.cg.get_script_setup("myprotocol", "2.9")
        self.assertIn("from opentrons import protocol_api", result)

    def test_contains_protocol_name(self):
        result = self.cg.get_script_setup("myprotocol", "2.9")
        joined = "\n".join(result)
        self.assertIn("myprotocol", joined)

    def test_contains_api_level(self):
        result = self.cg.get_script_setup("myprotocol", "2.9")
        joined = "\n".join(result)
        self.assertIn("2.9", joined)

    def test_contains_run_function(self):
        result = self.cg.get_script_setup("myprotocol", "2.9")
        self.assertTrue(
            any("def run(protocol" in line for line in result),
            "Expected 'def run(protocol...' in setup lines",
        )


class TestCommandGeneratorLabwareSetup(TestCase):
    """Tests for CommandGenerator.get_labware_setup()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_single_plate(self):
        plates = [{"name": "rxn_plate", "labware": "corning_96", "index": 4}]
        result = self.cg.get_labware_setup(plates)
        joined = "\n".join(result)
        self.assertIn("rxn_plate", joined)
        self.assertIn("load_labware", joined)
        self.assertIn("corning_96", joined)

    def test_multiple_plates(self):
        plates = [
            {"name": "p1", "labware": "lab1", "index": 1},
            {"name": "p2", "labware": "lab2", "index": 2},
        ]
        result = self.cg.get_labware_setup(plates)
        joined = "\n".join(result)
        self.assertIn("p1", joined)
        self.assertIn("p2", joined)

    def test_empty_plates(self):
        result = self.cg.get_labware_setup([])
        # Should still have the comment header
        self.assertTrue(len(result) >= 1)


class TestCommandGeneratorTiprackSetup(TestCase):
    """Tests for CommandGenerator.get_tiprack_setup()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_single_tiprack(self):
        tipracks = [{"name": "tr1", "labware": "opentrons_96_tiprack_300ul", "index": 1}]
        result = self.cg.get_tiprack_setup(tipracks)
        joined = "\n".join(result)
        self.assertIn("tr1", joined)
        self.assertIn("load_labware", joined)

    def test_empty_tipracks(self):
        result = self.cg.get_tiprack_setup([])
        self.assertTrue(len(result) >= 1)  # comment header


class TestCommandGeneratorPipetteSetup(TestCase):
    """Tests for CommandGenerator.get_pipette_setup()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_pipette_setup_content(self):
        result = self.cg.get_pipette_setup(
            pipette_name="p300_single",
            pipette_labware="p300_single_gen2",
            pipette_position="left",
            tiprack_names=["tr1", "tr2"],
        )
        joined = "\n".join(result)
        self.assertIn("load_instrument", joined)
        self.assertIn("p300_single_gen2", joined)
        self.assertIn("left", joined)
        self.assertIn("tr1", joined)
        self.assertIn("tr2", joined)


class TestCommandGeneratorTipTracking(TestCase):
    """Tests for get_number_tips_available_setup, pickup/drop tip functions."""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_tips_available_single_channel(self):
        result = self.cg.get_number_tips_available_setup(2, "single")
        joined = "\n".join(result)
        self.assertIn("192", joined)  # 2 * 96
        self.assertIn("single", joined)

    def test_tips_available_multi_channel(self):
        result = self.cg.get_number_tips_available_setup(1, "multi")
        joined = "\n".join(result)
        self.assertIn("96", joined)
        self.assertIn("multi", joined)

    def test_pickup_tip_function(self):
        result = self.cg.get_pickup_tip_function("p300_single")
        joined = "\n".join(result)
        self.assertIn("def pickUpTip()", joined)
        self.assertIn("pick_up_tip", joined)
        self.assertIn("p300_single", joined)

    def test_drop_tip_function(self):
        result = self.cg.get_drop_tip_function("p300_single")
        joined = "\n".join(result)
        self.assertIn("def dropTip()", joined)
        self.assertIn("drop_tip", joined)

    def test_pick_up_tip_call(self):
        result = self.cg.pick_up_tip()
        joined = "\n".join(result)
        self.assertIn("pickUpTip()", joined)

    def test_drop_tip_call(self):
        result = self.cg.drop_tip()
        joined = "\n".join(result)
        self.assertIn("dropTip()", joined)


class TestCommandGeneratorProtocolControl(TestCase):
    """Tests for pause_protocol, delay_protocol, speed setting."""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_pause_protocol(self):
        result = self.cg.pause_protocol("Wait for user")
        joined = "\n".join(result)
        self.assertIn("protocol.pause", joined)
        self.assertIn("Wait for user", joined)

    def test_delay_protocol(self):
        result = self.cg.delay_protocol(30)
        joined = "\n".join(result)
        self.assertIn("protocol.delay", joined)
        self.assertIn("30", joined)

    def test_set_aspirate_speed(self):
        result = self.cg.set_aspirate_speed(50)
        joined = "\n".join(result)
        self.assertIn("flow_rate.aspirate", joined)
        self.assertIn("50", joined)

    def test_set_dispense_speed(self):
        result = self.cg.set_dispense_speed(100)
        joined = "\n".join(result)
        self.assertIn("flow_rate.dispense", joined)
        self.assertIn("100", joined)


class TestCommandGeneratorTransferSingle(TestCase):
    """Tests for CommandGenerator.transfer_fluid_single()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_basic_transfer(self):
        result = self.cg.transfer_fluid_single(
            aspirateplatename="plate_a",
            dispenseplatename="plate_b",
            aspiratewellindex=0,
            dispensewellindex=5,
            transvolume=50.0,
        )
        joined = "\n".join(result)
        self.assertIn("plate_a", joined)
        self.assertIn("plate_b", joined)
        self.assertIn("50.0", joined)
        self.assertIn(".transfer(", joined)
        self.assertIn("wells()", joined)

    def test_transfer_with_custom_height(self):
        result = self.cg.transfer_fluid_single(
            aspirateplatename="p1",
            dispenseplatename="p2",
            aspiratewellindex=0,
            dispensewellindex=0,
            transvolume=25.0,
            aspirateheight=5.0,
        )
        joined = "\n".join(result)
        self.assertIn(".bottom(5.0)", joined)

    def test_transfer_type_in_comment(self):
        result = self.cg.transfer_fluid_single(
            aspirateplatename="p1",
            dispenseplatename="p2",
            aspiratewellindex=0,
            dispensewellindex=0,
            transvolume=10.0,
            transfertype="dilution",
        )
        joined = "\n".join(result)
        self.assertIn("dilution", joined)

    def test_returns_two_lines(self):
        result = self.cg.transfer_fluid_single(
            aspirateplatename="p1",
            dispenseplatename="p2",
            aspiratewellindex=0,
            dispensewellindex=0,
            transvolume=10.0,
        )
        self.assertEqual(len(result), 2)  # comment + command


class TestCommandGeneratorTransferMulti(TestCase):
    """Tests for CommandGenerator.transfer_fluid_multi()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_basic_multi_transfer(self):
        result = self.cg.transfer_fluid_multi(
            aspirateplatename="plate_a",
            dispenseplatename="plate_b",
            aspiratecolumnindex=0,
            dispensecolumnindex=3,
            transvolume=100.0,
        )
        joined = "\n".join(result)
        self.assertIn("plate_a", joined)
        self.assertIn("plate_b", joined)
        self.assertIn("100.0", joined)
        self.assertIn("columns()", joined)

    def test_multi_returns_two_lines(self):
        result = self.cg.transfer_fluid_multi(
            aspirateplatename="p1",
            dispenseplatename="p2",
            aspiratecolumnindex=0,
            dispensecolumnindex=0,
            transvolume=10.0,
        )
        self.assertEqual(len(result), 2)


class TestCommandGeneratorMixWell(TestCase):
    """Tests for CommandGenerator.mix_well()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_mix_well(self):
        result = self.cg.mix_well(wellindex=3, nomixes=5, plate="rxn_plate")
        joined = "\n".join(result)
        self.assertIn("rxn_plate", joined)
        self.assertIn(".mix(5", joined)
        self.assertIn("wells()", joined)

    def test_mix_well_returns_two_lines(self):
        result = self.cg.mix_well(wellindex=0, nomixes=3, plate="p1")
        self.assertEqual(len(result), 2)


class TestCommandGeneratorMixColumn(TestCase):
    """Tests for CommandGenerator.mix_column()"""

    def setUp(self):
        from backend.opentrons.otwriter.command_generator import CommandGenerator
        self.cg = CommandGenerator(_make_script_generator_stub())

    def test_mix_column(self):
        result = self.cg.mix_column(columnindex=2, nomixes=4, plate="rxn_plate")
        joined = "\n".join(result)
        self.assertIn("rxn_plate", joined)
        self.assertIn(".mix(4", joined)
        self.assertIn("columns()", joined)

    def test_mix_column_returns_two_lines(self):
        result = self.cg.mix_column(columnindex=0, nomixes=3, plate="p1")
        self.assertEqual(len(result), 2)


# ===================================================================
# VolumeManager tests
# ===================================================================
class TestVolumeManagerCheckVolumesClose(TestCase):
    """Tests for VolumeManager.check_volumes_close()"""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.vm = VolumeManager(_make_script_generator_stub())

    def test_identical_volumes(self):
        self.assertTrue(self.vm.check_volumes_close(100.0, 100.0))

    def test_close_volumes(self):
        self.assertTrue(self.vm.check_volumes_close(100.0, 100.05))

    def test_different_volumes(self):
        self.assertFalse(self.vm.check_volumes_close(100.0, 200.0))

    def test_both_zero(self):
        self.assertTrue(self.vm.check_volumes_close(0.0, 0.0))

    def test_near_zero(self):
        # math.isclose with rel_tol treats near-zero specially
        result = self.vm.check_volumes_close(0.0001, 0.0)
        # With rel_tol=0.001, this depends on abs_tol (default 0)
        self.assertIsInstance(result, bool)


class TestVolumeManagerDeadVolume(TestCase):
    """Tests for VolumeManager.get_dead_volume()"""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.vm = VolumeManager(_make_script_generator_stub())

    def test_standard_plate(self):
        result = self.vm.get_dead_volume(2500.0)
        self.assertAlmostEqual(result, 125.0)

    def test_small_plate(self):
        result = self.vm.get_dead_volume(100.0)
        self.assertAlmostEqual(result, 5.0)

    def test_zero_volume(self):
        result = self.vm.get_dead_volume(0.0)
        self.assertAlmostEqual(result, 0.0)


class TestVolumeManagerOptimalTransfer(TestCase):
    """Tests for VolumeManager.get_optimal_transfer_volume()"""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.vm = VolumeManager(_make_script_generator_stub())

    def test_sufficient_volume(self):
        result = self.vm.get_optimal_transfer_volume(200.0, 100.0)
        self.assertEqual(result, 100.0)

    def test_exact_volume(self):
        result = self.vm.get_optimal_transfer_volume(100.0, 100.0)
        self.assertEqual(result, 100.0)

    def test_insufficient_volume(self):
        result = self.vm.get_optimal_transfer_volume(50.0, 100.0)
        self.assertEqual(result, 50.0)

    def test_zero_available(self):
        result = self.vm.get_optimal_transfer_volume(0.0, 100.0)
        self.assertEqual(result, 0.0)

    def test_negative_available(self):
        result = self.vm.get_optimal_transfer_volume(-10.0, 100.0)
        self.assertEqual(result, 0.0)


class TestVolumeManagerAspirateHeight(TestCase):
    """Tests for VolumeManager.calculate_aspirate_height()"""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.vm = VolumeManager(_make_script_generator_stub())

    def test_known_labware(self):
        # plateone_96_wellplate_2500ul: m=0.01806, c=0.6015
        result = self.vm.calculate_aspirate_height(
            bottomlayervolume=100.0, labware="plateone_96_wellplate_2500ul"
        )
        raw = (0.01806 * 100.0) + 0.6015
        expected = raw * 1.15
        self.assertAlmostEqual(result, expected, places=4)

    def test_paradox_labware(self):
        # paradox_96_wellplate_961ul: m=0.02865, c=-0.4040
        result = self.vm.calculate_aspirate_height(
            bottomlayervolume=50.0, labware="paradox_96_wellplate_961ul"
        )
        raw = (0.02865 * 50.0) + (-0.4040)
        expected = raw * 1.15
        self.assertAlmostEqual(result, expected, places=4)

    def test_unknown_labware_fallback(self):
        result = self.vm.calculate_aspirate_height(
            bottomlayervolume=200.0, labware="unknown_plate"
        )
        expected = (200.0 / 200) + 0.5
        self.assertAlmostEqual(result, expected, places=4)

    def test_zero_volume(self):
        result = self.vm.calculate_aspirate_height(
            bottomlayervolume=0.0, labware="plateone_96_wellplate_2500ul"
        )
        raw = 0.6015
        expected = raw * 1.15
        self.assertAlmostEqual(result, expected, places=4)


class TestVolumeManagerGetWellVolumeAvailable(TestCase):
    """Tests for VolumeManager.get_well_volume_available()"""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.sg = _make_script_generator_stub()
        self.vm = VolumeManager(self.sg)

    def test_returns_volume_minus_dead(self):
        plate = _make_plate_stub(maxwellvolume=2500.0)
        self.sg.query_service.get_plate_by_id.return_value = plate

        well = _make_well_stub(volume=200.0, plate_id=10)
        result = self.vm.get_well_volume_available(well)
        # dead = 2500 * 0.05 = 125
        self.assertAlmostEqual(result, 200.0 - 125.0)

    def test_marks_unavailable_when_negative(self):
        plate = _make_plate_stub(maxwellvolume=2500.0)
        self.sg.query_service.get_plate_by_id.return_value = plate

        well = _make_well_stub(volume=50.0, plate_id=10)
        result = self.vm.get_well_volume_available(well)
        # dead = 125, available = 50 - 125 = -75
        self.assertLess(result, 0)
        well.save.assert_called()
        self.assertFalse(well.available)


class TestVolumeManagerUpdateWellVolume(TestCase):
    """Tests for VolumeManager.update_well_volume()"""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.vm = VolumeManager(_make_script_generator_stub())

    def test_subtracts_volume(self):
        well = _make_well_stub(volume=100.0)
        self.vm.update_well_volume(well, 30.0)
        self.assertAlmostEqual(well.volume, 70.0)
        well.save.assert_called_once()

    def test_full_depletion(self):
        well = _make_well_stub(volume=50.0)
        self.vm.update_well_volume(well, 50.0)
        self.assertAlmostEqual(well.volume, 0.0)


class TestVolumeManagerReactantStatus(TestCase):
    """Tests for update_well_reactant_status and update_column_reactant_status."""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.vm = VolumeManager(_make_script_generator_stub())

    def test_update_well_reactant_true(self):
        well = _make_well_stub()
        self.vm.update_well_reactant_status(well, True)
        self.assertTrue(well.reactantfornextstep)
        well.save.assert_called_once()

    def test_update_well_reactant_false(self):
        well = _make_well_stub()
        self.vm.update_well_reactant_status(well, False)
        self.assertFalse(well.reactantfornextstep)
        well.save.assert_called_once()

    @patch("backend.models.Well")
    def test_update_column_reactant_status(self, mock_Well):
        w1 = _make_well_stub()
        w2 = _make_well_stub()
        well_qs = MagicMock()
        well_qs.count.return_value = 2
        well_qs.__iter__ = lambda s: iter([w1, w2])
        well_qs.__bool__ = lambda s: True
        mock_Well.objects.filter.return_value = well_qs

        col = MagicMock()
        col.index = 0
        self.vm.update_column_reactant_status(col, True)
        self.assertTrue(w1.reactantfornextstep)
        self.assertTrue(w2.reactantfornextstep)
        w1.save.assert_called()
        w2.save.assert_called()


class TestVolumeManagerGetMaxWellVolume(TestCase):
    """Tests for VolumeManager.get_max_well_volume()"""

    def setUp(self):
        from backend.opentrons.otwriter.utils.volume_manager import VolumeManager
        self.sg = _make_script_generator_stub()
        self.vm = VolumeManager(self.sg)

    def test_returns_plate_max_volume(self):
        plate = _make_plate_stub(maxwellvolume=1500.0)
        self.sg.query_service.get_plate_by_id.return_value = plate
        result = self.vm.get_max_well_volume(10)
        self.assertEqual(result, 1500.0)


# ===================================================================
# FileManager tests
# ===================================================================
class TestFileManagerCreateFilePath(TestCase):
    """Tests for FileManager.create_file_path()"""

    @patch(f"{_FM_MOD}.settings")
    def test_file_path_format(self, mock_settings):
        mock_settings.MEDIA_ROOT = "/media"
        sg = _make_script_generator_stub(protocolname="reaction-bXYZ-r1-s1")
        from backend.opentrons.otwriter.file_manager import FileManager
        fm = FileManager(sg)
        filepath, filename = fm.create_file_path()
        self.assertEqual(filename, "reaction-bXYZ-r1-s1.txt")
        self.assertIn("/media", filepath)
        self.assertIn("tmp/", filepath)


class TestFileManagerWriteContent(TestCase):
    """Tests for FileManager.write_content()"""

    @patch(f"{_FM_MOD}.settings")
    def test_writes_lines_to_file(self, mock_settings):
        mock_settings.MEDIA_ROOT = "/media"
        sg = _make_script_generator_stub()
        from backend.opentrons.otwriter.file_manager import FileManager
        fm = FileManager(sg)

        m = mock_open()
        with patch("builtins.open", m):
            fm.write_content(["line1", "line2"])

        m.assert_called_once_with(fm.filepath, "w")
        handle = m()
        calls = handle.write.call_args_list
        self.assertEqual(len(calls), 2)
        self.assertIn("line1", calls[0][0][0])
        self.assertIn("line2", calls[1][0][0])


class TestFileManagerCreateOTScriptModel(TestCase):
    """Tests for FileManager.create_ot_script_model()"""

    @patch(f"{_FM_MOD}.default_storage")
    @patch(f"{_FM_MOD}.OTScript")
    @patch(f"{_FM_MOD}.settings")
    def test_creates_model(self, mock_settings, mock_OTScript, mock_storage):
        mock_settings.MEDIA_ROOT = "/media"
        sg = _make_script_generator_stub()
        from backend.opentrons.otwriter.file_manager import FileManager

        fm = FileManager(sg)
        otscript_instance = MagicMock()
        mock_OTScript.return_value = otscript_instance
        mock_storage.save.return_value = "otscripts/test.py"

        m = mock_open(read_data=b"content")
        with patch("builtins.open", m):
            result = fm.create_ot_script_model()

        otscript_instance.save.assert_called_once()
        self.assertEqual(result, otscript_instance)


# ===================================================================
# SessionHandler (base) tests
# ===================================================================
class TestSessionHandlerAddCommand(TestCase):
    """Tests for SessionHandler.add_command()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.base_handler import SessionHandler
        self.sg = _make_script_generator_stub()
        self.handler = SessionHandler(self.sg)

    def test_add_string_command(self):
        self.handler.add_command("test command")
        self.sg.add_command.assert_called_once_with("test command")

    def test_add_list_command(self):
        self.handler.add_command(["cmd1", "cmd2"])
        self.sg.add_command.assert_called_once_with(["cmd1", "cmd2"])


class TestSessionHandlerProcessSession(TestCase):
    """Tests for SessionHandler.process_session() abstract method."""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.base_handler import SessionHandler
        self.sg = _make_script_generator_stub()
        self.handler = SessionHandler(self.sg)

    def test_raises_not_implemented(self):
        with self.assertRaises(NotImplementedError):
            self.handler.process_session(MagicMock())


class TestSessionHandlerGetSessionNumber(TestCase):
    """Tests for SessionHandler.get_session_number()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.base_handler import SessionHandler
        self.sg = _make_script_generator_stub()
        self.handler = SessionHandler(self.sg)

    def test_returns_session_number(self):
        qs = MagicMock()
        qs.exists.return_value = True
        distinct_result = MagicMock()
        distinct_result.__len__ = lambda s: 1
        distinct_result.__getitem__ = lambda s, i: 3
        vl = MagicMock()
        vl.distinct.return_value = distinct_result
        qs.values_list.return_value = vl
        result = self.handler.get_session_number(qs)
        self.assertEqual(result, 3)

    def test_empty_queryset_returns_zero(self):
        qs = MagicMock()
        qs.exists.return_value = False
        result = self.handler.get_session_number(qs)
        self.assertEqual(result, 0)

    def test_multiple_session_numbers_uses_first(self):
        qs = MagicMock()
        qs.exists.return_value = True
        distinct_result = MagicMock()
        distinct_result.__len__ = lambda s: 2
        distinct_result.__getitem__ = lambda s, i: [1, 2][i]
        vl = MagicMock()
        vl.distinct.return_value = distinct_result
        qs.values_list.return_value = vl
        result = self.handler.get_session_number(qs)
        self.assertEqual(result, 1)


# ===================================================================
# WellFinder tests
# ===================================================================
class TestWellFinderFindReactionWell(TestCase):
    """Tests for WellFinder.find_reaction_well()"""

    @patch(f"{_WF_MOD}.get_product_smiles")
    @patch(f"{_WF_MOD}.Well")
    def test_finds_well_by_reaction_id(self, mock_Well, mock_get_smiles):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder

        mock_get_smiles.return_value = ["CCC"]
        expected_well = _make_well_stub()
        mock_Well.objects.get.return_value = expected_well

        sg = _make_script_generator_stub()
        wf = WellFinder(sg)
        result = wf.find_reaction_well(reaction_id=1, role="reaction", role_index=1)
        self.assertEqual(result, expected_well)
        mock_Well.objects.get.assert_called_once()

    @patch(f"{_WF_MOD}.get_product_smiles")
    @patch(f"{_WF_MOD}.Well")
    def test_raises_on_not_found(self, mock_Well, mock_get_smiles):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder
        from backend.models import Well as RealWell

        mock_get_smiles.return_value = ["CCC"]
        mock_Well.objects.get.side_effect = RealWell.DoesNotExist("not found")

        sg = _make_script_generator_stub()
        wf = WellFinder(sg)
        with self.assertRaises(Exception):
            wf.find_reaction_well(reaction_id=999, role="reaction")


class TestWellFinderFindSolventWells(TestCase):
    """Tests for WellFinder.find_solvent_wells()"""

    @patch(f"{_WF_MOD}.Plate")
    def test_finds_solvent_wells_single(self, mock_Plate):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder

        sg = _make_script_generator_stub()
        sg.volume_manager.check_volumes_close.side_effect = (
            lambda volume1, volume2: math.isclose(volume1, volume2, rel_tol=0.001)
        )
        sg.volume_manager.get_well_volume_available.return_value = 500.0

        well = _make_well_stub(volume=500.0, role="solvent", solvent="DMSO")
        plate = _make_plate_stub(role="solvent")
        well_qs = MagicMock()
        well_qs.exists.return_value = True
        well_qs.count.return_value = 1
        well_qs.__iter__ = lambda s: iter([well])
        plate.well_set.all.return_value.filter.return_value = well_qs

        plate_qs = MagicMock()
        plate_qs.exists.return_value = True
        plate_qs.count.return_value = 1
        plate_qs.__iter__ = lambda s: iter([plate])
        mock_Plate.objects.filter.return_value = plate_qs

        wf = WellFinder(sg)
        result = wf.find_solvent_wells("DMSO", 100.0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][0], well)
        self.assertAlmostEqual(result[0][1], 100.0)

    @patch(f"{_WF_MOD}.Plate")
    def test_no_solvent_plates(self, mock_Plate):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder

        sg = _make_script_generator_stub()
        plate_qs = MagicMock()
        plate_qs.exists.return_value = False
        mock_Plate.objects.filter.return_value = plate_qs

        wf = WellFinder(sg)
        result = wf.find_solvent_wells("DMSO", 100.0)
        self.assertEqual(len(result), 0)

    @patch(f"{_WF_MOD}.Plate")
    def test_partial_fill_across_wells(self, mock_Plate):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder

        sg = _make_script_generator_stub()
        sg.volume_manager.check_volumes_close.side_effect = (
            lambda volume1, volume2: math.isclose(volume1, volume2, rel_tol=0.001)
        )
        # First well has 60µL available, second has 500µL
        sg.volume_manager.get_well_volume_available.side_effect = [60.0, 500.0]

        w1 = _make_well_stub(index=0, volume=60.0)
        w2 = _make_well_stub(index=1, volume=500.0)

        plate = _make_plate_stub(role="solvent")
        well_qs = MagicMock()
        well_qs.exists.return_value = True
        well_qs.count.return_value = 2
        well_qs.__iter__ = lambda s: iter([w1, w2])
        plate.well_set.all.return_value.filter.return_value = well_qs

        plate_qs = MagicMock()
        plate_qs.exists.return_value = True
        plate_qs.count.return_value = 1
        plate_qs.__iter__ = lambda s: iter([plate])
        mock_Plate.objects.filter.return_value = plate_qs

        wf = WellFinder(sg)
        result = wf.find_solvent_wells("DMSO", 100.0)
        self.assertEqual(len(result), 2)
        # First well contributes 60, second contributes 40
        self.assertAlmostEqual(result[0][1], 60.0)
        self.assertAlmostEqual(result[1][1], 40.0)


class TestWellFinderFindStartingMaterialWells(TestCase):
    """Tests for WellFinder.find_starting_material_wells()"""

    @patch(f"{_WF_MOD}.get_previous_reaction_query_sets")
    @patch(f"{_WF_MOD}.Well")
    def test_step_1_with_solvent(self, mock_Well, mock_prev_rxn):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder

        sg = _make_script_generator_stub()
        sg.volume_manager.check_volumes_close.side_effect = (
            lambda volume1, volume2: math.isclose(volume1, volume2, rel_tol=0.001)
        )
        sg.volume_manager.get_well_volume_available.return_value = 200.0

        mock_prev_rxn.return_value = None

        well = _make_well_stub(volume=200.0, smiles="CCO", solvent="DMSO", concentration=0.5)
        well_qs = MagicMock()
        well_qs.exists.return_value = True
        well_qs.count.return_value = 1
        well_qs.__iter__ = lambda s: iter([well])
        mock_Well.objects.filter.return_value = well_qs
        mock_Well.objects.filter.return_value.order_by.return_value = well_qs

        wf = WellFinder(sg)
        result = wf.find_starting_material_wells(
            reaction_step_no=1,
            reaction_id=1,
            smiles="CCO",
            solvent="DMSO",
            concentration=0.5,
            transfer_volume=50.0,
        )
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 1)

    @patch(f"{_WF_MOD}.get_previous_reaction_query_sets")
    @patch(f"{_WF_MOD}.Well")
    def test_step_1_without_solvent(self, mock_Well, mock_prev_rxn):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder

        sg = _make_script_generator_stub()
        sg.volume_manager.check_volumes_close.side_effect = (
            lambda volume1, volume2: math.isclose(volume1, volume2, rel_tol=0.001)
        )
        sg.volume_manager.get_well_volume_available.return_value = 200.0

        mock_prev_rxn.return_value = None

        well = _make_well_stub(volume=200.0, smiles="CCO")
        well_qs = MagicMock()
        well_qs.exists.return_value = True
        well_qs.count.return_value = 1
        well_qs.__iter__ = lambda s: iter([well])
        mock_Well.objects.filter.return_value = well_qs
        mock_Well.objects.filter.return_value.order_by.return_value = well_qs

        wf = WellFinder(sg)
        result = wf.find_starting_material_wells(
            reaction_step_no=1,
            reaction_id=1,
            smiles="CCO",
            solvent=None,
            concentration=None,
            transfer_volume=50.0,
        )
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 1)

    @patch(f"{_WF_MOD}.get_previous_reaction_query_sets")
    @patch(f"{_WF_MOD}.Well")
    def test_step_2_fallback_to_reaction_plate(self, mock_Well, mock_prev_rxn):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder

        sg = _make_script_generator_stub()
        sg.volume_manager.check_volumes_close.side_effect = (
            lambda volume1, volume2: math.isclose(volume1, volume2, rel_tol=0.001)
        )

        prev_qs = MagicMock()
        mock_prev_rxn.return_value = prev_qs

        # No starting material wells
        empty_qs = MagicMock()
        empty_qs.exists.return_value = False
        empty_qs.count.return_value = 0
        empty_qs.__iter__ = lambda s: iter([])

        # Fallback well from reaction/workup plate
        fallback_well = _make_well_stub(role="reaction", reactantfornextstep=True)
        fallback_qs = MagicMock()
        fallback_qs.exists.return_value = True
        fallback_qs.__getitem__ = lambda s, i: fallback_well

        mock_Well.objects.filter.side_effect = [empty_qs, fallback_qs]
        mock_Well.objects.filter.return_value = empty_qs

        # Need to handle .order_by call on first filter
        empty_qs.order_by.return_value = empty_qs

        wf = WellFinder(sg)
        result = wf.find_starting_material_wells(
            reaction_step_no=2,
            reaction_id=1,
            smiles="CCO",
            solvent=None,
            concentration=None,
            transfer_volume=50.0,
        )
        self.assertIsNotNone(result)
        self.assertTrue(len(result) >= 1)


class TestWellFinderProcessWellObjects(TestCase):
    """Tests for WellFinder._process_well_objects helper."""

    def setUp(self):
        from backend.opentrons.otwriter.utils.well_finder import WellFinder
        self.sg = _make_script_generator_stub()
        self.sg.volume_manager.check_volumes_close.side_effect = (
            lambda volume1, volume2: math.isclose(volume1, volume2, rel_tol=0.001)
        )
        self.wf = WellFinder(self.sg)

    def test_single_well_sufficient(self):
        self.sg.volume_manager.get_well_volume_available.return_value = 200.0
        well = _make_well_stub(volume=200.0)
        well_info = []
        self.wf._process_well_objects(
            well_objects=[well],
            previous_reaction_queryset=None,
            transfer_volume=100.0,
            volume_manager=self.sg.volume_manager,
            well_info=well_info,
        )
        self.assertEqual(len(well_info), 1)
        self.assertAlmostEqual(well_info[0][2], 100.0)

    def test_two_wells_split_volume(self):
        self.sg.volume_manager.get_well_volume_available.side_effect = [60.0, 200.0]
        w1 = _make_well_stub(index=0, volume=60.0)
        w2 = _make_well_stub(index=1, volume=200.0)
        well_info = []
        self.wf._process_well_objects(
            well_objects=[w1, w2],
            previous_reaction_queryset=None,
            transfer_volume=100.0,
            volume_manager=self.sg.volume_manager,
            well_info=well_info,
        )
        self.assertEqual(len(well_info), 2)
        self.assertAlmostEqual(well_info[0][2], 60.0)
        self.assertAlmostEqual(well_info[1][2], 40.0)

    def test_empty_well_list(self):
        well_info = []
        self.wf._process_well_objects(
            well_objects=[],
            previous_reaction_queryset=None,
            transfer_volume=100.0,
            volume_manager=self.sg.volume_manager,
            well_info=well_info,
        )
        self.assertEqual(len(well_info), 0)


# ===================================================================
# QueryService tests
# ===================================================================
class TestQueryServiceInit(TestCase):
    """Tests for QueryService.__init__()"""

    def test_stores_references(self):
        from backend.opentrons.otwriter.utils.query_service import QueryService
        sg = _make_script_generator_stub()
        qs = QueryService(sg)
        self.assertEqual(qs.otsession_id, sg.otsession_id)
        self.assertEqual(qs.actionsession_ids, sg.actionsession_ids)


class TestQueryServiceGetActionSessions(TestCase):
    """Tests for QueryService.get_action_session_query_set()"""

    @patch(f"{_QS_MOD}.ActionSession")
    def test_queries_by_ids(self, mock_AS):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub(actionsession_ids=[5, 6])
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 2
        mock_AS.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_action_session_query_set()
        self.assertEqual(result, mock_result)


class TestQueryServiceGetAddActions(TestCase):
    """Tests for QueryService.get_add_action_query_set()"""

    @patch(f"{_QS_MOD}.AddAction")
    def test_no_filters(self, mock_AA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 3
        mock_AA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_add_action_query_set(reaction_ids=[1, 2])
        self.assertEqual(result, mock_result)

    @patch(f"{_QS_MOD}.AddAction")
    def test_with_session_type(self, mock_AA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 1
        mock_AA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_add_action_query_set(
            reaction_ids=[1], actionsessiontype="reaction"
        )
        self.assertEqual(result, mock_result)

    @patch(f"{_QS_MOD}.AddAction")
    def test_with_session_type_and_number(self, mock_AA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 1
        mock_AA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_add_action_query_set(
            reaction_ids=[1], actionsessiontype="reaction", actionnumber=2
        )
        self.assertEqual(result, mock_result)


class TestQueryServiceGetExtractActions(TestCase):
    """Tests for QueryService.get_extract_action_query_set()"""

    @patch(f"{_QS_MOD}.ExtractAction")
    def test_no_filters(self, mock_EA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 2
        mock_EA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_extract_action_query_set(reaction_ids=[1])
        self.assertEqual(result, mock_result)

    @patch(f"{_QS_MOD}.ExtractAction")
    def test_with_session_type(self, mock_EA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 1
        mock_EA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_extract_action_query_set(
            reaction_ids=[1], actionsessiontype="workup"
        )
        self.assertEqual(result, mock_result)


class TestQueryServiceGetMixActions(TestCase):
    """Tests for QueryService.get_mix_action_query_set()"""

    @patch(f"{_QS_MOD}.MixAction")
    def test_no_filters(self, mock_MA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 2
        mock_MA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_mix_action_query_set(reaction_ids=[1, 2])
        self.assertEqual(result, mock_result)

    @patch(f"{_QS_MOD}.MixAction")
    def test_with_action_session_ids(self, mock_MA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 1
        mock_MA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_mix_action_query_set(
            reaction_ids=[1], actionsession_ids=[10, 20]
        )
        self.assertEqual(result, mock_result)

    @patch(f"{_QS_MOD}.MixAction")
    def test_with_session_type_and_number(self, mock_MA):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 1
        mock_MA.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_mix_action_query_set(
            reaction_ids=[1], actionsessiontype="reaction", actionnumber=1
        )
        self.assertEqual(result, mock_result)


class TestQueryServiceGetPlates(TestCase):
    """Tests for QueryService.get_plates()"""

    @patch(f"{_QS_MOD}.Plate")
    def test_returns_queryset(self, mock_Plate):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        p1 = _make_plate_stub(name="plate_1")
        mock_result = MagicMock()
        mock_result.count.return_value = 1
        mock_result.__iter__ = lambda s: iter([p1])
        mock_Plate.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_plates()
        self.assertEqual(result, mock_result)


class TestQueryServiceGetPlateById(TestCase):
    """Tests for QueryService.get_plate_by_id()"""

    @patch(f"{_QS_MOD}.Plate")
    def test_returns_plate(self, mock_Plate):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        p = _make_plate_stub(id=42, name="test_plate")
        mock_Plate.objects.filter.return_value.first.return_value = p

        result = qs_obj.get_plate_by_id(42)
        self.assertEqual(result, p)

    @patch(f"{_QS_MOD}.Plate")
    def test_returns_none_when_missing(self, mock_Plate):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)
        mock_Plate.objects.filter.return_value.first.return_value = None

        result = qs_obj.get_plate_by_id(999)
        self.assertIsNone(result)


class TestQueryServiceGetTipRacks(TestCase):
    """Tests for QueryService.get_tip_racks()"""

    @patch(f"{_QS_MOD}.TipRack")
    def test_returns_queryset(self, mock_TR):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        tr = MagicMock()
        tr.name = "tiprack_1"
        tr.index = 1
        mock_result = MagicMock()
        mock_result.count.return_value = 1
        mock_result.__iter__ = lambda s: iter([tr])
        mock_TR.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_tip_racks()
        self.assertEqual(result, mock_result)


class TestQueryServiceGetPipette(TestCase):
    """Tests for QueryService.get_pipette()"""

    @patch(f"{_QS_MOD}.Pipette")
    def test_returns_pipette(self, mock_Pip):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        pip = MagicMock()
        pip.name = "p300_single"
        pip.type = "single"
        pip.position = "left"
        mock_Pip.objects.get.return_value = pip

        result = qs_obj.get_pipette()
        self.assertEqual(result, pip)

    @patch(f"{_QS_MOD}.Pipette")
    def test_raises_on_not_found(self, mock_Pip):
        from backend.opentrons.otwriter.utils.query_service import QueryService
        from backend.models import Pipette as RealPipette

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)
        mock_Pip.objects.get.side_effect = RealPipette.DoesNotExist("no pipette")

        with self.assertRaises(Exception):
            qs_obj.get_pipette()


class TestQueryServiceGetColumns(TestCase):
    """Tests for QueryService.get_column_query_set()"""

    @patch(f"{_QS_MOD}.Column")
    def test_returns_queryset(self, mock_Col):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_result = MagicMock()
        mock_result.count.return_value = 3
        mock_Col.objects.filter.return_value.order_by.return_value = mock_result

        result = qs_obj.get_column_query_set(
            role="reaction", role_index=1, reactionclass="amidation"
        )
        self.assertEqual(result, mock_result)


class TestQueryServiceGetUniqueClasses(TestCase):
    """Tests for QueryService.get_unique_reaction_classes()"""

    def test_returns_distinct_classes(self):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        rxn_qs = MagicMock()
        vl = MagicMock()
        vl.order_by.return_value.distinct.return_value = ["amidation", "suzuki"]
        rxn_qs.values_list.return_value = vl

        result = qs_obj.get_unique_reaction_classes(rxn_qs)
        self.assertEqual(list(result), ["amidation", "suzuki"])


class TestQueryServiceGetUniqueRecipes(TestCase):
    """Tests for QueryService.get_unique_reaction_recipes()"""

    def test_returns_recipes_for_class(self):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        rxn_qs = MagicMock()
        filtered = MagicMock()
        filtered.values_list.return_value.order_by.return_value.distinct.return_value = [
            "recipe_A",
            "recipe_B",
        ]
        rxn_qs.filter.return_value = filtered

        result = qs_obj.get_unique_reaction_recipes("amidation", rxn_qs)
        self.assertEqual(list(result), ["recipe_A", "recipe_B"])


class TestQueryServiceGetGroupedReactions(TestCase):
    """Tests for QueryService.get_grouped_reaction_by_class_recipe()"""

    def test_single_class_single_recipe(self):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        rxn_qs = MagicMock()
        rxn_qs.count.return_value = 2

        # Unique classes
        vl = MagicMock()
        vl.order_by.return_value.distinct.return_value = ["amidation"]
        rxn_qs.values_list.return_value = vl

        # Filter for class
        class_qs = MagicMock()
        class_qs.count.return_value = 2
        class_qs.__bool__ = lambda s: True
        class_qs.values_list.return_value.distinct.return_value.order_by.return_value = [
            "recipe_A"
        ]
        rxn_qs.filter.return_value = class_qs

        result = qs_obj.get_grouped_reaction_by_class_recipe(rxn_qs)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], class_qs)


class TestQueryServiceGetNextObjEntries(TestCase):
    """Tests for QueryService.get_next_obj_entries()"""

    def test_returns_proceeding_entries(self):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        obj = MagicMock()
        obj.pk = 5
        qs = MagicMock()
        next_qs = MagicMock()
        next_qs.count.return_value = 2
        qs.filter.return_value.order_by.return_value = next_qs

        result = qs_obj.get_next_obj_entries(qs, obj)
        self.assertEqual(result, next_qs)


# ===================================================================
# ScriptGenerator tests
# ===================================================================
class TestScriptGeneratorInit(TestCase):
    """Tests for ScriptGenerator.__init__() initialization."""

    @patch(f"{_SG_MOD}.WellFinder")
    @patch(f"{_SG_MOD}.VolumeManager")
    @patch(f"{_SG_MOD}.QueryService")
    @patch(f"{_SG_MOD}.FileManager")
    @patch(f"{_SG_MOD}.CommandGenerator")
    @patch(f"{_SG_MOD}.AnalysisSessionHandler")
    @patch(f"{_SG_MOD}.WorkupSessionHandler")
    @patch(f"{_SG_MOD}.ReactionSessionHandler")
    @patch(f"{_SG_MOD}.get_reaction_query_set")
    def test_sets_protocol_name(
        self, mock_rqs, mock_rh, mock_wh, mock_ah,
        mock_cg, mock_fm, mock_qs, mock_vm, mock_wf
    ):
        from backend.opentrons.otwriter.script_generator import ScriptGenerator

        otsession = MagicMock()
        otsession.id = 7
        otsession.sessiontype = "reaction"
        otsession.reactionstep = 2

        # Mock query results
        mock_qs_inst = MagicMock()
        mock_qs.return_value = mock_qs_inst
        as_qs = MagicMock()
        as_qs.exists.return_value = True
        as_qs.__iter__ = lambda s: iter([])
        mock_qs_inst.get_action_session_query_set.return_value = as_qs
        mock_qs_inst.get_tip_racks.return_value = MagicMock()
        mock_qs_inst.get_tip_racks.return_value.count.return_value = 1
        mock_qs_inst.get_tip_racks.return_value.__len__ = lambda s: 1
        pipette = MagicMock()
        pipette.name = "p300_single"
        mock_qs_inst.get_pipette.return_value = pipette
        mock_qs_inst.get_plates.return_value = MagicMock()
        mock_qs_inst.get_plates.return_value.count.return_value = 2

        fm_inst = MagicMock()
        fm_inst.create_file_path.return_value = ("/tmp/test.txt", "test.txt")
        mock_fm.return_value = fm_inst

        sg = ScriptGenerator(
            batchtag="BATCH1",
            otsessionobj=otsession,
            actionsession_ids=[1, 2],
        )
        self.assertEqual(sg.protocolname, "reaction-bBATCH1-r2-s7")
        self.assertEqual(sg.otsessiontype, "reaction")
        self.assertEqual(sg.batchtag, "BATCH1")


class TestScriptGeneratorAddCommand(TestCase):
    """Tests for ScriptGenerator.add_command()"""

    @patch(f"{_SG_MOD}.WellFinder")
    @patch(f"{_SG_MOD}.VolumeManager")
    @patch(f"{_SG_MOD}.QueryService")
    @patch(f"{_SG_MOD}.FileManager")
    @patch(f"{_SG_MOD}.CommandGenerator")
    @patch(f"{_SG_MOD}.AnalysisSessionHandler")
    @patch(f"{_SG_MOD}.WorkupSessionHandler")
    @patch(f"{_SG_MOD}.ReactionSessionHandler")
    @patch(f"{_SG_MOD}.get_reaction_query_set")
    def test_add_string(
        self, mock_rqs, mock_rh, mock_wh, mock_ah,
        mock_cg, mock_fm, mock_qs, mock_vm, mock_wf
    ):
        from backend.opentrons.otwriter.script_generator import ScriptGenerator

        otsession = MagicMock()
        otsession.id = 1
        otsession.sessiontype = "reaction"
        otsession.reactionstep = 1

        mock_qs_inst = MagicMock()
        mock_qs.return_value = mock_qs_inst
        as_qs = MagicMock()
        as_qs.exists.return_value = True
        as_qs.__iter__ = lambda s: iter([])
        mock_qs_inst.get_action_session_query_set.return_value = as_qs
        mock_qs_inst.get_tip_racks.return_value = MagicMock()
        mock_qs_inst.get_tip_racks.return_value.count.return_value = 1
        mock_qs_inst.get_tip_racks.return_value.__len__ = lambda s: 1
        pipette = MagicMock()
        pipette.name = "p300"
        mock_qs_inst.get_pipette.return_value = pipette
        mock_qs_inst.get_plates.return_value = MagicMock()
        mock_qs_inst.get_plates.return_value.count.return_value = 1

        fm_inst = MagicMock()
        fm_inst.create_file_path.return_value = ("/tmp/t.txt", "t.txt")
        mock_fm.return_value = fm_inst

        sg = ScriptGenerator("B", otsession, [1])
        sg.add_command("line1")
        self.assertIn("line1", sg.content)

        sg.add_command(["line2", "line3"])
        self.assertIn("line2", sg.content)
        self.assertIn("line3", sg.content)


class TestScriptGeneratorGenerateScript(TestCase):
    """Tests for ScriptGenerator.generate_script() routing."""

    @patch(f"{_SG_MOD}.WellFinder")
    @patch(f"{_SG_MOD}.VolumeManager")
    @patch(f"{_SG_MOD}.QueryService")
    @patch(f"{_SG_MOD}.FileManager")
    @patch(f"{_SG_MOD}.CommandGenerator")
    @patch(f"{_SG_MOD}.AnalysisSessionHandler")
    @patch(f"{_SG_MOD}.WorkupSessionHandler")
    @patch(f"{_SG_MOD}.ReactionSessionHandler")
    @patch(f"{_SG_MOD}.get_reaction_query_set")
    def _build_sg(
        self, session_type, mock_rqs, mock_rh, mock_wh, mock_ah,
        mock_cg, mock_fm, mock_qs, mock_vm, mock_wf
    ):
        from backend.opentrons.otwriter.script_generator import ScriptGenerator

        otsession = MagicMock()
        otsession.id = 1
        otsession.sessiontype = session_type
        otsession.reactionstep = 1

        mock_qs_inst = MagicMock()
        mock_qs.return_value = mock_qs_inst
        as_qs = MagicMock()
        as_qs.exists.return_value = True
        as_qs.__iter__ = lambda s: iter([])
        mock_qs_inst.get_action_session_query_set.return_value = as_qs
        mock_qs_inst.get_tip_racks.return_value = MagicMock()
        mock_qs_inst.get_tip_racks.return_value.count.return_value = 1
        mock_qs_inst.get_tip_racks.return_value.__len__ = lambda s: 1
        pipette = MagicMock()
        pipette.name = "p300"
        mock_qs_inst.get_pipette.return_value = pipette
        mock_qs_inst.get_plates.return_value = MagicMock()
        mock_qs_inst.get_plates.return_value.count.return_value = 1

        mock_cg_inst = MagicMock()
        mock_cg_inst.get_script_setup.return_value = ["# setup"]
        mock_cg_inst.get_labware_setup.return_value = ["# labware"]
        mock_cg_inst.get_tiprack_setup.return_value = ["# tipracks"]
        mock_cg_inst.get_pipette_setup.return_value = ["# pipette"]
        mock_cg_inst.get_number_tips_available_setup.return_value = ["# tips"]
        mock_cg_inst.get_pickup_tip_function.return_value = ["# pickup"]
        mock_cg_inst.get_drop_tip_function.return_value = ["# drop"]
        mock_cg.return_value = mock_cg_inst

        fm_inst = MagicMock()
        fm_inst.create_file_path.return_value = ("/tmp/t.txt", "t.txt")
        fm_inst.create_ot_script_model.return_value = MagicMock(id=42)
        mock_fm.return_value = fm_inst

        sg = ScriptGenerator("B", otsession, [1])
        return sg

    def test_reaction_routing(self):
        sg = self._build_sg("reaction")
        sg.generate_script()
        sg.reaction_handler.process_session.assert_called_once()
        sg.workup_handler.process_session.assert_not_called()
        sg.analysis_handler.process_session.assert_not_called()

    def test_workup_routing(self):
        sg = self._build_sg("workup")
        sg.generate_script()
        sg.workup_handler.process_session.assert_called_once()
        sg.reaction_handler.process_session.assert_not_called()

    def test_analyse_routing(self):
        sg = self._build_sg("analyse")
        sg.generate_script()
        sg.analysis_handler.process_session.assert_called_once()
        sg.reaction_handler.process_session.assert_not_called()


# ===================================================================
# ReactionSessionHandler tests
# ===================================================================
class TestReactionHandlerProcessSession(TestCase):
    """Tests for ReactionSessionHandler.process_session()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.reaction_handler import (
            ReactionSessionHandler,
        )
        self.sg = _make_script_generator_stub(reactionstep=1)
        self.handler = ReactionSessionHandler(self.sg)

    @patch(f"{_RH_MOD}.get_reaction")
    def test_iterates_action_sessions(self, mock_get_rxn):
        rxn_obj = MagicMock()
        rxn_obj.id = 1
        mock_get_rxn.return_value = rxn_obj

        qs = _make_action_session_qs(session_number=1)

        with patch.object(self.handler, "process_reaction_actions") as mock_pra:
            self.handler.process_session(qs)
            mock_pra.assert_called_once()

    @patch(f"{_RH_MOD}.get_reaction_query_set")
    @patch(f"{_RH_MOD}.get_reaction")
    def test_dilution_step_for_later_reactions(self, mock_get_rxn, mock_rqs):
        from backend.opentrons.otwriter.session_handlers.reaction_handler import (
            ReactionSessionHandler,
        )
        self.sg.reactionstep = 2
        handler = ReactionSessionHandler(self.sg)

        qs = _make_action_session_qs(session_number=1)

        rxn = MagicMock()
        rxn.id = 1
        mock_get_rxn.return_value = rxn

        with patch.object(handler, "process_dilution_step") as mock_dil, \
             patch.object(handler, "process_reaction_actions"):
            handler.process_session(qs)
            mock_dil.assert_called_once()


class TestReactionHandlerProcessAddAction(TestCase):
    """Tests for ReactionSessionHandler.process_add_action()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.reaction_handler import (
            ReactionSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = ReactionSessionHandler(self.sg)

    @patch(f"{_RH_MOD}.AddAction")
    def test_add_action_finds_wells_and_transfers(self, mock_AddAction):
        add_action_obj = MagicMock()
        add_action_obj.smiles = "CCCO"
        add_action_obj.solvent = "DMSO"
        add_action_obj.volume = 50.0
        add_action_obj.concentration = 0.5
        mock_AddAction.objects.get.return_value = add_action_obj

        # Mock well finder results
        from_well = _make_well_stub(index=0)
        to_well = _make_well_stub(index=5)
        self.sg.well_finder.find_starting_material_wells.return_value = [
            [None, from_well, 50.0]
        ]
        self.sg.well_finder.find_reaction_well.return_value = to_well

        from_plate = _make_plate_stub(name="sm_plate")
        to_plate = _make_plate_stub(name="rxn_plate")
        self.sg.query_service.get_plate_by_id.side_effect = [from_plate, to_plate]

        recipe_action = MagicMock()
        recipe_action.to_plate_role = "reaction"
        recipe_action.to_plate_role_index = 1
        recipe_action.from_plate_role = "startingmaterial"
        recipe_action.from_plate_role_index = 1

        # No mix action follows
        reaction_actions = [recipe_action]

        self.handler.process_add_action(
            reaction_action=recipe_action,
            action_number=1,
            index=0,
            reaction_actions=reaction_actions,
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        # Should have called transfer and tip commands
        self.sg.command_generator.pick_up_tip.assert_called()
        self.sg.command_generator.transfer_fluid_single.assert_called()
        self.sg.command_generator.drop_tip.assert_called()


class TestReactionHandlerProcessMixAction(TestCase):
    """Tests for ReactionSessionHandler.process_mix_action()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.reaction_handler import (
            ReactionSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = ReactionSessionHandler(self.sg)

    @patch(f"{_RH_MOD}.MixAction")
    def test_mix_action_generates_commands(self, mock_MixAction):
        mix_obj = MagicMock()
        mix_obj.plate_role = "reaction"
        mix_obj.plate_role_index = 1
        mix_obj.repetitions = 5
        mock_MixAction.objects.get.return_value = mix_obj

        mix_well = _make_well_stub(index=3)
        self.sg.well_finder.find_reaction_well.return_value = mix_well
        self.sg.query_service.get_plate_by_id.return_value = _make_plate_stub(
            name="rxn_plate"
        )

        self.handler.process_mix_action(
            action_number=1,
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        self.sg.command_generator.mix_well.assert_called_once_with(
            wellindex=3, nomixes=5, plate="rxn_plate"
        )
        self.sg.command_generator.drop_tip.assert_called()


class TestReactionHandlerProcessReactionActions(TestCase):
    """Tests for ReactionSessionHandler.process_reaction_actions()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.reaction_handler import (
            ReactionSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = ReactionSessionHandler(self.sg)

    @patch(f"{_RH_MOD}.get_session_recipe_actions")
    def test_routes_add_and_mix_actions(self, mock_recipe_actions):
        from backend.models import RecipeAddAction, RecipeMixAction

        add_action = MagicMock(spec=RecipeAddAction)
        add_action.action_number = 1
        mix_action = MagicMock(spec=RecipeMixAction)
        mix_action.action_number = 2
        mock_recipe_actions.return_value = [add_action, mix_action]

        rxn_obj = MagicMock()
        rxn_obj.id = 1
        rxn_obj.reactionclass = "amidation"
        rxn_obj.recipe = "standard"
        rxn_obj.intramolecular = False

        with patch.object(self.handler, "process_add_action") as mock_add, \
             patch.object(self.handler, "process_mix_action") as mock_mix:
            self.handler.process_reaction_actions(MagicMock(), rxn_obj, 1)
            mock_add.assert_called_once()
            mock_mix.assert_called_once()


# ===================================================================
# WorkupSessionHandler tests
# ===================================================================
class TestWorkupHandlerProcessSession(TestCase):
    """Tests for WorkupSessionHandler.process_session()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.workup_handler import (
            WorkupSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = WorkupSessionHandler(self.sg)

    @patch(f"{_WH_MOD}.get_reaction")
    def test_iterates_sessions(self, mock_get_rxn):
        rxn = MagicMock()
        rxn.id = 1
        mock_get_rxn.return_value = rxn

        qs = _make_action_session_qs(session_number=1)

        with patch.object(self.handler, "process_workup_actions") as mock_pwa:
            self.handler.process_session(qs)
            mock_pwa.assert_called_once()


class TestWorkupHandlerProcessWorkupActions(TestCase):
    """Tests for WorkupSessionHandler.process_workup_actions()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.workup_handler import (
            WorkupSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = WorkupSessionHandler(self.sg)

    @patch(f"{_WH_MOD}.get_session_recipe_actions")
    def test_routes_extract_and_mix(self, mock_recipe_actions):
        from backend.models import RecipeExtractAction, RecipeMixAction

        ext = MagicMock(spec=RecipeExtractAction)
        ext.action_number = 1
        mix = MagicMock(spec=RecipeMixAction)
        mix.action_number = 2
        mock_recipe_actions.return_value = [ext, mix]

        rxn = MagicMock()
        rxn.id = 1
        rxn.reactionclass = "amidation"
        rxn.recipe = "standard"
        rxn.intramolecular = False

        with patch.object(self.handler, "process_extract_action") as mock_ext, \
             patch.object(self.handler, "process_mix_action") as mock_mix:
            self.handler.process_workup_actions(MagicMock(), rxn, 1)
            mock_ext.assert_called_once()
            mock_mix.assert_called_once()

    @patch(f"{_WH_MOD}.get_session_recipe_actions")
    def test_no_actions_adds_comment(self, mock_recipe_actions):
        mock_recipe_actions.return_value = []
        rxn = MagicMock()
        rxn.id = 1
        rxn.reactionclass = "amidation"
        rxn.recipe = "standard"
        rxn.intramolecular = False

        self.handler.process_workup_actions(MagicMock(), rxn, 1)
        # Should add a comment about no actions found
        self.sg.add_command.assert_called()


class TestWorkupHandlerProcessExtractAction(TestCase):
    """Tests for WorkupSessionHandler.process_extract_action()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.workup_handler import (
            WorkupSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = WorkupSessionHandler(self.sg)

    @patch(f"{_WH_MOD}.ExtractAction")
    def test_extract_generates_transfer(self, mock_EA):
        ext_obj = MagicMock()
        ext_obj.from_plate_role = "reaction"
        ext_obj.from_plate_role_index = 1
        ext_obj.to_plate_role = "workup"
        ext_obj.to_plate_role_index = 1
        ext_obj.layer = "top"
        ext_obj.bottomlayervolume = 200.0
        ext_obj.volume = 100.0
        mock_EA.objects.get.return_value = ext_obj

        from_well = _make_well_stub(index=0)
        to_well = _make_well_stub(index=5)
        self.sg.well_finder.find_reaction_well.side_effect = [from_well, to_well]

        from_plate = _make_plate_stub(name="rxn_plate", labware="plateone_96_wellplate_2500ul")
        to_plate = _make_plate_stub(name="workup_plate")
        self.sg.query_service.get_plate_by_id.side_effect = [from_plate, to_plate]
        self.sg.volume_manager.calculate_aspirate_height.return_value = 5.0

        workup_action = MagicMock()
        workup_actions = [workup_action]  # Only one action, no mix follows

        self.handler.process_extract_action(
            workup_action=workup_action,
            action_number=1,
            index=0,
            workup_actions=workup_actions,
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        self.sg.command_generator.pick_up_tip.assert_called()
        self.sg.command_generator.transfer_fluid_single.assert_called_once()
        self.sg.command_generator.drop_tip.assert_called()


class TestWorkupHandlerProcessMixAction(TestCase):
    """Tests for WorkupSessionHandler.process_mix_action()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.workup_handler import (
            WorkupSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = WorkupSessionHandler(self.sg)

    @patch(f"{_WH_MOD}.MixAction")
    def test_mix_generates_commands(self, mock_MA):
        mix_obj = MagicMock()
        mix_obj.plate_role = "workup"
        mix_obj.plate_role_index = 1
        mix_obj.repetitions = 3
        mock_MA.objects.get.return_value = mix_obj

        mix_well = _make_well_stub(index=7)
        self.sg.well_finder.find_reaction_well.return_value = mix_well
        self.sg.query_service.get_plate_by_id.return_value = _make_plate_stub(
            name="workup_plate"
        )

        self.handler.process_mix_action(
            action_number=1,
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        self.sg.command_generator.mix_well.assert_called_once_with(
            wellindex=7, nomixes=3, plate="workup_plate"
        )
        self.sg.command_generator.drop_tip.assert_called()

    @patch(f"{_WH_MOD}.MixAction")
    def test_mix_not_found_adds_error_comment(self, mock_MA):
        mock_MA.objects.get.side_effect = mock_MA.DoesNotExist("not found")

        self.handler.process_mix_action(
            action_number=99,
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        # Should add error comment to content
        self.sg.add_command.assert_called()


# ===================================================================
# AnalysisSessionHandler tests
# ===================================================================
class TestAnalysisHandlerProcessSession(TestCase):
    """Tests for AnalysisSessionHandler.process_session()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.analysis_handler import (
            AnalysisSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = AnalysisSessionHandler(self.sg)

    @patch(f"{_AH_MOD}.get_reaction")
    def test_iterates_sessions(self, mock_get_rxn):
        rxn = MagicMock()
        rxn.id = 1
        mock_get_rxn.return_value = rxn

        qs = _make_action_session_qs(session_number=1)

        with patch.object(self.handler, "process_analysis_actions") as mock_paa:
            self.handler.process_session(qs)
            mock_paa.assert_called_once()


class TestAnalysisHandlerProcessAnalysisActions(TestCase):
    """Tests for AnalysisSessionHandler.process_analysis_actions()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.analysis_handler import (
            AnalysisSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = AnalysisSessionHandler(self.sg)

    @patch(f"{_AH_MOD}.get_session_recipe_actions")
    def test_routes_add_actions(self, mock_recipe_actions):
        from backend.models import RecipeAddAction

        add_action = MagicMock(spec=RecipeAddAction)
        add_action.action_number = 1
        mock_recipe_actions.return_value = [add_action]

        rxn = MagicMock()
        rxn.id = 1
        rxn.reactionclass = "amidation"
        rxn.recipe = "standard"
        rxn.intramolecular = False

        with patch.object(self.handler, "process_add_action") as mock_add:
            self.handler.process_analysis_actions(MagicMock(), rxn, 1)
            mock_add.assert_called_once()

    @patch(f"{_AH_MOD}.get_session_recipe_actions")
    def test_no_actions_adds_comment(self, mock_recipe_actions):
        mock_recipe_actions.return_value = []
        rxn = MagicMock()
        rxn.id = 1
        rxn.reactionclass = "amidation"
        rxn.recipe = "standard"
        rxn.intramolecular = False

        self.handler.process_analysis_actions(MagicMock(), rxn, 1)
        self.sg.add_command.assert_called()


class TestAnalysisHandlerProcessAddAction(TestCase):
    """Tests for AnalysisSessionHandler.process_add_action()"""

    def setUp(self):
        from backend.opentrons.otwriter.session_handlers.analysis_handler import (
            AnalysisSessionHandler,
        )
        self.sg = _make_script_generator_stub()
        self.handler = AnalysisSessionHandler(self.sg)

    @patch(f"{_AH_MOD}.AddAction")
    def test_add_action_generates_transfer(self, mock_AddAction):
        add_obj = MagicMock()
        add_obj.volume = 25.0
        add_obj.mix = False
        mock_AddAction.objects.get.return_value = add_obj

        from_well = _make_well_stub(index=0)
        to_well = _make_well_stub(index=10)
        self.sg.well_finder.find_reaction_well.side_effect = [from_well, to_well]

        from_plate = _make_plate_stub(name="rxn_plate")
        to_plate = _make_plate_stub(name="lcms_plate")
        self.sg.query_service.get_plate_by_id.side_effect = [from_plate, to_plate]

        analysis_action = MagicMock()
        analysis_action.from_plate_role = "reaction"
        analysis_action.from_plate_role_index = 1
        analysis_action.to_plate_role = "lcms"
        analysis_action.to_plate_role_index = 1

        self.handler.process_add_action(
            analysis_action=analysis_action,
            action_number=1,
            index=0,
            analysis_actions=[analysis_action],
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        self.sg.command_generator.pick_up_tip.assert_called()
        self.sg.command_generator.transfer_fluid_single.assert_called_once()
        self.sg.command_generator.drop_tip.assert_called()
        # mix should NOT be called since add_obj.mix is False
        self.sg.command_generator.mix_well.assert_not_called()

    @patch(f"{_AH_MOD}.AddAction")
    def test_add_action_with_mix(self, mock_AddAction):
        add_obj = MagicMock()
        add_obj.volume = 25.0
        add_obj.mix = True
        mock_AddAction.objects.get.return_value = add_obj

        from_well = _make_well_stub(index=0)
        to_well = _make_well_stub(index=10)
        self.sg.well_finder.find_reaction_well.side_effect = [from_well, to_well]

        from_plate = _make_plate_stub(name="rxn_plate")
        to_plate = _make_plate_stub(name="lcms_plate")
        self.sg.query_service.get_plate_by_id.side_effect = [from_plate, to_plate]

        analysis_action = MagicMock()
        analysis_action.from_plate_role = "reaction"
        analysis_action.from_plate_role_index = 1
        analysis_action.to_plate_role = "lcms"
        analysis_action.to_plate_role_index = 1

        self.handler.process_add_action(
            analysis_action=analysis_action,
            action_number=1,
            index=0,
            analysis_actions=[analysis_action],
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        self.sg.command_generator.mix_well.assert_called_once()

    @patch(f"{_AH_MOD}.AddAction")
    def test_add_action_not_found(self, mock_AddAction):
        mock_AddAction.objects.get.side_effect = mock_AddAction.DoesNotExist("nope")

        analysis_action = MagicMock()
        analysis_action.from_plate_role = "reaction"
        analysis_action.from_plate_role_index = 1
        analysis_action.to_plate_role = "lcms"
        analysis_action.to_plate_role_index = 1

        self.handler.process_add_action(
            analysis_action=analysis_action,
            action_number=1,
            index=0,
            analysis_actions=[analysis_action],
            actionsession_obj=MagicMock(),
            reaction_obj=MagicMock(id=1),
            reaction_id=1,
        )

        # Should add error comment
        self.sg.add_command.assert_called()


# ===================================================================
# QueryService - advanced query tests
# ===================================================================
class TestQueryServiceCheckNextReactionsAddActions(TestCase):
    """Tests for QueryService.check_next_reactions_add_actions()"""

    @patch(f"{_QS_MOD}.get_reaction_query_set")
    def test_finds_matching_add_actions(self, mock_rqs):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        rxn_obj = MagicMock()
        rxn_obj.id = 1
        rxn_obj.method_id.id = 10

        rxn_qs = MagicMock()
        next_rxn = MagicMock()
        next_rxn.id = 2
        next_qs = MagicMock()
        next_qs.count.return_value = 1
        next_qs.__iter__ = lambda s: iter([next_rxn])
        rxn_qs.filter.return_value.order_by.return_value = next_qs
        mock_rqs.return_value = rxn_qs

        add_match = MagicMock()
        add_match_qs = MagicMock()
        add_match_qs.__bool__ = lambda s: True
        add_match_qs.__getitem__ = lambda s, i: add_match
        with patch.object(qs_obj, "get_add_action_query_set") as mock_get_aa:
            mock_get_aa.return_value.filter.return_value = add_match_qs
            result = qs_obj.check_next_reactions_add_actions(rxn_obj, "CCC")

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0], add_match)

    @patch(f"{_QS_MOD}.get_reaction_query_set")
    def test_no_matching_add_actions(self, mock_rqs):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        rxn_obj = MagicMock()
        rxn_obj.id = 1
        rxn_obj.method_id.id = 10

        rxn_qs = MagicMock()
        next_qs = MagicMock()
        next_qs.count.return_value = 0
        next_qs.__iter__ = lambda s: iter([])
        rxn_qs.filter.return_value.order_by.return_value = next_qs
        mock_rqs.return_value = rxn_qs

        result = qs_obj.check_next_reactions_add_actions(rxn_obj, "CCC")
        self.assertEqual(len(result), 0)


class TestQueryServiceGetWellByReactionId(TestCase):
    """Tests for QueryService.get_well_by_reaction_id()"""

    @patch(f"{_QS_MOD}.get_product_smiles")
    @patch(f"{_QS_MOD}.Well")
    def test_finds_well(self, mock_Well, mock_get_smiles):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        qs_obj = QueryService(sg)

        mock_get_smiles.return_value = ["CCC"]
        expected_well = _make_well_stub()
        mock_Well.objects.get.return_value = expected_well

        result = qs_obj.get_well_by_reaction_id(
            reaction_id=1, role="reaction", role_index=1
        )
        self.assertEqual(result, expected_well)


class TestQueryServiceFindSolventPlateWellObj(TestCase):
    """Tests for QueryService.find_solvent_plate_well_obj()"""

    @patch(f"{_QS_MOD}.Plate")
    def test_finds_solvent_wells(self, mock_Plate):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        sg.volume_manager.check_volumes_close.side_effect = (
            lambda volume1, volume2: math.isclose(volume1, volume2, rel_tol=0.001)
        )
        sg.volume_manager.get_well_volume_available.return_value = 500.0

        well = _make_well_stub(volume=500.0, role="solvent", solvent="DMSO")

        plate = _make_plate_stub(role="solvent", name="solvent_plate")
        well_qs = MagicMock()
        well_qs.count.return_value = 1
        well_qs.__iter__ = lambda s: iter([well])
        plate.well_set.all.return_value.filter.return_value = well_qs

        plate_qs = MagicMock()
        plate_qs.count.return_value = 1
        plate_qs.__bool__ = lambda s: True
        plate_qs.__iter__ = lambda s: iter([plate])
        mock_Plate.objects.filter.return_value = plate_qs

        qs_obj = QueryService(sg)
        result = qs_obj.find_solvent_plate_well_obj("DMSO", 100.0)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0][0], well)

    @patch(f"{_QS_MOD}.Plate")
    def test_no_solvent_plates_returns_empty(self, mock_Plate):
        from backend.opentrons.otwriter.utils.query_service import QueryService

        sg = _make_script_generator_stub()
        plate_qs = MagicMock()
        plate_qs.count.return_value = 0
        plate_qs.__bool__ = lambda s: False
        mock_Plate.objects.filter.return_value = plate_qs

        qs_obj = QueryService(sg)
        result = qs_obj.find_solvent_plate_well_obj("DMSO", 100.0)
        self.assertEqual(len(result), 0)
