"""
Command Generator for OpenTrons protocols.

This module provides methods for generating protocol commands like
transfers, mixing, tip handling, etc.
"""

import logging
from typing import Dict, List

logger = logging.getLogger(__name__)


class CommandGenerator:
    """
    Generates OpenTrons protocol commands.

    This class is responsible for creating the actual Python commands
    that will be included in the generated script.
    """

    def __init__(self, script_generator):
        """
        Initialize the command generator.

        Parameters
        ----------
        script_generator : ScriptGenerator
            Parent script generator
        """
        self.script_generator = script_generator
        logger.info("CommandGenerator initialized")

    def comment(self, text: str, num_newlines: int = 1, level: int = 1) -> str:
        """Generate a properly indented comment.

        Parameters
        ----------
        text : str
            Comment text
        num_newlines : int, optional
            Number of newlines to add before the comment, by default 1
        level : int, optional
            Indentation level, by default 1

        Returns
        -------
        str
            Formatted comment with proper indentation
        """
        logger.info(f"Adding comment: {text}")
        prefix = "\n" * num_newlines if num_newlines > 0 else ""
        return f"{prefix}{self.indent(level)}# {text}"

    def indent(self, level: int = 1) -> str:
        """Return proper indentation string using tabs."""
        return "\t" * level

    def get_script_setup(self, protocolname: str, apiLevel: str) -> List[str]:
        """Generates header information for an OT script.

        Parameters
        ----------
        protocolname : str
            Name of the protocol
        apiLevel : str
            OpenTrons API level to use

        Returns
        -------
        List[str]
            Script header commands
        """
        logger.info(
            f"Generating script setup for protocol '{protocolname}' with API level {apiLevel}"
        )

        return [
            "from opentrons import protocol_api",
            f"# {protocolname} produced by XChem Car (https://.xchem.diamond.ac.uk)",
            "# metadata",
            f"metadata = {{'protocolName': '{protocolname}','apiLevel': '{apiLevel}'}}",
            "",
            "def run(protocol: protocol_api.ProtocolContext):",
        ]

    def get_labware_setup(self, plates: List[Dict]) -> List[str]:
        """Writes the plate setup instructions for an OT script.

        Parameters
        ----------
        plates : List[Dict]
            List of plate information dictionaries

        Returns
        -------
        List[str]
            Labware setup commands
        """
        logger.info(f"Generating labware setup for {len(plates)} plates")

        commands = [f"\n{self.indent()}# labware"]
        for plate in plates:
            platename = plate["name"]
            labware = plate["labware"]
            plateindex = plate["index"]
            logger.info(
                f"Adding plate '{platename}' ({labware}) at position {plateindex}"
            )
            commands.append(
                f"{self.indent()}{platename} = protocol.load_labware('{labware}', '{plateindex}')"
            )
        return commands

    def get_tiprack_setup(self, tipracks: List[Dict]) -> List[str]:
        """Writes the tipracks setup instructions for an OT script.

        Parameters
        ----------
        tipracks : List[Dict]
            List of tiprack information dictionaries

        Returns
        -------
        List[str]
            Tiprack setup commands
        """
        logger.info(f"Generating tiprack setup for {len(tipracks)} tipracks")

        commands = [f"\n{self.indent()}# tipracks"]
        for tiprack in tipracks:
            name = tiprack["name"]
            labware = tiprack["labware"]
            index = tiprack["index"]
            logger.info(f"Adding tiprack '{name}' ({labware}) at position {index}")
            commands.append(
                f"{self.indent()}{name} = protocol.load_labware('{labware}', '{index}')"
            )
        return commands

    def get_pipette_setup(
        self,
        pipette_name: str,
        pipette_labware: str,
        pipette_position: str,
        tiprack_names: List[str],
    ) -> List[str]:
        """Writes the pipette setup instructions for an OT script.

        Parameters
        ----------
        pipette_name : str
            Name of the pipette
        pipette_labware : str
            Type of pipette
        pipette_position : str
            Position on the robot (left/right)
        tiprack_names : List[str]
            List of tiprack names

        Returns
        -------
        List[str]
            Pipette setup commands
        """
        logger.info(
            f"Generating pipette setup for '{pipette_name}' ({pipette_labware}) at {pipette_position} position"
        )
        logger.info(
            f"Pipette will use {len(tiprack_names)} tiprack(s): {', '.join(tiprack_names)}"
        )

        return [
            self.comment("# pipettes"),
            f"{self.indent()}{pipette_name} = protocol.load_instrument('{pipette_labware}', '{pipette_position}', tip_racks=[{','.join(tiprack_names)}])",
        ]

    def get_number_tips_available_setup(
        self, num_tipracks: int, channel_type: str
    ) -> List[str]:
        """Captures number of tips available in tipracks.

        Parameters
        ----------
        num_tipracks : int
            Number of tip racks
        channel_type : str
            Type of channel ('single' or 'multi')

        Returns
        -------
        List[str]
            Tip tracking setup commands
        """
        numbertipsavailable = num_tipracks * 96
        logger.info(
            f"Setting up tip tracking for {channel_type} channel pipette with {num_tipracks} tiprack(s)"
        )
        logger.info(f"Total of {numbertipsavailable} tips available")

        return [
            f'{self.indent()}tipstate = {{"channeltype": "{channel_type}", "maxnumbertips": {numbertipsavailable}, "notipsavailable": {numbertipsavailable}}}'
        ]

    def get_pickup_tip_function(self, pipette_name: str) -> List[str]:
        """Function definition for picking up a tip.

        Parameters
        ----------
        pipette_name : str
            Name of the pipette

        Returns
        -------
        List[str]
            Pickup tip function definition
        """
        logger.info(f"Generating pickUpTip function for {pipette_name} pipette")

        return [
            f"\n\n{self.indent()}def pickUpTip():",
            f'{self.indent(2)}if tipstate["notipsavailable"] == 0:',
            f'{self.indent(3)}protocol.pause("Please replace tips")',
            f"{self.indent(3)}{pipette_name}.reset_tipracks()",
            f'{self.indent(3)}tipstate["notipsavailable"] = tipstate["maxnumbertips"]',
            f"{self.indent(2)}if not {pipette_name}.has_tip:",
            f"{self.indent(3)}{pipette_name}.pick_up_tip()",
            f'{self.indent(3)}if tipstate["channeltype"] == "multi":',
            f'{self.indent(4)}tipstate["notipsavailable"] = tipstate["notipsavailable"] - 8',
            f'{self.indent(3)}if tipstate["channeltype"] == "single":',
            f'{self.indent(4)}tipstate["notipsavailable"] = tipstate["notipsavailable"] - 1',
        ]

    def get_drop_tip_function(self, pipette_name: str) -> List[str]:
        """Function definition for dropping a tip.

        Parameters
        ----------
        pipette_name : str
            Name of the pipette

        Returns
        -------
        List[str]
            Drop tip function definition
        """
        logger.info(f"Generating dropTip function for {pipette_name} pipette")

        return [
            f"\n\n{self.indent()}def dropTip():",
            f"{self.indent(2)}if {pipette_name}.has_tip:",
            f"{self.indent(3)}{pipette_name}.drop_tip()",
        ]

    def pick_up_tip(self) -> List[str]:
        """Generate pick up tip command."""
        humanread = "Picking up a new tip"
        logger.info("Generating pick up tip command")

        return [
            self.comment(humanread),
            f"{self.indent()}pickUpTip()",  # Add self.indent() here
        ]

    def drop_tip(self) -> List[str]:
        """Generate drop tip command."""
        humanread = "Sending tip to waste"
        logger.info("Generating drop tip command")

        return [
            self.comment(humanread),
            f"{self.indent()}dropTip()",  # Add self.indent() here
        ]

    def pause_protocol(self, message: str) -> List[str]:
        """Generate pause protocol command."""
        humanread = "Pausing protocol operation"
        logger.info(f"Generating protocol pause with message: '{message}'")

        return [self.comment(humanread), f'{self.indent()}protocol.pause("{message}")']

    def delay_protocol(self, delay: int) -> List[str]:
        """Generate delay protocol command."""
        humanread = "Delaying protocol operation"
        logger.info(f"Generating protocol delay for {delay} seconds")

        return [
            self.comment(humanread),
            f"{self.indent()}protocol.delay(seconds={delay})",
        ]

    def set_aspirate_speed(self, speed: int) -> List[str]:
        """Generate set aspirate speed command."""
        humanread = "Setting aspirate speed"
        pipette_name = self.script_generator.pipettename
        logger.info(f"Setting aspirate speed for {pipette_name} to {speed} µL/s")

        return [
            self.comment(humanread),
            f"{self.indent()}{pipette_name}.flow_rate.aspirate={speed}",
        ]

    def set_dispense_speed(self, speed: int) -> List[str]:
        """Generate set dispense speed command."""
        humanread = "Setting dispense speed"
        pipette_name = self.script_generator.pipettename
        logger.info(f"Setting dispense speed for {pipette_name} to {speed} µL/s")

        return [
            self.comment(humanread),
            f"{self.indent()}{pipette_name}.flow_rate.dispense={speed}",
        ]

    def transfer_fluid_single(
        self,
        aspirateplatename: str,
        dispenseplatename: str,
        aspiratewellindex: int,
        dispensewellindex: int,
        transvolume: float,
        aspirateheight: int = 0.1,
        transfertype: str = "standard",
    ) -> List[str]:
        """Generate transfer fluid command for single channel.

        Parameters
        ----------
        aspirateplatename : str
            Name of the plate to aspirate from
        dispenseplatename : str
            Name of the plate to dispense to
        aspiratewellindex : int
            Well index to aspirate from
        dispensewellindex : int
            Well index to dispense to
        transvolume : float
            Volume to transfer in µL
        aspirateheight : int, optional
            Height to aspirate from in mm, by default 0.1
        transfertype : str, optional
            Type of transfer, by default "standard"

        Returns
        -------
        List[str]
            Transfer fluid commands
        """
        humanread = f"transfertype - {transfertype} - transfer - {transvolume:.1f}ul from {aspiratewellindex} to {dispensewellindex}"

        logger.info(f"Generating single-channel fluid transfer: {transfertype}")
        logger.info(
            f"Transfer: {transvolume:.1f} µL from {aspirateplatename}[{aspiratewellindex}] to {dispenseplatename}[{dispensewellindex}]"
        )
        logger.info(f"Aspirate height: {aspirateheight} mm")

        pipette_name = self.script_generator.pipettename
        max_volume = f"{pipette_name}.max_volume"
        air_gap_volume = f"{max_volume}*0.05"

        if transvolume <= 0:
            logger.warning(
                f"Transfer volume is {transvolume} µL which is <= 0, this may cause errors"
            )

        if transvolume > 1000:
            logger.warning(
                f"Transfer volume is {transvolume} µL which is unusually large"
            )

        return [
            self.comment(humanread),
            f"{self.indent()}{pipette_name}.transfer({transvolume}, {aspirateplatename}.wells()[{aspiratewellindex}].bottom({aspirateheight}), {dispenseplatename}.wells()[{dispensewellindex}].top({dispenseplatename}.highest_z*0.05), air_gap = {air_gap_volume}, touch_tip=True, new_tip='never', blow_out=True, blowout_location='destination well')",
        ]

    def transfer_fluid_multi(
        self,
        aspirateplatename: str,
        dispenseplatename: str,
        aspiratecolumnindex: int,
        dispensecolumnindex: int,
        transvolume: float,
        aspirateheight: int = 0.1,
        transfertype: str = "standard",
    ) -> List[str]:
        """Generate transfer fluid command for multi channel.

        Parameters
        ----------
        aspirateplatename : str
            Name of the plate to aspirate from
        dispenseplatename : str
            Name of the plate to dispense to
        aspiratecolumnindex : int
            Column index to aspirate from
        dispensecolumnindex : int
            Column index to dispense to
        transvolume : float
            Volume to transfer in µL
        aspirateheight : int, optional
            Height to aspirate from in mm, by default 0.1
        transfertype : str, optional
            Type of transfer, by default "standard"

        Returns
        -------
        List[str]
            Transfer fluid commands
        """
        humanread = f"transfertype - {transfertype} - transfer - {transvolume:.1f}ul from {aspiratecolumnindex} column to {dispensecolumnindex} column"

        logger.info(f"Generating multi-channel fluid transfer: {transfertype}")
        logger.info(
            f"Transfer: {transvolume:.1f} µL from {aspirateplatename} column {aspiratecolumnindex} to {dispenseplatename} column {dispensecolumnindex}"
        )
        logger.info(f"Aspirate height: {aspirateheight} mm")

        pipette_name = self.script_generator.pipettename
        max_volume = f"{pipette_name}.max_volume"
        air_gap_volume = f"{max_volume}*0.05"

        if transvolume <= 0:
            logger.warning(
                f"Transfer volume is {transvolume} µL which is <= 0, this may cause errors"
            )

        if transvolume > 1000:
            logger.warning(
                f"Transfer volume is {transvolume} µL which is unusually large"
            )

        return [
            self.comment(humanread),
            f"{self.indent()}{pipette_name}.transfer({transvolume}, {aspirateplatename}.columns()[{aspiratecolumnindex}][0].bottom({aspirateheight}), {dispenseplatename}.columns()[{dispensecolumnindex}][0].top({dispenseplatename}.highest_z*0.05), air_gap = {air_gap_volume}, touch_tip=True, new_tip='never', blow_out=True, blowout_location='destination well')",
        ]

    def mix_well(self, wellindex: int, nomixes: int, plate: str) -> List[str]:
        """Generate mix well command."""
        humanread = f"Mixing contents of plate: {plate} at well index: {wellindex}"

        logger.info(f"Generating mix well command for {plate}[{wellindex}]")
        logger.info(f"Mix repetitions: {nomixes}")

        pipette_name = self.script_generator.pipettename
        mix_volume = f"({pipette_name}.max_volume * 0.7)"

        if nomixes <= 0:
            logger.warning(
                f"Number of mixes is {nomixes}, which is <= 0. This may cause errors."
            )

        if nomixes > 10:
            logger.warning(f"Number of mixes is {nomixes}, which is unusually high.")

        return [
            self.comment(humanread),
            f"{self.indent()}{pipette_name}.mix({nomixes}, {mix_volume}, {plate}.wells()[{wellindex}])",
        ]

    def mix_column(self, columnindex: int, nomixes: int, plate: str) -> List[str]:
        """Generate mix column command."""
        humanread = f"Mixing contents of plate: {plate} at column index: {columnindex}"

        logger.info(f"Generating mix column command for {plate} column {columnindex}")
        logger.info(f"Mix repetitions: {nomixes}")

        pipette_name = self.script_generator.pipettename
        mix_volume = f"({pipette_name}.max_volume * 0.7)"

        if nomixes <= 0:
            logger.warning(
                f"Number of mixes is {nomixes}, which is <= 0. This may cause errors."
            )

        if nomixes > 10:
            logger.warning(f"Number of mixes is {nomixes}, which is unusually high.")

        if (
            columnindex < 0 or columnindex > 11
        ):  # Standard 96-well plates have 12 columns, indexed 0-11
            logger.warning(
                f"Column index {columnindex} may be out of range for standard plates"
            )

        return [
            self.comment(humanread),
            f"{self.indent()}{pipette_name}.mix({nomixes}, {mix_volume}, {plate}.columns()[{columnindex}][0])",
        ]
