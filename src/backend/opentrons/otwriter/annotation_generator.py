"""
Annotation Generator for OpenTrons protocol scripts.

This module centralises every human-readable comment that is inserted
into a generated OT script.  Keeping all annotation text here makes it
easy to adjust wording, avoids duplicated descriptions between the
session handlers and the command generator, and gives a single place to
extend when new action types are introduced.
"""

import logging
from typing import List, Optional

logger = logging.getLogger(__name__)

# Separator used for section banners
_SEP = "=" * 60


class AnnotationGenerator:
    """Single source for all human-readable comments in OT scripts."""

    def __init__(self, script_generator):
        """
        Initialise the annotation generator.

        Parameters
        ----------
        script_generator : ScriptGenerator
            Parent script generator (used only for reading state, never
            for writing commands).
        """
        self.script_generator = script_generator
        logger.info("AnnotationGenerator initialized")

    # ------------------------------------------------------------------
    # Low-level formatting helpers
    # ------------------------------------------------------------------

    def indent(self, level: int = 1) -> str:
        """Return proper indentation string using tabs."""
        return "\t" * level

    def comment(self, text: str, num_newlines: int = 1, level: int = 1) -> str:
        """Generate a properly indented comment line.

        Parameters
        ----------
        text : str
            Comment text
        num_newlines : int, optional
            Number of blank lines before the comment, by default 1
        level : int, optional
            Indentation level (tab count), by default 1

        Returns
        -------
        str
            Formatted comment with indentation and optional leading newlines
        """
        prefix = "\n" * num_newlines if num_newlines > 0 else ""
        return f"{prefix}{self.indent(level)}# {text}"

    # ------------------------------------------------------------------
    # Section / reaction headers
    # ------------------------------------------------------------------

    def reaction_header(
        self, reaction_obj, session_label: str = "Reaction"
    ) -> List[str]:
        """Build a banner block identifying a reaction in the script.

        Parameters
        ----------
        reaction_obj : Reaction
            Django Reaction model instance
        session_label : str, optional
            Session type label (``Reaction``, ``Workup``, ``Analysis``),
            by default ``"Reaction"``

        Returns
        -------
        List[str]
            Lines to append to the script content
        """
        reaction_id = reaction_obj.id
        reaction_class = reaction_obj.reactionclass

        reactant_smiles = list(
            reaction_obj.reactants.values_list("smiles", flat=True)
        )
        product_smiles = list(
            reaction_obj.products.values_list("smiles", flat=True)
        )

        lines: List[str] = [self.comment(_SEP, num_newlines=2)]

        lines.append(
            self.comment(
                f"{session_label}: {reaction_class} (Reaction ID: {reaction_id})",
                num_newlines=0,
            )
        )

        if reactant_smiles:
            bb_parts = [f"BB{i+1} ({s})" for i, s in enumerate(reactant_smiles)]
            lines.append(
                self.comment(f"Reactants: {', '.join(bb_parts)}", num_newlines=0)
            )

        if product_smiles:
            lines.append(
                self.comment(
                    f"Product: {', '.join(product_smiles)}", num_newlines=0
                )
            )

        lines.append(self.comment(_SEP, num_newlines=0))
        return lines

    def dilution_header(self) -> str:
        """Return a banner comment for the dilution section."""
        return self.comment(
            _SEP
            + "\n\t# Dilution step: adding solvent to previous reaction products"
            + "\n\t# "
            + _SEP,
            num_newlines=2,
        )

    def multichannel_header(self) -> str:
        """Return a banner comment for the multichannel transfer section."""
        return self.comment(
            _SEP
            + "\n\t# Multichannel transfers: column-to-column starting material additions"
            + "\n\t# "
            + _SEP,
            num_newlines=2,
        )

    # ------------------------------------------------------------------
    # Per-action summary comments (session-level)
    # ------------------------------------------------------------------

    def action_summary(
        self,
        action_label: str,
        action_number: int,
        total_actions: int,
        description: str,
    ) -> str:
        """Return a one-line action summary comment.

        Example output::

            # Add action 1 of 4: Transfer 10.0 uL of CCO …

        Parameters
        ----------
        action_label : str
            Short verb, e.g. ``"Add"``, ``"Extract"``, ``"Mix"``
        action_number : int
            1-based index of this action
        total_actions : int
            Total number of actions in the current session
        description : str
            Free-form text describing what this action does
        """
        return self.comment(
            f"{action_label} action {action_number} of {total_actions}: {description}",
            num_newlines=1,
        )

    # ------------------------------------------------------------------
    # Reaction-handler specific helpers
    # ------------------------------------------------------------------

    def add_action_description(
        self,
        volume: float,
        smiles: str,
        solvent: Optional[str],
        from_plate_role: str,
        to_plate_role: str,
    ) -> str:
        """Build the description text for an add-action summary.

        Parameters
        ----------
        volume : float
            Transfer volume in uL
        smiles : str
            SMILES of the material being transferred
        solvent : str or None
            Solvent name, if applicable
        from_plate_role : str
            Source plate role
        to_plate_role : str
            Destination plate role

        Returns
        -------
        str
            Descriptive text (without the ``#`` prefix)
        """
        solvent_info = f" in {solvent}" if solvent else ""
        return (
            f"Transfer {volume:.1f} uL of {smiles}{solvent_info} "
            f"from {from_plate_role} plate to {to_plate_role} plate"
        )

    def extract_action_description(
        self,
        layer: str,
        from_plate_role: str,
        to_plate_role: str,
    ) -> str:
        """Build description text for an extract-action summary."""
        return (
            f"Extract ({layer} layer) from "
            f"{from_plate_role} plate to {to_plate_role} plate"
        )

    def mix_action_description(self, plate_role: str) -> str:
        """Build description text for a mix-action summary."""
        return f"Mix contents in {plate_role} plate"

    def analysis_transfer_description(
        self, from_plate_role: str, to_plate_role: str
    ) -> str:
        """Build description text for an analysis transfer summary."""
        return (
            f"Transfer sample from {from_plate_role} plate "
            f"to {to_plate_role} plate"
        )

    # ------------------------------------------------------------------
    # Command-level annotations  (replace the old ``humanread`` strings)
    # ------------------------------------------------------------------

    def transfer_single(
        self,
        transfertype: str,
        transvolume: float,
        aspirateplatename: str,
        aspiratewellindex: int,
        dispenseplatename: str,
        dispensewellindex: int,
    ) -> str:
        """Human-readable comment for a single-channel transfer command."""
        return self.comment(
            f"Transfer ({transfertype}): {transvolume:.1f} uL from "
            f"{aspirateplatename} well {aspiratewellindex} to "
            f"{dispenseplatename} well {dispensewellindex}",
        )

    def transfer_multi(
        self,
        transfertype: str,
        transvolume: float,
        aspirateplatename: str,
        aspiratecolumnindex: int,
        dispenseplatename: str,
        dispensecolumnindex: int,
    ) -> str:
        """Human-readable comment for a multi-channel transfer command."""
        return self.comment(
            f"Transfer ({transfertype}): {transvolume:.1f} uL from "
            f"{aspirateplatename} column {aspiratecolumnindex} to "
            f"{dispenseplatename} column {dispensecolumnindex}",
        )

    def mix_well(self, nomixes: int, plate: str, wellindex: int) -> str:
        """Human-readable comment for a mix-well command."""
        return self.comment(f"Mix {nomixes} times in {plate} at well {wellindex}")

    def mix_column(self, nomixes: int, plate: str, columnindex: int) -> str:
        """Human-readable comment for a mix-column command."""
        return self.comment(f"Mix {nomixes} times in {plate} at column {columnindex}")

    def pick_up_tip(self) -> str:
        """Human-readable comment for picking up a tip."""
        return self.comment("Picking up a new tip")

    def drop_tip(self) -> str:
        """Human-readable comment for dropping a tip."""
        return self.comment("Sending tip to waste")

    def pause_protocol(self) -> str:
        """Human-readable comment for pausing the protocol."""
        return self.comment("Pausing protocol operation")

    def delay_protocol(self) -> str:
        """Human-readable comment for delaying the protocol."""
        return self.comment("Delaying protocol operation")

    def set_aspirate_speed(self) -> str:
        """Human-readable comment for setting aspirate speed."""
        return self.comment("Setting aspirate speed")

    def set_dispense_speed(self) -> str:
        """Human-readable comment for setting dispense speed."""
        return self.comment("Setting dispense speed")

    def dilution_transfer(
        self,
        volume: float,
        solvent: str,
        aspirateplatename: str,
        aspiratewellindex: int,
        dispenseplatename: str,
        dispensewellindex: int,
    ) -> str:
        """Human-readable comment for a dilution transfer."""
        return self.comment(
            f"Dilution: Transfer {volume:.1f} uL of {solvent} "
            f"from {aspirateplatename} well {aspiratewellindex} "
            f"to {dispenseplatename} well {dispensewellindex}",
        )

    def reaction_mix(
        self,
        repetitions: int,
        plate_role: str,
        plate_name: str,
        well_index: int,
    ) -> str:
        """Human-readable comment for a reaction-session mix action."""
        return self.comment(
            f"Mix {repetitions} times in {plate_role} plate "
            f"({plate_name}) at well {well_index}",
        )
