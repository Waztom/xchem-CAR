"""
Transfer Record for OpenTrons protocol scripts.

Lightweight dataclass that captures every liquid-handling event during
script generation.  The ``ScriptGenerator`` maintains a list of these
records (the *transfer ledger*) which downstream consumers — such as the
``SessionVisualizer`` — read to build visualizations, audit trails, or
validation reports.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import List, Optional

logger = logging.getLogger(__name__)


@dataclass
class TransferRecord:
    """One liquid-handling event captured during script generation.

    Attributes
    ----------
    action_type : str
        ``"add"``, ``"extract"``, ``"mix"``, ``"dilution"``, ``"analysis"``
    source_plate_name : str
        Human-readable plate name (e.g. ``"sm-plate-1"``)
    source_plate_role : str
        Plate role (e.g. ``"startingmaterial"``, ``"reaction"``)
    source_well_index : int
        0-based well index on the source plate
    source_well_name : str
        Human-readable well name (e.g. ``"A1"``)
    dest_plate_name : str
        Destination plate name
    dest_plate_role : str
        Destination plate role
    dest_well_index : int
        0-based well index on the destination plate
    dest_well_name : str
        Human-readable destination well name
    volume : float
        Volume transferred in µL
    smiles : str or None
        SMILES of the material being transferred
    solvent : str or None
        Solvent name, if applicable
    reaction_id : int or None
        Associated reaction ID
    reaction_class : str or None
        Reaction class label (e.g. ``"amidation"``)
    recipe : str or None
        Recipe name
    transfer_mode : str
        ``"single"`` (default) or ``"multichannel"``
    """

    action_type: str
    source_plate_name: str
    source_plate_role: str
    source_well_index: int
    source_well_name: str
    dest_plate_name: str
    dest_plate_role: str
    dest_well_index: int
    dest_well_name: str
    volume: float
    smiles: Optional[str] = None
    solvent: Optional[str] = None
    reaction_id: Optional[int] = None
    reaction_class: Optional[str] = None
    recipe: Optional[str] = None
    transfer_mode: str = "single"


class TransferLedger:
    """Ordered collection of :class:`TransferRecord` instances.

    Attached to ``ScriptGenerator`` so that session handlers can call
    ``self.script_generator.transfer_ledger.record(...)`` while generating
    commands.
    """

    def __init__(self) -> None:
        self.records: List[TransferRecord] = []

    # ------------------------------------------------------------------
    # Recording helpers
    # ------------------------------------------------------------------

    def record(self, **kwargs) -> TransferRecord:
        """Create a :class:`TransferRecord` and append it to the ledger.

        All keyword arguments are forwarded to the dataclass constructor.
        Returns the newly created record.
        """
        rec = TransferRecord(**kwargs)
        self.records.append(rec)
        logger.debug(
            f"TransferLedger: {rec.action_type} "
            f"{rec.source_plate_name}:{rec.source_well_name} → "
            f"{rec.dest_plate_name}:{rec.dest_well_name} "
            f"({rec.volume} µL)"
        )
        return rec

    # ------------------------------------------------------------------
    # Query helpers (used by SessionVisualizer)
    # ------------------------------------------------------------------

    def get_records_for_plate(self, plate_name: str) -> List[TransferRecord]:
        """Return all records where *plate_name* is source or dest."""
        return [
            r
            for r in self.records
            if r.source_plate_name == plate_name or r.dest_plate_name == plate_name
        ]

    def get_incoming(self, plate_name: str) -> List[TransferRecord]:
        """Return records where material arrives at *plate_name*."""
        return [r for r in self.records if r.dest_plate_name == plate_name]

    def get_outgoing(self, plate_name: str) -> List[TransferRecord]:
        """Return records where material leaves *plate_name*."""
        return [r for r in self.records if r.source_plate_name == plate_name]

    def count_bb_usage(self) -> dict:
        """Count how many distinct reactions each source-well building
        block participates in.

        Returns a dict keyed by ``(plate_name, well_index)`` whose values
        are the number of unique ``reaction_id`` values that received
        material from that well.
        """
        usage: dict = {}
        for r in self.records:
            if r.action_type in ("add", "dilution") and r.reaction_id is not None:
                key = (r.source_plate_name, r.source_well_index)
                if key not in usage:
                    usage[key] = set()
                usage[key].add(r.reaction_id)
        # Convert sets → counts
        return {k: len(v) for k, v in usage.items()}

    # ------------------------------------------------------------------
    # Efficiency statistics
    # ------------------------------------------------------------------

    def compute_pipette_stats(self) -> dict:
        """Compute pipette operation statistics from the ledger.

        Groups consecutive transfer records into *tip cycles* — sequences
        of transfers that share a single tip pick-up / drop cycle on the
        robot — then classifies each cycle as either a *tip-change*
        operation or a *no-tip-change* operation.

        Tip-cycle grouping rules (mirror the handler logic):

        * **Multichannel (MC)** – each physical MC operation picks up
          a tip, does one ``transfer_fluid_multi`` call that moves
          liquid through all 8 channels simultaneously, then drops the
          tip.  In the ledger this appears as a contiguous block of
          8 MC records (one per channel / well) that share the same
          ``source_plate_name``, ``dest_plate_name`` and ``volume``.
          Multiple back-to-back MC operations with the same plate/volume
          are merged into one group; the number of physical MC
          operations is ``ceil(len(group) / 8)``.
        * **Dilution** – consecutive ``action_type="dilution"`` records
          with the same ``reaction_class`` share a tip.  Each record is a
          separate physical aspirate/dispense within that tip cycle.
        * **All other SC transfers** (add, extract, analysis) – each
          record is its own 1-record tip cycle (one tip per transfer).

        Time estimates
        ~~~~~~~~~~~~~~
        * Tip-change operation  : **30 s** (pick-up → aspirate → dispense
          → drop)
        * No-tip-change operation: **10 s** (aspirate → dispense, tip
          already held)

        Returns
        -------
        dict
            Keys: ``total_transfers``, ``tip_change_ops``,
            ``no_tip_change_ops``, ``mc_operations``, ``sc_operations``,
            ``tip_change_time_s``, ``no_tip_change_time_s``,
            ``estimated_time_s``.
        """
        TIP_CHANGE_TIME_S = 30.0
        NO_TIP_CHANGE_TIME_S = 10.0

        empty = {
            "total_transfers": 0,
            "tip_change_ops": 0,
            "no_tip_change_ops": 0,
            "mc_operations": 0,
            "sc_operations": 0,
            "tip_change_time_s": 0.0,
            "no_tip_change_time_s": 0.0,
            "estimated_time_s": 0.0,
        }
        if not self.records:
            return empty

        # ---- build tip-cycle groups ----
        tip_groups: List[List[TransferRecord]] = []
        current_group: List[TransferRecord] = [self.records[0]]

        for rec in self.records[1:]:
            prev = current_group[-1]
            same_group = False

            # Consecutive MC records for the same src/dest plate & volume
            if (
                rec.transfer_mode == "multichannel"
                and prev.transfer_mode == "multichannel"
                and rec.source_plate_name == prev.source_plate_name
                and rec.dest_plate_name == prev.dest_plate_name
                and abs(rec.volume - prev.volume) < 0.01
            ):
                same_group = True

            # Consecutive dilution records within the same reaction-class
            # group (one tip per class group in the handler)
            elif (
                rec.action_type == "dilution"
                and prev.action_type == "dilution"
                and (rec.reaction_class or "") == (prev.reaction_class or "")
            ):
                same_group = True

            if same_group:
                current_group.append(rec)
            else:
                tip_groups.append(current_group)
                current_group = [rec]

        tip_groups.append(current_group)

        # ---- tally operations ----
        tip_change_ops = 0
        no_tip_change_ops = 0
        mc_operations = 0
        sc_operations = 0

        # Number of channels on the multichannel pipette.  Each physical
        # MC operation records exactly this many ledger entries.
        _MC_CHANNELS = 8

        for group in tip_groups:
            is_mc = group[0].transfer_mode == "multichannel"

            if is_mc:
                # An MC "group" may contain multiple physical MC
                # operations that were merged because they share the
                # same src plate, dest plate, and volume.  Each
                # physical operation uses one pick-up / transfer / drop
                # cycle and records exactly _MC_CHANNELS ledger entries
                # (one per pipette channel).
                ops_in_group = math.ceil(len(group) / _MC_CHANNELS)
                tip_change_ops += ops_in_group
                mc_operations += ops_in_group
            else:
                # SC group: 1 tip change + (N-1) no-tip-change operations
                tip_change_ops += 1
                no_tip_change_ops += len(group) - 1
                sc_operations += len(group)

        tip_change_time = tip_change_ops * TIP_CHANGE_TIME_S
        no_tip_change_time = no_tip_change_ops * NO_TIP_CHANGE_TIME_S

        return {
            "total_transfers": len(self.records),
            "tip_change_ops": tip_change_ops,
            "no_tip_change_ops": no_tip_change_ops,
            "mc_operations": mc_operations,
            "sc_operations": sc_operations,
            "tip_change_time_s": tip_change_time,
            "no_tip_change_time_s": no_tip_change_time,
            "estimated_time_s": tip_change_time + no_tip_change_time,
        }

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self):
        return iter(self.records)
