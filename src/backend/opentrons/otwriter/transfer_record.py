"""
Transfer Record for OpenTrons protocol scripts.

Lightweight dataclass that captures every liquid-handling event during
script generation.  The ``ScriptGenerator`` maintains a list of these
records (the *transfer ledger*) which downstream consumers — such as the
``PlateVisualizer`` — read to build visualizations, audit trails, or
validation reports.
"""

from __future__ import annotations

import logging
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
    # Query helpers (used by PlateVisualizer)
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

    def __len__(self) -> int:
        return len(self.records)

    def __iter__(self):
        return iter(self.records)
