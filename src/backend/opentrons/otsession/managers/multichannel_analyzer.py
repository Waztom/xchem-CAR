"""
Multichannel compatibility analyzer for OpenTrons plate layouts.

Analyzes starter plates and materials to determine which transfers can
use multi-channel pipetting (full columns of the same reagent at the
same volume) vs single-channel cherry-picks.

This module is pure analysis — it does not modify the database.
"""

import logging
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


@dataclass
class WellInfo:
    """Represents a well's contents for multichannel analysis.

    Attributes
    ----------
    index : int
        0-based well index on the plate.
    smiles : str or None
        Chemical SMILES string in this well.
    volume : float
        Volume in µL.
    concentration : float or None
        Concentration in mM.
    solvent : str or None
        Solvent name.
    """

    index: int
    smiles: Optional[str] = None
    volume: float = 0.0
    concentration: Optional[float] = None
    solvent: Optional[str] = None


@dataclass
class MultichannelGroup:
    """A group of wells forming a full column eligible for multi-channel transfer.

    Attributes
    ----------
    column_index : int
        0-based column index on the plate.
    well_indices : list of int
        The well indices that make up this column.
    smiles : str
        The shared SMILES for all wells in the column.
    volume : float
        The shared transfer volume.
    concentration : float or None
        Shared concentration.
    solvent : str or None
        Shared solvent.
    """

    column_index: int
    well_indices: List[int] = field(default_factory=list)
    smiles: str = ""
    volume: float = 0.0
    concentration: Optional[float] = None
    solvent: Optional[str] = None


@dataclass
class CherryPickWell:
    """A well that requires single-channel cherry-pick transfer.

    Attributes
    ----------
    well_index : int
        0-based well index on the plate.
    smiles : str
        Chemical SMILES.
    volume : float
        Transfer volume in µL.
    concentration : float or None
        Concentration.
    solvent : str or None
        Solvent.
    reason : str
        Why this well cannot be multichannel (e.g. 'incomplete_column',
        'mixed_volumes', 'mixed_reagents').
    """

    well_index: int
    smiles: str = ""
    volume: float = 0.0
    concentration: Optional[float] = None
    solvent: Optional[str] = None
    reason: str = ""


@dataclass
class ColumnAnalysis:
    """Analysis result for a single column.

    Attributes
    ----------
    column_index : int
        0-based column index.
    wells : list of WellInfo
        Wells in this column.
    is_multichannel_eligible : bool
        True if all wells share the same identity key and volume.
    reason : str
        Explanation if not eligible.
    identity_key : str or None
        The shared '{smiles}-{solvent}-{concentration}' key if uniform.
    volume : float or None
        The shared volume if uniform.
    """

    column_index: int
    wells: List[WellInfo] = field(default_factory=list)
    is_multichannel_eligible: bool = False
    reason: str = ""
    identity_key: Optional[str] = None
    volume: Optional[float] = None


@dataclass
class PlateAnalysisResult:
    """Complete analysis result for a plate.

    Attributes
    ----------
    multichannel_groups : list of MultichannelGroup
        Columns eligible for multi-channel transfer.
    cherry_pick_wells : list of CherryPickWell
        Wells requiring single-channel transfer.
    column_analyses : list of ColumnAnalysis
        Per-column breakdown.
    total_wells : int
        Total occupied wells.
    multichannel_well_count : int
        Wells covered by multichannel.
    cherry_pick_well_count : int
        Wells requiring cherry-pick.
    efficiency : float
        Fraction of wells that can use multichannel (0.0 – 1.0).
    """

    multichannel_groups: List[MultichannelGroup] = field(default_factory=list)
    cherry_pick_wells: List[CherryPickWell] = field(default_factory=list)
    column_analyses: List[ColumnAnalysis] = field(default_factory=list)
    total_wells: int = 0
    multichannel_well_count: int = 0
    cherry_pick_well_count: int = 0
    efficiency: float = 0.0


@dataclass
class MaterialGroup:
    """A group of materials (reactions) sharing the same reagent and volume.

    Used by the prioritisation algorithm to decide which reagents get
    full-column placement on the starter plate.

    Attributes
    ----------
    identity_key : str
        '{smiles}-{solvent}-{concentration}' identifier.
    smiles : str
        Chemical SMILES.
    volume : float
        Per-well transfer volume in µL.
    concentration : float or None
        Concentration in mM.
    solvent : str or None
        Solvent name.
    reaction_count : int
        Number of reactions that use this material at this volume.
    columns_needed : int
        Number of full columns needed (ceil(reaction_count / wells_per_column)).
    leftover_wells : int
        Reactions that don't fill a complete column.
    """

    identity_key: str
    smiles: str
    volume: float
    concentration: Optional[float] = None
    solvent: Optional[str] = None
    reaction_count: int = 0
    columns_needed: int = 0
    leftover_wells: int = 0


class MultichannelAnalyzer:
    """Analyzes plates and materials for multi-channel pipette compatibility.

    A multi-channel pipette operates on a full column (all rows) simultaneously.
    For a transfer to be multichannel-eligible:
    1. Every well in the column must contain the same reagent (SMILES + solvent +
       concentration).
    2. Every well in the column must transfer the same volume.
    3. The column must be fully occupied (no empty wells).

    Parameters
    ----------
    wells_per_column : int
        Number of wells in each column (8 for 96-well, 16 for 384-well, 4 for
        24-well). Default 8.
    num_columns : int
        Total number of columns on the plate. Default 12.
    volume_tolerance : float
        Relative tolerance for comparing volumes (default 0.001, i.e. 0.1%).
    """

    def __init__(
        self,
        wells_per_column: int = 8,
        num_columns: int = 12,
        volume_tolerance: float = 0.001,
    ):
        self.wells_per_column = wells_per_column
        self.num_columns = num_columns
        self.volume_tolerance = volume_tolerance

    # ------------------------------------------------------------------
    # Identity key helper
    # ------------------------------------------------------------------

    @staticmethod
    def make_identity_key(
        smiles: str,
        solvent: Optional[str] = None,
        concentration: Optional[float] = None,
    ) -> str:
        """Create a unique identity key for a reagent solution.

        Parameters
        ----------
        smiles : str
        solvent : str or None
        concentration : float or None

        Returns
        -------
        key : str
            Format: '{smiles}-{solvent}-{concentration}'
        """
        return f"{smiles or ''}-{solvent or ''}-{concentration or ''}"

    # ------------------------------------------------------------------
    # Column helpers
    # ------------------------------------------------------------------

    def get_column_index_for_well(self, well_index: int) -> int:
        """Return the 0-based column index that a well belongs to.

        Wells are stored in column-major order: indices 0..(wells_per_column-1)
        are column 0, etc.

        Parameters
        ----------
        well_index : int

        Returns
        -------
        column_index : int
        """
        return well_index // self.wells_per_column

    def get_well_indices_for_column(self, column_index: int) -> List[int]:
        """Return the well indices that make up a given column.

        Parameters
        ----------
        column_index : int

        Returns
        -------
        indices : list of int
        """
        start = column_index * self.wells_per_column
        return list(range(start, start + self.wells_per_column))

    # ------------------------------------------------------------------
    # 1. Analyze an existing plate (custom starter plate CSV or DB wells)
    # ------------------------------------------------------------------

    def analyze_plate(self, wells: List[WellInfo]) -> PlateAnalysisResult:
        """Analyze an existing plate layout for multichannel compatibility.

        Groups wells by column, then checks each column for uniformity of
        reagent identity and transfer volume.

        Parameters
        ----------
        wells : list of WellInfo
            The occupied wells on the plate. Empty wells should not be
            included.

        Returns
        -------
        result : PlateAnalysisResult
        """
        # Build a lookup: well_index -> WellInfo
        well_map: Dict[int, WellInfo] = {w.index: w for w in wells}

        column_analyses = []
        multichannel_groups = []
        cherry_pick_wells = []

        for col_idx in range(self.num_columns):
            col_well_indices = self.get_well_indices_for_column(col_idx)
            col_wells = [
                well_map[idx] for idx in col_well_indices if idx in well_map
            ]

            analysis = self._analyze_column(col_idx, col_well_indices, col_wells)
            column_analyses.append(analysis)

            if analysis.is_multichannel_eligible:
                group = MultichannelGroup(
                    column_index=col_idx,
                    well_indices=[w.index for w in col_wells],
                    smiles=col_wells[0].smiles,
                    volume=col_wells[0].volume,
                    concentration=col_wells[0].concentration,
                    solvent=col_wells[0].solvent,
                )
                multichannel_groups.append(group)
            else:
                # Every occupied well in this column is a cherry-pick
                for w in col_wells:
                    cherry_pick_wells.append(
                        CherryPickWell(
                            well_index=w.index,
                            smiles=w.smiles or "",
                            volume=w.volume,
                            concentration=w.concentration,
                            solvent=w.solvent,
                            reason=analysis.reason,
                        )
                    )

        mc_count = sum(len(g.well_indices) for g in multichannel_groups)
        cp_count = len(cherry_pick_wells)
        total = mc_count + cp_count

        return PlateAnalysisResult(
            multichannel_groups=multichannel_groups,
            cherry_pick_wells=cherry_pick_wells,
            column_analyses=column_analyses,
            total_wells=total,
            multichannel_well_count=mc_count,
            cherry_pick_well_count=cp_count,
            efficiency=mc_count / total if total > 0 else 0.0,
        )

    def _analyze_column(
        self,
        col_idx: int,
        expected_indices: List[int],
        occupied_wells: List[WellInfo],
    ) -> ColumnAnalysis:
        """Analyze a single column for multichannel eligibility.

        Parameters
        ----------
        col_idx : int
            Column index.
        expected_indices : list of int
            All well indices that should be filled for a full column.
        occupied_wells : list of WellInfo
            Wells actually present in this column.

        Returns
        -------
        analysis : ColumnAnalysis
        """
        analysis = ColumnAnalysis(column_index=col_idx, wells=occupied_wells)

        # Column must be fully occupied
        if len(occupied_wells) == 0:
            analysis.reason = "empty_column"
            return analysis

        if len(occupied_wells) < self.wells_per_column:
            analysis.reason = "incomplete_column"
            return analysis

        # All wells must contain a reagent
        if any(w.smiles is None or w.smiles == "" for w in occupied_wells):
            analysis.reason = "empty_well_contents"
            return analysis

        # Check identity uniformity
        first_key = self.make_identity_key(
            occupied_wells[0].smiles,
            occupied_wells[0].solvent,
            occupied_wells[0].concentration,
        )
        for w in occupied_wells[1:]:
            key = self.make_identity_key(w.smiles, w.solvent, w.concentration)
            if key != first_key:
                analysis.reason = "mixed_reagents"
                return analysis

        # Check volume uniformity
        ref_vol = occupied_wells[0].volume
        for w in occupied_wells[1:]:
            if not math.isclose(w.volume, ref_vol, rel_tol=self.volume_tolerance):
                analysis.reason = "mixed_volumes"
                return analysis

        # All checks passed
        analysis.is_multichannel_eligible = True
        analysis.identity_key = first_key
        analysis.volume = ref_vol
        return analysis

    # ------------------------------------------------------------------
    # 2. Prioritise materials for multichannel placement
    # ------------------------------------------------------------------

    def prioritize_materials(
        self,
        materials: List[Dict],
        wells_per_column: Optional[int] = None,
    ) -> Tuple[List[MaterialGroup], List[MaterialGroup]]:
        """Prioritise materials for multichannel column placement.

        Groups materials by identity key (SMILES + solvent + concentration)
        AND volume. Materials with enough reactions to fill at least one
        full column are prioritised for multichannel placement; the rest
        are marked for sequential single-channel filling.

        Priority ordering (highest first):
        1. Number of full columns (most reactions = most time saved)
        2. Total reaction count (tiebreaker — maximise utilisation)

        Parameters
        ----------
        materials : list of dict
            Each dict must have keys: 'smiles', 'volume', 'reaction_count'.
            Optional keys: 'solvent', 'concentration'.
        wells_per_column : int or None
            Override instance default if needed.

        Returns
        -------
        multichannel_materials : list of MaterialGroup
            Materials that can fill at least one full column, sorted by
            priority (descending).
        single_channel_materials : list of MaterialGroup
            Materials that cannot fill a full column.
        """
        wpc = wells_per_column or self.wells_per_column

        # Group by identity key + volume
        groups: Dict[str, MaterialGroup] = {}
        for mat in materials:
            key = self.make_identity_key(
                mat["smiles"],
                mat.get("solvent"),
                mat.get("concentration"),
            )
            # Include volume in the group key since multichannel requires
            # uniform volume across the column
            group_key = f"{key}|{mat['volume']}"

            if group_key not in groups:
                groups[group_key] = MaterialGroup(
                    identity_key=key,
                    smiles=mat["smiles"],
                    volume=mat["volume"],
                    concentration=mat.get("concentration"),
                    solvent=mat.get("solvent"),
                    reaction_count=0,
                )

            groups[group_key].reaction_count += mat["reaction_count"]

        # Calculate columns needed and leftovers
        for g in groups.values():
            g.columns_needed = g.reaction_count // wpc
            g.leftover_wells = g.reaction_count % wpc

        # Split into multichannel-eligible (≥1 full column) vs single-channel
        multichannel = []
        single_channel = []

        for g in groups.values():
            if g.columns_needed >= 1:
                multichannel.append(g)
            else:
                single_channel.append(g)

        # Sort multichannel by priority: most full columns first, then
        # total reaction count as tiebreaker
        multichannel.sort(
            key=lambda g: (g.columns_needed, g.reaction_count), reverse=True
        )

        # Single channel sorted by reaction count descending (for deterministic
        # output and to place higher-use materials first)
        single_channel.sort(key=lambda g: g.reaction_count, reverse=True)

        logger.info(
            f"Multichannel prioritisation: {len(multichannel)} groups eligible "
            f"for multichannel ({sum(g.columns_needed for g in multichannel)} "
            f"full columns), {len(single_channel)} groups for single-channel"
        )

        return multichannel, single_channel

    # ------------------------------------------------------------------
    # 3. Plan starter plate layout for multichannel + single-channel combo
    # ------------------------------------------------------------------

    def plan_plate_layout(
        self,
        multichannel_materials: List[MaterialGroup],
        single_channel_materials: List[MaterialGroup],
        num_columns: Optional[int] = None,
        wells_per_column: Optional[int] = None,
    ) -> Tuple[List[MultichannelGroup], List[CherryPickWell]]:
        """Plan well assignments for a starter plate combining multichannel
        columns with single-channel cherry-pick wells.

        Multichannel materials get full columns first (priority order).
        Single-channel materials (and leftover wells from multichannel
        groups) fill remaining wells sequentially.

        Parameters
        ----------
        multichannel_materials : list of MaterialGroup
            From ``prioritize_materials()``, sorted by priority.
        single_channel_materials : list of MaterialGroup
            Materials that don't fill a full column.
        num_columns : int or None
            Override for number of columns.
        wells_per_column : int or None
            Override for wells per column.

        Returns
        -------
        multichannel_groups : list of MultichannelGroup
            Column assignments for multichannel transfers.
        cherry_pick_wells : list of CherryPickWell
            Well assignments for single-channel transfers.
        """
        nc = num_columns or self.num_columns
        wpc = wells_per_column or self.wells_per_column

        next_column = 0  # Next available column for multichannel
        next_sequential_index = None  # Will be set after multichannel columns

        multichannel_groups = []
        cherry_pick_wells = []

        # --- Phase A: Assign full columns for multichannel materials ---
        leftover_single_channel = []

        for mat in multichannel_materials:
            for _ in range(mat.columns_needed):
                if next_column >= nc:
                    logger.warning(
                        "Ran out of columns for multichannel placement; "
                        "remaining reactions will use single-channel"
                    )
                    # Remaining full columns become leftovers
                    leftover_single_channel.append(
                        MaterialGroup(
                            identity_key=mat.identity_key,
                            smiles=mat.smiles,
                            volume=mat.volume,
                            concentration=mat.concentration,
                            solvent=mat.solvent,
                            reaction_count=wpc,
                        )
                    )
                    continue

                well_indices = self.get_well_indices_for_column(next_column)
                multichannel_groups.append(
                    MultichannelGroup(
                        column_index=next_column,
                        well_indices=well_indices,
                        smiles=mat.smiles,
                        volume=mat.volume,
                        concentration=mat.concentration,
                        solvent=mat.solvent,
                    )
                )
                next_column += 1

            # Leftover reactions from this material that don't fill a column
            if mat.leftover_wells > 0:
                leftover_single_channel.append(
                    MaterialGroup(
                        identity_key=mat.identity_key,
                        smiles=mat.smiles,
                        volume=mat.volume,
                        concentration=mat.concentration,
                        solvent=mat.solvent,
                        reaction_count=mat.leftover_wells,
                    )
                )

        # --- Phase B: Fill remaining wells sequentially for single-channel ---
        # Start sequential filling after the last multichannel column
        next_sequential_index = next_column * wpc

        all_single = leftover_single_channel + list(single_channel_materials)

        total_plate_wells = nc * wpc

        for mat in all_single:
            for _ in range(mat.reaction_count):
                if next_sequential_index >= total_plate_wells:
                    logger.warning(
                        "Plate is full; additional wells would require "
                        "overflow plate (not handled in layout planning)"
                    )
                    break

                cherry_pick_wells.append(
                    CherryPickWell(
                        well_index=next_sequential_index,
                        smiles=mat.smiles,
                        volume=mat.volume,
                        concentration=mat.concentration,
                        solvent=mat.solvent,
                        reason="single_channel",
                    )
                )
                next_sequential_index += 1

        logger.info(
            f"Plate layout planned: {len(multichannel_groups)} multichannel columns, "
            f"{len(cherry_pick_wells)} cherry-pick wells"
        )

        return multichannel_groups, cherry_pick_wells

    # ------------------------------------------------------------------
    # 4. Analyze a custom starter plate CSV for multichannel compatibility
    # ------------------------------------------------------------------

    def analyze_custom_plate_csv(
        self,
        csv_rows: List[Dict],
        wells_per_column: Optional[int] = None,
        num_columns: Optional[int] = None,
    ) -> PlateAnalysisResult:
        """Analyze a parsed custom starter plate CSV for multichannel
        compatibility.

        Parameters
        ----------
        csv_rows : list of dict
            Each dict should have keys matching the CSV columns:
            'well-index', 'SMILES', 'amount-uL'. Optional: 'concentration',
            'solvent'.
        wells_per_column : int or None
            Override for wells per column based on the plate labware.
        num_columns : int or None
            Override for column count.

        Returns
        -------
        result : PlateAnalysisResult
        """
        wpc = wells_per_column or self.wells_per_column
        nc = num_columns or self.num_columns

        # Temporarily override instance values for this analysis
        original_wpc = self.wells_per_column
        original_nc = self.num_columns
        self.wells_per_column = wpc
        self.num_columns = nc

        try:
            well_infos = []
            for row in csv_rows:
                well_infos.append(
                    WellInfo(
                        index=int(row["well-index"]),
                        smiles=row.get("SMILES", ""),
                        volume=float(row.get("amount-uL", 0)),
                        concentration=(
                            float(row["concentration"])
                            if row.get("concentration")
                            else None
                        ),
                        solvent=row.get("solvent"),
                    )
                )

            return self.analyze_plate(well_infos)
        finally:
            self.wells_per_column = original_wpc
            self.num_columns = original_nc
