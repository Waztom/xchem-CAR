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
    """A group of wells eligible for multi-channel transfer.

    For 96-well plates this is a full column (8 wells).  For 384-well
    plates the 8-channel pipette accesses every other row, so each
    physical column contains two sub-columns of 8 wells each.

    Attributes
    ----------
    column_index : int
        0-based physical column index on the plate.
    sub_column_index : int
        Sub-column within the physical column (0 for 96/24-well;
        0 or 1 for 384-well).
    well_indices : list of int
        The well indices that make up this group.
    smiles : str
        The shared SMILES for all wells in the group.
    volume : float
        The shared transfer volume.
    concentration : float or None
        Shared concentration.
    solvent : str or None
        Shared solvent.
    plate_number : int
        0-based plate index (0 = first plate, 1 = overflow plate, etc.).
    """

    column_index: int
    sub_column_index: int = 0
    well_indices: List[int] = field(default_factory=list)
    smiles: str = ""
    volume: float = 0.0
    concentration: Optional[float] = None
    solvent: Optional[str] = None
    plate_number: int = 0


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
    plate_number : int
        0-based plate index (0 = first plate, 1 = overflow plate, etc.).
    """

    well_index: int
    smiles: str = ""
    volume: float = 0.0
    concentration: Optional[float] = None
    solvent: Optional[str] = None
    reason: str = ""
    plate_number: int = 0


@dataclass
class ColumnAnalysis:
    """Analysis result for a sub-column (or full column on 96/24-well plates).

    On 384-well plates each physical column produces two
    ``ColumnAnalysis`` entries (one per sub-column the 8-channel
    pipette can reach).

    Attributes
    ----------
    column_index : int
        0-based physical column index.
    sub_column_index : int
        Sub-column within the physical column (0 for 96/24-well;
        0 or 1 for 384-well).
    wells : list of WellInfo
        Wells in this sub-column.
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
    sub_column_index: int = 0
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
    plates_needed: int = 1


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
        Number of full MC groups needed (ceil(reaction_count /
        wells_per_group)).  For 384-well plates each physical column
        holds two groups.
    leftover_wells : int
        Reactions that don't fill a complete MC group.
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

    An 8-channel pipette has tips spaced at 9 mm (96-well pitch).  On a
    384-well plate (4.5 mm pitch, 16 rows) the pipette reaches every
    other row, producing two *sub-columns* of 8 wells per physical
    column.  For a transfer to be multichannel-eligible:

    1. Every well in the sub-column must contain the same reagent
       (SMILES + solvent + concentration).
    2. Every well in the sub-column must transfer the same volume.
    3. The sub-column must be fully occupied (no empty wells).

    Parameters
    ----------
    wells_per_column : int
        Physical rows per column (8 for 96-well, 16 for 384-well, 4 for
        24-well).  Default 8.
    num_columns : int
        Total physical columns on the plate.  Default 12.
    volume_tolerance : float
        Relative tolerance for comparing volumes (default 0.001, i.e. 0.1%).
    pipette_channels : int
        Number of channels on the multi-channel pipette.  Default 8.
    """

    def __init__(
        self,
        wells_per_column: int = 8,
        num_columns: int = 12,
        volume_tolerance: float = 0.001,
        pipette_channels: int = 8,
    ):
        self.wells_per_column = wells_per_column
        self.num_columns = num_columns
        self.volume_tolerance = volume_tolerance
        self.pipette_channels = pipette_channels

        #: Wells the pipette can access in one pass within a column.
        #: Equals ``min(wells_per_column, pipette_channels)``.
        self.wells_per_group: int = min(wells_per_column, pipette_channels)

        #: Number of independent sub-columns per physical column.
        #: 1 for 96/24-well, 2 for 384-well with an 8-channel pipette.
        self.sub_columns_per_column: int = max(
            1, wells_per_column // pipette_channels
        )

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

    def get_sub_column_well_indices(
        self, column_index: int, sub_column_index: int = 0
    ) -> List[int]:
        """Return the well indices for one sub-column within a physical column.

        For plates where ``wells_per_column <= pipette_channels`` (e.g. 96-well)
        this returns the same result as :meth:`get_well_indices_for_column`.

        For 384-well plates (``wells_per_column=16, pipette_channels=8``) each
        column is split into two interleaved sub-columns:

        * sub-column 0 → even row indices (A, C, E, G, I, K, M, O)
        * sub-column 1 → odd row indices  (B, D, F, H, J, L, N, P)

        Parameters
        ----------
        column_index : int
            0-based physical column.
        sub_column_index : int
            0-based sub-column within the physical column.

        Returns
        -------
        indices : list of int
        """
        col_start = column_index * self.wells_per_column
        if self.sub_columns_per_column <= 1:
            return list(range(col_start, col_start + self.wells_per_column))
        return [
            col_start + sub_column_index + i * self.sub_columns_per_column
            for i in range(self.wells_per_group)
        ]

    # ------------------------------------------------------------------
    # 1. Analyze an existing plate (custom starter plate CSV or DB wells)
    # ------------------------------------------------------------------

    def analyze_plate(self, wells: List[WellInfo]) -> PlateAnalysisResult:
        """Analyze an existing plate layout for multichannel compatibility.

        Groups wells by column (and sub-column on 384-well plates), then
        checks each group for uniformity of reagent identity and transfer
        volume.

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
            # Track which wells in this column are claimed by an MC group
            # so the remainder become cherry-picks.
            mc_claimed_indices: set = set()

            for sub_col_idx in range(self.sub_columns_per_column):
                sub_col_well_indices = self.get_sub_column_well_indices(
                    col_idx, sub_col_idx
                )
                sub_col_wells = [
                    well_map[idx]
                    for idx in sub_col_well_indices
                    if idx in well_map
                ]

                analysis = self._analyze_column(
                    col_idx, sub_col_well_indices, sub_col_wells,
                    sub_column_index=sub_col_idx,
                )
                column_analyses.append(analysis)

                if analysis.is_multichannel_eligible:
                    group = MultichannelGroup(
                        column_index=col_idx,
                        sub_column_index=sub_col_idx,
                        well_indices=[w.index for w in sub_col_wells],
                        smiles=sub_col_wells[0].smiles,
                        volume=sub_col_wells[0].volume,
                        concentration=sub_col_wells[0].concentration,
                        solvent=sub_col_wells[0].solvent,
                    )
                    multichannel_groups.append(group)
                    mc_claimed_indices.update(w.index for w in sub_col_wells)
                else:
                    for w in sub_col_wells:
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
        sub_column_index: int = 0,
    ) -> ColumnAnalysis:
        """Analyze a single sub-column for multichannel eligibility.

        Parameters
        ----------
        col_idx : int
            Physical column index.
        expected_indices : list of int
            All well indices that should be filled for a full group.
        occupied_wells : list of WellInfo
            Wells actually present in this sub-column.
        sub_column_index : int
            Sub-column within the physical column.

        Returns
        -------
        analysis : ColumnAnalysis
        """
        analysis = ColumnAnalysis(
            column_index=col_idx,
            sub_column_index=sub_column_index,
            wells=occupied_wells,
        )

        # Sub-column must be fully occupied
        if len(occupied_wells) == 0:
            analysis.reason = "empty_column"
            return analysis

        if len(occupied_wells) < self.wells_per_group:
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
            Override instance default if needed.  Ignored when
            ``wells_per_column > pipette_channels``; in that case the
            effective group size (``wells_per_group``) is used instead.

        Returns
        -------
        multichannel_materials : list of MaterialGroup
            Materials that can fill at least one full MC group, sorted by
            priority (descending).
        single_channel_materials : list of MaterialGroup
            Materials that cannot fill a full MC group.
        """
        # The group size is how many wells the pipette can transfer in
        # one pass — equals min(wells_per_column, pipette_channels).
        if wells_per_column is not None:
            wpg = min(wells_per_column, self.pipette_channels)
        else:
            wpg = self.wells_per_group

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

        # Calculate groups needed and leftovers using the effective
        # pipette group size (8 for both 96- and 384-well plates).
        for g in groups.values():
            g.columns_needed = g.reaction_count // wpg
            g.leftover_wells = g.reaction_count % wpg

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
        wpg = min(wpc, self.pipette_channels)
        spc = max(1, wpc // self.pipette_channels)  # sub-columns per column

        # Track placement by (column, sub_column) slots.
        next_column = 0
        next_sub_column = 0  # within current column
        current_plate = 0

        multichannel_groups = []
        cherry_pick_wells = []

        total_mc_slots = nc * spc  # total sub-column slots per plate
        used_mc_slots = 0

        # --- Phase A: Assign sub-column slots for multichannel materials ---
        leftover_single_channel = []

        for mat in multichannel_materials:
            for _ in range(mat.columns_needed):
                if used_mc_slots >= total_mc_slots:
                    # Current plate is full — start a new overflow plate
                    # for the remaining MC groups.
                    current_plate += 1
                    next_column = 0
                    next_sub_column = 0
                    used_mc_slots = 0
                    logger.info(
                        f"MC plate {current_plate - 1} full; creating "
                        f"overflow plate {current_plate} for remaining "
                        f"multichannel groups"
                    )

                well_indices = self.get_sub_column_well_indices(
                    next_column, next_sub_column
                )
                multichannel_groups.append(
                    MultichannelGroup(
                        column_index=next_column,
                        sub_column_index=next_sub_column,
                        well_indices=well_indices,
                        smiles=mat.smiles,
                        volume=mat.volume,
                        concentration=mat.concentration,
                        solvent=mat.solvent,
                        plate_number=current_plate,
                    )
                )
                used_mc_slots += 1
                next_sub_column += 1
                if next_sub_column >= spc:
                    next_sub_column = 0
                    next_column += 1

            # Leftover reactions from this material that don't fill a group
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
        # Derive the first SC column from the actually-placed MC groups
        # on the current (last) plate, rather than relying on counter
        # algebra.  This avoids subtle bugs when next_sub_column resets
        # to 0 after completing all sub-columns of a physical column.
        groups_on_current_plate = [
            g for g in multichannel_groups if g.plate_number == current_plate
        ]
        if groups_on_current_plate:
            first_sc_column = (
                max(g.column_index for g in groups_on_current_plate) + 1
            )
        else:
            first_sc_column = 0
        next_sequential_index = first_sc_column * wpc

        all_single = leftover_single_channel + list(single_channel_materials)

        total_plate_wells = nc * wpc

        for mat in all_single:
            for _ in range(mat.reaction_count):
                if next_sequential_index >= total_plate_wells:
                    # Current plate is full — start a new overflow plate.
                    current_plate += 1
                    next_sequential_index = 0
                    logger.info(
                        f"SC plate {current_plate - 1} full; creating "
                        f"overflow plate {current_plate} for remaining "
                        f"single-channel wells"
                    )

                cherry_pick_wells.append(
                    CherryPickWell(
                        well_index=next_sequential_index,
                        smiles=mat.smiles,
                        volume=mat.volume,
                        concentration=mat.concentration,
                        solvent=mat.solvent,
                        reason="single_channel",
                        plate_number=current_plate,
                    )
                )
                next_sequential_index += 1

        plates_needed = current_plate + 1

        logger.info(
            f"Plate layout planned: {len(multichannel_groups)} multichannel groups, "
            f"{len(cherry_pick_wells)} cherry-pick wells, "
            f"{plates_needed} plate(s) needed"
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
        original_wpg = self.wells_per_group
        original_spc = self.sub_columns_per_column

        self.wells_per_column = wpc
        self.num_columns = nc
        self.wells_per_group = min(wpc, self.pipette_channels)
        self.sub_columns_per_column = max(1, wpc // self.pipette_channels)

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
            self.wells_per_group = original_wpg
            self.sub_columns_per_column = original_spc
