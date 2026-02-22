"""
Plate Visualizer for OpenTrons protocol scripts.

Generates self-contained interactive HTML platemaps that show what is in
each well and how reagents move between plates during a session.  The
HTML files are packaged alongside the OT scripts so that chemists can
click any well to see its contents, molecule structure (SVG), and the
transfers that touch it.

Colour conventions
------------------
* **Starting-material plates** – wells on a green-to-red scale reflecting
  how many reactions the building block participates in (green = 1,
  red = many).
* **Reaction / workup / analysis plates** – wells coloured by
  ``(reaction_class, recipe)`` combination so that different chemistries
  are visually distinct.
* **Solvent plates** – uniform light-blue fill.
* **Empty wells** – white / light grey.
"""

from __future__ import annotations

import base64
import html as html_mod
import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from django.conf import settings

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# SVG helper – wrapped so a missing RDKit won't crash the module
# ---------------------------------------------------------------------------
_RDKIT_AVAILABLE = False
try:
    from rdkit import Chem
    from rdkit.Chem import Draw

    _RDKIT_AVAILABLE = True
except ImportError:
    logger.warning("RDKit not available – molecule SVGs will be omitted")


def _mol_svg(smiles: str, width: int = 200, height: int = 120) -> str:
    """Return an inline SVG string for *smiles*, or an empty string on failure."""
    if not _RDKIT_AVAILABLE or not smiles:
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return ""
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.SetFontSize(10)
        drawer.SetLineWidth(1)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------

# Distinct hues for reaction-class/recipe combos (up to 12, then recycles)
_PALETTE = [
    "#6baed6",  # blue
    "#fd8d3c",  # orange
    "#74c476",  # green
    "#9e9ac8",  # purple
    "#fb6a4a",  # red-ish
    "#fdd0a2",  # peach
    "#a1d99b",  # light green
    "#bcbddc",  # lavender
    "#fdae6b",  # light orange
    "#c7e9c0",  # pale green
    "#dadaeb",  # pale purple
    "#fcbba1",  # salmon
]

_SOLVENT_COLOUR = "#bde0fe"
_EMPTY_COLOUR = "#f5f5f5"


def _usage_colour(usage_count: int, max_usage: int) -> str:
    """Green-to-red gradient based on how many reactions use this BB.

    ``usage_count == 1`` → pure green, ``usage_count == max_usage`` → pure
    red; intermediate values are linearly interpolated via an
    orange midpoint.
    """
    if max_usage <= 1:
        return "#2ca02c"  # green
    t = (usage_count - 1) / (max_usage - 1)
    # Green → Yellow → Red
    if t <= 0.5:
        # green(44,160,44) → yellow(255,255,0)
        s = t * 2
        r = int(44 + (255 - 44) * s)
        g = int(160 + (255 - 160) * s)
        b = int(44 + (0 - 44) * s)
    else:
        # yellow(255,255,0) → red(214,39,40)
        s = (t - 0.5) * 2
        r = int(255 + (214 - 255) * s)
        g = int(255 + (39 - 255) * s)
        b = int(0 + (40 - 0) * s)
    return f"#{r:02x}{g:02x}{b:02x}"


# ---------------------------------------------------------------------------
# Well-name helper (index → "A1" etc.)
# ---------------------------------------------------------------------------

def _well_name_from_index(index: int, rows: int = 8) -> str:
    """Convert a 0-based well index to a human-readable name like ``A1``.

    Wells are numbered column-major: A1, B1, …, H1, A2, B2, …
    """
    row = index % rows
    col = index // rows + 1
    return f"{chr(65 + row)}{col}"


# ===================================================================
# Main class
# ===================================================================

class PlateVisualizer:
    """Generates interactive HTML plate visualizations for an OT session.

    Parameters
    ----------
    script_generator : ScriptGenerator
        The parent ``ScriptGenerator`` whose ``platequeryset``,
        ``transfer_ledger`` and session metadata are used.
    """

    # Role display order – solvent first, then starting materials, then
    # reaction / analysis plates.
    _ROLE_PRIORITY = {
        "solvent": 0,
        "startingmaterial": 1,
        "reaction": 2,
        "workup": 3,
        "spefilter": 4,
        "analyse": 5,
        "lcms": 6,
        "xchem": 7,
        "nmr": 8,
    }

    def __init__(self, script_generator):
        self.sg = script_generator
        self.plates = sorted(
            script_generator.platequeryset,
            key=lambda p: self._ROLE_PRIORITY.get(getattr(p, "role", "") or "", 99),
        )
        self.ledger = script_generator.transfer_ledger
        self.session_type = script_generator.otsessiontype
        self.protocol_name = script_generator.protocolname

        # Pre-compute data structures used by the renderer
        self._bb_usage = self.ledger.count_bb_usage()
        self._max_bb_usage = max(self._bb_usage.values()) if self._bb_usage else 1
        self._reaction_colour_map = self._build_reaction_colour_map()

        # Gather well data keyed by plate id
        self._plate_wells: Dict[int, list] = {}
        for plate in self.plates:
            from backend.models import Well
            wells = list(
                Well.objects.filter(plate_id=plate).order_by("index")
            )
            self._plate_wells[plate.id] = wells

    # ------------------------------------------------------------------
    # Colour map for reaction plates
    # ------------------------------------------------------------------

    def _build_reaction_colour_map(self) -> Dict[Tuple[str, str], str]:
        """Assign a distinct colour to each unique ``(reaction_class, recipe)``
        pair found in the transfer ledger."""
        combos: List[Tuple[str, str]] = []
        seen = set()
        for rec in self.ledger:
            key = (rec.reaction_class or "", rec.recipe or "")
            if key not in seen:
                seen.add(key)
                combos.append(key)
        return {
            combo: _PALETTE[i % len(_PALETTE)] for i, combo in enumerate(combos)
        }

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def generate_html(self) -> str:
        """Return a complete self-contained HTML document."""
        plate_sections = "\n".join(
            self._render_plate(plate, idx)
            for idx, plate in enumerate(self.plates)
        )
        transfer_table = self._render_transfer_table()
        legend = self._render_legend()
        stats_section = self._render_stats()

        # Build <option> list for the plate selector
        plate_options = "\n".join(
            f'<option value="{i}">{html_mod.escape(p.name)} '
            f'({p.role or "unknown"})</option>'
            for i, p in enumerate(self.plates)
        )

        return _HTML_TEMPLATE.format(
            title=html_mod.escape(self.protocol_name),
            session_type=html_mod.escape(self.session_type),
            legend=legend,
            plates=plate_sections,
            plate_options=plate_options,
            plate_count=len(self.plates),
            transfer_table=transfer_table,
            record_count=len(self.ledger),
            stats_section=stats_section,
        )

    def write_html(self, directory: Optional[str] = None) -> str:
        """Write the visualization to an HTML file.

        Parameters
        ----------
        directory : str, optional
            Target directory.  Defaults to ``MEDIA_ROOT/plate_visualizations/``.

        Returns
        -------
        str
            Absolute path to the written file.
        """
        if directory is None:
            directory = os.path.join(settings.MEDIA_ROOT, "plate_visualizations")
        os.makedirs(directory, exist_ok=True)

        filename = f"{self.protocol_name}-plates.html"
        filepath = os.path.join(directory, filename)

        content = self.generate_html()
        with open(filepath, "w", encoding="utf-8") as fh:
            fh.write(content)

        logger.info(f"Plate visualization written to {filepath}")
        return filepath

    # ------------------------------------------------------------------
    # Plate rendering
    # ------------------------------------------------------------------

    def _render_plate(self, plate, plate_idx: int = 0) -> str:
        """Build the HTML for a single plate grid."""
        wells = self._plate_wells.get(plate.id, [])
        rows = plate.numberwellsincolumn or 8
        cols = plate.numbercolumns or 12

        # Build lookup: well_index → well object
        well_map = {w.index: w for w in wells}

        # Build lookup: (plate_name, well_index) → list of transfer summaries
        incoming = defaultdict(list)
        outgoing = defaultdict(list)
        for rec in self.ledger:
            if rec.dest_plate_name == plate.name:
                incoming[rec.dest_well_index].append(rec)
            if rec.source_plate_name == plate.name:
                outgoing[rec.source_well_index].append(rec)

        cells_html = []
        for idx in range(rows * cols):
            well = well_map.get(idx)
            well_label = _well_name_from_index(idx, rows)
            colour = self._well_colour(plate, well, idx)

            # Build tooltip data
            smiles = getattr(well, "smiles", None) if well else None
            volume = getattr(well, "volume", None) if well else None
            solvent = getattr(well, "solvent", None) if well else None
            conc = getattr(well, "concentration", None) if well else None
            db_name = getattr(well, "name", None) if well else None

            # SVG for molecule
            svg_b64 = ""
            if smiles:
                raw_svg = _mol_svg(smiles)
                if raw_svg:
                    svg_b64 = base64.b64encode(raw_svg.encode("utf-8")).decode("ascii")

            # Transfer summary strings
            in_summaries = self._transfer_summaries(incoming.get(idx, []), "in")
            out_summaries = self._transfer_summaries(outgoing.get(idx, []), "out")
            transfer_html = in_summaries + out_summaries

            data_attrs = (
                f'data-well-label="{html_mod.escape(well_label)}"'
                f' data-well-index="{idx}"'
                f' data-db-name="{html_mod.escape(db_name or well_label)}"'
                f' data-smiles="{html_mod.escape(smiles or "")}"'
                f' data-volume="{volume if volume is not None else ""}"'
                f' data-solvent="{html_mod.escape(solvent or "")}"'
                f' data-concentration="{conc if conc is not None else ""}"'
                f' data-svg="{svg_b64}"'
                f' data-transfers="{html_mod.escape(transfer_html)}"'
            )

            # Multichannel indicator — a well is MC if:
            #   1. The Well model has transfer_type="multichannel"
            #      (set on starting-material wells during plate creation), OR
            #   2. Any incoming/outgoing transfer record has
            #      transfer_mode="multichannel" (reaction plate destinations).
            transfer_type = (
                getattr(well, "transfer_type", None) if well else None
            )
            has_mc_transfer = any(
                getattr(r, "transfer_mode", "single") == "multichannel"
                for r in incoming.get(idx, []) + outgoing.get(idx, [])
            )
            is_mc = transfer_type == "multichannel" or has_mc_transfer
            mc_cls = " mc-well" if is_mc else ""

            occupied_cls = " occupied" if (well and smiles) else ""
            mc_mode = "multichannel" if is_mc else (transfer_type or "")
            cells_html.append(
                f'<div class="well{occupied_cls}{mc_cls}" style="background:{colour};" '
                f'data-transfer-type="{mc_mode}" '
                f'{data_attrs} onclick="showWell(this)">'
                f'<span class="well-label">{well_label}</span>'
                f'{"<span class=mc-badge>MC</span>" if is_mc else ""}'
                f"</div>"
            )

        # Column headers
        col_headers = "".join(
            f'<div class="col-header">{c + 1}</div>' for c in range(cols)
        )
        # Row headers
        row_labels = [chr(65 + r) for r in range(rows)]

        grid_with_headers = self._build_grid_html(
            cells_html, rows, cols, row_labels, col_headers
        )

        role_display = f"{plate.role or 'unknown'}"
        if hasattr(plate, "role_index") and plate.role_index and plate.role_index > 1:
            role_display += f" #{plate.role_index}"

        hidden_style = ' style="display:none;"' if plate_idx > 0 else ''
        return (
            f'<div class="plate-section" data-plate-index="{plate_idx}"{hidden_style}>'
            f'<h2>{html_mod.escape(plate.name)} '
            f'<span class="plate-role">({role_display})</span></h2>'
            f'<div class="plate-info">Labware: {html_mod.escape(plate.labware)} '
            f"| {plate.numberwells} wells | Deck slot {plate.index}</div>"
            f"{grid_with_headers}"
            f"</div>"
        )

    def _build_grid_html(
        self, cells: List[str], rows: int, cols: int,
        row_labels: List[str], col_headers: str,
    ) -> str:
        """Assemble the CSS-grid HTML with row/column headers."""
        # Top-left corner blank + column headers
        header_row = '<div class="corner"></div>' + col_headers

        body = ""
        for r in range(rows):
            body += f'<div class="row-header">{row_labels[r]}</div>'
            for c in range(cols):
                idx = c * rows + r  # column-major: A1,B1,...,H1,A2,...
                body += cells[idx]

        return (
            f'<div class="plate-grid" style="'
            f"grid-template-columns: 30px repeat({cols}, 1fr); "
            f'grid-template-rows: 24px repeat({rows}, 1fr);">'
            f"{header_row}{body}</div>"
        )

    # ------------------------------------------------------------------
    # Well colouring
    # ------------------------------------------------------------------

    def _well_colour(self, plate, well, well_index: int) -> str:
        """Return the CSS background colour for a well."""
        role = getattr(plate, "role", None) or ""
        smiles = getattr(well, "smiles", None) if well else None

        if not well or not smiles:
            return _EMPTY_COLOUR

        if role == "startingmaterial":
            usage = self._bb_usage.get((plate.name, well_index), 0)
            if usage == 0:
                return _EMPTY_COLOUR
            return _usage_colour(usage, self._max_bb_usage)

        if role == "solvent":
            return _SOLVENT_COLOUR

        # Reaction / workup / analysis – colour by (class, recipe)
        rxn_class = ""
        recipe = ""
        # Try to determine from incoming transfers
        for rec in self.ledger:
            if rec.dest_plate_name == plate.name and rec.dest_well_index == well_index:
                rxn_class = rec.reaction_class or ""
                recipe = rec.recipe or ""
                break
        key = (rxn_class, recipe)
        return self._reaction_colour_map.get(key, _PALETTE[0])

    # ------------------------------------------------------------------
    # Transfer summaries
    # ------------------------------------------------------------------

    @staticmethod
    def _transfer_summaries(records, direction: str) -> str:
        """Build HTML list items for transfer records."""
        if not records:
            return ""
        items = []
        label = "From" if direction == "in" else "To"
        for rec in records:
            other_plate = rec.source_plate_name if direction == "in" else rec.dest_plate_name
            other_well = rec.source_well_name if direction == "in" else rec.dest_well_name
            other_index = rec.source_well_index if direction == "in" else rec.dest_well_index
            mat = rec.smiles or rec.solvent or rec.action_type
            mode_tag = (
                ' <span class="mc-tag">MC</span>'
                if getattr(rec, "transfer_mode", "single") == "multichannel"
                else ""
            )
            items.append(
                f"<li><b>{label}</b> {html_mod.escape(other_plate)} "
                f"{html_mod.escape(other_well)} (idx {other_index}): "
                f"{rec.volume:.1f} µL "
                f"<i>{html_mod.escape(str(mat))}</i>{mode_tag}</li>"
            )
        return "\n".join(items)

    # ------------------------------------------------------------------
    # Transfer table
    # ------------------------------------------------------------------

    def _render_transfer_table(self) -> str:
        """Build an HTML summary table of all transfers."""
        if not self.ledger.records:
            return "<p>No transfers recorded.</p>"

        rows_html = []
        for i, rec in enumerate(self.ledger, 1):
            smiles_cell = html_mod.escape(rec.smiles or "")
            if rec.smiles:
                svg = _mol_svg(rec.smiles, 120, 60)
                if svg:
                    b64 = base64.b64encode(svg.encode("utf-8")).decode("ascii")
                    smiles_cell += (
                        f'<br><img src="data:image/svg+xml;base64,{b64}" '
                        f'alt="molecule" class="table-mol">'
                    )

            mode = getattr(rec, "transfer_mode", "single")
            mode_cell = (
                '<span class="mc-tag">MC</span>' if mode == "multichannel"
                else "SC"
            )

            rows_html.append(
                f"<tr>"
                f"<td>{i}</td>"
                f"<td>{html_mod.escape(rec.action_type)}</td>"
                f"<td>{mode_cell}</td>"
                f"<td>{html_mod.escape(rec.source_plate_name)}<br>"
                f"<small>{html_mod.escape(rec.source_well_name)} (idx {rec.source_well_index})</small></td>"
                f"<td>{html_mod.escape(rec.dest_plate_name)}<br>"
                f"<small>{html_mod.escape(rec.dest_well_name)} (idx {rec.dest_well_index})</small></td>"
                f"<td>{rec.volume:.1f}</td>"
                f"<td>{smiles_cell}</td>"
                f"<td>{html_mod.escape(rec.solvent or '')}</td>"
                f"<td>{html_mod.escape(rec.reaction_class or '')}</td>"
                f"</tr>"
            )

        return (
            '<table class="transfer-table">'
            "<thead><tr>"
            "<th>#</th><th>Type</th><th>Mode</th><th>Source</th><th>Dest</th>"
            "<th>Vol (µL)</th><th>SMILES / Structure</th><th>Solvent</th><th>Rxn Class</th>"
            "</tr></thead>"
            "<tbody>" + "\n".join(rows_html) + "</tbody></table>"
        )

    # ------------------------------------------------------------------
    # Efficiency statistics
    # ------------------------------------------------------------------

    def _render_stats(self) -> str:
        """Build the HTML for the efficiency-metrics panel."""
        stats = self.ledger.compute_pipette_stats()

        total_transfers = stats["total_transfers"]
        tip_ops = stats["tip_change_ops"]
        no_tip_ops = stats["no_tip_change_ops"]
        mc_ops = stats["mc_operations"]
        sc_ops = stats["sc_operations"]
        tip_time = stats["tip_change_time_s"]
        no_tip_time = stats["no_tip_change_time_s"]
        total_time = stats["estimated_time_s"]

        # Format time as min:sec or just seconds
        def _fmt_time(seconds: float) -> str:
            if seconds >= 60:
                mins = int(seconds // 60)
                secs = seconds % 60
                return f"{mins}m {secs:.0f}s"
            return f"{seconds:.0f}s"

        return (
            '<div class="stats-grid">'
            # Row 1 – transfer counts
            '<div class="stat-card">'
            f'<div class="stat-value">{total_transfers}</div>'
            '<div class="stat-label">Total Transfers</div>'
            "</div>"
            '<div class="stat-card">'
            f'<div class="stat-value">{mc_ops}</div>'
            '<div class="stat-label">MC Operations</div>'
            "</div>"
            '<div class="stat-card">'
            f'<div class="stat-value">{sc_ops}</div>'
            '<div class="stat-label">SC Operations</div>'
            "</div>"
            # Row 2 – tip-change breakdown
            '<div class="stat-card">'
            f'<div class="stat-value">{tip_ops}</div>'
            '<div class="stat-label">Tip Changes</div>'
            f'<div class="stat-detail">{_fmt_time(tip_time)} @ 10 s ea.</div>'
            "</div>"
            '<div class="stat-card">'
            f'<div class="stat-value">{no_tip_ops}</div>'
            '<div class="stat-label">No-Tip-Change Ops</div>'
            f'<div class="stat-detail">{_fmt_time(no_tip_time)} @ 2.5 s ea.</div>'
            "</div>"
            # Row 2, last card – total time
            '<div class="stat-card stat-card-accent">'
            f'<div class="stat-value">{_fmt_time(total_time)}</div>'
            '<div class="stat-label">Est. Run Time</div>'
            "</div>"
            "</div>"
        )

    # ------------------------------------------------------------------
    # Legend
    # ------------------------------------------------------------------

    def _render_legend(self) -> str:
        """Build the colour legend HTML."""
        items = []

        # Starting-material gradient
        items.append(
            '<div class="legend-group">'
            "<h4>Starting Material (BB usage count)</h4>"
            '<div class="gradient-bar">'
            '<span class="gradient" style="background:linear-gradient(to right, #2ca02c, #ffff00, #d6272b);"></span>'
            "<span>1 reaction</span><span>many reactions</span>"
            "</div></div>"
        )

        # Solvent
        items.append(
            '<div class="legend-group">'
            f'<div class="legend-swatch" style="background:{_SOLVENT_COLOUR};"></div>'
            " Solvent</div>"
        )

        # Reaction class / recipe combos
        if self._reaction_colour_map:
            items.append('<div class="legend-group"><h4>Reaction Class / Recipe</h4>')
            for (rxn_cls, recipe), colour in self._reaction_colour_map.items():
                label = f"{rxn_cls}" + (f" – {recipe}" if recipe else "")
                items.append(
                    f'<div><div class="legend-swatch" style="background:{colour};"></div>'
                    f" {html_mod.escape(label)}</div>"
                )
            items.append("</div>")

        # Multichannel indicator
        items.append(
            '<div class="legend-group">'
            '<div class="legend-swatch mc-legend-swatch"></div>'
            ' Multichannel transfer (dashed border + <b>MC</b> badge)</div>'
        )

        # Empty
        items.append(
            '<div class="legend-group">'
            f'<div class="legend-swatch" style="background:{_EMPTY_COLOUR}; '
            f'border:1px solid #ccc;"></div> Empty well</div>'
        )

        return "\n".join(items)


# ===================================================================
# HTML template – everything inlined (CSS, JS) for portability
# ===================================================================

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title} – Session Visualization</title>
<style>
  :root {{
    --bg: #fafafa;
    --card-bg: #fff;
    --text: #333;
    --border: #ddd;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         background: var(--bg); color: var(--text); padding: 20px; }}

  /* Sticky header band – compact */
  .sticky-header {{ position: sticky; top: 0; z-index: 50; background: var(--bg);
                    padding: 4px 0 2px; }}
  h1 {{ margin-bottom: 0; font-size: 1.15rem; }}
  .subtitle {{ color: #666; margin-bottom: 0; font-size: 0.8rem; }}

  /* Plate navigator – inside legend, sticks with it */
  .plate-nav {{ display: flex; align-items: center; gap: 10px; margin-top: 6px;
                padding-top: 6px; border-top: 1px solid var(--border); }}
  .plate-nav select {{ font-size: 0.9rem; padding: 4px 8px; border: 1px solid var(--border);
                       border-radius: 6px; background: var(--card-bg); cursor: pointer;
                       min-width: 180px; }}
  .plate-nav button {{ background: var(--card-bg); border: 1px solid var(--border);
                       border-radius: 6px; padding: 4px 12px; cursor: pointer;
                       font-size: 1rem; line-height: 1; transition: background 0.15s; }}
  .plate-nav button:hover:not(:disabled) {{ background: #e8e8e8; }}
  .plate-nav button:disabled {{ opacity: 0.35; cursor: default; }}
  .plate-nav .plate-counter {{ font-size: 0.82rem; color: #888; }}

  .plate-section {{ background: var(--card-bg); border: 1px solid var(--border);
                    border-radius: 8px; padding: 16px; margin-bottom: 24px; }}
  .plate-section h2 {{ margin-bottom: 4px; font-size: 1.1rem; }}
  .plate-role {{ color: #888; font-weight: normal; font-size: 0.9rem; }}
  .plate-info {{ color: #666; font-size: 0.82rem; margin-bottom: 10px; }}

  /* Grid */
  .plate-grid {{ display: grid; gap: 2px; }}
  .corner {{ }}
  .col-header, .row-header {{ display: flex; align-items: center; justify-content: center;
                              font-size: 0.7rem; color: #888; font-weight: 600; }}
  .well {{ aspect-ratio: 1; border-radius: 50%; border: 1px solid #ccc; display: flex;
           align-items: center; justify-content: center; cursor: pointer;
           transition: transform 0.1s, box-shadow 0.1s; position: relative;
           min-width: 0; min-height: 0; }}
  .well:hover {{ transform: scale(1.25); box-shadow: 0 2px 8px rgba(0,0,0,.25); z-index: 2; }}
  .well-label {{ font-size: 0.55rem; color: rgba(0,0,0,.45); pointer-events: none;
                 user-select: none; }}
  .well.occupied .well-label {{ color: rgba(0,0,0,.7); font-weight: 600; }}

  /* Multichannel well indicator */
  .well.mc-well {{ border: 2px dashed #e65100; }}
  .mc-badge {{ position: absolute; bottom: -2px; right: -2px; font-size: 0.4rem;
               background: #e65100; color: #fff; border-radius: 3px;
               padding: 0 2px; line-height: 1.3; font-weight: 700;
               pointer-events: none; }}
  .mc-tag {{ display: inline-block; background: #e65100; color: #fff; font-size: 0.7rem;
             border-radius: 3px; padding: 0 4px; margin-left: 4px; font-weight: 600;
             vertical-align: middle; }}
  .mc-legend-swatch {{ display: inline-block; width: 13px; height: 13px;
                       border-radius: 50%; vertical-align: middle; margin-right: 4px;
                       border: 2px dashed #e65100; background: #fff; }}

  /* Detail panel */
  #detail-panel {{ position: fixed; top: 20px; right: 20px; width: 340px;
                   background: var(--card-bg); border: 1px solid var(--border);
                   border-radius: 10px; padding: 18px; box-shadow: 0 4px 20px rgba(0,0,0,.12);
                   display: none; z-index: 100; max-height: 90vh; overflow-y: auto; }}
  #detail-panel h3 {{ margin-bottom: 8px; }}
  #detail-panel .field {{ margin-bottom: 6px; font-size: 0.88rem; }}
  #detail-panel .field b {{ display: inline-block; width: 100px; }}
  #detail-panel .mol-img {{ max-width: 100%; margin: 8px 0; }}
  #detail-panel .close-btn {{ position: absolute; top: 8px; right: 12px; cursor: pointer;
                               font-size: 1.2rem; color: #888; background: none; border: none; }}
  #detail-panel .close-btn:hover {{ color: #333; }}
  #detail-panel ul {{ font-size: 0.82rem; padding-left: 18px; margin-top: 4px; }}

  /* Legend – compact, sticky below header */
  .legend {{ position: sticky; top: 34px; z-index: 40;
             background: var(--card-bg); border: 1px solid var(--border);
             border-radius: 6px; padding: 8px 12px; margin-bottom: 8px; }}
  .legend h3 {{ margin-bottom: 4px; font-size: 0.9rem; }}
  .legend-group {{ margin-bottom: 4px; }}
  .legend-group h4 {{ font-size: 0.8rem; margin-bottom: 2px; }}
  .legend-swatch {{ display: inline-block; width: 13px; height: 13px;
                    border-radius: 50%; vertical-align: middle; margin-right: 4px; }}
  .gradient-bar {{ display: flex; align-items: center; gap: 6px; font-size: 0.74rem; }}
  .gradient {{ display: inline-block; width: 90px; height: 12px; border-radius: 6px; }}

  /* Transfer table */
  .transfer-section {{ background: var(--card-bg); border: 1px solid var(--border);
                       border-radius: 8px; padding: 16px; margin-top: 24px; }}
  .transfer-table {{ width: 100%; border-collapse: collapse; font-size: 0.82rem; }}
  .transfer-table th {{ background: #f0f0f0; text-align: left; padding: 6px 8px;
                        border-bottom: 2px solid var(--border);
                        position: sticky; top: 0; z-index: 5; }}
  .transfer-table td {{ padding: 6px 8px; border-bottom: 1px solid #eee; vertical-align: top; }}
  .transfer-table tr:hover {{ background: #fafafa; }}
  .table-mol {{ height: 40px; }}

  /* Stats panel */
  .stats-section {{ background: var(--card-bg); border: 1px solid var(--border);
                    border-radius: 8px; padding: 16px; margin-bottom: 16px; }}
  .stats-section h3 {{ margin-bottom: 10px; font-size: 1rem; }}
  .stats-grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; }}
  .stat-card {{ background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 8px;
                padding: 12px; text-align: center; }}
  .stat-card-accent {{ background: #e8f5e9; border-color: #a5d6a7; }}
  .stat-value {{ font-size: 1.5rem; font-weight: 700; color: var(--text); }}
  .stat-card-accent .stat-value {{ color: #2e7d32; }}
  .stat-label {{ font-size: 0.78rem; color: #666; margin-top: 2px; }}
  .stat-detail {{ font-size: 0.72rem; color: #999; margin-top: 2px; }}
</style>
</head>
<body>
<div class="sticky-header">
<h1>{title}</h1>
<p class="subtitle">Session type: <b>{session_type}</b></p>
</div>

<div class="stats-section">
<h3>Efficiency Metrics</h3>
{stats_section}
</div>

<div class="legend">
<h3>Legend</h3>
{legend}
<div class="plate-nav">
  <button id="prev-plate" onclick="navigatePlate(-1)" disabled>&larr;</button>
  <select id="plate-select" onchange="showPlate(parseInt(this.value))">
    {plate_options}
  </select>
  <button id="next-plate" onclick="navigatePlate(1)">&rarr;</button>
  <span class="plate-counter" id="plate-counter">1 / {plate_count}</span>
</div>
</div>

{plates}

<div class="transfer-section">
<h3>Transfer Log ({record_count} transfers)</h3>
{transfer_table}
</div>

<!-- detail panel -->
<div id="detail-panel">
  <button class="close-btn" onclick="document.getElementById('detail-panel').style.display='none'">&times;</button>
  <h3 id="dp-title"></h3>
  <div id="dp-body"></div>
</div>

<script>
function showWell(el) {{
  const panel = document.getElementById('detail-panel');
  const title = document.getElementById('dp-title');
  const body  = document.getElementById('dp-body');

  const label = el.dataset.wellLabel;
  const idx   = el.dataset.wellIndex;
  const dbName = el.dataset.dbName;
  const smiles = el.dataset.smiles;
  const volume = el.dataset.volume;
  const solvent = el.dataset.solvent;
  const conc   = el.dataset.concentration;
  const svgB64 = el.dataset.svg;
  const transfers = el.dataset.transfers;
  const transferType = el.dataset.transferType;

  title.textContent = 'Well ' + label;

  let html = '';
  html += '<div class="field"><b>Well name:</b> ' + dbName + '</div>';
  html += '<div class="field"><b>Well index:</b> ' + idx + '</div>';
  if (transferType === 'multichannel') {{
    html += '<div class="field"><b>Mode:</b> <span style="background:#e65100;color:#fff;border-radius:3px;padding:1px 6px;font-size:0.82rem;font-weight:600;">Multichannel</span></div>';
  }}

  if (smiles) {{
    html += '<div class="field"><b>SMILES:</b> ' + smiles + '</div>';
  }}
  if (volume) {{
    html += '<div class="field"><b>Volume:</b> ' + volume + ' µL</div>';
  }}
  if (solvent) {{
    html += '<div class="field"><b>Solvent:</b> ' + solvent + '</div>';
  }}
  if (conc) {{
    html += '<div class="field"><b>Conc:</b> ' + conc + '</div>';
  }}
  if (svgB64) {{
    html += '<img class="mol-img" src="data:image/svg+xml;base64,' + svgB64 + '" alt="molecule">';
  }}
  if (transfers) {{
    html += '<h4 style="margin-top:10px;">Transfers</h4><ul>' + transfers + '</ul>';
  }}

  body.innerHTML = html;
  panel.style.display = 'block';
}}

/* ---------- plate navigator ---------- */
var currentPlate = 0;
var plateCount = {plate_count};

function showPlate(idx) {{
  var sections = document.querySelectorAll('.plate-section');
  sections.forEach(function(s) {{ s.style.display = 'none'; }});
  var target = document.querySelector('.plate-section[data-plate-index="' + idx + '"]');
  if (target) target.style.display = '';
  currentPlate = idx;
  document.getElementById('plate-select').value = idx;
  document.getElementById('plate-counter').textContent = (idx + 1) + ' / ' + plateCount;
  document.getElementById('prev-plate').disabled = (idx === 0);
  document.getElementById('next-plate').disabled = (idx === plateCount - 1);
}}

function navigatePlate(delta) {{
  var next = currentPlate + delta;
  if (next >= 0 && next < plateCount) showPlate(next);
}}

document.addEventListener('keydown', function(e) {{
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'TEXTAREA') return;
  if (e.key === 'ArrowLeft')  navigatePlate(-1);
  if (e.key === 'ArrowRight') navigatePlate(1);
}});
</script>
</body>
</html>
"""
