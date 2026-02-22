====================
Session Visualizer
====================

Every OT session produces an interactive **HTML visualization** that
lets you inspect the deck layout, well contents, and transfer log
*before* you run the physical protocol.

Opening a Visualization
-----------------------

Session visualizations are bundled inside the protocol ``.zip`` at::

    session_visualizations/step_<N>_session_<M>.html

Open the HTML file in any modern browser.  No internet connection is
required — the file is fully self-contained.

Plate Grid Maps
---------------

Each plate on the deck is rendered as a CSS grid matching the physical
format (96-well or 384-well).  Wells are colour-coded by role:

.. list-table::
   :header-rows: 1

   * - Colour
     - Meaning
   * - **Green → Red gradient**
     - Starting-material well; intensity reflects how many reactions use
       that building block (green = few, red = many)
   * - **Per-class colours**
     - Reaction / workup / analysis wells, coloured by
       ``(reaction_class, recipe)`` pair
   * - **Light blue**
     - Solvent wells
   * - **White / grey**
     - Empty or unavailable wells

Click any well to open a popup showing:

- SMILES and 2-D molecule rendering (via RDKit SVG)
- Well name (e.g. ``A03``)
- Volume (µL)
- Solvent and concentration
- Transfer type (single-channel or multichannel)
- Incoming and outgoing transfer summaries

Multichannel Annotations
^^^^^^^^^^^^^^^^^^^^^^^^

Wells involved in multichannel transfers display an **"MC"** badge.
Multichannel groups span 8 contiguous rows within a single sub-column
(see :doc:`/technical/multichannel`).

Transfer Log Table
------------------

Below the plate grids, a numbered table lists every liquid-handling
transfer in execution order:

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - #
     - Transfer index
   * - Source plate / well
     - Aspirate location
   * - Dest plate / well
     - Dispense location
   * - Volume (µL)
     - Transfer volume
   * - Molecule
     - 2-D SVG of the transferred compound

Legend and Stats
----------------

- **Colour legend** — maps each reaction class / recipe combination to
  its grid colour.
- **Tip-change count** — how many tip changes are required.
- **Estimated run time** — approximate wall-clock time for the session.
- **Plate selector dropdown** — switch between multiple plates without
  scrolling.

Interpreting Volumes
--------------------

Starting-material wells show the **initial provisioned volume**
(i.e. the total amount loaded before any transfers).  This is the
volume the compound-order CSV instructs you to dispense.

The provisioned volume includes:

- Raw transfer volume × number of reactions served
- 15 % error margin (×1.15)
- Dead volume for the labware (e.g. 7.56 µL for a Labcyte 384LDV)

After all transfers, any remaining liquid in the well is waste.

.. seealso::

   :doc:`/technical/multichannel` for details on how multichannel vs
   single-channel provisioning differs.
