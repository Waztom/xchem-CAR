========================
Multichannel Pipetting
========================

CAR supports **8-channel (multichannel) pipette optimisation** for
high-throughput synthesis.  When enabled, identical starting materials
are grouped into blocks of 8 wells so the multichannel pipette can
transfer them in a single aspirate-dispense cycle.

.. contents:: On this page
   :local:
   :depth: 2

Enabling Multichannel Mode
--------------------------

Pass ``use_multichannel=True`` when creating an OT protocol:

- **Frontend:** tick the "Use Multichannel" checkbox
- **API:** include ``"use_multichannel": true`` in the ``createotproject``
  request body

384-Well Plate Layout
---------------------

CAR uses **384-well plates** (16 rows × 24 columns) for starting
materials.  The 24 physical columns are split into two interleaved
**sub-columns**:

.. code-block:: text

   Sub-column 0 (odd):  columns 1, 3, 5, 7, …, 23
   Sub-column 1 (even): columns 2, 4, 6, 8, …, 24

Each sub-column contains **8 wells** spaced at **9 mm** — exactly
matching the OT multichannel pipette tip spacing.

A multichannel group is therefore 8 wells in the same sub-column,
all containing the same compound at the same concentration.

MC vs SC Wells
--------------

When multichannel mode is enabled, starting-material wells are split
into two categories:

**Multichannel (MC) wells:**

- Provision volume for **multichannel transfers only**
- Grouped in blocks of 8 within a sub-column
- ``transfer_type = "multichannel"``

**Single-channel (SC) wells:**

- Provision volume for **single-channel transfers only**
- Handle leftover reactions not covered by MC groups
- ``transfer_type = "single"``

This clean separation ensures that:

1. MC wells hold exactly enough volume for MC transfers
2. SC wells hold exactly enough volume for SC transfers
3. The ``well_finder`` queries the correct wells at script-generation
   time by filtering on ``transfer_type``

Volume Calculation
------------------

MC well volume
^^^^^^^^^^^^^^

.. code-block:: text

   transfers_per_channel = reaction_count // wells_per_group
   mc_volume = raw_volume × transfers_per_channel × 1.15 + dead_volume

- ``wells_per_group`` = 8 (multichannel tip count)
- ``1.15`` = 15 % error margin
- ``dead_volume`` = labware-dependent (e.g. 7.56 µL for Labcyte
  384LDV 100 µL plate)
- Floor division (``//``) is used so that leftover reactions go to SC
  wells

SC well volume
^^^^^^^^^^^^^^

.. code-block:: text

   sc_volume = raw_volume × sc_reaction_count × 1.15 + dead_volume

``sc_reaction_count`` is the number of reactions not covered by MC
groups, plus any materials that are inherently single-channel only.

Example
^^^^^^^

20 reactions need NaBH4 (2 µL each), on a plate with 8-channel MC:

.. code-block:: text

   transfers_per_channel = 20 // 8 = 2
   mc_reactions_covered  = 2 × 8 = 16
   sc_leftover           = 20 - 16 = 4

   MC well volume = 2 × 2 × 1.15 + 7.56 = 12.16 µL   (× 8 wells)
   SC well volume = 2 × 4 × 1.15 + 7.56 = 16.76 µL   (× 1 well)

Plate Factory Phases
--------------------

The ``PlateFactory`` builds starting-material plates in two phases:

**Phase A — Multichannel layout:**

1. Group materials by ``(smiles, concentration)``
2. For each material, compute ``transfers_per_channel = n // 8``
3. Create MC wells across sub-columns with transfer-type
   ``"multichannel"``
4. Track uncovered reactions as ``sc_leftover``

**Phase B — Single-channel layout:**

1. Collect all SC-only materials + Phase A leftovers
2. Create individual SC wells with transfer-type ``"single"``
3. Handle volume overflow by splitting across multiple wells

Script Generation
-----------------

During OT script generation, the ``well_finder`` respects the
``transfer_type`` field:

- ``_process_step_multichannel()`` looks up wells with
  ``transfer_type="multichannel"``
- ``_process_step_single_channel()`` looks up wells with
  ``transfer_type="single"``

This prevents cross-contamination between MC and SC volume pools.

.. seealso::

   - :doc:`/user/session_visualiser` — MC wells are annotated with an
     "MC" badge in visualizations
   - :doc:`/user/ot_protocols` — overall protocol creation workflow
