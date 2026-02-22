============================
Creating OT Protocol Bundles
============================

CAR generates ready-to-run `Opentrons <https://opentrons.com/>`_
protocol bundles.  A bundle contains Python protocol scripts, compound
ordering sheets, solvent-preparation guides, and interactive session
visualizations — all packaged in a single **.zip** file.

Requesting a Protocol
---------------------

1. Navigate to your project and open the batch you want to synthesise.
2. Click **Create OT Protocol**.
3. Fill in the form:

   - **Protocol name** — a human-readable label
   - **Batch IDs** — select one or more batches
   - **Use Multichannel** — enable 8-channel pipette optimisation
     (see :doc:`/technical/multichannel`)
   - **Custom Starting-Material CSV** *(optional)* — override the
     compound ordering sheet with pre-plated starting materials

4. Click **Submit** — a Celery task begins processing.
5. Poll the status with **Get Task Status** until ``SUCCESS``.
6. Download the ``.zip`` from the OT Batch Protocol detail page.

.. note::

   Protocol generation can take 30 seconds to several minutes for large
   batches.  The Celery worker must be running.

What Is Inside the .zip
------------------------

::

    <protocol-name>.zip
    ├── compound_orders/
    │   └── step_<N>_compound_order.csv   # One per reaction step
    ├── ot_scripts/
    │   └── step_<N>_session_<M>.py       # Opentrons Python protocols
    ├── solvent_prep/
    │   └── step_<N>_solvent_prep.csv      # Solvent preparation guide
    └── session_visualizations/
        └── step_<N>_session_<M>.html      # Interactive HTML visualizations

Compound Order CSV
^^^^^^^^^^^^^^^^^^

Lists every starting material you need to purchase and plate.  Columns
include SMILES, vendor, catalogue ID, mass/volume required, and the
destination well on the source plate.

OT Script (Python)
^^^^^^^^^^^^^^^^^^

Each ``.py`` file is a standalone Opentrons protocol.  Upload it to the
Opentrons App or run it on a connected robot.  The script:

- Declares the deck layout (plates, tipracks, pipettes)
- Performs aspirate/dispense or multichannel transfer commands
- Includes human-readable comments annotating each transfer

Solvent Prep CSV
^^^^^^^^^^^^^^^^

Lists solvents to pre-load into solvent plates or reservoirs before the
run, with well positions and volumes.

Session Visualization
^^^^^^^^^^^^^^^^^^^^^

Self-contained HTML pages for reviewing the session before running.
See :doc:`session_visualiser` for a full walkthrough.

Protocol-Creation Pipeline (Under the Hood)
--------------------------------------------

.. code-block:: text

   Frontend POST → api/otprojects/createotproject
     → protocol_service.initiate_ot_project()
       → Celery task: create_ot_script()
         → For each batch:
             → For each reaction step:
                 → create_multiple_ot_sessions()
                   → SessionOrchestrator.execute()
                     → PlateFactory (wells, columns, MC/SC layout)
                     → ScriptGenerator (OT Python script)
                     → SessionVisualizer (HTML)
                 → ZipOTBatchProtocol (package .zip)

Sessions are split automatically if more than ~600 reactions would
overflow the deck.

Deck Overflow Handling
^^^^^^^^^^^^^^^^^^^^^^

If the deck cannot accommodate all plates for a group of reactions, CAR
automatically splits the group in half and retries recursively until
each sub-group fits on the 11-slot OT deck.
