==========================
Uploading Chemistry to CAR
==========================

CAR accepts target molecules via CSV upload.  Each row defines a target
compound; CAR will plan retrosynthetic routes, identify building blocks,
and prepare everything needed for automated synthesis.

CSV Format
----------

The upload CSV must contain the following columns:

.. list-table::
   :header-rows: 1

   * - Column
     - Required
     - Description
   * - ``targets``
     - Yes
     - SMILES string of the target molecule
   * - ``amount-required-mg``
     - Yes
     - Amount of target compound required (mg)
   * - ``batch-tag``
     - Yes
     - Tag to group targets into batches (e.g. ``"batch-1"``, ``"multi-step"``)

Example CSV::

    targets,amount-required-mg,batch-tag
    O=C(N1CCCCC1)c1ccnn1C,10,single-step
    CNCC(=O)Nc1cc(CNCCC(NC(=O)C2CCCC2)C(=O)O)ccn1,10,multi-step

Upload Process
--------------

1. Navigate to the CAR web interface at ``http://localhost:3000/``
2. Click **Upload** or **Create Project**
3. Provide a **project name** and select your CSV file
4. Click **Submit** — a Celery background task will process the upload

CAR will:

- Parse and canonicalise the SMILES
- Group targets by ``batch-tag`` into separate batches
- Query retrosynthesis services (Manifold) for synthetic routes
- Extract building blocks with vendor availability (eMolecules, MCule)
- Create ``Project → Batch → Target → Method → Reaction`` records

.. note::

   The Celery worker must be running for uploads to process.  Start it with::

       cd /container/src/
       celery -A CAR worker -l info

Custom Starting Materials
-------------------------

If you have pre-made starting-material plates, you can upload a custom
CSV when creating an OT protocol (see :doc:`ot_protocols`).  The custom
CSV overrides the auto-generated compound ordering sheet for the
specified batch.

SMILES Canonicalisation
-----------------------

CAR canonicalises all uploaded SMILES using RDKit to ensure consistent
matching across the system.  You can also use the **Canonicalize SMILES**
endpoint on the Batch detail page to validate SMILES before upload.
