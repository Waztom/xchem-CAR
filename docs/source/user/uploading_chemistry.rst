==========================
Uploading Chemistry to CAR
==========================




CAR accepts target molecules via CSV upload.  There are three types of upload possible for CAR.

The first type of upload follows a system where a target compound is defined per row, and run through postera; CAR will plan retrosynthetic routes, identify building blocks,
and prepare everything needed for automated synthesis for this type of upload.

The second type of upload follows a system where one target compound is defined per row and uses custom chemistry following a recipe that has been pre-defined in the recipe database (see :doc:`recipe_generator`).

The third type of upload follows a system where a column of reactant 1 SMILES, and a second column of reactant 2 SMILES and CAR will produce the Combinatorial chemistry of all possible combinations of reactants 1 and reactants 2, following a recipe that has been pre-defined in the recipe database (see :doc:`recipe_generator`).

Postera CSV Format
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

Postera Upload Process
--------------

1. Navigate to the CAR web interface at ``http://localhost:3000/``
2. Click **Upload** or **Create Project**
3. Provide a **project name** and select your CSV file
4. Click the **upload** button in validate choice section
5. Select the **Postera** option from the API choice selection
6. Click **Submit** — a Celery background task will process the upload

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


Custom Chemistry CSV Format
----------

The upload CSV must contain the following columns:

.. list-table::
   :header-rows: 1

   * - Column
     - Required
     - Description
   * - ``target-name``
     - Yes
     - Name of the target compound
   * - ``no-steps``
     - Yes
     - Number of steps of the reaction
   * - ``concentration-required-mM``
     - Yes
     - Concentration of target compound desired (mM)
   * - ``amount-required-uL``
     - Yes
     - Volume of target compound required (uL)
   * - ``batch-tag``
     - Yes
     - Tag to group targets into batches (e.g. ``"batch-1"``, ``"multi-step"``)
   * - ``reactant-1-1``
     - Yes
     - SMILES of reactant 1 of ``reaction-recipe-1``
   * - ``reactant-2-1``
     - If the recipe requires 2 reactants, then yes.
     - SMILES of reactant 2 of ``reaction-recipe-1``
   * - ``reaction-product-smiles-1``
     - Yes
     - SMILES of the product formed by the reactants following the ``reaction-recipe-1``
   * - ``reaction-name-1``
     - Yes
     - Name of reaction type defined in the reaction recipe database
   * - ``reaction-recipe-1``
     - Yes
     - Name of specific recipe in the ``reaction-name`` category
   * - ``reaction-groupby-column-1``
     - Yes
     - Whether reactions start on a new column or not when recipe and/or reaction type changes
   
The upload CSV can contain the following columns (where x is for every integer up to the number defined in the ``no-steps`` column):

.. list-table::
   :header-rows: 1

   * - ``reactant-1-x``
     - No
     - SMILES of reactant 1 of ``reaction-recipe-x``
   * - ``reactant-2-x``
     - If the recipe requires 2 reactants, then yes.
     - SMILES of reactant 2 of ``reaction-recipe-x``
   * - ``reaction-product-smiles-x``
     - No
     - SMILES of the product formed by the reactants following the ``reaction-recipe-x``
   * - ``reaction-name-x``
     - No
     - Name of reaction type defined in the reaction recipe database
   * - ``reaction-recipe-x``
     - No
     - Name of specific recipe in the ``reaction-name`` category
   * - ``reaction-groupby-column-x``
     - No
     - Whether reactions start on a new column or not when recipe and/or reaction type changes

Custom Chemistry Upload Process
--------------

1. Navigate to the CAR web interface at ``http://localhost:3000/``
2. Click **Upload** or **Create Project**
3. Provide a **project name** and select your CSV file
4. Click the **upload** button in validate choice section
5. Select the **Custom chemistry** option from the API choice selection
6. Click **Submit** — a Celery background task will process the upload

CAR will:

- Parse and canonicalise the SMILES
- Group targets by ``batch-tag`` into separate batches
- Create ``Project → Batch → Target → Method → Reaction`` records

.. note::

   The Celery worker must be running for uploads to process.  Start it with::

       cd /container/src/
       celery -A CAR worker -l info

Combi custom Chemistry CSV Format
----------

The upload CSV must contain the following columns:

.. list-table::
   :header-rows: 1

   * - Column
     - Required
     - Description
   * - ``combi-group``
     - Yes
     - What group of reactants are in the same Combinatorial chemistry grouping (e.g. reactants defined in ``"group"`` 1, will only mix with other reactants defined in ``"group 1"``, but not ones defined in ``"group 2"``)
   * - ``no-steps``
     - Yes
     - Number of steps of the reaction
   * - ``concentration-required-mM``
     - Yes
     - Concentration of target compound desired (mM)
   * - ``amount-required-uL``
     - Yes
     - Volume of target compound required (uL)
   * - ``batch-tag``
     - Yes
     - Tag to group targets into batches (e.g. ``"batch-1"``, ``"multi-step"``)
   * - ``reactant-1-1``
     - Yes
     - SMILES of reactant 1 of ``reaction-recipe-1``
   * - ``reactant-2-1``
     - If the recipe requires 2 reactants, then yes.
     - SMILES of reactant 2 of ``reaction-recipe-1``
   * - ``reaction-name-1``
     - Yes
     - Name of reaction type defined in the reaction recipe database
   * - ``reaction-recipe-1``
     - Yes
     - Name of specific recipe in the ``reaction-name`` category
   
The upload CSV can contain the following columns (where x is for every integer up to the number defined in the ``no-steps`` column):

.. list-table::
   :header-rows: 1

   * - ``reactant-1-x``
     - No
     - SMILES of reactant 1 of ``reaction-recipe-x``
   * - ``reactant-2-x``
     - If the recipe requires 2 reactants, then yes.
     - SMILES of reactant 2 of ``reaction-recipe-x``
   * - ``reaction-name-x``
     - No
     - Name of reaction type defined in the reaction recipe database
   * - ``reaction-recipe-x``
     - No
     - Name of specific recipe in the ``reaction-name`` category


Combi custom Chemistry Upload Process
--------------

1. Navigate to the CAR web interface at ``http://localhost:3000/``
2. Click **Upload** or **Create Project**
3. Provide a **project name** and select your CSV file
4. Click the **upload** button in validate choice section
5. Select the **Custom combi chemistry** option from the API choice selection
6. Click **Submit** — a Celery background task will process the upload

CAR will:

- Parse and canonicalise the SMILES
- Group targets by ``batch-tag`` into separate batches
- Create all possible Combinatorial chemistry reaction products recipes
- Create the product SMILES based on the SMARTS pattern in the reaction recipe, and the reactants defined
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

The template must contain the following columns:

.. list-table::
   :header-rows: 1

   * - Column
     - Required
     - Description
   * - ``plate-ID``
     - Yes
     - The plate ID of the starting material
   * - ``labware-type``
     - Yes
     - The type of labware the plate is
   * - ``well-index``
     - Yes
     - Either the well index (0,1,2, etc.) or well ID (A1, B1, C1, etc. respectively) of the starting material
   * - ``SMILES``
     - Yes
     - SMILES pattern of the starting material
   * - ``amount-uL``
     - Yes
     - Volume of starting material in the plate well
   * - ``concentration``
     - Yes
     - Concentration of the starting material (M)
   * - ``solvent``
     - Yes
     - Solvent the starting material is dissolved in
SMILES Canonicalisation
-----------------------

CAR canonicalises all uploaded SMILES using RDKit to ensure consistent
matching across the system.  You can also use the **Canonicalize SMILES**
endpoint on the Batch detail page to validate SMILES before upload.
