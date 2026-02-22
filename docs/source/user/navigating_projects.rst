===================
Navigating Projects
===================

CAR organises chemistry into a hierarchy of objects.  Understanding this
hierarchy is key to navigating the interface and the API.

Data Hierarchy
--------------

::

    Project
    └── Batch (batch-tag from CSV)
        └── Target (one SMILES)
            └── Method (retrosynthetic route, N steps)
                └── Reaction (one step)
                    ├── Product
                    ├── Reactant → CatalogEntry[]
                    ├── ActionSession[]
                    │   ├── AddAction
                    │   ├── StirAction
                    │   ├── ExtractAction
                    │   └── MixAction
                    └── Recipe (versioned)

Project List
------------

The landing page shows all projects sorted by date.  Click a project to
expand it.

Batch View
^^^^^^^^^^

Each project contains one or more **batches**, grouped by the
``batch-tag`` column from the upload CSV.

.. list-table::
   :header-rows: 1

   * - Feature
     - Description
   * - **Clone / Sub-batch**
     - Select specific methods and create a new child batch
   * - **Canonicalize SMILES**
     - Validate and canonicalize a CSV of SMILES strings
   * - **Mark Reactions**
     - Flag reactions as unsuccessful

Target & Method
^^^^^^^^^^^^^^^

Click a target to see its retrosynthetic methods.  Each method shows
the number of reaction steps (``nosteps``) and whether it is compatible
with OpenTrons chemistry (``otchem`` flag).

Reaction Detail Dialog
^^^^^^^^^^^^^^^^^^^^^^

Click a reaction to open the detail dialog.  This shows:

- Reaction SMILES and image
- Reaction class and recipe
- Reactants with vendor catalogue entries and pricing
- Products
- Action sessions with individual add / stir / extract / mix actions

Filtering
---------

Most list views support filtering via query parameters:

.. list-table::
   :header-rows: 1

   * - Resource
     - Filter
   * - Batches
     - ``project_id``
   * - Targets
     - ``batch_id``
   * - Methods
     - ``target_id``, ``nosteps``
   * - Reactions
     - ``method_id``
   * - Products / Reactants
     - ``reaction_id``
   * - Action Sessions
     - ``reaction_id``

``fetchall`` Mode
^^^^^^^^^^^^^^^^^

Appending ``?fetchall=yes`` to any Project, Batch, Target, Method, or
Reaction endpoint returns a deeply nested response that includes all
child objects.  This is convenient for one-shot data retrieval but can
be slow for large projects.
