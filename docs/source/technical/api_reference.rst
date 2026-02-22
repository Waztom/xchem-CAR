=============
API Reference
=============

CAR exposes a RESTful API built with Django REST Framework.  All
endpoints return JSON and live under the ``/api/`` prefix.

.. contents:: On this page
   :local:
   :depth: 2

Base URL
--------

::

    http://127.0.0.1:8000/api/

The DRF browsable API is accessible from a web browser at the same URL.

Core Resources
--------------

These endpoints follow standard REST conventions (``GET``, ``POST``,
``PUT``, ``PATCH``, ``DELETE``).

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Endpoint
     - ViewSet
     - Description
   * - ``/api/projects/``
     - ``ProjectViewSet``
     - Synthesis projects
   * - ``/api/batches/``
     - ``BatchViewSet``
     - Compound batches (filter: ``project_id``)
   * - ``/api/targets/``
     - ``TargetViewSet``
     - Target molecules (filter: ``batch_id``)
   * - ``/api/methods/``
     - ``MethodViewSet``
     - Retrosynthetic methods (filter: ``target_id``, ``nosteps``)
   * - ``/api/reactions/``
     - ``ReactionViewSet``
     - Reaction steps (filter: ``method_id``)
   * - ``/api/products/``
     - ``ProductViewSet``
     - Reaction products (filter: ``reaction_id``)
   * - ``/api/reactants/``
     - ``ReactantViewSet``
     - Reactants (filter: ``reaction_id``)
   * - ``/api/catalogentries/``
     - ``CatalogEntryViewSet``
     - Vendor catalog entries
   * - ``/api/pubcheminfo/``
     - ``PubChemInfoViewSet``
     - PubChem compound data

Action Resources
----------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Endpoint
     - ViewSet
     - Description
   * - ``/api/actionsessions/``
     - ``ActionSessionViewSet``
     - Action session groups (filter: ``reaction_id``)
   * - ``/api/addactions/``
     - ``AddActionViewSet``
     - Liquid-handling add steps
   * - ``/api/extractactions/``
     - ``ExtractActionViewSet``
     - Liquid-liquid extraction steps
   * - ``/api/mixactions/``
     - ``MixActionViewSet``
     - Pipette mixing steps
   * - ``/api/stiractions/``
     - ``StirActionViewSet``
     - Stir / incubation steps

OT Resources
-------------

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Endpoint
     - ViewSet
     - Description
   * - ``/api/otprojects/``
     - ``OTProjectViewSet``
     - OT project containers
   * - ``/api/otbatchprotocols/``
     - ``OTBatchProtocolViewSet``
     - Batch protocols (.zip) (filter: ``otproject_id``, ``batch_id``)
   * - ``/api/otsessions/``
     - ``OTSessionViewSet``
     - Individual sessions (filter: ``otbatchprotocol_id``)
   * - ``/api/decks/``
     - ``DeckViewSet``
     - Deck models
   * - ``/api/pipettes/``
     - ``PipetteViewSet``
     - Pipette models
   * - ``/api/tipracks/``
     - ``TipRackViewSet``
     - Tip rack models
   * - ``/api/plates/``
     - ``PlateViewSet``
     - Plate models (filter: ``otbatchprotocol_id``)
   * - ``/api/columns/``
     - ``ColumnViewSet``
     - Plate columns
   * - ``/api/wells/``
     - ``WellViewSet``
     - Individual wells
   * - ``/api/compoundorders/``
     - ``CompoundOrderViewSet``
     - Compound ordering CSVs
   * - ``/api/otscripts/``
     - ``OTScriptViewSet``
     - Generated OT Python scripts

Custom Actions
--------------

createproject
^^^^^^^^^^^^^

.. code-block:: text

   POST /api/projects/createproject/

Create a new project from a CSV upload.

**Request** (``multipart/form-data``):

- ``project_name`` — project name
- ``csv_file`` — CSV with ``targets``, ``amount-required-mg``,
  ``batch-tag`` columns

**Response**: ``{ "taskid": "<celery-task-id>" }``

gettaskstatus (Project)
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   GET /api/projects/gettaskstatus/?taskid=<id>

Poll the Celery task status for project creation.

**Response**: ``{ "state": "PENDING|STARTED|SUCCESS|FAILURE", ... }``

createotproject
^^^^^^^^^^^^^^^

.. code-block:: text

   POST /api/otprojects/createotproject/

Create OT protocols for selected batches.

**Request** (``multipart/form-data``):

- ``protocol_name`` — human-readable label
- ``batch_ids`` — JSON array of batch IDs
- ``use_multichannel`` — boolean (default ``false``)
- ``custom_sm_files`` *(optional)* — custom starting-material CSVs

**Response**: ``{ "taskid": "<celery-task-id>" }``

canonicalizesmiles
^^^^^^^^^^^^^^^^^^

.. code-block:: text

   POST /api/batches/<id>/canonicalizesmiles/

Canonicalize a list of SMILES strings or a CSV.

updatereactionsuccess
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   POST /api/batches/<id>/updatereactionsuccess/

Mark reactions as failed/unsuccessful.

Nested Serialisation (fetchall)
-------------------------------

Append ``?fetchall=yes`` to any core resource endpoint to get a deeply
nested response with all child objects inlined:

.. code-block:: text

   GET /api/projects/<id>/?fetchall=yes

Returns the project with all batches, targets, methods, reactions,
products, reactants, and catalog entries fully expanded.

.. warning::

   ``fetchall`` responses can be very large for projects with many
   targets.  Use filtering when you only need a subset.

Filtering
---------

All ViewSets that support filtering use ``django_filters``.  Pass
filter values as query parameters:

.. code-block:: text

   GET /api/targets/?batch_id=42
   GET /api/methods/?target_id=7&nosteps=2
   GET /api/reactions/?method_id=15
