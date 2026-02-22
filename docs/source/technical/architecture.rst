============
Architecture
============

CAR (Chemist Assisted Robotics) is a full-stack web application for
automated retrosynthesis and liquid-handling protocol generation.

.. contents:: On this page
   :local:
   :depth: 2

High-Level Overview
-------------------

.. code-block:: text

   ┌─────────────┐     REST/JSON     ┌───────────────────┐
   │  React 18   │ ◄──────────────►  │  Django 3.1 (DRF) │
   │  (Vite)     │    :3000/:8000    │                   │
   └─────────────┘                   │  ├── api.py       │
                                     │  ├── models/      │
                                     │  ├── serializers   │
                                     │  ├── services/     │
                                     │  └── tasks/        │
                                     └────────┬──────────┘
                                              │  Celery task
                                              ▼
                                     ┌────────────────────┐
                                     │  RabbitMQ  │ Redis  │
                                     │  (broker)  │(result)│
                                     └────────┬───────────┘
                                              │
                                              ▼
                                     ┌────────────────────┐
                                     │    PostgreSQL 10   │
                                     └────────────────────┘

Technology Stack
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Layer
     - Technology
     - Version
   * - Frontend
     - React, Vite
     - 18, 4.x
   * - Backend
     - Django, Django REST Framework
     - 3.1.7, 3.x
   * - Task queue
     - Celery
     - 5.x
   * - Broker
     - RabbitMQ
     - alpine
   * - Result backend
     - Redis
     - alpine
   * - Database
     - PostgreSQL
     - 10
   * - Chemistry
     - RDKit
     - 2022+
   * - Containerisation
     - Docker, Docker Compose
     - —

Backend Structure
-----------------

The Django backend lives in ``src/backend/``:

.. code-block:: text

   backend/
   ├── models/
   │   ├── core.py        # Project → Batch → Target → Method → Reaction
   │   ├── actions.py     # ActionSession, AddAction, StirAction, …
   │   ├── ot.py          # OTProject, OTSession, Plate, Well, …
   │   └── recipes.py     # Recipe, RecipeActionSession, RecipeAddAction, …
   ├── api.py              # DRF ViewSets and custom actions
   ├── serializers.py      # DRF serializers (flat and nested)
   ├── services/           # Business logic (protocol_service)
   ├── tasks/              # Celery tasks (ot_protocol, validation)
   ├── opentrons/          # OT protocol generation engine
   │   ├── otsession/      # Session orchestration, plate/deck/pipette managers
   │   └── otwriter/       # Script generation, session visualization
   ├── recipebuilder/      # Legacy encoded recipes
   ├── recipegenerator/    # Template-driven recipe variant generation
   ├── enamine/            # Enamine REAL Tools API integration
   ├── manifold/           # Manifold retrosynthesis API
   ├── mcule/              # MCule vendor API
   └── IBM/                # IBM RXN API

Model Hierarchy
^^^^^^^^^^^^^^^

The core domain model follows a strict hierarchy:

.. code-block:: text

   Project
   └── Batch
       └── Target
           └── Method
               └── Reaction
                   ├── Product
                   ├── Reactant → CatalogEntry[]
                   └── ActionSession[]
                       ├── AddAction
                       ├── StirAction
                       ├── ExtractAction
                       └── MixAction

The OT (OpenTrons) domain mirrors this with:

.. code-block:: text

   OTProject
   └── OTBatchProtocol → .zipfile
       └── OTSession (per reaction step)
           ├── Deck → Plate[] → Column[] → Well[]
           ├── Pipette
           ├── TipRack
           ├── CompoundOrder
           ├── OTScript
           └── SolventPrep

Frontend Structure
------------------

The React frontend lives in ``src/frontend/src/``:

.. code-block:: text

   src/
   ├── App/              # Application shell
   ├── Layout/           # Navigation, header, sidebar
   ├── Body/             # Route-level page components
   ├── Project/          # Project view (batch table, reactions, sub-batches)
   ├── Actions/          # IBM RXN action display components (19 types)
   ├── SetActionInputs/  # Form controls for action parameters
   ├── MolDrawer/        # JSME molecular editor modal
   ├── RDKit/            # Client-side RDKit (WASM) for 2-D rendering
   └── common/           # Shared API hooks, stores, constants, utilities

The frontend communicates with the backend exclusively through the DRF
REST API.

Request Lifecycle
-----------------

A typical project creation request flows through the following layers:

1. **React frontend** — user fills in the form and clicks Submit
2. **DRF ViewSet** (``api.py``) — validates input, calls service layer
3. **Service layer** (``services/``) — orchestrates business logic,
   dispatches Celery task
4. **Celery task** (``tasks/``) — performs heavy processing
   asynchronously (retrosynthesis calls, OT protocol generation)
5. **Models** — persist results to PostgreSQL
6. **Frontend polling** — periodically calls ``gettaskstatus`` until
   the task completes, then fetches results

OT Protocol Generation Pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :doc:`/user/ot_protocols` for the user-facing description.
Internally:

1. ``create_ot_script`` (Celery) iterates over batches and reaction
   steps
2. ``SessionOrchestrator`` sets up deck, pipettes, and plates via
   manager classes
3. ``PlateFactory`` creates well, column, and MC/SC plate layouts
4. ``ScriptGenerator`` emits the Python OT protocol and triggers the
   ``SessionVisualizer``
5. ``ZipOTBatchProtocol`` packages everything into a ``.zip``
