================
Local Deployment
================

This guide walks through setting up CAR for local development using
the VS Code devcontainer.

.. contents:: On this page
   :local:
   :depth: 2

Prerequisites
-------------

- **Docker Desktop** (Windows / macOS) or **Docker Engine** (Linux)
- **Docker Compose** (Linux only — bundled with Docker Desktop on
  Windows / macOS)
- **Visual Studio Code** with the **Remote - Containers** extension
- **git-crypt** — for decrypting secrets
- *Optional:* **WSL 2** (Windows users)

Windows / WSL 2 Setup
^^^^^^^^^^^^^^^^^^^^^

If you are on Windows, install WSL 2 first:

- `WSL 2 installation guide <https://docs.microsoft.com/en-gb/windows/wsl/install-win10>`_
- Enable WSL integration in Docker Desktop → Settings → Resources →
  WSL Integration

Install the **Remote - WSL** extension in VS Code.

Clone the Repository
--------------------

.. code-block:: bash

   sudo apt install git          # if git is not installed
   cd ~
   git clone https://github.com/Waztom/xchem-CAR.git

Unlock Secrets
--------------

CAR uses `git-crypt <https://github.com/AGWA/git-crypt>`_ to encrypt
environment variables and API keys.  Obtain the ``crypt-key`` file from
the project maintainer.

.. code-block:: bash

   # Install git-crypt
   sudo apt-get update && sudo apt-get install git-crypt   # Debian/Ubuntu
   brew install git-crypt                                  # macOS

   # Unlock
   cd xchem-CAR
   git-crypt unlock /path/to/crypt-key

Open in Devcontainer
--------------------

1. Open the repository folder in VS Code.
2. Press **Ctrl+Shift+P** → **Remote-Containers: Open Folder in
   Container**.
3. Select the repository root and click **Open**.
4. Wait for the container to build (first time takes a few minutes).

Launch Services
---------------

Open **three** integrated terminals inside the devcontainer:

**Terminal 1 — Backend:**

.. code-block:: bash

   cd /container/.devcontainer/
   chown -R root launch-backend.sh launch-frontend.sh
   ./launch-backend.sh
   ./launch-frontend.sh

**Terminal 2 — Django server:**

.. code-block:: bash

   cd /container/src/
   python3 manage.py runserver

**Terminal 3 — Celery worker** (required for file uploads and protocol
generation):

.. code-block:: bash

   cd /container/src/
   celery -A CAR worker -l info

**Terminal 4** *(optional)* **— Vite dev server:**

.. code-block:: bash

   cd /container/src/frontend/
   npm run dev


Utils available on local deployment backend
---------------
To make use of the various utils available in the codebase, 
you need to first run the command below to start a python shell:
.. code-block:: bash

   cd /container/src/
   python3 manage.py shell

Then you can import the utils you need from the relevant modules. For example:

**chem utils:**
first run the command below to import chem_utils into your python shell:
.. code-block:: python

   from backend.chem_utils import *

Then you can use the various functions available in chem_utils, such as:
.. list-table::
   :header-rows: 1

   * - Utility tool name
     - Description
   * - canon_smiles
     - Canonicalise a SMILES string
   * - get_mws
     - Get molecular weights for a list of SMILES strings
   * - get_inchi_keys
     - Get inchi key for a SMILES string
   * - get_molecular_formula
     - Get molecular formula for a SMILES string
   * - strip_salts
     - Remove salts from a SMILES string, returning the largest fragment by
       molecular weight
   * - match_smarts
     - Matches a SMILES pattern to a SMARTS pattern
   
and more!

**db utils:**
first run the command below to import chem_utils into your python shell:
.. code-block:: python

   from backend.db_utils import *

Then you can use the various functions available in chem_utils, such as:
.. list-table::
   :header-rows: 1

   * - Utility tool name
     - Description
   * - get_plate_map
     - Generates a Plate Map for a list of plate IDs
   * - get_add_action_query_set
     - Get add actions queryset for a given reaction_id
   * - get_ot_batch_protocol_query_set
     - Get OT-2 protocol query set for a given batch_id
   * - update_target_mols
     - Updates the Target mols in a Batch - using batch_id, concentration, and volume

and more!

**recipe utils:**
first run the command below to import chem_utils into your python shell:
.. code-block:: python

   from backend.recipe_utils import *

Then you can use the various functions available in chem_utils, such as:
.. list-table::
   :header-rows: 1

   * - Utility tool name
     - Description
   * - get_latest_recipe
     - Returns the latest (or a specific) Recipe for a reaction class + name.
   * - recipe_exists
     - Checks whether *any* recipe exists for the given reaction class
   * - get_session_recipe_actions
     - Returns a sorted list of recipe action **model instances** for a session

and more!

Access Points
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Service
     - URL
   * - Frontend
     - http://127.0.0.1:3000/
   * - Django REST API
     - http://127.0.0.1:8000/

On subsequent launches only the Django server and Celery worker are
needed — migrations and ``npm install`` do not need to be re-run unless
models or packages have changed.

Model Migrations
^^^^^^^^^^^^^^^^

If Django models change:

.. code-block:: bash

   cd /container/src/
   python3 manage.py makemigrations
   python3 manage.py migrate

Or re-run ``launch-backend.sh``.

Docker Compose (Production)
---------------------------

For standalone deployment without a devcontainer:

.. code-block:: bash

   cd /container/src/
   docker-compose up --build

This starts five services:

.. list-table::
   :header-rows: 1

   * - Service
     - Image
     - Purpose
   * - **app**
     - ``Dockerfile.prod``
     - Django + Vite production build; port 8000
   * - **db**
     - ``postgres:10-alpine``
     - PostgreSQL database
   * - **redis**
     - ``redis:alpine``
     - Celery result backend
   * - **celery**
     - ``Dockerfile.prod``
     - Celery worker
   * - **rabbit**
     - ``rabbitmq:alpine``
     - Celery message broker (AMQP)

Persistent data is stored in a ``pgdata`` Docker volume.

Troubleshooting
^^^^^^^^^^^^^^^

**Permission denied on WSL 2:**

If ``npm run dev`` or environment variable access fails after switching
branches:

.. code-block:: bash

   # Run this OUTSIDE the devcontainer, in your WSL terminal
   sudo chown -R $(whoami) /path/to/xchem-CAR
