===============
Troubleshooting
===============

Common issues and solutions when running CAR.

.. contents:: On this page
   :local:
   :depth: 2

Uploads & Celery
----------------

**"Upload stuck at PENDING"**
   The Celery worker is not running.  Start it:

   .. code-block:: bash

      cd /container/src/
      celery -A CAR worker -l info

**"Task failed: connection refused to RabbitMQ"**
   RabbitMQ is not running.  If using Docker Compose,
   ``docker-compose up rabbit``.  In the devcontainer the broker is
   started automatically.

**"No module named 'backend'"**
   You are not in the correct working directory.  Run commands from
   ``/container/src/``.

OT Protocol Generation
-----------------------

**"ValueError: deck overflow"**
   Too many plates for the 11-slot OT deck.  CAR automatically
   splits large reaction sets, but if you still see this, try
   reducing the batch size or splitting into sub-batches.

**"No starting-material wells found"**
   Ensure the Celery task completed successfully and that the batch
   has reactions with valid reactants and catalog entries.

**Empty compound-order CSV**
   Likely caused by missing reactant data or catalogue entries.  Check
   that retrosynthesis ran successfully for the batch.

Database
--------

**"Stale test_postgres database"**
   If test runs are interrupted, a leftover ``test_postgres`` database
   can block future test runs.  Use the ``--noinput`` flag to
   auto-confirm deletion:

   .. code-block:: bash

      python3 manage.py test backend --noinput

**"Relation does not exist"**
   Run migrations:

   .. code-block:: bash

      python3 manage.py makemigrations
      python3 manage.py migrate

Frontend
--------

**"EACCES: permission denied" on npm**
   On WSL 2, file permissions can get corrupted when switching
   branches.  Fix from *outside* the devcontainer:

   .. code-block:: bash

      sudo chown -R $(whoami) /path/to/xchem-CAR

**Vite dev server not hot-reloading**
   Ensure the Vite server is running in a separate terminal:

   .. code-block:: bash

      cd /container/src/frontend/
      npm run dev

**CORS errors**
   The Django settings whitelist ``localhost:3000`` and
   ``127.0.0.1:3000``.  If you are using a different port, update
   ``CORS_ALLOWED_ORIGINS`` in ``CAR/settings.py``.

git-crypt
---------

**"Error: no key found"**
   You need the ``crypt-key`` file from the project maintainer.  See
   :doc:`/technical/local_deployment`.

**"Already unlocked"**
   Safe to ignore — the repository secrets are already decrypted.

Docker
------

**Container build fails on M1 / ARM Mac**
   The ``Dockerfile.prod`` targets ``linux/amd64``.  Use Docker
   Desktop's Rosetta emulation or build with
   ``--platform linux/amd64``.

**PostgreSQL data lost after rebuild**
   Ensure the ``pgdata`` Docker volume is not being removed.
   ``docker-compose down`` preserves volumes;
   ``docker-compose down -v`` **deletes** them.
