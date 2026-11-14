============
Contributing
============

Thank you for your interest in contributing to CAR.  This guide covers
the development workflow, coding conventions, and testing procedures.

.. contents:: On this page
   :local:
   :depth: 2

Development Setup
-----------------

Follow the :doc:`/technical/local_deployment` guide to set up a local
devcontainer.  Once running, you will have:

- Django dev server on ``http://127.0.0.1:8000/``
- Vite dev server on ``http://127.0.0.1:3000/``
- Celery worker processing background tasks
- PostgreSQL database

Python Dependency Migration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The application requirements now target **Python 3.15** and **Django 6.1**;
compatibility with older Python releases is no longer a requirement.
``src/requirements-dev.txt`` includes the application dependencies plus
development and documentation tools. Keep its application pins aligned with
``src/requirements.txt``.

REALwrap is temporarily excluded until a PyPI release is available. Pycule
remains included: its updated release installs on Python 3.15, although live
vendor API compatibility still requires verification.

Development, production and CI containers use Python 3.15 and install Python
packages with pip into an isolated ``/opt/venv``. RDKit, pandas and psycopg2
are no longer installed through apt. The virtual environment does not expose
system site-packages; ``PATH`` selects its ``python``, ``python3``, ``pip``
and ``celery`` without requiring shell activation. In CI, the venv's bin
directory is added through ``GITHUB_PATH`` for subsequent steps.

Psycopg2 is built from source against Debian's ``libpq``. Development and CI
images include the compiler and ``libpq-dev``; the production image copies
the venv from a build stage with the same Python/base image and includes only
the native ``libpq5`` runtime dependency. Never copy a host-created venv into
these images or reuse it across Python versions or CPU architectures.

Rebuild existing devcontainers after changing requirements. VS Code's default
interpreter is ``/opt/venv/bin/python``; use **Python: Select Interpreter** if
an existing workspace selection still points elsewhere.

Branch Strategy
^^^^^^^^^^^^^^^

- ``master`` — stable release branch
- ``staging`` — integration testing
- Feature branches — branch off ``staging``, merge back via pull request

Running Tests
-------------

All tests live in ``src/backend/tests/``.

.. code-block:: bash

   cd /container/src/

   # Run the full test suite
   python3 manage.py test backend --noinput

   # Run a specific test module
   python3 manage.py test backend.tests.test_api --noinput

   # Run a specific test class
   python3 manage.py test backend.tests.test_api.TestOTProjectCreateAction --noinput

   # Run a specific test method
   python3 manage.py test backend.tests.test_api.TestOTProjectCreateAction.test_create_action --noinput

.. tip::

   Always use ``--noinput`` to avoid interactive prompts when a stale
   ``test_postgres`` database exists from interrupted runs.

The test suite includes **1080+ tests** covering:

- API endpoints and serialisation
- Model creation and validation
- OT session orchestration and plate layout
- Script generation and session visualization
- Recipe loading, generation, and validation
- Chemistry utilities (SMILES canonicalisation, conversions)
- Multichannel pipette grouping and volume calculations

Code Conventions
----------------

Backend (Python)
^^^^^^^^^^^^^^^^

- **Python 3.15** — use type hints where practical
- **Django 6.1** coding style
- **PEP 8** — enforced by linting
- Use ``snake_case`` for functions and variables
- Use ``CamelCase`` for classes
- Keep business logic in ``services/`` or dedicated modules, not in
  ViewSets
- Use descriptive variable names — avoid single-letter names except in
  tight loops

Frontend (JavaScript / React)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **React 18** with functional components and hooks
- **Vite** for build tooling
- Components in ``PascalCase.jsx``
- Keep API calls in ``common/api/``
- Use Zustand stores (``common/stores/``) for global state

Models & Migrations
^^^^^^^^^^^^^^^^^^^

When modifying Django models:

1. Make your changes in the appropriate ``models/*.py`` file
2. Generate migrations:

   .. code-block:: bash

      python3 manage.py makemigrations

3. Apply migrations:

   .. code-block:: bash

      python3 manage.py migrate

4. Run the test suite to verify nothing is broken
5. Commit the migration file with your code changes

Adding a New API Endpoint
^^^^^^^^^^^^^^^^^^^^^^^^^

1. Add or modify the ViewSet in ``api.py``
2. Register the route in ``urls.py`` (if new)
3. Add/update serializers in ``serializers.py``
4. Write tests in ``tests/``

Documentation
-------------

Documentation is built with Sphinx and hosted on Read the Docs.

.. code-block:: bash

   cd /container/docs/
   sphinx-build -b html source build/html

   # Or use make (if available)
   make html

Preview the output by opening ``build/html/index.html``.

When adding a new documentation page:

1. Create a ``.rst`` file in the appropriate subdirectory
   (``user/``, ``technical/``, or ``reference/``)
2. Add the file to the ``toctree`` in ``source/index.rst``
3. Build and verify

Filing Issues
-------------

When reporting a bug, include:

- Steps to reproduce
- Expected vs actual behaviour
- Relevant log output (``logs/logfile.log`` or Celery terminal)
- Browser console errors (for frontend issues)
