================
Recipe Generator
================

CAR uses **recipes** to encode the exact sequence of actions (add,
stir, extract, mix) needed to execute a reaction class.  The
**Recipe Generator** lets you explore and create variants of these
recipes for condition screening.

.. contents:: On this page
   :local:
   :depth: 2

What Is a Recipe?
-----------------

A recipe is a versioned, immutable definition of a reaction procedure.
It consists of:

- **Reaction class** — e.g. ``Amidation``, ``Suzuki``, ``Mitsunobu``
- **Name** — e.g. ``standard``, ``high-temp``
- **Version** — auto-incremented
- One or more **action sessions**, each containing action steps

Each action session has:

- ``session_type``: ``reaction``, ``stir``, ``workup``, or ``analyse``
- ``driver``: ``robot`` or ``human``
- Ordered list of **actions**: add, stir, extract, or mix

Recipe JSON Format
------------------

Recipes can be defined as JSON fixtures and loaded with the
``load_recipe_json`` management command.

Minimal example:

.. code-block:: json

   {
     "reaction_class": "Amidation",
     "name": "standard",
     "action_sessions": [
       {
         "session_type": "reaction",
         "driver": "robot",
         "actions": [
           {
             "type": "add",
             "material_smarts": "[#6:1](=[#8:2])-[#8;H1]",
             "amount": 1.0,
             "solvent": "DMA",
             "concentration": 0.5,
             "from_plate_role": "startingmaterial",
             "to_plate_role": "reaction"
           }
         ]
       },
       {
         "session_type": "stir",
         "driver": "human",
         "actions": [
           {
             "type": "stir",
             "temperature": 25,
             "duration": 12.0,
             "plate_role": "reaction"
           }
         ]
       }
     ]
   }

Loading Recipes
^^^^^^^^^^^^^^^

.. code-block:: bash

   # Validate without writing to the database
   python3 manage.py load_recipe_json --dry-run backend/tests/fixtures/amidation_standard.json

   # Load into the database
   python3 manage.py load_recipe_json backend/tests/fixtures/amidation_standard.json

   # Create a new version of an existing recipe
   python3 manage.py load_recipe_json --update backend/tests/fixtures/amidation_standard.json

Units Reference
^^^^^^^^^^^^^^^

**Quantity units** (``quantity_unit`` field on add actions):

.. list-table::
   :header-rows: 1
   :widths: 10 30 30

   * - Value
     - Description
     - Example
   * - ``moleq``
     - Molar equivalents (default)
     - ``1.5`` → 1.5 eq
   * - ``masseq``
     - Mass equivalents
     - ``2.0`` → 2× mass
   * - ``uL``
     - Microlitres
     - ``50`` → 50 µL
   * - ``mL``
     - Millilitres
     - ``0.5`` → 500 µL
   * - ``mg``
     - Milligrams
     - ``10`` → 10 mg
   * - ``g``
     - Grams
     - ``0.1`` → 100 mg
   * - ``M``
     - Molarity
     - ``0.1`` → 0.1 M
   * - ``uM``
     - Micromolarity
     - ``100`` → 100 µM

**Duration units** (``duration_unit`` field on stir actions):

.. list-table::
   :header-rows: 1

   * - Value
     - Description
   * - ``s``
     - Seconds
   * - ``m``
     - Minutes
   * - ``h``
     - Hours (default)

**Plate roles** (valid for ``plate_role`` / ``from_plate_role`` /
``to_plate_role``):

``reaction`` · ``workup`` · ``spefilter`` · ``lcms`` · ``xchem`` ·
``nmr`` · ``startingmaterial`` · ``solvent`` · ``analyse``

Generating Recipe Variants
--------------------------

The ``RecipeGenerator`` module lets you systematically vary reaction
conditions (solvent, temperature, equivalents, addition order) across
a design matrix of experiments.

Step 1: Create a Template
^^^^^^^^^^^^^^^^^^^^^^^^^^

Auto-generate a template from an existing recipe:

.. code-block:: python

   from backend.recipegenerator.parsers import create_template_from_recipe

   template = create_template_from_recipe(
       reaction_class="Amidation",
       recipe_name="standard",
       template_name="my-amidation-screen",
   )

   import json
   with open("my_template.json", "w") as f:
       json.dump(template.to_dict(), f, indent=2)

Or load a pre-built template:

.. code-block:: python

   from backend.recipegenerator.parsers import TemplateParser

   parser = TemplateParser()
   template = parser.parse_json("amidation_template.json")

A template maps **semantic action IDs** to positions in the base
recipe and defines **variables** pointing to the parameters you want to
screen (e.g. ``solvent``, ``concentration``, ``temperature``).

Step 2: Prepare a Design Matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The design matrix is a CSV (or Excel) file where each row defines one
experiment and each column maps to a template variable.

For example:

in the amidation_template.json, this code block defines the variables.
``"variables": {
        "acid_solvent": {
            "action_id": "add_acid",
            "path": "content.material.solvent"
        },
        "acid_concentration": {
            "action_id": "add_acid",
            "path": "content.material.concentration"
        },
        "coupling_agent_equiv": {
            "action_id": "add_coupling_agent",
            "path": "content.material.quantity.value"
        },
        "base_equiv": {
            "action_id": "add_base",
            "path": "content.material.quantity.value"
        },
        "amine_equiv": {
            "action_id": "add_amine",
            "path": "content.material.quantity.value"
        },
        "reaction_temperature": {
            "action_id": "stir_reaction",
            "path": "content.temperature.value"
        },
        "reaction_duration": {
            "action_id": "stir_reaction",
            "path": "content.duration.value"
        }
    },``

In the design csv, these values can be asigned values, and each row is a unique and new recipe:
.. list-table::
   :header-rows: 1
   * - acid_solvent
     - DMF
     - DMSO
   * - acid_concentration
     - 0.25
     - 0.5
   * - coupling_agent_equiv
     - 1
     - 2
   * - base_equiv
     - 3
     - 5
     - 10
   * - amine_equiv
     - 1
     - 1.2
   * - reaction_temperature
     - 25
     - 50
   * - reaction_duration
     - 12
     - 48


 Below is a list of columns not related to variables, but are still useful to know. Please note that ``recipe_name`` is a mandatory column.

Special columns:

.. list-table::
   :header-rows: 1

   * - Column
     - Purpose
   * - ``experiment_id``
     - Optional identifier for the row
   * - ``recipe_name``
     - Name for the generated recipe variant
   * - ``order_preset``
     - Named ordering preset to apply
   * - ``{session_id}_action_order``
     - Comma-separated action IDs for reordering within a session
   * - ``session_order``
     - Comma-separated session IDs for session reordering
   * - ``notes``
     - Ignored — for human documentation

Columns that do not match any template variable are ignored (with a
warning).

Step 3: Generate Recipes
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

   from backend.recipegenerator import RecipeGenerator

   generator = RecipeGenerator(template)
   recipes = generator.from_csv(
       "amidation_design_matrix.csv",
       recipe_name_column="recipe_name",
   )

   for recipe_data in recipes:
       print(f"Generated: {recipe_data['name']}")

Each generated ``recipe_data`` dict contains the full encoded recipe
ready for use with CAR's recipe loading infrastructure.

.. seealso::

   - Example files in ``src/backend/recipegenerator/examples/``
   - Fixture documentation in ``src/backend/tests/fixtures/README.md``