==============
Django Models
==============

CAR's data layer is organised into four model modules.

.. contents:: On this page
   :local:
   :depth: 2

Core Models (``models/core.py``)
---------------------------------

Project
^^^^^^^

Top-level entity representing a synthesis campaign.

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``name``
     - SlugField
     - URL-safe project name (unique)
   * - ``init_date``
     - DateTimeField
     - Creation timestamp
   * - ``submittername``
     - CharField
     - Name of the submitter
   * - ``submitterorganisation``
     - CharField
     - Organisation
   * - ``proteintarget``
     - CharField
     - Protein target of interest
   * - ``quotedcost``
     - CharField
     - Estimated cost
   * - ``quoteurl``
     - URLField
     - Vendor quote URL

Batch
^^^^^

A group of target compounds, derived from the ``batch-tag`` CSV column.
Supports sub-batches via a self-referential foreign key.

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``project_id``
     - FK → Project
     - Parent project
   * - ``batch_id``
     - FK → Batch (nullable)
     - Parent batch (for sub-batches)
   * - ``batchtag``
     - CharField
     - Tag from the upload CSV

Target
^^^^^^

A single target molecule for synthesis.

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``batch_id``
     - FK → Batch
     - Parent batch
   * - ``smiles``
     - CharField
     - Canonical SMILES string
   * - ``image``
     - ImageField
     - 2-D molecule rendering
   * - ``name``
     - CharField
     - Human-readable name
   * - ``concentration``
     - FloatField
     - Target concentration
   * - ``volume``
     - FloatField
     - Target volume
   * - ``mass``
     - FloatField
     - Target mass (mg)
   * - ``mols``
     - FloatField
     - Target amount (mol)

Method
^^^^^^

A retrosynthetic route to reach a target.

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``target_id``
     - FK → Target
     - Parent target
   * - ``nosteps``
     - IntegerField
     - Number of reaction steps
   * - ``otchem``
     - BooleanField
     - Compatible with OT chemistry

Reaction
^^^^^^^^

A single reaction step within a method.

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``method_id``
     - FK → Method
     - Parent method
   * - ``reactionclass``
     - CharField
     - Reaction class name (e.g. ``Amidation``)
   * - ``recipe``
     - CharField
     - Recipe name variant
   * - ``recipe_id``
     - FK → Recipe (nullable)
     - Linked versioned recipe
   * - ``groupbycolumn``
     - CharField
     - Column grouping key
   * - ``number``
     - IntegerField
     - Step number within the method
   * - ``success``
     - BooleanField
     - Whether synthesis succeeded
   * - ``temperature``
     - FloatField
     - Reaction temperature
   * - ``intramolecular``
     - BooleanField
     - Intramolecular reaction flag
   * - ``image``
     - ImageField
     - Reaction diagram

Product / Reactant / CatalogEntry / PubChemInfo
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- **Product** — output compound of a reaction (FK → Reaction,
  FK → PubChemInfo)
- **Reactant** — input compound (FK → Reaction, FK → PubChemInfo);
  ``previousreactionproduct`` flag indicates it comes from a prior step
- **CatalogEntry** — vendor listing (FK → Reactant or Target); stores
  ``vendor``, ``catalogid``, ``priceinfo``, ``upperprice``,
  ``leadtime``
- **PubChemInfo** — cached PubChem data: ``compoundid``,
  ``summaryurl``, ``lcssurl``, ``smiles``, ``cas``

Action Models (``models/actions.py``)
--------------------------------------

ActionSession
^^^^^^^^^^^^^

Groups related actions into a logical execution unit.

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``reaction_id``
     - FK → Reaction
     - Parent reaction
   * - ``sessionnumber``
     - IntegerField
     - Ordering within the reaction
   * - ``type``
     - CharField
     - ``reaction`` / ``stir`` / ``workup`` / ``analyse``
   * - ``driver``
     - CharField
     - ``robot`` or ``human``
   * - ``continuation``
     - BooleanField
     - Continues from previous session

AddAction
^^^^^^^^^

A liquid-handling transfer.

Key fields: ``smiles``, ``volume``, ``molecularweight``, ``solvent``,
``concentration``, ``calcunit``, ``from_plate_role``,
``to_plate_role``, ``from_plate_role_index``, ``to_plate_role_index``.

StirAction
^^^^^^^^^^

A stir/incubation step.

Key fields: ``duration``, ``durationunit``, ``temperature``,
``temperatureunit``, ``stirringspeed``, ``plate_role``.

ExtractAction
^^^^^^^^^^^^^

A liquid-liquid extraction step.

Key fields: ``layer`` (top/bottom), ``smiles``, ``volume``,
``bottomlayervolume``, ``solvent``.

MixAction
^^^^^^^^^

A pipette mixing step.

Key fields: ``repetitions``, ``plate_role``.

OT Models (``models/ot.py``)
------------------------------

OTProject → OTBatchProtocol → OTSession
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - Model
     - Key Fields
     - Description
   * - **OTProject**
     - ``project_id``, ``name``
     - Container for a protocol run
   * - **OTBatchProtocol**
     - ``otproject_id``, ``batch_id``, ``zipfile``
     - Per-batch protocol bundle
   * - **OTSession**
     - ``otbatchprotocol_id``, ``reactionstep``, ``sessiontype``
     - One robot session (reaction/workup/lcmsprep)

Deck, Plate, Column, Well
^^^^^^^^^^^^^^^^^^^^^^^^^^

Physical deck layout representation.

- **Deck** — 11-slot OT deck
- **Plate** — a plate on the deck, with ``role`` (PlateRole enum),
  ``labware``, ``maxwellvolume``, ``numberwells``
- **Column** — groups wells by column index, ``reactionclass``
- **Well** — individual well with ``name`` (e.g. ``A03``), ``volume``,
  ``smiles``, ``concentration``, ``solvent``,
  ``transfer_type`` (``single`` or ``multichannel``), ``available``

Pipette, TipRack, CompoundOrder, OTScript, SolventPrep
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Supporting OT session models for hardware configuration and generated
artefacts.

Recipe Models (``models/recipes.py``)
--------------------------------------

Recipe
^^^^^^

Immutable, versioned recipe definition.

.. list-table::
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``reaction_class``
     - CharField
     - Reaction class (e.g. ``Amidation``)
   * - ``name``
     - CharField
     - Recipe variant name
   * - ``version``
     - IntegerField
     - Auto-incremented version
   * - ``parent``
     - FK → Recipe (nullable)
     - Previous version
   * - ``description``
     - TextField
     - Human-readable description
   * - ``reaction_smarts``
     - JSONField
     - SMARTS patterns for reaction
   * - ``estimated_yield``
     - FloatField
     - Expected yield percentage

Uniqueness constraint: ``(reaction_class, name, version)``.

RecipeActionSession / RecipeAddAction / RecipeStirAction / RecipeExtractAction / RecipeMixAction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mirror the Action models but define the *template* actions that will be
instantiated for each reaction.  Key difference: recipe actions use
``material_smarts`` (SMARTS pattern) instead of a concrete ``smiles``.
