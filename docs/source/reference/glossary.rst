========
Glossary
========

.. glossary::
   :sorted:

   Action Session
      A group of actions (add, stir, extract, mix) executed together as
      one logical step.  Each action session has a type (``reaction``,
      ``stir``, ``workup``, ``analyse``) and a driver (``robot`` or
      ``human``).

   Batch
      A group of target compounds within a project, derived from the
      ``batch-tag`` column in the upload CSV.  Batches can have child
      sub-batches.

   Catalog Entry
      A vendor listing for a reactant or target, including vendor name,
      catalog ID, price, and lead time.

   Compound Order
      A CSV file listing the starting materials to purchase and the
      wells to load them into on the source plate.

   Dead Volume
      The minimum volume that must remain in a well after all transfers,
      due to the labware geometry.  For a Labcyte 384LDV 100 µL plate
      this is approximately **7.56 µL**.

   Deck
      The Opentrons robot deck with 11 slots for plates and tipracks.

   Design Matrix
      A CSV or Excel file where each row defines one experiment in a
      recipe-generator screen.  Columns map to template variables.

   Error Margin
      A 15 % multiplicative safety factor (×1.15) applied to calculated
      transfer volumes to account for pipetting imprecision.

   Identity Key
      The combination of ``(smiles, concentration)`` that uniquely
      identifies a starting material for well assignment and grouping.

   Method
      A retrosynthetic route for a target molecule, consisting of one
      or more reaction steps.

   Multichannel (MC)
      An 8-channel pipette transfer mode.  All 8 tips aspirate and
      dispense simultaneously, requiring 8 identically loaded wells in
      a sub-column.

   OT Batch Protocol
      The protocol bundle for one batch, containing OT scripts,
      compound orders, solvent preps, and session visualizations,
      packaged as a ``.zip`` file.

   OT Session
      One execution session on the Opentrons robot, corresponding to a
      single reaction step or workup step.

   Plate Factory
      The backend component that creates plates, columns, and wells for
      an OT session, handling both multichannel and single-channel
      layouts.

   PlateRole
      An enum identifying a plate's function: ``reaction``, ``workup``,
      ``startingmaterial``, ``solvent``, ``spefilter``, ``lcms``,
      ``xchem``, ``nmr``, ``analyse``.

   Project
      The top-level entity in CAR.  A project groups one or more
      batches of target compounds.

   Reaction
      A single synthetic step within a method, defined by a reaction
      class and recipe.

   Reaction Class
      A named category of reaction (e.g. Amidation, Suzuki coupling,
      Mitsunobu).  Each class has one or more recipe variants.

   Recipe
      An immutable, versioned definition of a reaction procedure,
      specifying action sessions and their actions (add, stir, extract,
      mix).

   Recipe Generator
      A tool for creating recipe variants by substituting parameters
      across a design matrix of experiments, using a recipe template.

   Retrosynthesis
      Working backwards from a target molecule to identify commercially
      available building blocks and a synthetic route.

   Session Visualizer
      A component that generates self-contained HTML pages showing
      plate layouts, well contents, and transfer logs for an OT
      session.

   Single-Channel (SC)
      A single-tip pipette transfer.  Used for reactions not covered by
      multichannel groups or when multichannel mode is disabled.

   SMARTS
      SMILES Arbitrary Target Specification — a pattern language for
      substructure matching used in recipes to identify which reactant
      plays which role.

   SMILES
      Simplified Molecular-Input Line-Entry System — a string notation
      for chemical structures.  CAR canonicalises all SMILES with RDKit.

   Sub-Column
      One of two interleaved well groups within a 384-well plate.
      Sub-column 0 contains odd-numbered physical columns (1, 3, 5, …);
      sub-column 1 contains even-numbered columns (2, 4, 6, …).  Each
      sub-column has 8 wells spaced at 9 mm — matching the multichannel
      pipette tip spacing.

   Target
      A molecule to be synthesised within a batch, specified by a
      SMILES string and a required quantity.

   Transfer Type
      A field on the Well model indicating whether the well is used for
      single-channel (``single``) or multichannel (``multichannel``)
      transfers.

   Well
      A single position on a plate, identified by name (e.g. ``A03``).
      Stores volume, SMILES, concentration, solvent, and transfer type.

   Well Finder
      The backend utility that locates source wells for a given
      starting material during OT script generation, optionally
      filtering by transfer type.
