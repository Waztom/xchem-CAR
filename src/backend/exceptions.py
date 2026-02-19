"""Domain exception hierarchy for the CAR backend.

Provides structured, catchable error types so callers can distinguish
between *expected* failures (invalid SMILES, no reaction products, …)
and *unexpected* internal errors.  Every exception carries the original
cause via standard ``raise … from e`` chaining.

Hierarchy
---------
CARError
├── ChemistryError            – RDKit / molecular-processing failures
│   ├── SMILESParsingError    – invalid or un-parseable SMILES input
│   ├── SMARTSReactionError   – reaction-SMARTS application failures
│   ├── MolecularPropertyError– MW, formula, InChI, etc. failures
│   └── SVGGenerationError    – molecule/reaction drawing failures
├── ModelCreationError        – Django model build / save failures
│   └── ActionModelError      – action-model (Add/Extract/Mix/Stir) failures
└── ServiceError              – service-layer / API-level failures
"""


class CARError(Exception):
    """Base exception for all CAR backend errors."""


# ── Chemistry layer ─────────────────────────────────────────────────


class ChemistryError(CARError):
    """Any RDKit / molecular-processing failure."""


class SMILESParsingError(ChemistryError):
    """A SMILES string could not be parsed into a valid molecule."""


class SMARTSReactionError(ChemistryError):
    """A reaction-SMARTS pattern could not be applied to the given reactants."""


class MolecularPropertyError(ChemistryError):
    """A molecular property (MW, formula, InChI key, …) could not be computed."""


class SVGGenerationError(ChemistryError):
    """An SVG image could not be generated for a molecule or reaction."""


# ── Model-creation layer ────────────────────────────────────────────


class ModelCreationError(CARError):
    """A Django model instance could not be built or saved."""


class ActionModelError(ModelCreationError):
    """An action model (AddAction, ExtractAction, …) could not be created."""


# ── Service / API layer ─────────────────────────────────────────────


class ServiceError(CARError):
    """A service-layer operation failed."""
