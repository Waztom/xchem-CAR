"""Backend services package.

Thin service layer that encapsulates business logic previously
scattered across api.py ViewSets, createmodels.py, and db_utils.py.

Modules
-------
batch_service
    Batch lifecycle: cloning targets/methods, creating batches,
    marking reaction outcomes.
chemistry_service
    SMILES operations, MW calculations, PubChem look-ups,
    product SMILES retrieval.
protocol_service
    OT protocol creation orchestration: parsing request data,
    dispatching Celery tasks.
"""
