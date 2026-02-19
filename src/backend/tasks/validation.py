"""Validation and SMILES canonicalization tasks."""
from __future__ import annotations

import os

import pandas as pd
from celery import shared_task
from rdkit import Chem

from ..chem_utils import canon_smiles
from ..validate import ValidateFile


def delete_tmp_file(filepath):
    os.remove(filepath)


@shared_task
def validate_file_upload(
    csv_fp, validate_type: str = None, project_info=None, validate_only=True
):
    """Celery task to process validate the uploaded files for retrosynthesis planning.

    Parameters
    ----------
    csv_fp: str
        filepath of the uploaded csv file, which is saved to temporary storage by `viewer.views.UploadCSV`
    validate_type: validate different types of upload files
    project_info: dict
        dictionary of project details (name, email and project_name) that will be used to create the project model
        if the csv file is validated
    validate_only: boolean
        set to True to delete tmp file saved for validation and set to False to save tmp file and keep
        for uploading and creating model objects

    Returns
    -------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the validation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded csv file
    """

    validation = ValidateFile(csv_to_validate=csv_fp, validate_type=validate_type)

    if validate_only:
        delete_tmp_file(csv_fp)
        csv_fp = None

    uploaded_dict = validation.df.to_dict("list")
    validated = validation.validated
    validate_dict = validation.validate_dict

    return (
        validate_type,
        validate_dict,
        validated,
        project_info,
        csv_fp,
        uploaded_dict,
    )


@shared_task
def canonicalize_smiles(csvfile: str = None, smiles: list = None):
    """ "
    Canonicalizes smiles from csv file uploaded from frontend
    """
    validated = True

    if csvfile:
        csvdf = pd.read_csv(csvfile, encoding="utf8")
        smiles = [smi for smi in csvdf["SMILES"]]
        delete_tmp_file(csvfile)
    if smiles:
        smiles = smiles

    molcheck = [Chem.MolFromSmiles(smi) for smi in smiles]

    if None not in molcheck:
        canonicalizedsmiles = [canon_smiles(smi) for smi in smiles]
        return validated, canonicalizedsmiles
    else:
        validated = False
        indexerrors = [i for i, v in enumerate(molcheck) if v is None]
        errorsummary = "There was an error with the smiles csv at index: {}".format(
            indexerrors
        )
        return validated, errorsummary
