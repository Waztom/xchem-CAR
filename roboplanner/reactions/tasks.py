from celery import shared_task
from django.core.files.storage import default_storage

from rdkit import Chem
from rxn4chemistry import RXN4ChemistryWrapper
import pandas as pd

import time
from .validate import (
    add_warning,
    checkColumnName,
    checkNumberColumns,
    checkSMILES)

# Setup IBM RxN API
# api_key=os.environ['IBM_API_KEY'] 
# rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=api_key)
# rxn4chemistry_wrapper.create_project('Test actions')

def delete_tmp_file(filepath):
    default_storage.delete(filepath)

# Celery validate task
@shared_task
def validateFileUpload(csv_fp, validate_only=True):
    """ Celery task to process validate the uploaded files for retrosynthesis planning.
    
    Parameters
    ----------
    csv_fp: str
        filepath of the uploaded csv file, which is saved to temporary storage by `viewer.views.UploadCSV`
    project_details: dict
        dictionary of project details (name, email and project_name) that will be used to create the project model
        if the csv file is validated
    
    Returns
    -------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the validation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded csv file
    """
    
    validated = True
    smiles_list = []
    
    validate_dict = {'target_name': [],
                     'field': [],
                     'warning_string': []}
    
    # Open csv file as Pandas df
    uploaded_csv_df = pd.read_csv(csv_fp)
    
    # Check no of column headings and name of column heading
    columns = uploaded_csv_df.columns
    validate_dict = checkNumberColumns(columns, validate_dict)
    
    if len(validate_dict['target_name']) != 0:
        validated = False
    else:
        # Check SMILES
        indexes = [i for i,smi in enumerate(uploaded_csv_df['Targets'])]
        smiles_list = [smi for i, smi in enumerate(uploaded_csv_df['Targets'])]
        
        for smi, index in zip(smiles_list, indexes): 
            validate_dict = checkSMILES(smi, index, validate_dict)

    if len(validate_dict['target_name']) !=0:
        validated = False
    
    # Delete tempory file if only validate selected
    if validate_only:
        default_storage.delete(csv_fp)
        csv_fp = None
      
    return (validate_dict, validated, csv_fp, smiles_list)
