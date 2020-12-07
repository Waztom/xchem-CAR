from celery import shared_task
from django.core.files.storage import default_storage

import pandas as pd

import time
from .validate import (
    add_warning,
    checkColumnNames,
    checkNumberColumns,
    checkSMILES,
    checkIsNumber)

from .createmodels import (
    createProjectModel,
    createTargetModel
)

def delete_tmp_file(filepath):
    default_storage.delete(filepath)

# Celery validate task
@shared_task
def validateFileUpload(csv_fp, project_info=None, validate_only=True):
    """ Celery task to process validate the uploaded files for retrosynthesis planning.
    
    Parameters
    ----------
    csv_fp: str
        filepath of the uploaded csv file, which is saved to temporary storage by `viewer.views.UploadCSV`
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
    
    validated = True
    smiles_list = []
    
    validate_dict = {'field': [],
                     'warning_string': []}
    
    # Open csv file as Pandas df
    uploaded_csv_df = pd.read_csv(csv_fp)
    
    # Check no of column headings and name of column headings
    columns = uploaded_csv_df.columns
    validate_dict = checkNumberColumns(columns, validate_dict)
    
    if len(validate_dict['warning_string']) != 0:
        validated = False
    else:
        # Check SMILES 
        indexes = [i for i,smi in enumerate(uploaded_csv_df['Targets'])]
        smiles_list = [smi for smi in uploaded_csv_df['Targets']]
        ammounts_list = [amount for amount in uploaded_csv_df['Ammount_required (mg)']]
        
        for index, smi, ammount in zip(indexes, smiles_list, ammounts_list): 
            validate_dict = checkSMILES(smi, index, validate_dict)
            validate_dict = checkIsNumber(ammount, index, validate_dict)

    if len(validate_dict['warning_string']) !=0:
        validated = False
    
    # Delete tempory file if only validate selected
    if validate_only:
        default_storage.delete(csv_fp)
        csv_fp = None
      
    return (validate_dict, validated, csv_fp, project_info, smiles_list)

    # Functions to request info from IBM API
    
@shared_task
def uploadIBMReaction(validate_output):
    # Validate output is a list - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    validate_dict, validated, csv_fp, project_info, smiles_list = validate_output
    
    if not validated:
        # Delete tempory file if only validate selected
        default_storage.delete(csv_fp)
        return (validate_dict, validated)
    
    if validated:
        # Create project model and return project id
        project_id  = createProjectModel(project_info)
        
        slug_no = 1
        for smiles in smiles_list:

            # Number of steps/name of reactants from IBM API???????
            # Create a Target model and return target id
            createTargetModel(project_id = project_id, smiles=smiles, slug_no=slug_no)
                # Create a Method model

                    # Create a Reaction model

                        # Create Reactant models

                            # Create Action models
            slug_no += 1
        
        # Delete tempory file if only validate selected
        default_storage.delete(csv_fp)
    
    return validate_dict, validated
 
            
