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
    createTargetModel,
    createMethodModel,
    createReactionModel,
    createProductModel,
    createReactantModel,
    createActionModel
)

from.IBMapicalls import (
    createIBMProject,
    getIBMRetroSyn,
    collectIBMReactionInfo    
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
        smiles_list = [smi.strip() for smi in uploaded_csv_df['Targets']]
        uploaded_csv_df['Targets'] = smiles_list   
        print(uploaded_csv_df['Targets'])
        amounts_list = [amount for amount in uploaded_csv_df['Ammount_required (mg)']]
        
        for index, smi, ammount in zip(indexes, smiles_list, amounts_list): 
            validate_dict = checkSMILES(smi, index, validate_dict)
            validate_dict = checkIsNumber(ammount, index, validate_dict)

    if len(validate_dict['warning_string']) !=0:
        validated = False
    
    # Delete tempory file if only validate selected
    if validate_only:
        default_storage.delete(csv_fp)
        csv_fp = None

    # Convert dataframe to dictionary to make it JSON serializable
    uploaded_dict = uploaded_csv_df.to_dict('list') 
      
    return (validate_dict, validated, csv_fp, project_info, uploaded_dict)

    # Functions to request info from IBM API
    
@shared_task
def uploadIBMReaction(validate_output):
    # Validate output is a list - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    validate_dict, validated, csv_fp, project_info, uploaded_dict = validate_output
    
    if not validated:
        # Delete tempory file if only validate selected
        default_storage.delete(csv_fp)
        return (validate_dict, validated)
    
    if validated:
        # Create project model and return project id
        project_id, project_name  = createProjectModel(project_info)
        
        # Create an IBM project
        IBM_project_id, rxn4chemistry_wrapper  = createIBMProject(project_name)

        target_no = 1
        for smiles, amount in zip(uploaded_dict['Targets'], uploaded_dict['Ammount_required (mg)']):

            # Number of steps/name of reactants from IBM API???????
            # Create a Target model and return target id
            target_id = createTargetModel(project_id = project_id, smiles=smiles, target_no=target_no)
            
            # Create IBM 
            # Run IBM API cal to get retrosyn info
            max_steps = 3
            results = getIBMRetroSyn(rxn4chemistry_wrapper=rxn4chemistry_wrapper, smiles=smiles, max_steps=max_steps)

            # Set maximum number of methods/pathways to collect
            no_pathways_found = len(results['retrosynthetic_paths'])

            if no_pathways_found < 3:
                max_pathways = no_pathways_found
            if no_pathways_found > 3:
                max_pathways = 3
            
            pathway_no = 1
            
            while pathway_no < max_pathways:
                for pathway in results['retrosynthetic_paths']:
                    # Get pathways with confidence above threshold
                    if pathway['confidence'] > 0.90:
                        # Create a Method model
                        method_id = createMethodModel(target_id=target_id, smiles=smiles, max_steps=max_steps, amount=amount)
                        
                        # Get reaction info about pathway
                        reaction_info = collectIBMReactionInfo(rxn4chemistry_wrapper=rxn4chemistry_wrapper, pathway=pathway)
                        
                        product_no = 1
                        for product_smiles, reaction_class, reactants, actions in zip(reaction_info['product_smiles'],reaction_info['rclass'],reaction_info['reactants'],reaction_info['actions']):
                            # Product_smiles and reaction class is a list of individual elements
                            # Reactants and actions is a list of lists    
                            
                            # Create a Reaction model
                            reaction_id = createReactionModel(method_id=method_id,reaction_class=reaction_class)

                            # Create a Product model
                            createProductModel(reaction_id=reaction_id, project_name=project_name, target_no=target_no, pathway_no=pathway_no, product_no=product_no, product_smiles=product_smiles)
                            
                            # Create Reactant models
                            reactant_no = 1
                            for reactant_smiles in reactants:
                                createReactantModel(reaction_id=reaction_id, project_name=project_name, target_no=target_no, pathway_no=pathway_no, product_no=product_no, reactant_no=reactant_no, reactant_smiles=reactant_smiles)
                                reactant_no += 1

                            # Create robotic actions for each reaction - link robo actions to known working methods????????
                            action_no = 1
                            for action in actions:
                                createActionModel(reaction_id=reaction_id, action_no=action_no, action=action)
                                action_no += 1
                            
                            product_no += 1                          

                    pathway_no += 1     

            target_no += 1
        
    # Delete tempory file
    default_storage.delete(csv_fp)
    
    return validate_dict, validated  