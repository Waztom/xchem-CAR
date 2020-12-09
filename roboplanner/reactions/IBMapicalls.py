from rxn4chemistry import RXN4ChemistryWrapper
import os
import time
from rdkit.Chem import AllChem

def createIBMProject(project_name):
    # Setup IBM RxN API
    api_key=os.environ['IBM_API_KEY'] 
    rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=api_key)
    rxn4chemistry_wrapper.create_project(project_name)
    IBM_project_id = rxn4chemistry_wrapper.project_id

    return IBM_project_id, rxn4chemistry_wrapper

def getIBMRetroSyn(rxn4chemistry_wrapper, smiles, max_steps):
        """
        Use the IBM API to get some possible retrosynthesis routes
        """
        # Create dummy dictionary to create while loop to catch when status is a SUCCESS
        results = {}
        results['status'] = None
        results['results'] = None

        while results['results'] is None:  
            try:
                time.sleep(30)
                response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(product=smiles, max_steps=max_steps)
                while results['status'] != 'SUCCESS': 
                    time.sleep(130)
                    results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(response['prediction_id'])
                    results['results'] = results
            except Exception as e:
                print(e)
        return results['results']

def collectIBMReactionInfo(rxn4chemistry_wrapper, pathway):
    reaction_info = {}
    reaction_info['rclass'] = []
    reaction_info['product_smiles'] = []
    reaction_info['reactants'] = []
    reaction_info['actions'] = []

    time.sleep(30)
    pathway_synthesis_response = rxn4chemistry_wrapper.create_synthesis_from_sequence(sequence_id=pathway['sequenceId'])
    pathway_synthesis_id = pathway_synthesis_response['synthesis_id']
    synthesis_tree, reactions, actions = rxn4chemistry_wrapper.get_synthesis_plan(synthesis_id=pathway_synthesis_id)

    for node in reactions:
        # NB gives list of actions for each reaction
        if len(node['actions']) > 0:
            reaction_info['actions'].append(node['actions'])

    if 'children' in pathway and len(pathway['children']):
        reaction_info['rclass'] .append(pathway['rclass'])
        reaction_info['product_smiles'] .append(pathway['smiles'])
        reaction_info['reactants'].append([node['smiles'] for node in pathway['children']])

    for node in pathway['children']:
        if 'children' in node and len(node['children']):
            reaction_info['rclass'].append(node['rclass'])
            reaction_info['product_smiles'].append(node['smiles'])
            reaction_info['reactants'].append([node['smiles'] for node in node['children']])

    return reaction_info
    
    
    
