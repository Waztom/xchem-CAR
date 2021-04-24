from rxn4chemistry import RXN4ChemistryWrapper
import os
import time
from rdkit.Chem import AllChem
from rdkit import Chem
import requests
import json
from .common_solvents import common_solvents


def canonSmiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    canon_smiles = Chem.MolToSmiles(mol)
    return canon_smiles


def convertIBMNameToSmiles(chemical_name):
    try:
        api_key = os.environ["IBM_API_KEY"]
        data = [chemical_name]
        headers = {
            "Authorization": api_key,
            "Content-Type": "application/json",
            "Accept": "application/json",
        }
        url = "https://rxn.res.ibm.com/rxn/api/api/v1/actions/convert-material-to-smiles"
        r = requests.post(url=url, data=json.dumps(data), headers=headers, cookies={})
        response_dict = r.json()
        smiles = response_dict["payload"][chemical_name]
        return smiles
    except:
        return False


def createIBMProject(project_name):
    # Setup IBM RxN API
    api_key = os.environ["IBM_API_KEY"]
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
    results["status"] = None
    results["results"] = None

    while results["results"] is None:
        try:
            time.sleep(10)
            response = rxn4chemistry_wrapper.predict_automatic_retrosynthesis(
                product=smiles, max_steps=max_steps
            )
            while results["status"] != "SUCCESS":
                time.sleep(130)
                results = rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(
                    response["prediction_id"]
                )
                results["results"] = results
        except Exception as e:
            return None
    return results["results"]


def collect_actions(tree):
    actions = []

    if "children" in tree and len(tree["children"]):
        actions.append(tree["actions"])

    for node in tree["children"]:
        actions.extend(collect_actions(node))

    return actions


def collect_rclass(tree):
    rclass = []
    if "children" in tree and len(tree["children"]):
        rclass.append(tree["rclass"])

    for node in tree["children"]:
        rclass.extend(collect_rclass(node))
    return rclass


def collect_products(tree):
    products = []
    if "children" in tree and len(tree["children"]):
        products.append(canonSmiles(tree["smiles"]))

    for node in tree["children"]:
        products.extend(collect_products(node))
    return products


def collect_reactants(tree):
    reactants = []
    if "children" in tree and len(tree["children"]):
        reactants.append([canonSmiles(node["smiles"]) for node in tree["children"]])

    for node in tree["children"]:
        reactants.extend(collect_reactants(node))
    return reactants


def collect_reactions(tree):
    reactions = []
    if "children" in tree and len(tree["children"]):
        reactions.append(
            AllChem.ReactionFromSmarts(
                "{}>>{}".format(
                    ".".join(
                        [
                            canonSmiles(node["smiles"])
                            for node in tree["children"]
                            if canonSmiles(node["smiles"]) not in common_solvents
                        ]
                    ),
                    tree["smiles"],
                ),
                useSmiles=True,
            )
        )

    for node in tree["children"]:
        reactions.extend(collect_reactions(node))
    return reactions


def collectIBMReactionInfo(rxn4chemistry_wrapper, pathway):
    reaction_info = {}

    try:
        time.sleep(10)
        pathway_synthesis_response = rxn4chemistry_wrapper.create_synthesis_from_sequence(
            sequence_id=pathway["sequenceId"]
        )
        pathway_synthesis_id = pathway_synthesis_response["synthesis_id"]
        time.sleep(10)
        synthesis_tree, reactions, actions = rxn4chemistry_wrapper.get_synthesis_plan(
            synthesis_id=pathway_synthesis_id
        )

        reaction_info["actions"] = collect_actions(synthesis_tree)
        # Can we join these into one loop?
        reaction_info["rclass"] = collect_rclass(pathway)
        reaction_info["product_smiles"] = collect_products(pathway)
        reaction_info["reactants"] = collect_reactants(pathway)
        reaction_info["reactions"] = collect_reactions(pathway)

        return reaction_info
    except:
        return None
