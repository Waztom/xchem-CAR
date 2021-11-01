import json
import requests
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import pubchempy as pcp
import itertools


def calculateproductmols(target_mass, target_SMILES):
    target_MW = Descriptors.ExactMolWt(Chem.MolFromSmiles(target_SMILES))
    target_mass = target_mass / 1e3
    product_moles = target_mass / target_MW
    return product_moles


def canonSmiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    canon_smiles = Chem.MolToSmiles(mol)
    return canon_smiles


def combichem(reactant_1_SMILES: list, reactant_2_SMILES: list):
    """ "
    Gets all possible combinations between two uneven lists of
    reactants
    Args:
        reactant_1_SMILES (list): List of reactant one smiles
        reactant_2_SMILES (list): List of reactant two smiles
    Returns:
        all_possible_combinations (list): All possible combinations possible
                           between reactat 1 and reactant two lists
                           as a list of tuples
    """
    all_possible_combinations = list(itertools.product(reactant_1_SMILES, reactant_2_SMILES))

    return all_possible_combinations


def convertIBMNameToSmiles(chemical_name):
    try:
        data = [chemical_name]
        headers = {
            "Authorization": self.api_key,
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


def createSVGString(smiles):
    """
    Function that creates a SVG image string from smiles string

    target_name: string
        unique name of target
    smiles: string
        a valid smiles
    """
    mol = Chem.MolFromSmiles(smiles)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(100, 50)
    drawer.SetFontSize(8)
    drawer.SetLineWidth(1)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg_string = drawer.GetDrawingText()

    return svg_string


def createReactionSVGString(smarts):
    """
    Function that creates a SVG image string from smarts string

    target_name: string
        unique name of target
    smiles: string
        a valid smiles
    """
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(900, 200)
    drawer.DrawReaction(smarts)
    drawer.FinishDrawing()
    svg_string = drawer.GetDrawingText()
    return svg_string


def convertNameToSmiles(chemical_name):
    try:
        smiles = pcp.get_compounds(chemical_name, "name")[0].isomeric_smiles
        return smiles
    except:
        try:
            smiles = pcp.get_compounds(chemical_name, "formula")[0].isomeric_smiles
            return smiles
        except:
            try:
                smiles = convertIBMNameToSmiles(chemical_name)
                return smiles
            except:
                print("PubChemPy/IBM could not convert {}".format(chemical_name))
                return False


def checkSMILES(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return smiles
    if not mol:
        converted_smiles = convertNameToSmiles(smiles)
        return converted_smiles


def checkSMARTSPattern(SMILES, SMARTS_pattern):
    """function which checks whether the SMILES contains SMARTS"""
    pattern = Chem.MolFromSmarts(SMARTS_pattern)
    mol = Chem.MolFromSmiles(SMILES)
    if mol.HasSubstructMatch(pattern):
        return True
    else:
        return False


def getAddtionOrder(product_smi: str, reactant_SMILES: tuple, reaction_SMARTS: str):
    """
    Gets reactant pair addition order from reaction_smarts

    Args:
        product_smi (str): product SMILES
        reactant_SMILES_pair (tuple): Tuple of reactant smiles
        reaction_SMARTS (str): reaction SMARTS pattern

    Returns:
        reactant_SMILES_pair (list): List of reactant smiles in correct order
        None: If no match is found between the reactants and the reaction smarts
    """
    # Need to check if reaction works and then get corect order.
    rxn = AllChem.ReactionFromSmarts(reaction_SMARTS)
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactant_SMILES]

    for reactant_permutation in list(itertools.permutations(reactant_mols)):
        try:
            products = rxn.RunReactants(reactant_permutation)
            product_mols = [product[0] for product in products]
            if not product_mols:
                continue  # reactants were in wrong order so no product
        except Exception as e:
            print(e)
            print(reactant_permutation)
            continue
        product_smis = [Chem.MolToSmiles(m) for m in product_mols if m is not None]
        if product_smi in product_smis:
            ordered_smis = [Chem.MolToSmiles(m) for m in reactant_permutation]
    if "ordered_smis" in locals():
        return ordered_smis
    else:
        print(reaction_SMARTS)
        print(reactant_SMILES)
        return None


def checkReactantSMARTS(reactant_SMILES: tuple, reaction_SMARTS: str):
    """
    Checks if reactant pair can produce a product

    Args:
        reactant_SMILES_pair (tuple): Tuple of reactant smiles
        reaction_SMARTS (str): reaction SMARTS pattern

    Returns:
        products (rdkit obj): Rdkit object of products of reaction
        None: If no product mols formed
    """
    rxn = AllChem.ReactionFromSmarts(reaction_SMARTS)
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactant_SMILES]

    for reactant_permutation in list(itertools.permutations(reactant_mols)):
        try:
            products = rxn.RunReactants(reactant_permutation)
            product_mols = [product[0] for product in products]
            if product_mols:
                break
            if not product_mols:
                continue  # reactants were in wrong order so no product
        except Exception as e:
            print(e)
            print(reactant_permutation)
            continue
    if "product_mols" in locals():
        return product_mols
    else:
        print(reaction_SMARTS)
        print(reactant_SMILES)
        return None


def getChemicalName(smiles):
    try:
        name = pcp.get_compounds(smiles, "smiles")[0].iupac_name
        if not name:
            return None
        else:
            return name
    except:
        print("Pubchempy could not convert SMILES to a IUPAC name")
        return None