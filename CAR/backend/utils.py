import json
import requests
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp


def calculateproductmols(target_mass, target_SMILES):
    target_MW = Descriptors.ExactMolWt(Chem.MolFromSmiles(target_SMILES))
    target_mass = target_mass / 1e3
    product_moles = target_mass / target_MW
    return product_moles


def canonSmiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    canon_smiles = Chem.MolToSmiles(mol)
    return canon_smiles


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


def getChemicalName(smiles):
    try:
        name = pcp.get_compounds(smiles, "smiles")[0].iupac_name
        return name
    except:
        print("Pubchempy could not convert SMILES to a IUPAC name")
        return smiles