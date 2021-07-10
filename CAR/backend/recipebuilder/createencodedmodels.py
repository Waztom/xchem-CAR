from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

import pubchempy as pcp

import sys

sys.path.append("..")

# Import standard models
from ..models import Project, Target, Method, Reaction, Product, AnalyseAction

# Import IBM models
from ..models import (
    IBMAddAction,
    IBMCollectLayerAction,
    IBMConcentrateAction,
    IBMDegasAction,
    IBMDrySolidAction,
    IBMDrySolutionAction,
    IBMExtractAction,
    IBMFilterAction,
    IBMMakeSolutionAction,
    IBMPartitionAction,
    IBMpHAction,
    IBMPhaseSeparationAction,
    IBMQuenchAction,
    IBMRefluxAction,
    IBMSetTemperatureAction,
    IBMStirAction,
    IBMStoreAction,
    IBMWaitAction,
    IBMWashAction,
)

from ..IBM.createibmmodels import createSVGString

# import the logging library
import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)


# Do Harry's stuff here!
# 1. Calc expected mols product from target_mass
# NB need to write func that can do single plus multiple
# step calcs (Do we use SA score to estimate yield?)
# 2. Get recipe
# 2. Loop over recipe to populate relevant action models -> do
#    we need extra createmodels functions or is it possible to mix w
#    exisitng modelcreator?

# Collection of util function for creating models
def checkSMARTSPattern(SMILES, SMARTS_pattern):
    """function which checks whether the SMILES contains SMARTS"""
    pattern = Chem.MolFromSmarts(SMARTS_pattern)
    mol = Chem.MolFromSmiles(SMILES)
    if mol.HasSubstructMatch(pattern):
        return True
    else:
        return False


def calculateVolume(molar_eqv, product_moles, conc_reagents):
    # NB need addition_order added to ncoded recipes - can we rather use SA score?
    # estimated_loss = (
    #     0.7 ** addition_order
    # )  # Assume 70% conversion for each step//This needs to come from
    # reaction number NOT addtion order.....
    mol_material = molar_eqv * product_moles
    vol_material = (mol_material / conc_reagents) * 1e6  # in uL
    return vol_material


def getChemicalName(smiles):
    try:
        name = pcp.get_compounds(smiles, "smiles")[0].iupac_name
        return name
    except:
        print("Pubchempy could not convert SMILES to a IUPAC name")
        return smiles


def createEncodedActionModel(reaction_id, action, reactants, target_id):
    # Create a dictionary of key (action name from encoded recipe) and
    # funtion name to create the appropriate model
    actionMethods = {
        "add": createEncodedAddAction,
        "stir": createEncodedStirAction,
    }

    action_type = action["name"]

    if action_type in actionMethods:
        if action_type == "add":
            actionMethods[action_type](action_type, reaction_id, action, reactants, target_id)
        else:
            actionMethods[action_type](action_type, reaction_id, action)
        return True
    else:
        logger.info(action_type)
        print(action)


def createEncodedAddAction(action_type, reaction_id, action, reactants, target_id):
    try:
        # Does the action have a SMARTS pattern?
        if action["content"]["material"]["SMARTS"]:
            SMARTS_pattern = action["content"]["material"]["SMARTS"]
            for pattern in SMARTS_pattern:
                for reactant in reactants:
                    pattern_check = checkSMARTSPattern(reactant, pattern)
                    if pattern_check:
                        reactant_SMILES = reactant
        if action["content"]["material"]["SMILES"]:
            reactant_SMILES = action["content"]["material"]["SMILES"]

        action_no = action["content"]["action_no"]
        molar_eqv = action["content"]["material"]["quantity"]["value"]
        conc_reagents = action["content"]["material"]["concentration"]

        add = IBMAddAction()
        reaction_obj = Reaction.objects.get(id=reaction_id)
        add.reaction_id = reaction_obj
        add.actiontype = action_type
        add.actionno = action_no
        # NB function to convert reactant SMILES to IUPAC/common name
        # For now use material smiles
        material = getChemicalName(reactant_SMILES)
        mol = Chem.MolFromSmiles(reactant_SMILES)
        molecular_weight = Descriptors.ExactMolWt(mol)
        add.materialsmiles = reactant_SMILES
        add.molecularweight = molecular_weight
        add_svg_string = createSVGString(reactant_SMILES)
        add_svg_fn = default_storage.save(
            "addactionimages/{}-{}-{}.svg".format(reaction_id, action_no, material),
            ContentFile(add_svg_string),
        )
        add.materialimage = add_svg_fn  # need material (common name)
        target_obj = Target.objects.get(id=target_id)
        target_mols = target_obj.targetmols
        add.materialquantity = calculateVolume(molar_eqv, target_mols, conc_reagents)
        add.atmosphere = "air"
        add.save()

    except Exception as error:
        print(error)
        print(action)


def createEncodedStirAction(action_type, reaction_id, action):
    try:
        action_no = action["content"]["action_no"]
        duration = action["content"]["duration"]["value"]
        durationunit = action["content"]["duration"]["unit"]
        temperature = action["content"]["temperature"]

        stir = IBMStirAction()
        reaction_obj = Reaction.objects.get(id=reaction_id)
        stir.reaction_id = reaction_obj
        stir.actiontype = action_type
        stir.actionno = action_no
        stir.duration = duration
        stir.durationunit = durationunit
        stir.temperature = temperature
        stir.save()

    except Exception as error:
        print(action_type)
        print(error)
        print(action)