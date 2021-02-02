from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
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

# import the logging library
import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)

# Check if reactant_smiles is a common solvent
# Get smiles from csv file!!!!!
common_solvents = {
    "CC(C)=O": "acetone",
    "CC#N": "acetonitrile",
    "c1ccccc1": "benzene",
    "CCCCO": "1-butanol",
    "CCC(C)O": "2-butanol",
    "CCC(C)=O": "2-butanone",
    "CC(C)(C)O": "t-butyl alcohol",
    "ClC(Cl)(Cl)Cl": "carbon tetrachloride",
    "Clc1ccccc1": "chlorobenzene",
    "ClC(Cl)Cl": "chloroform",
    "C1CCCCC1": "cyclohexane",
    "ClCCCl": "1,2-dichloroethane",
    "OCCOCCO": "diethylene glycol",
    "CCOCC": "diethyl ether",
    "COC": "dimethyl ether",
    "COCCOC": "1,2-dimethoxy-ethane",
    "CN(C)C=O": "dimethyl-formamide ",
    "C[S](C)=O": "dimethyl sulfoxide",
    "C1COCCO1": "1,4-dioxane",
    "CCO": "ethanol",
    "CCOC(C)=O": "ethyl acetate",
    "OCCO": "ethylene glycol",
    "OCC(O)CO": "glycerin",
    "CCCCCCC": "heptane",
    "CN(C)[P](=O)(N(C)C)N(C)C": "hexamethylphosphoramide",
    "CN(C)P(N(C)C)N(C)C": "hexamethylphosphoroustriamide",
    "CCCCCC": "hexane",
    "CO": "methanol",
    "COC(C)(C)C": "methyl t-butyl ether",
    "ClCCl": "methylene chloride",
    "CN1CCCC1=O": "N-methyl-2-pyrrolidinone",
    "C[N+]([O-])=O": "nitromethane",
    "CCCCC": "pentane",
    "CCCC(C)C": "isohexane",
    "CCCO": "1-propanol",
    "CC(C)O": "2-propanol",
    "C1CCOC1": "tetrahydrofuran",
    "Cc1ccccc1": "toluene",
    "O": "water",
    "Cc1ccccc1C": "o-xylene",
    "Cc1cccc(C)c1": "m-xylene",
    "Cc1ccc(C)cc1": "p-xylene",
}


def createSVGString(smiles):
    """
    Function that creates a SVG image string from smiles string

    target_name: string
        unique name of target
    smiles: string
        a valid smiles
    """
    mol = Chem.MolFromSmiles(smiles)

    # Initiate drawer and set size/font size
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(900, 200)
    drawer.SetFontSize(12)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg_string = drawer.GetDrawingText()

    return svg_string


def convertNameToSmiles(chemical_name):
    try:
        smiles = pcp.get_compounds("Pd/C", "name")[0].isomeric_smiles
        return smiles
    except:
        return False


def checkSMILES(smiles):
    # Sometimes IBM yields chemical name and not smiles
    # check if this is the case using rdkit mol and if not
    # use NIH Cactus Resolver tool to convert name to smiles
    mol = Chem.MolFromSmiles(smiles)

    if mol:
        return smiles
    if not mol:
        converted_smiles = convertNameToSmiles(smiles)
        return converted_smiles


def createProjectModel(project_info):
    # Function that creates a project object
    # if the csv file uploaded is validated and
    # the user wants to upload the data

    # project_info is a dictionary of info passed into the validate/upload
    # Celery tasks

    # Create Project object
    project = Project()
    project.submittername = project_info["submittername"]
    project.submitterorganisation = project_info["submitterorganisation"]
    project.submitteremail = project_info["submitteremail"]
    project.save()

    return project.id, project.name


def createTargetModel(project_id, smiles, target_no, target_mass):
    """
    Function that creates a Target object
    if the csv file uploaded is validated and
    the user wants to upload the data
    
    project_id: string
        id of project created for upload
    smiles: string
        a valid smiles
    """
    target = Target()
    project_obj = Project.objects.get(id=project_id)
    project_name = project_obj.name
    target.project_id = project_obj
    target.smiles = smiles
    target.name = "{}-{}".format(project_name, target_no)
    # Create and Write svg to file
    target_svg_string = createSVGString(smiles)
    target_svg_fn = default_storage.save(
        "targetimages/" + target.name + ".svg", ContentFile(target_svg_string)
    )
    target.image = target_svg_fn
    target.targetmass = target_mass
    target.unit = "mg"
    target.save()

    return target.id


def createMethodModel(target_id, smiles, max_steps):
    method = Method()
    target_obj = Target.objects.get(id=target_id)
    method.target_id = target_obj
    method.nosteps = max_steps
    method.save()

    return method.id


def createReactionModel(method_id, reaction_class):
    # Function that takes in all the info from IBM API call
    # and creates a reaction object

    # Create Reaction object
    reaction = Reaction()
    method_obj = Method.objects.get(id=method_id)
    reaction.method_id = method_obj
    reaction.reactionclass = reaction_class
    reaction.save()

    return reaction.id


def createProductModel(
    reaction_id, project_name, target_no, pathway_no, product_no, product_smiles
):
    product = Product()
    product.name = "{}-{}-{}-{}".format(project_name, target_no, pathway_no, product_no)
    reaction_obj = Reaction.objects.get(id=reaction_id)
    product.reaction_id = reaction_obj
    product.smiles = product_smiles
    product_svg_string = createSVGString(product_smiles)
    product_svg_fn = default_storage.save(
        "productimages/" + product.name + ".svg", ContentFile(product_svg_string)
    )
    product.image = product_svg_fn
    product.save()


def createAnalyseActionModel(reaction_id, action_no):
    analyse = AnalyseAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    analyse.reaction_id = reaction_obj
    analyse.actiontype = "analyse"
    analyse.actiontype = action_no
    analyse.save()


def createActionModel(reaction_id, action_no, action):
    # Create a dictionary of key (action name from IBM API) and
    # funtion name to create the appropriate model
    actionMethods = {
        "add": createIBMAddAction,
        "collect-layer": createIBMCollectLayerAction,
        "concentrate": createIBMConcentrateAction,
        "degas": createIBMDegasAction,
        "dry-solid": createIBMDrySolidAction,
        "dry-solution": createIBMDrySolutionAction,
        "extract": createIBMExtractAction,
        "filter": createIBMFilterAction,
        "make-solution": createIBMMakeSolutionAction,
        "partition": createIBMPartitionAction,
        "ph": createIBMpHAction,
        "phase-separation": createIBMPhaseSeparationAction,
        "quench": createIBMQuenchAction,
        "reflux": createIBMRefluxAction,
        "set-temperature": createIBMSetTemperatureAction,
        "stir": createIBMStirAction,
        "store": createIBMStoreAction,
        "wait": createIBMWaitAction,
        "wash": createIBMWashAction,
    }

    action_type = action["name"]

    if action_type == "add":
        actionMethods[action_name](action_type, reaction_id, action_no, action)
    else:
        logger.info(action_type)


def createIBMAddAction(action_type, reaction_id, action_no, action):
    material = action["content"]["material"]["value"]
    materialquantity = action["content"]["material"]["quantity"]["value"]
    materialquantityunit = action["content"]["material"]["quantity"]["unit"]
    dropwise = action["content"]["dropwise"]["value"]
    atmosphere = action["content"]["atmosphere"]["value"]

    add = IBMAddAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    add.reaction_id = reaction_obj
    add.actiontype = action_type
    add.actionno = action_no
    add.material = material
    # Check if smiles
    smiles = checkSMILES(material)
    if smiles:
        # Get MW
        mol = Chem.MolFromSmiles(smiles)
        molecular_weight = Descriptors.ExactMolWt(mol)

        add.materialsmiles = smiles
        add.molecularweight = molecular_weight
        add_svg_string = createSVGString(smiles)
        add_svg_fn = default_storage.save(
            "addactionimages/" + reaction_id + "_" + action_no + "_" + material + ".svg",
            ContentFile(add_svg_string),
        )
        add.materialimage = add_svg_fn
    # Check if solvent then use ml as quantity

    add.materialquantity = materialquantity
    add.materialquantityunit = materialquantityunit
    add.dropwise = dropwise

    if atmosphere:
        add.atmosphere = atmosphere
    else:
        add.atmosphere = "air"

    add.save()

