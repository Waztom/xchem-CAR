from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

import sys

sys.path.append("..")

from ..models import (
    Project,
    Target,
    Method,
    Reaction,
    Product,
    AddAction,
    MakeSolutionAction,
    StirAction,
    WashAction,
    DrySolutionAction,
    ConcentrateAction,
    AnalyseAction,
)

from urllib.request import urlopen
from urllib.parse import quote

# Certificate to NIH Cactus stopped working (10/12/2020)
# use this for devlopment purposes only by setting use
# of unverified certs - defo not for prod. This is not a
# cert issue on my machine, did update certs and have
# pointed Conda to fresh certs. Same thing when using on Binder
# at 'https://github.com/xchem/strucbio_practical' and
# 'https://stackoverflow.com/questions/33699577/conda-update-failed-ssl-error-ssl-certificate-verify-failed-certificate-ver'

import ssl

ssl._create_default_https_context = ssl._create_unverified_context

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


def createActionModel(reaction_id, action_no, action):
    # Create a dictionary of key (action name from IBM API) and
    # funtion name to create the appropriate model
    actionMethods = {
        "stir": createStirActionModel,
        "wash": createWashActionModel,
        "dry-solution": createDrySolutionActionModel,
        "concentrate": createConcentrateActionModel,
    }

    action_name = action["name"]

    if action_name in actionMethods:
        actionMethods[action_name](reaction_id, action_no, action)


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


def convertNameToSmiles(chemical_name):
    try:
        name_converted = quote(chemical_name)
        url = (
            "https://cactus.nci.nih.gov/chemical/structure/" + name_converted + "/smiles"
        )
        ans = urlopen(url).read().decode("utf8")
        smiles = ans.split(" ")[0]
        return smiles
    except:
        return False


def createAddActionModel(
    reaction_id,
    project_name,
    target_no,
    pathway_no,
    product_no,
    action_no,
    reactant_smiles,
):

    if reactant_smiles not in common_solvents:
        # Get MW
        mol = Chem.MolFromSmiles(reactant_smiles)
        molecular_weight = Descriptors.ExactMolWt(mol)

        add = AddAction()
        add.name = "{}-{}-{}-{}-{}".format(
            project_name, target_no, pathway_no, product_no, action_no
        )
        reaction_obj = Reaction.objects.get(id=reaction_id)
        add.reaction_id = reaction_obj
        add.actionno = action_no
        add.smiles = reactant_smiles
        add_svg_string = createSVGString(reactant_smiles)
        add_svg_fn = default_storage.save(
            "addactionimages/" + add.name + ".svg", ContentFile(add_svg_string)
        )
        add.image = add_svg_fn
        add.molecularweight = molecular_weight
        add.save()


def createStirActionModel(reaction_id, action_no, action):
    # Get info from action JSON
    duration = action["content"]["duration"]["value"]
    unit = action["content"]["duration"]["unit"]
    try:
        temperature = action["content"]["temperature"]["value"]
    except Exception as e:
        temperature = 25
        # Add error log or something

    stir = StirAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    stir.reaction_id = reaction_obj
    stir.actiontype = action["name"]
    stir.actionno = action_no
    stir.duration = duration
    stir.unit = unit
    stir.temperature = temperature
    stir.save()


def createWashActionModel(reaction_id, action_no, action):
    # Get info from action JSON
    material = action["content"]["material"]["value"]
    no_repetitions = action["content"]["repetitions"]["value"]
    amount = action["content"]["material"]["quantity"]["value"]
    unit = action["content"]["material"]["quantity"]["unit"]

    # Update as we find more
    common_wash_materials = {
        "brine": ["[Na+].[Cl-]", "O"],
        "water": ["O"],
    }

    wash = WashAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    wash.reaction_id = reaction_obj
    wash.actiontype = action["name"]
    wash.actionno = action_no
    wash.material = material
    wash.norepetitions = no_repetitions
    wash.unit = unit
    wash.save()


def createDrySolutionActionModel(reaction_id, action_no, action):
    # Get info from action JSON
    material = action["content"]["material"]["value"]
    dry = DrySolutionAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    dry.reaction_id = reaction_obj
    dry.actiontype = action["name"]
    dry.actionno = action_no
    dry.dryingagent = material
    dry.save()


def createConcentrateActionModel(reaction_id, action_no, action):
    concentrate = ConcentrateAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    concentrate.reaction_id = reaction_obj
    concentrate.actiontype = action["name"]
    concentrate.actionno = action_no
    concentrate.concentrate = True
    concentrate.save()


# Check if need these models?
def createMakeSolutionActionModel(reaction_id, action_no, action):
    # Get info from action JSON
    materials = action["content"]["materials"]["value"]
    materials_smiles_check = [checkSMILES(material["value"]) for material in materials]

    if all(materials_smiles_check):
        solute_smiles = materials_smiles_check[0]
        solvent_smiles = materials_smiles_check[1]

        quantities = [material["quantity"]["value"] for material in materials]
        solute_quantity = quantities[0]
        solvent_quantity = quantities[1]

        units = [material["quantity"]["unit"] for material in materials]
        solute_unit = units[0]
        solvent_unit = units[1]

        makesoln = MakeSolutionAction()
        reaction_obj = Reaction.objects.get(id=reaction_id)
        makesoln.reaction_id = reaction_obj
        makesoln.actionno = action_no
        makesoln.solute = solute_smiles
        makesoln.solutequantity = solute_quantity
        makesoln.soluteunit = solute_unit
        makesoln.solvent = solvent_smiles
        makesoln.solventequantity = solvent_quantity
        makesoln.solventunit = solvent_unit
        makesoln.save()

