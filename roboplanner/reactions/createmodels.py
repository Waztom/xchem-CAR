from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from .models import (
    Project, 
    Target, 
    Method, 
    Reaction,
    Product, 
    Reactant,
    AddAction,
    MakeSolutionAction,
    StirAction,
    WashAction,
    DrySolutionAction,
    ConcentrateAction,
    AnalyseAction)

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
    project.submittername = project_info['submittername']
    project.submitterorganisation =  project_info['submitterorganisation'] 
    project.submitteremail = project_info['submitteremail']
    project.save()

    return project.id, project.name


def createTargetModel(project_id, smiles, target_no):
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
    target.name = '{}-{}'.format(project_name, target_no)
    # Create and Write svg to file
    target_svg_string = createSVGString(smiles)
    target_svg_fn = default_storage.save('targetimages/' + target.name, ContentFile(target_svg_string)) 
    target.image = target_svg_fn     
    target.save()

    return target.id

def createMethodModel(target_id, smiles, max_steps, amount):
    method = Method()
    target_obj = Target.objects.get(id=target_id)
    method.target_id = target_obj
    method.nosteps = max_steps
    method.targetmass = amount
    method.unit = 'mg'
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


def createProductModel(reaction_id, project_name, target_no, pathway_no, product_no, product_smiles):    
    product = Product()
    product.name = '{}-{}-{}-{}'.format(project_name, target_no, pathway_no, product_no)
    reaction_obj = Reaction.objects.get(id=reaction_id)
    product.reaction_id = reaction_obj
    product.smiles = product_smiles
    product_svg_string = createSVGString(product_smiles)
    product_svg_fn = default_storage.save('productimages/' + product.name, ContentFile(product_svg_string))
    product.image = product_svg_fn
    product.save()


def createReactantModel(reaction_id, project_name, target_no, pathway_no, product_no, reactant_no, reactant_smiles):    
    reactant = Reactant()
    reactant.name = '{}-{}-{}-{}-{}'.format(project_name, target_no, pathway_no, product_no, reactant_no)
    reaction_obj = Reaction.objects.get(id=reaction_id)
    reactant.reaction_id = reaction_obj
    reactant.smiles = reactant_smiles
    reactant_svg_string = createSVGString(reactant_smiles)
    reactant_svg_fn = default_storage.save('reactantimages/' + reactant.name, ContentFile(reactant_svg_string))
    reactant.image = reactant_svg_fn
    reactant.save()

def createActionModel(reaction_id, action_no, action):
    # action is a JSON

    # Create a dictionary of key (action name from IBM API) and
    # funtion name to create the appropriate model
    actionMethods = {
        "add"           : createAddActionModel,
        "make-solution" : createMakeSolutionActionModel,
        "stir"          : createStirActionModel,
        "wash"          : createWashActionModel,
        "dry-solution"  : createDrySolutionActionModel,
        "concentrate"   : createConcentrateActionModel
    }

    action_name = action['name']

    if action_name in actionMethods:
        actionMethods[action_name](reaction_id, action_no, action)
    # else:
    #     raise Exception("Method {} not implemented".format(action_name))


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
        url= 'https://cactus.nci.nih.gov/chemical/structure/' + name_converted + '/smiles'
        ans = urlopen(url).read().decode('utf8')
        smiles = ans.split(' ')[0]
        return smiles
    except:
        return False
      

def createAddActionModel(reaction_id, action_no, action):    
    # Get info from action JSON
    # Check if smiles and convert chem name to smiles if not
    material = action['content']['material']['value']
    material_smiles_check = checkSMILES(material)

    if material_smiles_check: 
        quantity = action['content']['material']['quantity']['value']
        unit = action['content']['material']['quantity']['unit']
        dropwise_bool = action['content']['dropwise']['value'] 
        
        add = AddAction()
        reaction_obj = Reaction.objects.get(id=reaction_id)
        add.reaction_id = reaction_obj
        add.actionno = action_no
        add.material = material_smiles_check
        add.quantity = quantity
        add.unit = unit
        add.dropwise = dropwise_bool   
        add.save()

def createMakeSolutionActionModel(reaction_id, action_no, action):    
    # Get info from action JSON
    materials = action['content']['materials']['value']
    materials_smiles_check = [checkSMILES(material['value']) for material in materials]

    if all(materials_smiles_check):
        solute_smiles = materials_smiles_check[0]
        solvent_smiles = materials_smiles_check[1]

        quantities = [material['quantity']['value'] for material in materials]
        solute_quantity = quantities[0]
        solvent_quantity = quantities[1]

        units = [material['quantity']['unit'] for material in materials]
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

def createStirActionModel(reaction_id, action_no, action):
    # Get info from action JSON
    duration = action['content']['duration']['value']
    unit = action['content']['duration']['unit']

    stir = StirAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    stir.reaction_id = reaction_obj
    stir.actionno = action_no
    stir.duration = duration
    stir.unit = unit
    stir.save()


def createWashActionModel(reaction_id, action_no, action):
    # Get info from action JSON
    material = action['content']['material']['value']
    no_repetitions = action['content']['repetitions']['value']
    amount = action['content']['material']['quantity']['value']
    unit = action['content']['material']['quantity']['unit']

    # Update as we find more
    common_wash_materials = {
        "brine" : ['[Na+].[Cl-]', 'O'],
        "water" : ['O'],

    }

    wash = WashAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    wash.reaction_id = reaction_obj
    wash.actionno = action_no
    wash.material = material
    wash.norepetitions = no_repetitions
    wash.unit = unit
    wash.save()


def createDrySolutionActionModel(reaction_id, action_no, action):
    # Get info from action JSON
    material = action['content']['material']['value']
    dry = DrySolutionAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    dry.reaction_id = reaction_obj
    dry.actionno = action_no
    dry.dryingagent = material
    dry.save()


def createConcentrateActionModel(reaction_id, action_no, action):
    concentrate = ConcentrateAction()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    concentrate.reaction_id = reaction_obj
    concentrate.actionno = action_no
    concentrate.concentrate = True
    concentrate.save()



    
