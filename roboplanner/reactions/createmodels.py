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
    reactant_svg_fn = default_storage.save('reactantimages/', ContentFile(reactant_svg_string))
    reactant.image = reactant_svg_fn
    reactant.save()

def createActionModel(action):
    # action is a JSON

    # Create a dictionary of key (action name from IBM API) and
    # funtion name to create the appropriate model
    actionMethods = {
        "add"           : createAddActionModel,
        "make-solution" : createMakeSolutionActionModel,
        "stir"          : createStirActionModel,
        "wash"          : createWashActionModel,
        "dry-solution"  : createDrySolutionActionModel,
        "concentrate"   : createConcentrateActionModel,
        "analyze"       : createAnalyseActionModel
    }

    action_name = action['name']

    if action_name in actionMethods:
        actionMethods[action_name](action)
    else:
        raise Exception("Method {} not implemented".format(action_name))

def createAddActionModel(action):    
    add = 
    reactant.name = '{}-{}-{}-{}-{}'.format(project_name, target_no, pathway_no, product_no, reactant_no)
    reaction_obj = Reaction.objects.get(id=reaction_id)
    reactant.reaction_id = reaction_obj
    reactant.smiles = reactant_smiles
    reactant_svg_string = createSVGString(reactant_smiles)
    reactant_svg_fn = default_storage.save('reactantimages/', ContentFile(reactant_svg_string))
    reactant.image = reactant_svg_fn
    reactant.save()

# Need to add action models and create them - how to create depending on action name? 

    
