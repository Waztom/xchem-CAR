from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from reactions.models import (
    Project, 
    Target, 
    Method, 
    Reaction, 
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
    svg_string =  drawer.GetDrawingText()
    
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

    return project.id


def createTargetModel(project_id, smiles, slug_no):
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
    target.name = project_name + '-{}'.format(slug_no)
    # Create and Write svg to file
    svg_string = createSVGString(smiles)
    svg_fn = default_storage.save('targetimages/' + target.name, ContentFile(svg_string)) 
    target.image = svg_fn       
    target.save()

def createReactionModel(smiles):
    # Function that takes in all the info from IBM API call
    # and creates a reaction object

    # Create Reaction object
    reaction = Reaction()
    # Need the method_id which is passed in as obj in argument
    # Need to have uploaded svg image to productimages folder
    method_obj = Method.objects.get(id=method_id)
    reaction.method_id = method_obj
    reaction.name = name
    reaction.productsmiles = product_smiles
    reaction.productimage = 'productimages/' + product_image_fn
    reaction.productname = product_name
    reaction.save()


    
