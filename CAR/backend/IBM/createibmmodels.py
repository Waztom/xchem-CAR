"""Create Django models from IBM API output"""
from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
import pubchempy as pcp
import os
import json
import requests

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


def createProjectModel(project_info):
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
    target.targetmols = calculateproductmols(target_mass, smiles)
    target.name = "{}-{}".format(project_name, target_no)
    target_svg_string = createSVGString(smiles)
    target_svg_fn = default_storage.save(
        "targetimages/" + target.name + ".svg", ContentFile(target_svg_string)
    )
    target.image = target_svg_fn
    target.targetmass = target_mass
    target.unit = "mg"
    target.save()

    return target.id


def createMethodModel(target_id, nosteps):
    method = Method()
    target_obj = Target.objects.get(id=target_id)
    method.target_id = target_obj
    method.nosteps = nosteps
    method.save()

    return method.id


def createReactionModel(method_id, reaction_class, reaction_smarts):
    reaction = Reaction()
    method_obj = Method.objects.get(id=method_id)
    reaction.method_id = method_obj
    reaction.reactionclass = reaction_class
    reaction_svg_string = createReactionSVGString(reaction_smarts)
    reaction_svg_fn = default_storage.save(
        "reactionimages/" + reaction_class + ".svg", ContentFile(reaction_svg_string)
    )
    reaction.reactionimage = reaction_svg_fn
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


class CreateIBMActionModels(object):
    """
    Creates a createIBMActions object for creating action models
    for a reaction
    """

    def __init__(self, actions: list, reaction_id: int):
        """
        ValidateFile constructor
        Args:
            actions (list): List of actions
            reaction_id (int): Reaction model id for actions
        """
        self.actions = actions
        self.api_key = os.environ["IBM_API_KEY"]
        self.reaction_id = reaction_id
        self.reaction_obj = Reaction.objects.get(id=reaction_id)
        self.action_numbers_list = range(1, len(actions), 1)

        for action, action_no in zip(self.actions, self.action_numbers_list):
            self.createActionModel(action_no=action_no, action=action)

    def createActionModel(self, action_no, action):
        actionMethods = {
            "add": self.createIBMAddAction,
            "collect-layer": self.createIBMCollectLayerAction,
            "concentrate": self.createIBMConcentrateAction,
            "degas": self.createIBMDegasAction,
            "dry-solid": self.createIBMDrySolidAction,
            "dry-solution": self.createIBMDrySolutionAction,
            "extract": self.createIBMExtractAction,
            "filter": self.createIBMFilterAction,
            "make-solution": self.createIBMMakeSolutionAction,
            "partition": self.createIBMPartitionAction,
            "ph": self.createIBMpHAction,
            "phase-separation": self.createIBMPhaseSeparationAction,
            "quench": self.createIBMQuenchAction,
            "reflux": self.createIBMRefluxAction,
            "set-temperature": self.createIBMSetTemperatureAction,
            "stir": self.createIBMStirAction,
            "store": self.createIBMStoreAction,
            "wait": self.createIBMWaitAction,
            "wash": self.createIBMWashAction,
        }

        action_type = action["name"]

        if action_type in actionMethods:
            actionMethods[action_type](action_type, action_no, action)
            return True
        else:
            logger.info(action_type)
            print(action)

    def convertIBMNameToSmiles(self, chemical_name):
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

    def calculateproductmols(self, target_mass, target_SMILES):
        target_MW = Descriptors.ExactMolWt(Chem.MolFromSmiles(target_SMILES))
        target_mass = target_mass / 1e3
        product_moles = target_mass / target_MW
        return product_moles

    def createSVGString(self, smiles):
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

    def createReactionSVGString(self, smarts):
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

    def convertNameToSmiles(self, chemical_name):
        try:
            smiles = pcp.get_compounds(chemical_name, "name")[0].isomeric_smiles
            return smiles
        except:
            try:
                smiles = pcp.get_compounds(chemical_name, "formula")[0].isomeric_smiles
                return smiles
            except:
                try:
                    smiles = self.convertIBMNameToSmiles(chemical_name)
                    return smiles
                except:
                    print("PubChemPy/IBM could not convert {}".format(chemical_name))
                    return False

    def checkSMILES(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return smiles
        if not mol:
            converted_smiles = self.convertNameToSmiles(smiles)
            return converted_smiles

    def createAnalyseActionModel(self, action_no):
        analyse = AnalyseAction()
        analyse.reaction_id = self.reaction_obj
        analyse.actiontype = "analyse"
        analyse.actiontype = action_no
        analyse.save()

    def createIBMAddAction(self, action_type, action_no, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            dropwise = action["content"]["dropwise"]["value"]
            atmosphere = action["content"]["atmosphere"]

            add = IBMAddAction()
            add.reaction_id = self.reaction_obj
            add.actiontype = action_type
            add.actionno = action_no
            add.material = material
            if material == "SLN":
                add.materialquantity = materialquantity
                add.materialquantityunit = materialquantityunit
            if material != "SLN":
                smiles = self.checkSMILES(material)
                if not smiles:
                    add.materialquantity = materialquantity
                    add.materialquantityunit = materialquantityunit
                if smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    molecular_weight = Descriptors.ExactMolWt(mol)

                    add.materialsmiles = smiles
                    add.molecularweight = molecular_weight
                    add_svg_string = self.createSVGString(smiles)
                    add_svg_fn = default_storage.save(
                        "addactionimages/{}-{}-{}.svg".format(
                            self.reaction_id, action_no, material
                        ),
                        ContentFile(add_svg_string),
                    )
                    add.materialimage = add_svg_fn
                    add.materialquantity = materialquantity
                    add.materialquantityunit = materialquantityunit
            if atmosphere:
                add.atmosphere = atmosphere["value"]
            else:
                add.atmosphere = "air"

            add.save()

        except Exception as error:
            print(error)
            print(action)

    def createIBMCollectLayerAction(self, action_type, action_no, action):
        try:
            layer = action["content"]["layer"]["value"]

            collect = IBMCollectLayerAction()
            collect.reaction_id = self.reaction_obj
            collect.actiontype = action_type
            collect.actionno = action_no
            collect.layer = layer
            collect.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMConcentrateAction(self, action_type, action_no, action):
        try:
            concentrate = IBMConcentrateAction()
            concentrate.reaction_id = self.reaction_obj
            concentrate.actiontype = action_type
            concentrate.actionno = action_no
            concentrate.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMDegasAction(self, action_type, action_no, action):
        try:
            gas = action["content"]["gas"]["value"]
            duration = action["content"]["duration"]["value"]
            unit = action["content"]["duration"]["unit"]

            degas = IBMDegasAction()
            degas.reaction_id = self.reaction_obj
            degas.actiontype = action_type
            degas.actionno = action_no
            degas.gas = gas
            degas.duration = duration
            degas.durationunit = unit
            degas.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMDrySolidAction(self, action_type, action_no, action):
        try:
            temperature = action["content"]["gas"]["value"]
            duration = action["content"]["duration"]["value"]
            unit = action["content"]["duration"]["unit"]
            atmosphere = action["content"]["atmosphere"]

            drysolid = IBMDrySolidAction()
            drysolid.reaction_id = self.reaction_obj
            drysolid.actiontype = action_type
            drysolid.actionno = action_no
            drysolid.temperature = temperature
            drysolid.duration = duration
            drysolid.durationunit = unit
            if atmosphere:
                drysolid.atmosphere = atmosphere["value"]
            drysolid.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMDrySolutionAction(self, action_type, action_no, action):
        try:
            dryingagent = action["content"]["material"]["value"]

            drysolution = IBMDrySolutionAction()
            drysolution.reaction_id = self.reaction_obj
            drysolution.actiontype = action_type
            drysolution.actionno = action_no
            drysolution.dryingagent = dryingagent
            drysolution.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMExtractAction(self, action_type, action_no, action):
        try:
            solvent = action["content"]["solvent"]["value"]
            quantity = action["content"]["solvent"]["quantity"]["value"]
            unit = action["content"]["solvent"]["quantity"]["unit"]
            repetitions = action["content"]["repetitions"]["value"]

            extract = IBMExtractAction()
            extract.reaction_id = self.reaction_obj
            extract.actiontype = action_type
            extract.actionno = action_no
            extract.solvent = solvent
            extract.solventquantity = quantity
            extract.solventquantityunit = unit
            extract.numberofrepetitions = repetitions
            extract.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMFilterAction(self, action_type, action_no, action):
        try:
            phasetokeep = action["content"]["phase_to_keep"]["value"]
            rinsingsolvent = action["content"]["rinsing_solvent"]["value"]
            rinsingsolventquantity = action["content"]["rinsing_solvent"]["quantity"]["value"]
            rinsingsolventquantityunit = action["content"]["rinsing_solvent"]["quantity"]["unit"]
            extractionforprecipitatesolvent = action["content"]["extraction_solvent"]["value"]
            extractionforprecipitatesolventquantity = action["content"]["extraction_solvent"][
                "quantity"
            ]["value"]
            extractionforprecipitatesolventquantityunit = action["content"]["extraction_solvent"][
                "quantity"
            ]["unit"]

            filteraction = IBMFilterAction()
            filteraction.reaction_id = self.reaction_obj
            filteraction.actiontype = action_type
            filteraction.actionno = action_no
            filteraction.phasetokeep = phasetokeep
            filteraction.rinsingsolvent = rinsingsolvent
            filteraction.rinsingsolventquantity = rinsingsolventquantity
            filteraction.rinsingsolventquantityunit = rinsingsolventquantityunit
            filteraction.extractionforprecipitatesolvent = extractionforprecipitatesolvent
            filteraction.extractionforprecipitatesolventquantity = (
                extractionforprecipitatesolventquantity
            )
            filteraction.extractionforprecipitatesolventquantityunit = (
                extractionforprecipitatesolventquantityunit
            )
            filteraction.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMMakeSolutionAction(self, action_type, action_no, action):
        try:
            materials = action["content"]["materials"]["value"]
            solute = materials[0]["value"]
            solutesmiles = self.checkSMILES(solute)
            solvent = materials[1]["value"]
            solventsmiles = self.checkSMILES(solvent)
            solutequantity = materials[0]["quantity"]["value"]
            solutequantityunit = materials[0]["quantity"]["unit"]
            solventquantity = materials[1]["quantity"]["value"]
            solventquantityunit = materials[1]["quantity"]["unit"]

            soln = IBMMakeSolutionAction()
            soln.reaction_id = self.reaction_obj
            soln.actiontype = action_type
            soln.actionno = action_no
            soln.solute = solute
            if solutesmiles:
                soln.solutesmiles = solutesmiles
                soln_svg_string = self.createSVGString(solutesmiles)
                soln_svg_fn = default_storage.save(
                    "IBMmakesolnimages/{}-{}-{}.svg".format(self.reaction_id, action_no, solute),
                    ContentFile(soln_svg_string),
                )
                soln.soluteimage = soln_svg_fn
            soln.solvent = solvent
            if solventsmiles:
                soln.solventsmiles = solventsmiles
                soln_svg_string = self.createSVGString(solventsmiles)
                soln_svg_fn = default_storage.save(
                    "IBMmakesolnimages/{}-{}-{}.svg".format(self.reaction_id, action_no, solute),
                    ContentFile(soln_svg_string),
                )
                soln.solventimage = soln_svg_fn

            soln.solutequantity = solutequantity
            soln.solutequantityunit = solutequantityunit
            soln.solventquantity = solventquantity
            soln.solventquantityunit = solventquantityunit
            soln.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMPartitionAction(self, action_type, action_no, action):
        try:
            firstpartitionsolvent = action["content"]["material_1"]["value"]
            secondpartitionsolvent = action["content"]["material_2"]["value"]
            firstpartitionsolventquantity = action["content"]["material_1"]["quantity"]["value"]
            secondpartitionsolventquantity = action["content"]["material_2"]["quantity"]["value"]
            firstpartitionsolventquantityunit = action["content"]["material_1"]["qauntity"]["unit"]
            secondpartitionsolventquantityunit = action["content"]["material_2"]["qauntity"]["unit"]

            parition = IBMPartitionAction()
            parition.reaction_id = self.reaction_obj
            parition.actiontype = action_type
            parition.actionno = action_no
            parition.firstpartitionsolvent = firstpartitionsolvent
            parition.firstpartitionsolventquantity = firstpartitionsolventquantity
            parition.firstpartitionsolventquantityunit = firstpartitionsolventquantityunit
            parition.secondpartitionsolvent = secondpartitionsolvent
            parition.secondpartitionsolventquantity = secondpartitionsolventquantity
            parition.secondpartitionsolventquantityunit = secondpartitionsolventquantityunit
            parition.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMpHAction(self, action_type, action_no, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            dropwise = action["content"]["dropwise"]["value"]
            temperature = action["content"]["temperature"]["value"]

            pH = IBMpHAction()
            pH.reaction_id = self.reaction_obj
            pH.actiontype = action_type
            pH.actionno = action_no
            pH.material = material
            pH.materialquantity = materialquantity
            pH.materialquantityunit = materialquantityunit
            pH.dropwise = dropwise
            pH.temperature = temperature
            pH.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMPhaseSeparationAction(self, action_type, action_no, action):
        try:
            phase = IBMPhaseSeparationAction()
            phase.reaction_id = self.reaction_obj
            phase.actiontype = action_type
            phase.actionno = action_no
            phase.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMQuenchAction(self, action_type, action_no, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            dropwise = action["content"]["dropwise"]["value"]
            temperature = action["content"]["temperature"]

            quench = IBMQuenchAction()
            quench.reaction_id = self.reaction_obj
            quench.actiontype = action_type
            quench.actionno = action_no
            quench.material = material
            quench.materialquantity = materialquantity
            quench.materialquantityunit = materialquantityunit
            if temperature:
                quench.temperature = temperature["value"]
            quench.dropwise = dropwise
            quench.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMRefluxAction(self, action_type, action_no, action):
        try:
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            stirringspeed = action["content"]["stirring_speed"]["value"]
            deanstarkapparatus = action["content"]["dean_stark"]["value"]
            atmosphere = action["content"]["atmosphere"]

            reflux = IBMRefluxAction()
            reflux.reaction_id = self.reaction_obj
            reflux.actiontype = action_type
            reflux.actionno = action_no
            reflux.duration = duration
            reflux.durationunit = durationunit
            reflux.stirringspeed = stirringspeed
            reflux.deanstarkapparatus = deanstarkapparatus
            if atmosphere:
                reflux.atmosphere = atmosphere["value"]
            reflux.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMSetTemperatureAction(self, action_type, action_no, action):
        try:
            temperature = action["content"]["temperature"]["value"]

            temp = IBMSetTemperatureAction()
            temp.reaction_id = self.reaction_obj
            temp.actiontype = action_type
            temp.actionno = action_no
            temp.temperature = temperature
            temp.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMStirAction(self, action_type, action_no, action):
        try:
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]
            stirringspeed = action["content"]["stirring_speed"]["value"]
            atmosphere = action["content"]["atmosphere"]

            stir = IBMStirAction()
            stir.reaction_id = self.reaction_obj
            stir.actiontype = action_type
            stir.actionno = action_no
            stir.duration = duration
            stir.durationunit = durationunit
            if temperature:
                stir.temperature = temperature["value"]
            else:
                stir.temperature = 25
            stir.stirringspeed = stirringspeed
            if atmosphere:
                stir.atmosphere = atmosphere["value"]
            else:
                stir.atmosphere = "air"
            stir.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMStoreAction(self, action_type, action_no, action):
        try:
            material = action["content"]["sample_name"]["value"]

            store = IBMStoreAction()
            store.reaction_id = self.reaction_obj
            store.actiontype = action_type
            store.actionno = action_no
            store.material = material
            smiles = self.checkSMILES(material)
            if smiles:
                store.materialsmiles = smiles
                store_svg_string = self.createSVGString(smiles)
                store_svg_fn = default_storage.save(
                    "storeimages/{}-{}-{}.svg".format(self.reaction_id, action_no, material),
                    ContentFile(store_svg_string),
                )
                store.materialimage = store_svg_fn
            store.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMWaitAction(self, action_type, action_no, action):
        try:
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]["value"]

            wait = IBMWaitAction()
            wait.reaction_id = self.reaction_obj
            wait.actiontype = action_type
            wait.actionno = action_no
            wait.duration = duration
            wait.durationunit = durationunit
            wait.temperature = temperature
            wait.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMWashAction(self, action_type, action_no, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            numberofrepetitions = action["content"]["repetitions"]["value"]

            wash = IBMWashAction()
            wash.reaction_id = self.reaction_obj
            wash.actiontype = action_type
            wash.actionno = action_no
            wash.material = material
            wash.materialquantity = materialquantity
            wash.materialquantityunit = materialquantityunit
            wash.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)
