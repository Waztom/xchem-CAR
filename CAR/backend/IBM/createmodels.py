"""Create Django models from IBM API output"""
from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
import pubchempy as pcp
import os

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

from ..utils import calculateproductmols, createSVGString, createReactionSVGString, checkSMILES

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

    project_id: string id of project created for upload
    smiles: string a valid smiles
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


def createReactionModel(method_id, reaction_class, reaction_temperature, reaction_smarts):
    reaction = Reaction()
    method_obj = Method.objects.get(id=method_id)
    reaction.method_id = method_obj
    reaction.reactionclass = reaction_class
    if reaction_temperature:
        reaction.reactiontemperature = reaction_temperature
    reaction_svg_string = createReactionSVGString(reaction_smarts)
    reaction_svg_fn = default_storage.save(
        "reactionimages/" + reaction_class + ".svg", ContentFile(reaction_svg_string)
    )
    reaction.reactionimage = reaction_svg_fn
    reaction.save()

    return reaction.id


def createProductModel(reaction_id, project_name, target_no, method_no, product_no, product_smiles):
    product = Product()
    product.name = "{}-{}-{}-{}".format(project_name, target_no, method_no, product_no)
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
        self.action_no = 1

        for action in self.actions:
            created_action = self.createActionModel(action=action)
            if created_action:
                self.action_no += 1

    def createActionModel(self, action):
        actionMethods = {
            "add": self.createIBMAddActionModel,
            "collect-layer": self.createIBMCollectLayerActionModel,
            "concentrate": self.createIBMConcentrateActionModel,
            "degas": self.createIBMDegasActionModel,
            "dry-solid": self.createIBMDrySolidActionModel,
            "dry-solution": self.createIBMDrySolutionActionModel,
            "extract": self.createIBMExtractActionModel,
            "filter": self.createIBMFilterActionModel,
            "make-solution": self.createIBMMakeSolutionActionModel,
            "partition": self.createIBMPartitionActionModel,
            "ph": self.createIBMpHActionModel,
            "phase-separation": self.createIBMPhaseSeparationActionModel,
            "quench": self.createIBMQuenchActionModel,
            "reflux": self.createIBMRefluxActionModel,
            "set-temperature": self.createIBMSetTemperatureActionModel,
            "stir": self.createIBMStirActionModel,
            "store": self.createIBMStoreActionModel,
            "wait": self.createIBMWaitActionModel,
            "wash": self.createIBMWashActionModel,
        }

        action_type = action["name"]

        if action_type in actionMethods:
            actionMethods[action_type](action_type, action)
            return True
        else:
            logger.info(action_type)
            print(action)

    def calculateproductmols(self, target_mass, target_SMILES):
        target_MW = Descriptors.ExactMolWt(Chem.MolFromSmiles(target_SMILES))
        target_mass = target_mass / 1e3
        product_moles = target_mass / target_MW
        return product_moles

    def createAnalyseActionModel(self):
        analyse = AnalyseAction()
        analyse.reaction_id = self.reaction_obj
        analyse.actiontype = "analyse"
        analyse.actionno = self.action_no
        analyse.save()

    def createIBMAddActionModel(self, action_type, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            dropwise = action["content"]["dropwise"]["value"]
            atmosphere = action["content"]["atmosphere"]

            add = IBMAddAction()
            add.reaction_id = self.reaction_obj
            add.actiontype = action_type
            add.actionno = self.action_no
            add.material = material
            if material == "SLN":
                add.materialquantity = materialquantity
                add.materialquantityunit = materialquantityunit
            if material != "SLN":
                smiles = checkSMILES(material)
                if not smiles:
                    add.materialquantity = materialquantity
                    add.materialquantityunit = materialquantityunit
                if smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    molecular_weight = Descriptors.ExactMolWt(mol)

                    add.materialsmiles = smiles
                    add.molecularweight = molecular_weight
                    add_svg_string = createSVGString(smiles)
                    add_svg_fn = default_storage.save(
                        "addactionimages/{}-{}-{}.svg".format(
                            self.reaction_id, self.action_no, material
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

    def createIBMCollectLayerActionModel(self, action_type, action):
        try:
            layer = action["content"]["layer"]["value"]

            collect = IBMCollectLayerAction()
            collect.reaction_id = self.reaction_obj
            collect.actiontype = action_type
            collect.actionno = self.action_no
            collect.layer = layer
            collect.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMConcentrateActionModel(self, action_type, action):
        try:
            concentrate = IBMConcentrateAction()
            concentrate.reaction_id = self.reaction_obj
            concentrate.actiontype = action_type
            concentrate.actionno = self.action_no
            concentrate.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMDegasActionModel(self, action_type, action):
        try:
            gas = action["content"]["gas"]["value"]
            duration = action["content"]["duration"]["value"]
            unit = action["content"]["duration"]["unit"]

            degas = IBMDegasAction()
            degas.reaction_id = self.reaction_obj
            degas.actiontype = action_type
            degas.actionno = self.action_no
            degas.gas = gas
            degas.duration = duration
            degas.durationunit = unit
            degas.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMDrySolidActionModel(self, action_type, action):
        try:
            temperature = action["content"]["gas"]["value"]
            duration = action["content"]["duration"]["value"]
            unit = action["content"]["duration"]["unit"]
            atmosphere = action["content"]["atmosphere"]

            drysolid = IBMDrySolidAction()
            drysolid.reaction_id = self.reaction_obj
            drysolid.actiontype = action_type
            drysolid.actionno = self.action_no
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

    def createIBMDrySolutionActionModel(self, action_type, action):
        try:
            dryingagent = action["content"]["material"]["value"]

            drysolution = IBMDrySolutionAction()
            drysolution.reaction_id = self.reaction_obj
            drysolution.actiontype = action_type
            drysolution.actionno = self.action_no
            drysolution.dryingagent = dryingagent
            drysolution.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMExtractActionModel(self, action_type, action):
        try:
            solvent = action["content"]["solvent"]["value"]
            quantity = action["content"]["solvent"]["quantity"]["value"]
            unit = action["content"]["solvent"]["quantity"]["unit"]
            repetitions = action["content"]["repetitions"]["value"]

            extract = IBMExtractAction()
            extract.reaction_id = self.reaction_obj
            extract.actiontype = action_type
            extract.actionno = self.action_no
            extract.solvent = solvent
            extract.solventquantity = quantity
            extract.solventquantityunit = unit
            extract.numberofrepetitions = repetitions
            extract.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMFilterActionModel(self, action_type, action):
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
            filteraction.actionno = self.action_no
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

    def createIBMMakeSolutionActionModel(self, action_type, action):
        try:
            materials = action["content"]["materials"]["value"]
            solute = materials[0]["value"]
            solutesmiles = checkSMILES(solute)
            solvent = materials[1]["value"]
            solventsmiles = checkSMILES(solvent)
            solutequantity = materials[0]["quantity"]["value"]
            solutequantityunit = materials[0]["quantity"]["unit"]
            solventquantity = materials[1]["quantity"]["value"]
            solventquantityunit = materials[1]["quantity"]["unit"]

            soln = IBMMakeSolutionAction()
            soln.reaction_id = self.reaction_obj
            soln.actiontype = action_type
            soln.actionno = self.action_no
            soln.solute = solute
            if solutesmiles:
                soln.solutesmiles = solutesmiles
                soln_svg_string = createSVGString(solutesmiles)
                soln_svg_fn = default_storage.save(
                    "IBMmakesolnimages/{}-{}-{}.svg".format(
                        self.reaction_id, self.action_no, solute
                    ),
                    ContentFile(soln_svg_string),
                )
                soln.soluteimage = soln_svg_fn
            soln.solvent = solvent
            if solventsmiles:
                soln.solventsmiles = solventsmiles
                soln_svg_string = createSVGString(solventsmiles)
                soln_svg_fn = default_storage.save(
                    "IBMmakesolnimages/{}-{}-{}.svg".format(
                        self.reaction_id, self.action_no, solute
                    ),
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

    def createIBMPartitionActionModel(self, action_type, action):
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
            parition.actionno = self.action_no
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

    def createIBMpHActionModel(self, action_type, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            dropwise = action["content"]["dropwise"]["value"]
            temperature = action["content"]["temperature"]["value"]

            pH = IBMpHAction()
            pH.reaction_id = self.reaction_obj
            pH.actiontype = action_type
            pH.actionno = self.action_no
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

    def createIBMPhaseSeparationActionModel(self, action_type, action):
        try:
            phase = IBMPhaseSeparationAction()
            phase.reaction_id = self.reaction_obj
            phase.actiontype = action_type
            phase.actionno = self.action_no
            phase.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMQuenchActionModel(self, action_type, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            dropwise = action["content"]["dropwise"]["value"]
            temperature = action["content"]["temperature"]

            quench = IBMQuenchAction()
            quench.reaction_id = self.reaction_obj
            quench.actiontype = action_type
            quench.actionno = self.action_no
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

    def createIBMRefluxActionModel(self, action_type, action):
        try:
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            stirringspeed = action["content"]["stirring_speed"]["value"]
            deanstarkapparatus = action["content"]["dean_stark"]["value"]
            atmosphere = action["content"]["atmosphere"]

            reflux = IBMRefluxAction()
            reflux.reaction_id = self.reaction_obj
            reflux.actiontype = action_type
            reflux.actionno = self.action_no
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

    def createIBMSetTemperatureActionModel(self, action_type, action):
        try:
            temperature = action["content"]["temperature"]["value"]

            temp = IBMSetTemperatureAction()
            temp.reaction_id = self.reaction_obj
            temp.actiontype = action_type
            temp.actionno = self.action_no
            temp.temperature = temperature
            temp.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMStirActionModel(self, action_type, action):
        try:
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]
            stirringspeed = action["content"]["stirring_speed"]["value"]
            atmosphere = action["content"]["atmosphere"]

            stir = IBMStirAction()
            stir.reaction_id = self.reaction_obj
            stir.actiontype = action_type
            stir.actionno = self.action_no
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

    def createIBMStoreActionModel(self, action_type, action):
        try:
            material = action["content"]["sample_name"]["value"]

            store = IBMStoreAction()
            store.reaction_id = self.reaction_obj
            store.actiontype = action_type
            store.actionno = self.action_no
            store.material = material
            smiles = checkSMILES(material)
            if smiles:
                store.materialsmiles = smiles
                store_svg_string = createSVGString(smiles)
                store_svg_fn = default_storage.save(
                    "storeimages/{}-{}-{}.svg".format(self.reaction_id, self.action_no, material),
                    ContentFile(store_svg_string),
                )
                store.materialimage = store_svg_fn
            store.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMWaitActionModel(self, action_type, action):
        try:
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]["value"]

            wait = IBMWaitAction()
            wait.reaction_id = self.reaction_obj
            wait.actiontype = action_type
            wait.actionno = self.action_no
            wait.duration = duration
            wait.durationunit = durationunit
            wait.temperature = temperature
            wait.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)

    def createIBMWashActionModel(self, action_type, action):
        try:
            material = action["content"]["material"]["value"]
            materialquantity = action["content"]["material"]["quantity"]["value"]
            materialquantityunit = action["content"]["material"]["quantity"]["unit"]
            numberofrepetitions = action["content"]["repetitions"]["value"]

            wash = IBMWashAction()
            wash.reaction_id = self.reaction_obj
            wash.actiontype = action_type
            wash.actionno = self.action_no
            wash.material = material
            wash.materialquantity = materialquantity
            wash.materialquantityunit = materialquantityunit
            wash.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)
