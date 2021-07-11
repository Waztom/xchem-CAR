from rdkit import Chem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

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

from ..utils import checkSMARTSPattern, getChemicalName, createSVGString

# import the logging library
import logging

# Get an instance of a logger
logger = logging.getLogger(__name__)


class CreateEncodedActionModels(object):
    """
    Creates a createEncodedActionModels object for creating action models
    for a reaction
    """

    def __init__(self, actions: list, target_id: int, reaction_id: int, reactant_pair_smiles: list):
        """
        ValidateFile constructor
        Args:
            actions (list): List of actions
            reaction_id (int): Reaction model id for actions
            target_id (int): Target model id for reaction
            reactant_pair_smiles (list): List of reactant smiles
        """
        self.actions = actions
        self.reaction_id = reaction_id
        self.reaction_obj = Reaction.objects.get(id=reaction_id)
        self.reactant_pair_smiles = reactant_pair_smiles
        self.target_mols = Target.objects.get(id=target_id).targetmols

        for action in self.actions:
            self.createEncodedActionModel(action)

    def createEncodedActionModel(self, action):
        actionMethods = {
            "add": self.createEncodedAddAction,
            "stir": self.createEncodedStirAction,
        }

        action_type = action["name"]

        if action_type in actionMethods:
            actionMethods[action_type](action_type, action)
            return True
        else:
            logger.info(action_type)
            print(action)

    def calculateVolume(
        self, molar_eqv, conc_reagents=None, reactant_density=None, reactant_MW=None
    ):
        # NB need addition_order added to ncoded recipes - can we rather use SA score?
        mol_material = molar_eqv * self.target_mols
        if reactant_density:
            vol_material = ((mol_material * reactant_MW) / reactant_density) * 1e3
        else:
            vol_material = (mol_material / conc_reagents) * 1e6  # in uL
        return vol_material

    def createEncodedAddAction(self, action_type, action):
        try:
            if action["content"]["material"]["SMARTS"]:
                SMARTS_pattern = action["content"]["material"]["SMARTS"]
                for pattern in SMARTS_pattern:
                    for reactant in self.reactant_pair_smiles:
                        pattern_check = checkSMARTSPattern(reactant, pattern)
                        if pattern_check:
                            reactant_SMILES = reactant
            if action["content"]["material"]["SMILES"]:
                reactant_SMILES = action["content"]["material"]["SMILES"]

            action_no = action["content"]["action_no"]
            molar_eqv = action["content"]["material"]["quantity"]["value"]
            conc_reagents = action["content"]["material"]["concentration"]
            solvent = action["content"]["material"]["solvent"]

            add = IBMAddAction()
            add.reaction_id = self.reaction_obj
            add.actiontype = action_type
            add.actionno = action_no
            material = getChemicalName(reactant_SMILES)
            mol = Chem.MolFromSmiles(reactant_SMILES)
            molecular_weight = Descriptors.ExactMolWt(mol)
            add.materialsmiles = reactant_SMILES
            add.molecularweight = molecular_weight
            add_svg_string = createSVGString(reactant_SMILES)
            add_svg_fn = default_storage.save(
                "addactionimages/{}-{}-{}.svg".format(self.reaction_id, action_no, material),
                ContentFile(add_svg_string),
            )
            add.materialimage = add_svg_fn  # need material (common name)
            add.atmosphere = "air"

            if solvent is None:
                reactant_density = action["content"]["material"]["density"]
                mol = Chem.MolFromSmiles(reactant_SMILES)
                reactant_MW = Descriptors.ExactMolWt(mol)
                add.materialquantity = self.calculateVolume(
                    molar_eqv=molar_eqv, reactant_density=reactant_density, reactant_MW=reactant_MW
                )

            if solvent:
                add.materialquantity = self.calculateVolume(
                    molar_eqv=molar_eqv, conc_reagents=conc_reagents
                )

            add.save()

        except Exception as error:
            print(error)
            print(action)

    def createEncodedStirAction(self, action_type, action):
        try:
            action_no = action["content"]["action_no"]
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]

            stir = IBMStirAction()
            stir.reaction_id = self.reaction_obj
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