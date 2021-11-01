from rdkit import Chem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from statistics import mean

import sys
from ..mcule.apicalls import MCuleAPI

sys.path.append("..")


# Import standard models
from ..models import Project, MculeQuote, Target, Method, Reaction, Product, AnalyseAction

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

    def __init__(
        self,
        actions: list,
        target_id: int,
        reaction_id: int,
        reactant_pair_smiles: list,
        reaction_name: str,
    ):
        """
        ValidateFile constructor
        Args:
            actions (list): List of actions
            project_id (int): Project model id
            reaction_id (int): Reaction model id for actions
            target_id (int): Target model id for reaction
            reactant_pair_smiles (list): List of reactant smiles
            reaction_name (str): Reaction name
        """
        self.mculeapi = MCuleAPI()
        self.actions = actions
        self.reaction_id = reaction_id
        self.reaction_obj = Reaction.objects.get(id=reaction_id)
        self.reactant_pair_smiles = reactant_pair_smiles
        self.reaction_name = reaction_name
        self.target_mols = Target.objects.get(id=target_id).targetmols
        self.mculeidlist = []
        self.amountslist = []

        for action in self.actions:
            self.createEncodedActionModel(action)

    def createEncodedActionModel(self, action):
        actionMethods = {
            "add": self.createEncodedAddActionModel,
            "stir": self.createEncodedStirActionModel,
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

    def calculateAmount(self, molar_eqv: float, reactant_MW: float):
        """ "
        Calculates amount of compound needed in mg

        Args:
            molar_eq (float): Molar equivalents required
            reactant_MW (float): Molecular weight of compound
        Returns:
            amount (float): Amount required in mg
        """
        mol_material = molar_eqv * self.target_mols
        mass_material = (mol_material / reactant_MW) * 1e6
        return mass_material

    def createEncodedAddActionModel(self, action_type, action):
        try:
            if action["content"]["material"]["SMARTS"]:
                reactant_SMILES = self.reactant_pair_smiles[0]
                del self.reactant_pair_smiles[0]
                print(self.reactant_pair_smiles)
            if action["content"]["material"]["SMILES"]:
                reactant_SMILES = action["content"]["material"]["SMILES"]
            action_no = action["content"]["action_no"]
            molar_eqv = action["content"]["material"]["quantity"]["value"]
            concentration = action["content"]["material"]["concentration"]
            if not concentration:
                concentration = 0
            solvent = action["content"]["material"]["solvent"]
            mol = Chem.MolFromSmiles(reactant_SMILES)
            reactant_MW = Descriptors.ExactMolWt(mol)
            amount = self.calculateAmount(molar_eqv=molar_eqv, reactant_MW=reactant_MW)

            add = IBMAddAction()
            add.reaction_id = self.reaction_obj
            add.actiontype = action_type
            add.actionno = action_no
            # material = getChemicalName(reactant_SMILES)
            material = reactant_SMILES
            add.material = material
            # if not material:
            #     add.material = str(self.reaction_id) + str(action_no) + "-" + reactant_SMILES
            # else:
            #     add.material = material
            mol = Chem.MolFromSmiles(reactant_SMILES)
            molecular_weight = Descriptors.ExactMolWt(mol)
            add.materialsmiles = reactant_SMILES
            # mculeinfo = self.mculeapi.getMCuleInfo(smiles=reactant_SMILES)
            # if mculeinfo:
            #     mculeid = mculeinfo[0]
            #     self.mculeidlist.append(mculeid)
            #     self.amountslist.append(amount)
            #     add.mculeid = mculeid
            #     add.mculeurl = mculeinfo[1]
            #     priceinfo = self.mculeapi.getMCulePrice(mculeid=mculeid, amount=amount)
            #     if priceinfo:
            #         add.mculeprice = priceinfo[0]
            #         add.mculedeliverytime = priceinfo[1]
            add.molecularweight = molecular_weight
            add_svg_string = createSVGString(reactant_SMILES)
            add_svg_fn = default_storage.save(
                "addactionimages/{}-{}-{}.svg".format(self.reaction_id, action_no, material),
                ContentFile(add_svg_string),
            )
            add.materialimage = add_svg_fn  # need material (common name)
            add.atmosphere = "air"

            if not solvent:
                reactant_density = action["content"]["material"]["density"]
                add.materialquantity = self.calculateVolume(
                    molar_eqv=molar_eqv, reactant_density=reactant_density, reactant_MW=reactant_MW
                )

            if solvent:
                add.materialquantity = self.calculateVolume(
                    molar_eqv=molar_eqv,
                    conc_reagents=concentration,
                )
                add.solvent = solvent
            add.concentration = concentration
            add.save()

        except Exception as error:
            print(error)
            print(action)

    def createEncodedStirActionModel(self, action_type, action):
        try:
            action_no = action["content"]["action_no"]
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]["value"]
            temperatureunit = action["content"]["temperature"]["unit"]

            stir = IBMStirAction()
            stir.reaction_id = self.reaction_obj
            stir.actiontype = action_type
            stir.actionno = action_no
            stir.duration = duration
            stir.durationunit = durationunit
            stir.temperature = temperature
            stir.temperatureunit = temperatureunit
            stir.save()

        except Exception as error:
            print(action_type)
            print(error)
            print(action)


class CreateMculeQuoteModel(object):
    """
    Creates a CreateMculeQuoteModel object for creating a Mcule quote
    for a project
    """

    def __init__(
        self,
        mculeids: list,
        amounts: list,
        project_id: int,
    ):
        """
        ValidateFile constructor
        Args:
            mculeids (list): List of mcule ids
            project_id (int): Project model id
        """
        self.mculeidlist = [item for sublist in mculeids for item in sublist]
        self.amounts = amounts
        self.amountaverage = self.getAmountAverage()
        self.project_obj = Project.objects.get(id=project_id)
        self.mculeapi = MCuleAPI()
        self.createMculeQuoteModel()

    def getAmountAverage(self):
        if len(self.amounts) == 0:
            return 0
        if len(self.amounts) == 1:
            return self.amounts[0]
        else:
            mean([item for sublist in self.amounts for item in sublist])
            return mean

    def createMculeQuoteModel(self):
        quote_info = self.mculeapi.getTotalQuote(
            mculeids=self.mculeidlist, amount=self.amountaverage
        )

        if quote_info:
            try:
                quote = MculeQuote()
                quote.project_id = self.project_obj
                quote.quoteid = quote_info["quoteid"]
                quote.quoteurl = quote_info["quoteurl"]
                quote.quotecost = quote_info["quotecost"]
                quote.quotevaliduntil = quote_info["quotevaliduntil"]
                quote.save()

            except Exception as error:
                print(error)