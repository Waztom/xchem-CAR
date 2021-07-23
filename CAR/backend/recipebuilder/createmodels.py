from rdkit import Chem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

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
        project_id: int,
        target_id: int,
        reaction_id: int,
        reactant_pair_smiles: list,
    ):
        """
        ValidateFile constructor
        Args:
            actions (list): List of actions
            project_id (int): Project model id
            reaction_id (int): Reaction model id for actions
            target_id (int): Target model id for reaction
            reactant_pair_smiles (list): List of reactant smiles
        """
        self.actions = actions
        self.project_obj = Project.objects.get(id=project_id)
        self.reaction_id = reaction_id
        self.reaction_obj = Reaction.objects.get(id=reaction_id)
        self.reactant_pair_smiles = reactant_pair_smiles
        self.target_mols = Target.objects.get(id=target_id).targetmols
        self.mculeapi = MCuleAPI()
        self.mculeidlist = []

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

    def createEncodedAddActionModel(self, action_type, action):
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
            concentration = action["content"]["material"]["concentration"]
            solvent = action["content"]["material"]["solvent"]

            add = IBMAddAction()
            add.reaction_id = self.reaction_obj
            add.actiontype = action_type
            add.actionno = action_no
            material = getChemicalName(reactant_SMILES)
            add.material = material
            mol = Chem.MolFromSmiles(reactant_SMILES)
            molecular_weight = Descriptors.ExactMolWt(mol)
            add.materialsmiles = reactant_SMILES
            mculeinfo = self.mculeapi.getMCuleInfo(smiles=reactant_SMILES)
            if mculeinfo:
                mculeid = mculeinfo[0]
                self.mculeidlist.append(mculeid)
                add.mculeid = mculeid
                add.mculeurl = mculeinfo[1]
                priceinfo = self.mculeapi.getMCulePrice(mculeid=mculeid, amount=10)
                if priceinfo:
                    add.mculeprice = priceinfo[0]
                    add.mculedeliverytime = priceinfo[1]
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
                    molar_eqv=molar_eqv, conc_reagents=concentration
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

    def createMculeQuoteModel(self):
        quote_info = self.mculeapi.getTotalQuote(mculeids=self.mculeidlist)

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


{
    "response": {
        "id": 36301,
        "state_display": "Done",
        "api_url": "https://mcule.com/api/v1/iquote-queries/36301/",
        "site_url": "https://mcule.com/quote/query/36301/",
        "group": {
            "id": 28234,
            "quotes": [
                {
                    "id": 42992,
                    "reference_id_full": "I-42992",
                    "valid_until": "2021-08-20T16:32:40.336",
                    "po_number": None,
                    "is_saved": False,
                    "is_expired": False,
                    "payment_due_days": 30,
                    "avg_product_price": 50.830000000000005,
                    "total_cost": 218.9,
                    "total_cost_without_discount": 218.9,
                    "type_display": "Best price",
                    "state_display": "Displayed",
                    "api_url": "https://mcule.com/api/v1/iquotes/42992/",
                    "site_url": "https://mcule.com/quote/I-42992/",
                    "created": "2021-07-21T16:32:40.336440",
                    "updated": "2021-07-21T16:32:40.336456",
                    "name": "",
                    "description": "",
                    "type": 10,
                    "coverage_percent": 100,
                    "price": "218.90",
                    "products_price": "152.49",
                    "delivery_price": "39.00",
                    "shipping_price": "27.41",
                    "reformatting_price": "0.00",
                    "special_formatting_price": None,
                    "analytical_services_price": "0.00",
                    "total_discount_price": "0.00",
                    "product_discount_price": "0.00",
                    "extra_handling_price": "0.00",
                    "discount": 0,
                    "delivery_days": 9,
                    "suppliers_count": 1,
                    "state": 10,
                    "stock_with_amount_ratio": 0,
                    "duplicate_count": 0,
                    "additional_document_notes": "",
                    "group": 28234,
                    "order_data": None,
                }
            ],
            "site_url": "https://mcule.com/quote/group/28234/",
            "created": "2021-07-21T16:32:40.330731",
            "updated": "2021-07-21T16:32:40.330750",
            "query": 36301,
        },
        "item_filters": {},
        "customer_email": "warren.thompson@diamond.ac.uk",
        "amount": 1,
        "min_amount": 1,
        "target_volume": None,
        "target_cc": None,
        "extra_amount": None,
        "min_extra_amount": None,
        "purity": None,
        "delivery_time": 21,
        "higher_amounts": False,
        "keep_original_salt_form": False,
        "keep_original_tautomer_form": False,
        "keep_original_stereo_form": False,
        "deliver_multiple_salt_forms": False,
        "additional_document_notes": "",
        "notes": "",
        "thoroughness": 10,
        "customer_first_name": "Warren",
        "customer_last_name": "Thompson",
        "delivery_country": "GB",
        "delivery_contact_person_name": "",
        "delivery_contact_person_email": "",
        "delivery_contact_person_phone": "",
        "delivery_post_code": "",
        "delivery_city": "",
        "delivery_address": "",
        "created": "2021-07-21T16:32:36.549086",
        "promo_code": "",
        "state": 30,
        "start_date": "2021-07-21T16:32:36.615390",
        "end_date": "2021-07-21T16:32:40.346162",
        "user": 20829,
        "scheme": None,
    }
}
