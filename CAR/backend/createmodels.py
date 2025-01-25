"""Create Django models from Manifold API/custom chemistry outputs"""
from __future__ import annotations
import math
from typing import Tuple
import inspect
from rdkit import Chem
from rdkit.Chem import Descriptors
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from .mcule.apicalls import MCuleAPI

# Import standard models
from .models import (
    Project,
    Batch,
    PubChemInfo,
    ActionSession,
    Target,
    Method,
    Reaction,
    Product,
    Reactant,
    CatalogEntry,
)

# Import action models
from .models import (
    ActionSession,
    AddAction,
    ExtractAction,
    MixAction,
    StirAction,
)

from .utils import (
    getProduct,
    getProductSmiles,
    canonSmiles,
    calculateMassFromMols,
    checkProceedingReactions,
    calculateMolsFromConc,
    calculateMassFromMols,
    createSVGString,
    createReactionSVGString,
    getPubChemCAS,
    getPubChemCompound,
    getReactionYields,
    getInchiKey,
)

import logging

logger = logging.getLogger(__name__)


def createProjectModel(project_info: dict) -> Tuple[int, str]:
    """Creates a Django project object - a project model entry

    Parameters
    ----------
    project_info: dict
        The project info captured from the frontend project upload page

    Returns
    -------
    project.id: int
        The id of the project model object created
    project.name: str
        The name of the project model object created
    """
    project = Project()
    project.name = project_info["projectname"]
    project.submittername = project_info["submittername"]
    project.submitterorganisation = project_info["submitterorganisation"]
    project.proteintarget = project_info["proteintarget"]
    project.save()
    return project.id


def createBatchModel(project_id: int, batchtag: str, batch_id: int = None) -> int:
    """Creates a Django batch object - a batch of target compounds model entry

    Parameters
    ----------
    project_id: int
        The project model object id
    batchtag: str
        The tag or name used to create a batch
    batch_id: int
        Optional batch id, batch id included if batch is created from a previous parent batch.
        New parent batches will not have a batch id.

    Returns
    -------
    batch.id: int
        The id of the batch model object created
    """
    batch = Batch()
    project_obj = Project.objects.get(id=project_id)
    batch.project_id = project_obj
    if batch_id:
        fk_batch_obj = Batch.objects.get(pk=batch_id)
        batch.batch_id = fk_batch_obj
    batch.batchtag = batchtag
    batch.save()
    return batch.id


def createTargetModel(
    batch_id: int, name: str, smiles: str, concentration: float, volume: float
) -> int:
    """Creates a Django target object - a target compound model entry

    Parameters
    ----------
    batch_id: int
        The batch model object id that the target is linked to
    name: str
        The name of the Target compound
    smiles: str
        The SMILES of the target compound
    concentration: float
        The concentration (mM) to be made for the target compound
    volume: float
        The volume (uL) to be made for the target compound

    Returns
    -------
    target.id: int
        The id of the target model object created
    """
    target = Target()
    batch_obj = Batch.objects.get(id=batch_id)
    target.batch_id = batch_obj
    target.smiles = canonSmiles(smiles=smiles)
    mols = calculateMolsFromConc(
        target_concentration=concentration, target_volume=volume
    )
    mass = calculateMassFromMols(mols=mols, SMILES=smiles)
    target.mols = mols
    target.concentration = concentration
    target.volume = volume
    target.mass = mass
    target.name = name
    target_svg_string = createSVGString(smiles)
    target_svg_fn = default_storage.save(
        "targetimages/" + target.name + ".svg", ContentFile(target_svg_string)
    )
    target.image = target_svg_fn
    target.save()
    return target.id


def createMethodModel(target_id: int, nosteps: int, otchem: bool) -> int:
    """Creates a Django method object - a method is a collection of reactions
       for a target compound

    Parameters
    ----------
    target_id: int
        The target model object id the method is linked to
    nosteps: str
        The number of reaction steps in a method
    otchem: bool
        Set to True if all the reactions in a method can be executed on the
        OpenTrons

    Returns
    -------
    method.id: int
        The id of the method model object created
    """
    method = Method()
    target_obj = Target.objects.get(id=target_id)
    method.target_id = target_obj
    method.nosteps = nosteps
    method.otchem = otchem
    method.save()
    return method.id


def createReactionModel(
    method_id: int,
    reaction_class: str,
    reaction_number: int,
    intramolecular: bool,
    reaction_smarts: str,
    reaction_temperature: float = None,
    reaction_recipe: str = None,
    groupby_column: bool = True,
) -> int:
    """Creates a Django reaction object - a chemical reaction

    Parameters
    ----------
    method_id: int
        The method model object id the reaction is linked to
    reaction_class: str
        The name of the reaction eg. Buchwald-Hartwig amination
    reaction_number: int
        The number of the reaction
    intramolecular: bool
        Set to True if the reaction is intramolecular
    reaction_SMARTS: str
        The SMARTS for the reaction
    reaction_temperature: float
        The opotional reaction temperature
    reaction_recipe: str
        The optional (if found in encoded recipes) type of encoded recipe used to execute the reaction
    grouby_column: bool
        Set to True (default) to group reactions in columns by recipe type

    Returns
    -------
    reaction.id: int
        The id of the reaction model object created
    """
    reaction = Reaction()
    method_obj = Method.objects.get(id=method_id)
    reaction.method_id = method_obj
    reaction.reactionclass = reaction_class
    reaction.number = reaction_number
    reaction.intramolecular = intramolecular
    if reaction_temperature:
        reaction.temperature = reaction_temperature
    if reaction_recipe:
        reaction.recipe = reaction_recipe
    reaction.groupbycolumn = groupby_column
    reaction_svg_string = createReactionSVGString(reaction_smarts)
    reaction_svg_fn = default_storage.save(
        "reactionimages/" + reaction_class + ".svg", ContentFile(reaction_svg_string)
    )
    reaction.image = reaction_svg_fn
    reaction.save()

    return reaction.id


def createPubChemInfoModel(compoundid: int, smiles: str, cas: str = None) -> object:
    """Creates a Django pubcheminfo object - the PubChem info captured for a
    compound

    Parameters
    ----------
    compoundid: str
        The PubChem DB compound id
    smiles: str
        The SMILES of the compound
    cas: str
        The optional CAS number for the compound

    Returns
    -------
    pubcheminfo: object
        The PubChem model object created
    """
    pubcheminfo = PubChemInfo()
    pubcheminfo.smiles = canonSmiles(smiles=smiles)
    if cas:
        pubcheminfo.cas = cas
    pubcheminfo.compoundid = compoundid
    pubcheminfo.summaryurl = "https://pubchem.ncbi.nlm.nih.gov/compound/{}".format(
        compoundid
    )
    pubcheminfo.lcssurl = (
        "https://pubchem.ncbi.nlm.nih.gov/compound/{}#datasheet=LCSS".format(compoundid)
    )
    pubcheminfo.save()
    return pubcheminfo


def getPubChemInfo(smiles: str) -> object:
    """Searches if Django PubChemInfo object exists for smiles. If not checks if
    an entry exists on PuBChem and if it does, creates a PubChemInfo model object

    Parameters
    ----------
    smiles: str
        The SMILES of the compound

    Returns
    -------
    pubcheminfo: object
        The PubChem model object found or created
    status: bool
        Returns False if no Django PubChemInfo model instance found or
        PubChem DB entry found for the compound
    """
    smiles = canonSmiles(smiles=smiles)
    pubcheminfoqueryset = PubChemInfo.objects.filter(smiles=smiles)
    if pubcheminfoqueryset:
        pubcheminfo = pubcheminfoqueryset[0]
        return pubcheminfo
    else:
        inchikey = getInchiKey(smiles=smiles)
        compound = getPubChemCompound(inchikey=inchikey)
        if compound:
            compoundid = compound.cid
            cas = getPubChemCAS(compound=compound)
            if cas:
                pubcheminfo = createPubChemInfoModel(
                    compoundid=compoundid, smiles=smiles, cas=cas
                )
            if not cas:
                pubcheminfo = createPubChemInfoModel(
                    compoundid=compoundid, smiles=smiles
                )
            return pubcheminfo
        if not compound:
            return False


def createProductModel(reaction_id: int, product_smiles: str):
    """Creates a Django product object - the product of a reaction

    Parameters
    ----------
    reaction_id: int
        The reaction id the product is linked to
    product_smiles: str
        The SMILES of the product
    """
    product_smiles = canonSmiles(smiles=product_smiles)
    # pubcheminfoobj = getPubChemInfo(smiles=product_smiles)
    product = Product()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    product.reaction_id = reaction_obj
    product.smiles = product_smiles
    # if pubcheminfoobj:
    #     product.pubcheminfo_id = pubcheminfoobj
    product_svg_string = createSVGString(product_smiles)
    product_svg_fn = default_storage.save(
        "productimages/.svg", ContentFile(product_svg_string)
    )
    product.image = product_svg_fn
    product.save()


def createReactantModel(
    reaction_id: int, reactant_smiles: str, previous_reaction_product: bool
) -> int:
    """Creates a Django reactant object - the reactant in a reaction

    Parameters
    ----------
    reaction_id: int
        The PubChem DB compound id
    reactant_smiles: str
        The SMILES of the reactant
    previous_reaction_product: bool
        If the reactant is the product from a previous reaction in the method. Used
        for later determining if reactant needs to be purchased vs is made in the
        synthesis

    Returns
    -------
    reactant_id: int
        The id of the reactant model object created
    """
    # pubcheminfoobj = getPubChemInfo(smiles=reactant_smiles)
    reactant = Reactant()
    reaction_obj = Reaction.objects.get(id=reaction_id)
    reactant.reaction_id = reaction_obj
    reactant.smiles = reactant_smiles
    # if pubcheminfoobj:
    #     reactant.pubcheminfo_id = pubcheminfoobj
    reactant.previousreactionproduct = previous_reaction_product
    reactant.save()
    return reactant.id


def createCatalogEntryModel(
    catalog_entry: dict = None,
    target_id: int = None,
    reactant_id: int = None,
    previous_reaction_product: bool = False,
    lab_inventory: bool = False,
):
    """Creates a Django catalogentry object - the catalog details
       for a reactant or target

    Parameters
    ----------
    catalog_entry: dict
        The optional Manifold catalog entry information
    target_id: int
        The id of the Django target model object
    reactant_id: int
        The optional id of the Django reactant model object
    previous_reaction_product: bool
        If the reactant is a previous product from a reaction, therfor does not need to be purchased
    lab_inventory: bool
        Set to True if the compound is in the lab inventory
    """
    catalogentry = CatalogEntry()
    if target_id:
        target_obj = Target.objects.get(id=target_id)
        catalogentry.target_id = target_obj
    if reactant_id:
        reactant_obj = Reactant.objects.get(id=reactant_id)
        catalogentry.reactant_id = reactant_obj

    if previous_reaction_product or lab_inventory:
        catalogentry.vendor = "reaction product"
        catalogentry.catalogid = "NA"
        catalogentry.priceinfo = "< $100 / g"
        catalogentry.upperprice = 0
        catalogentry.leadtime = 0

    if not previous_reaction_product and not lab_inventory:
        catalogentry.vendor = catalog_entry["catalogName"]
        catalogentry.catalogid = catalog_entry["catalogId"]

        if catalog_entry["catalogName"] == "generic":
            catalogentry.upperprice = None
            catalogentry.leadtime = None

        if catalog_entry["purchaseInfo"]["isScreening"]:
            if catalog_entry["purchaseInfo"]["scrLeadTimeWeeks"] != "unknown":
                catalogentry.leadtime = catalog_entry["purchaseInfo"][
                    "scrLeadTimeWeeks"
                ]
            else:
                catalogentry.leadtime = None
            if catalog_entry["purchaseInfo"]["scrPriceRange"] != "unknown":
                priceinfo = catalog_entry["purchaseInfo"]["scrPriceRange"]
                catalogentry.priceinfo = priceinfo
                priceinfo = priceinfo.replace(" ", "")
                if priceinfo[0] == "<" or priceinfo[0] == ">":
                    if "k" in priceinfo:
                        upperprice = int("".join(filter(str.isdigit, priceinfo))) * 1000
                    else:
                        upperprice = int("".join(filter(str.isdigit, priceinfo)))
                if priceinfo[0] == "$":
                    if "k" in priceinfo:
                        upperprice = (
                            int("".join(filter(str.isdigit, priceinfo.split("-")[1])))
                            * 1000
                        )
                    else:
                        upperprice = int(
                            "".join(filter(str.isdigit, priceinfo.split("-")[1]))
                        )
                catalogentry.upperprice = upperprice
            else:
                catalogentry.upperprice = None

        if not catalog_entry["purchaseInfo"]["isScreening"]:
            if catalog_entry["purchaseInfo"]["bbLeadTimeWeeks"] != "unknown":
                catalogentry.leadtime = catalog_entry["purchaseInfo"]["bbLeadTimeWeeks"]
            else:
                catalogentry.leadtime = None
            if catalog_entry["purchaseInfo"]["bbPriceRange"] != "unknown":
                priceinfo = catalog_entry["purchaseInfo"]["bbPriceRange"]
                catalogentry.priceinfo = priceinfo
                priceinfo = priceinfo.replace(" ", "")
                if priceinfo[0] == "<" or priceinfo[0] == ">":
                    if "k" in priceinfo:
                        upperprice = int("".join(filter(str.isdigit, priceinfo))) * 1000
                    else:
                        upperprice = int("".join(filter(str.isdigit, priceinfo)))
                if priceinfo[0] == "$":
                    if "k" in priceinfo:
                        upperprice = (
                            int("".join(filter(str.isdigit, priceinfo.split("-")[1])))
                            * 1000
                        )
                    else:
                        upperprice = int(
                            "".join(filter(str.isdigit, priceinfo.split("-")[1]))
                        )
                catalogentry.upperprice = upperprice
            else:
                catalogentry.upperprice = None
    catalogentry.save()


class CreateEncodedActionModels(object):
    """
    Creates a createEncodedActionModels object for creating action models
    for a reaction
    """

    def __init__(
        self,
        intramolecular: bool,
        actionsessions: dict,
        target_id: int,
        reaction_id: int,
        reactant_pair_smiles: list,
        reaction_name: str,
    ):
        """
        CreateEncodedActionModels

        Parameters
        ----------
        intramolecular: bool
            Set to true if the reaction is intramolecular
        actionsessions: dict
            The action (add, stir, analyse) actionsessions  for the execution of a synthesis
        project_id: int
            The Django project model object id
        reaction_id: int
            The Django reaction model object id
        target_id: int
            The Django target model object id that the reaction is related
            to
        reactant_pair_smiles: list
            The reactant smiles used in the reatcion
        reaction_name: str
            The name of the reaction type
        """
        self.intramolecular = intramolecular
        self.mculeapi = MCuleAPI()
        self.actionsessions = actionsessions
        self.reaction_id = reaction_id
        self.reaction_obj = Reaction.objects.get(id=reaction_id)
        self.reactant_pair_smiles = reactant_pair_smiles
        self.reaction_name = reaction_name
        self.product_obj = getProduct(reaction_id=reaction_id)
        self.target_obj = Target.objects.get(id=target_id)
        self.productmols = self.getProductMols()
        self.productmass = calculateMassFromMols(
            mols=self.productmols, SMILES=self.product_obj.smiles
        )
        self.productsmiles = getProductSmiles(reaction_ids=[reaction_id])[0]
        self.mculeidlist = []
        self.amountslist = []

        for actionsession in actionsessions:
            actionsessiontype = actionsession["type"]
            if "continuation" in actionsession:
                continuation = actionsession["continuation"]
            else:
                continuation = False
            actionsession_obj = self.createActionSessionModel(
                actionsessiontype=actionsessiontype,
                driver=actionsession["driver"],
                sessionnumber=actionsession["sessionnumber"],
                continuation=continuation,
            )
            if actionsessiontype == "reaction" and self.intramolecular:
                reaction_actions = actionsession["intramolecular"]["actions"]
                [
                    self.createEncodedActionModel(
                        actionsession_obj=actionsession_obj, action=action
                    )
                    for action in reaction_actions
                ]
            if actionsessiontype == "reaction" and not self.intramolecular:
                reaction_actions = actionsession["intermolecular"]["actions"]
                [
                    self.createEncodedActionModel(
                        actionsession_obj=actionsession_obj, action=action
                    )
                    for action in reaction_actions
                ]
            if actionsessiontype != "reaction":
                actions = actionsession["actions"]
                [
                    self.createEncodedActionModel(
                        actionsession_obj=actionsession_obj, action=action
                    )
                    for action in actions
                ]

    def getProductMols(self):
        """Calculates mols fo product required based on reaction
        success and amount final target required

        Returns
        -------
        product_mols: float
            The product mols required
        """
        proceedingreactionqueryset = checkProceedingReactions(
            reaction_id=self.reaction_obj.id
        )

        if proceedingreactionqueryset:
            reactionclasslist = [
                reactionobj.reactionclass for reactionobj in proceedingreactionqueryset
            ] + [self.reaction_name]
            recipelist = [
                reactionobj.recipe for reactionobj in proceedingreactionqueryset
            ] + [self.reaction_obj.recipe]
            reactionyields = getReactionYields(
                reactionclasslist=reactionclasslist, recipelist=recipelist
            )
            yieldcorrection = math.prod(reactionyields)

        if not proceedingreactionqueryset:
            reactionyields = getReactionYields(
                reactionclasslist=[self.reaction_name],
                recipelist=[self.reaction_obj.recipe],
            )
            yieldcorrection = math.prod(reactionyields)

        product_mols = self.target_obj.mols / yieldcorrection
        return product_mols

    def createEncodedActionModel(self, actionsession_obj: ActionSession, action: dict):
        """Calls the functions to create the appropriate Django action
           model object for the action

        Parameters
        ----------
        actionsession_obj: ActionSession
            The type of action session being excuted
        action: dict
            The action that needs to be executed
        """
        actionMethods = {
            "add": self.createAddActionModel,
            "extract": self.createExtractActionModel,
            "mix": self.createMixActionModel,
            "stir": self.createStirActionModel,
        }

        action_type = action["type"]

        if action_type in actionMethods:
            actionMethods[action_type](
                actionsession_obj=actionsession_obj, action=action
            )
            return True
        else:
            logger.info(action_type)

    def calculateMass(
        self,
        calcunit: str,
        calcvalue: float,
        reactant_MW: float,
    ) -> float:
        """Calculates the reactant mass (mg) required for an add action step

        Parameters
        ----------
        calcunit: str
            The unit used for the calculation eg. mass equivalents
        calcvalue: int
            The the quivalents to use in the calculation
        reactant_MW: float
            The molecular weight (g/mol) of the reactant

        Returns
        -------
        mass_material: float
            The mass (mg) of the material required for the add action step
        """
        try:
            if calcunit == "moleq":
                mol_material = float(calcvalue) * self.productmols
                mass_material = mol_material * reactant_MW * 1000
                return mass_material
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))

    def calculateVolume(
        self,
        calcunit: str,
        calcvalue: float,
        conc_reagents: float = None,
        reactant_density: float = None,
        reactant_MW: float = None,
    ) -> float:
        """Calculates the reactant volume (ul) required for an add action step

        Parameters
        ----------
        calcunit: str
            The unit used for the calculation eg. mass equivalents
        calcvalue: int
            The the quivalents to use in the calculation
        conc_reagents: float
            The optional concentration of the reactant
        reactant_density: float
            The optional density (g/ml) of the reactant
        reactant_MW: float
            The optional molecular weight (g/mol) of the reactant

        Returns
        -------
        vol_material: float
            The volume (ul) of the material required for the add action step
        """
        try:
            if calcunit == "masseq":
                vol_material = float(calcvalue) * self.productmass
                return vol_material
            if calcunit == "moleq":
                mol_material = float(calcvalue) * self.productmols
                if reactant_density:
                    vol_material = (
                        (mol_material * reactant_MW) / reactant_density
                    ) * 1e3
                if not reactant_density:
                    vol_material = (mol_material / conc_reagents) * 1e6  # in uL
                return vol_material
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))

    def createActionSessionModel(
        self,
        actionsessiontype: str,
        driver: str,
        sessionnumber: int,
        continuation: bool,
    ):
        """Creates a Django action session object - a session are colelctions of
           actions that can be collectively exceuted. Eg. a reaction will
           inlcude a series of add actions

        Parameters
        ----------
        type: str
            The session type eg. reaction, stir etc
        driver: str
            The main driver operating the session. If any operations can be
            automated then the driver is a robot else human
        sessionnumber: int
            The session sequence number eg. reaction session one happens first
        continuation: bool
            Does the session continue from a previous session eg. reaction sessions broken
            up by stir action session where another reagent/reactant is added after stirring
        """
        try:
            actionsession = ActionSession()
            actionsession.reaction_id = self.reaction_obj
            actionsession.type = actionsessiontype
            actionsession.driver = driver
            actionsession.sessionnumber = sessionnumber
            actionsession.continuation = continuation
            actionsession.save()
            return actionsession
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))

    def createAddActionModel(self, actionsession_obj: ActionSession, action: dict):
        """Creates a Django add action object - an add action

        Parameters
        ----------
        actionsession_obj: ActionSession
            The type of action session being excuted
        action: dict
            The analyse action that needs to be executed
        """
        try:
            actionnumber = action["actionnumber"]
            fromplatetype = action["content"]["plates"]["fromplatetype"]
            toplatetype = action["content"]["plates"]["toplatetype"]
            materialinfo = action["content"]["material"]
            if materialinfo["SMARTS"]:
                # MUst fix this to match SMARTS with molecule vs just taking first reactant SMILES!!!!!!
                # FIX!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # Do not use reactant pair smiles here
                smiles = self.reactant_pair_smiles[0]
                del self.reactant_pair_smiles[0]
            if materialinfo["SMILES"]:
                smiles = materialinfo["SMILES"]
            if not materialinfo["SMILES"] and not materialinfo["SMARTS"]:
                smiles = self.productsmiles
            calcvalue = materialinfo["quantity"]["value"]
            calcunit = materialinfo["quantity"]["unit"]
            concentration = materialinfo["concentration"]
            solvent = materialinfo["solvent"]
            mol = Chem.MolFromSmiles(smiles)
            molecular_weight = Descriptors.MolWt(mol)
            add = AddAction()
            add.reaction_id = self.reaction_obj
            add.actionsession_id = actionsession_obj
            add.number = actionnumber
            add.fromplatetype = fromplatetype
            add.toplatetype = toplatetype
            add.smiles = smiles
            add.molecularweight = molecular_weight
            if calcunit == "ul":
                add.volume = calcvalue
                add.solvent = solvent
            if calcunit == "masseq":
                add.volume = self.calculateVolume(
                    calcunit=calcunit,
                    calcvalue=calcvalue,
                )
                add.solvent = solvent
            if calcunit == "moleq":
                if not solvent and not concentration:

                    add.mass = self.calculateMass(
                        calcunit=calcunit,
                        calcvalue=calcvalue,
                        reactant_MW=molecular_weight,
                    )
                if "density" in materialinfo:
                    reactant_density = materialinfo["density"]
                    add.volume = self.calculateVolume(
                        calcunit=calcunit,
                        calcvalue=calcvalue,
                        reactant_density=reactant_density,
                        reactant_MW=molecular_weight,
                    )
                if solvent:
                    add.volume = self.calculateVolume(
                        calcunit=calcunit,
                        calcvalue=calcvalue,
                        conc_reagents=concentration,
                    )
                    add.solvent = solvent
            if concentration:
                add.concentration = concentration
            add.save()

        except Exception as e:
            print(smiles)
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))

    def createExtractActionModel(self, actionsession_obj: ActionSession, action: dict):
        """Creates a Django extract action object - an add action

        Parameters
        ----------
        actionsession_obj: ActionSession
            The type of action session being excuted
        action: dict
            The analyse action that needs to be executed
        """
        try:
            actionnumber = action["actionnumber"]
            fromplatetype = action["content"]["plates"]["fromplatetype"]
            toplatetype = action["content"]["plates"]["toplatetype"]
            materialinfo = action["content"]["material"]
            calcvalue = materialinfo["quantity"]["value"]
            calcunit = materialinfo["quantity"]["unit"]
            concentration = materialinfo["concentration"]
            if not concentration:
                concentration = 0
            solvent = materialinfo["solvent"]
            smiles = self.productsmiles
            mol = Chem.MolFromSmiles(smiles)
            molecular_weight = Descriptors.MolWt(mol)
            extract = ExtractAction()
            extract.reaction_id = self.reaction_obj
            extract.actionsession_id = actionsession_obj
            extract.number = actionnumber
            extract.fromplatetype = fromplatetype
            extract.toplatetype = toplatetype
            extract.smiles = smiles
            extract.molecularweight = molecular_weight
            if calcunit == "ul":
                extract.volume = calcvalue
                extract.solvent = solvent
            if calcunit == "masseq":
                extract.volume = self.calculateVolume(
                    calcunit=calcunit,
                    calcvalue=calcvalue,
                )
                extract.solvent = solvent
            if calcunit == "moleq":
                if not solvent:
                    reactant_density = materialinfo["density"]
                    extract.volume = self.calculateVolume(
                        calcunit=calcunit,
                        calcvalue=calcvalue,
                        reactant_density=reactant_density,
                        reactant_MW=molecular_weight,
                    )
                if solvent:
                    extract.volume = self.calculateVolume(
                        calcunit=calcunit,
                        calcvalue=calcvalue,
                        conc_reagents=concentration,
                    )
                    extract.solvent = solvent
            if "bottomlayerquantity" in materialinfo:
                bottomlayercalcvalue = materialinfo["bottomlayerquantity"]["value"]
                bottomlayercalcunit = materialinfo["bottomlayerquantity"]["unit"]
                extract.bottomlayervolume = self.calculateVolume(
                    calcunit=bottomlayercalcunit,
                    calcvalue=bottomlayercalcvalue,
                )
            extract.concentration = concentration
            extract.save()

        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))

    def createMixActionModel(self, actionsession_obj: ActionSession, action: dict):
        """Creates a Django mix action object - a mix action

        Parameters
        ----------
        actionsession_obj: ActionSession
            The type of action session being excuted
        action: dict
            The analyse action that needs to be executed
        """
        try:
            actionnumber = action["actionnumber"]
            platetype = action["content"]["platetype"]
            repetitions = action["content"]["repetitions"]["value"]

            mix = MixAction()
            mix.reaction_id = self.reaction_obj
            mix.actionsession_id = actionsession_obj
            mix.number = actionnumber
            mix.platetype = platetype
            mix.repetitions = repetitions
            mix.save()

        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))

    def createStirActionModel(self, actionsession_obj: ActionSession, action: dict):
        """Creates a Django stir action object - a stir action

        Parameters
        ----------
        actionsession_obj: ActionSession
            The type of action session being excuted
        action: dict
            The analyse action that needs to be executed
        """
        try:
            actionnumber = action["actionnumber"]
            platetype = action["content"]["platetype"]
            duration = action["content"]["duration"]["value"]
            durationunit = action["content"]["duration"]["unit"]
            temperature = action["content"]["temperature"]["value"]
            temperatureunit = action["content"]["temperature"]["unit"]

            stir = StirAction()
            stir.reaction_id = self.reaction_obj
            stir.actionsession_id = actionsession_obj
            stir.number = actionnumber
            stir.platetype = platetype
            stir.duration = duration
            stir.durationunit = durationunit
            stir.temperature = temperature
            stir.temperatureunit = temperatureunit
            stir.save()

        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
