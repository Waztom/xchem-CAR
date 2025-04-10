from django.db.models import QuerySet
from django.db.models import Q, Max
from rdkit.Chem import Descriptors
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

import pubchempy as pcp
import itertools
import re
import inspect
import logging
import pandas as pd
import csv
import datetime

from .models import (
    ActionSession,
    Batch,
    Method,
    AddAction,
    ExtractAction,
    OTBatchProtocol,
    Product,
    Reaction,
    Target,
    Plate,
)

from .recipebuilder.encodedrecipes import encoded_recipes

logger = logging.getLogger(__name__)


def getAddActionQuerySet(
    reaction_ids: list,
    actionsession_ids: list = None,
    actionsessiontype: str = None,
) -> QuerySet[AddAction]:
    """Get add actions queryset for reaction_id

    Parameters
    ----------
    reaction_ids: list
        The reactions to search for related add actions
    actionsession_ids: list
        Optional action session ids to match add actions with
    actionsessiontype: str
        The optional action session type to look for add actions for

    Returns
    -------
    addactionqueryset: QuerySet[AddAction]
        The add actions related to the reaction
    """

    if actionsession_ids:
        criterion1 = Q(reaction_id__in=reaction_ids)
        criterion2 = Q(actionsession_id__in=actionsession_ids)
        addactionqueryset = AddAction.objects.filter(criterion1 & criterion2).order_by(
            "id"
        )
        return addactionqueryset
    if actionsessiontype:
        criterion1 = Q(reaction_id__in=reaction_ids)
        criterion2 = Q(actionsession_id__type=actionsessiontype)
        addactionqueryset = AddAction.objects.filter(criterion1 & criterion2).order_by(
            "id"
        )
        return addactionqueryset


def getExtractActionQuerySet(
    reaction_ids: list,
    actionsession_ids: list = None,
    actionsessiontype: str = None,
) -> QuerySet[ExtractAction]:
    """Get extract actions queryset for reaction_id

    Parameters
    ----------
    reaction_ids: list
        The reactions to search for related add actions
    actionsession_ids: list
        Optional action session ids to match add actions with
    actionsessiontype: str
        The optional action session type to look for add actions for

    Returns
    -------
    extractactionqueryset: QuerySet[ExtractAction]
        The extract actions related to the reaction
    """

    if actionsession_ids:
        criterion1 = Q(reaction_id__in=reaction_ids)
        criterion2 = Q(actionsession_id__in=actionsession_ids)
        extractactionqueryset = ExtractAction.objects.filter(
            criterion1 & criterion2
        ).order_by("id")
        return extractactionqueryset
    if actionsessiontype:
        criterion1 = Q(reaction_id__in=reaction_ids)
        extractactionqueryset = ExtractAction.objects.filter(
            criterion1, actionsession_id__type=actionsessiontype
        ).order_by("id")
        return extractactionqueryset


def getOTBatchProtocolQuerySet(batch_id: int) -> QuerySet[OTBatchProtocol]:
    """Gets the OT batch protocol queryset

    Parameters
    ----------
    batch_id: int
        The batch id to saerch for an associated OT protocol

    Returns
    -------
        The OT Batch protocol queryset
    """
    otbatchprotocolqueryset = OTBatchProtocol.objects.filter(batch_id=batch_id)
    return otbatchprotocolqueryset


def getActionSessionSequenceNumbers(
    actionsessionqueryset: QuerySet[ActionSession],
) -> list:
    """Set of action session sequence numbers

    Returns
    ------
    sessionnumbers: list
        The set of session numbers in an action session
        queryset eg. [1,2,3,4....n]
    """
    maxsessionnumber = actionsessionqueryset.aggregate(Max("sessionnumber"))[
        "sessionnumber__max"
    ]
    sessionnumbers = list(range(1, maxsessionnumber + 1))
    return sessionnumbers


def getActionSessionTypes(actionsessionqueryset: QuerySet[ActionSession]) -> QuerySet:
    """Set of action session types

    Returns
    ------
    actionsessiontypes: QuerySet
        The set of action session types in a queryset
        eg. ["reaction", "workup", "stir"]
    """
    actionsessiontypes = set(list(actionsessionqueryset.values_list("type", flat=True)))
    return actionsessiontypes


def getGroupedActionSessionSequences(
    sessionnumbers: list, actionsessionqueryset: QuerySet[ActionSession]
) -> list:
    """Group action sessions by sequence number

    Parameters
    ----------
    sessionnumbers: list
        The list of action session sequence numbers
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset to group by sequence number

    Returns
    -------
    groupedactionsessionsequences: list
        List of sub-lists of action sessions grouped by sequence number
    """
    groupedactionsessionsequences = []
    for sessionnumber in sessionnumbers:
        actionsessiongroup = actionsessionqueryset.filter(
            sessionnumber=sessionnumber
        ).order_by("-pk")
        groupedactionsessionsequences.append(actionsessiongroup)
    return groupedactionsessionsequences


def getGroupedActionSessionTypes(
    actionsessiontypes: QuerySet, actionsessionqueryset: QuerySet[ActionSession]
) -> list:
    """Group action sessions by type

    Parameters
    ----------
    actionsessiontypes: QuerySet
        The list of action session sequence numbers
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset to group by sequence number

    Returns
    -------
    groupedactionsessionquerysettypes: list
        List of sub-lists of action sessions grouped by types
    """
    groupedactionsessiontypes = []
    for actionsessiontype in actionsessiontypes:
        actionsessiongrouptype = actionsessionqueryset.filter(
            type=actionsessiontype
        ).order_by("-pk")
        if actionsessiongrouptype:
            groupedactionsessiontypes.append(actionsessiongrouptype)
    return groupedactionsessiontypes


def getActionSessionQuerySet(
    reaction_ids: QuerySet[Reaction],
    driver: str = None,
) -> QuerySet[ActionSession]:
    """Returns the action session queryset for a type of driver
       (human or robot)

    Parameters
    ----------
    reactions_ids: QuerySet[Reaction]
        The reactions that the action session will execute
    driver: str
        The optional main driver of the action session

    Returns
    -------
    actionsessionqueryset: QuerySet[ActionSession]
        The action session queryset for a given driver
    """
    if driver:
        criterion1 = Q(reaction_id__in=reaction_ids)
        criterion2 = Q(driver=driver)
        actionsessionqueryset = ActionSession.objects.filter(criterion1 & criterion2)
        if actionsessionqueryset:
            return actionsessionqueryset
    if not driver:
        criterion1 = Q(reaction_id__in=reaction_ids)
        actionsessionqueryset = ActionSession.objects.filter(criterion1)
        if actionsessionqueryset:
            return actionsessionqueryset


def getPreviousObjEntries(queryset: list, obj: object) -> QuerySet:
    """Finds all previous objects relative to obj of queryset"""
    previousqueryset = queryset.filter(pk__lt=obj.pk).order_by("-pk")
    return previousqueryset


def checkPreviousReactionFailures(reactionobj: Reaction) -> bool:
    """Check if any previous reaction failures for a method"""
    reactionqueryset = getReactions(method_ids=[reactionobj.method_id.id])
    previousreactionqueryset = getPreviousObjEntries(
        queryset=reactionqueryset, obj=reactionobj
    )
    failedreactions = previousreactionqueryset.filter(success=False)
    if failedreactions.exists():
        return True
    else:
        return False


def checkNoMethodSteps(reactionobj: Reaction) -> bool:
    """Check no reaction steps in method is > 1"""
    methodobj = reactionobj.method_id
    noreactionsteps = methodobj.nosteps
    if noreactionsteps > 1:
        return True
    else:
        return False


def getReactionsToDo(groupreactionqueryset: QuerySet[Reaction]) -> QuerySet[Reaction]:
    """Get reactions that need to be done. Exclude those in methods that had
    failed previous reaction step

    Parameters
    ---------
    groupreactionqueryset: QuerySet[Reaction]
        The group of reactions to find to do based on if the previous reaction was successful

    Returns
    -------
    groupreactiontodoqueryset: QuerySet[Reaction]
        The reactions that need to be done
    """
    reactionstodo = []
    for reactionobj in groupreactionqueryset:
        if checkNoMethodSteps(reactionobj=reactionobj):
            if not checkPreviousReactionFailures(reactionobj=reactionobj):
                reactionstodo.append(reactionobj.id)
    groupreactiontodoqueryset = groupreactionqueryset.filter(id__in=reactionstodo)
    return groupreactiontodoqueryset


def getTargets(batch_ids: QuerySet[Batch]) -> QuerySet[Target]:
    targetqueryset = Target.objects.filter(batch_id__in=batch_ids).order_by("id")
    return targetqueryset


def getMethods(target_ids: QuerySet[Target]) -> QuerySet[Method]:
    methodqueryset = (
        Method.objects.filter(target_id__in=target_ids)
        .filter(otchem=True)
        .order_by("id")
    )
    return methodqueryset


def getReactions(method_ids: QuerySet[Method]) -> QuerySet[Reaction]:
    reactionqueryset = Reaction.objects.filter(method_id__in=method_ids).order_by("id")
    return reactionqueryset


def getBatchTag(batchid):
    batch_obj = Batch.objects.get(id=batchid)
    batchtag = batch_obj.batchtag
    return batchtag


def getBatchReactions(batchid: int) -> QuerySet[Reaction]:
    targetqueryset = getTargets(batch_ids=[batchid])
    if targetqueryset:
        methodqueryset = getMethods(target_ids=targetqueryset)
        if methodqueryset:
            reactionqueryset = getReactions(method_ids=methodqueryset)
            if reactionqueryset:
                return reactionqueryset


def getMaxReactionNumber(reactionqueryset: QuerySet[Reaction]) -> int:
    """Get the maximum number of reaction steps in a reaction queryset

    Parameters
    ----------
    reactionqueryset: QuerySet[Reaction]
        The reaction queryset to get the max number of reaction steps for

    Returns
    -------
    maxreactionnumber: int
        The maximum reaction number in a set of reactions

    """

    maxreactionnumber = reactionqueryset.aggregate(Max("number"))["number__max"]
    return maxreactionnumber


def groupReactions(reactionqueryset: QuerySet[Reaction], maxreactionnumber: int):
    """
    Groups reactionqueries into first reactions, second reactions and so on
    """

    groupedreactionquerysets = []
    for i in range(1, maxreactionnumber + 1):
        reactionnumberqueryset = (
            reactionqueryset.filter(number=i).distinct().order_by("id")
        )
        if reactionnumberqueryset:
            groupedreactionquerysets.append(reactionnumberqueryset)
    return groupedreactionquerysets


def getReactantsToBuy(batch_ids: list[int]) -> list:
    """Finds the reactnats that need to be bought to execute a batch/batches
    synthesis. Finds recatants that are not made in previous method's reactions

    Parameters
    ----------
    batch_ids: list[int]
        The batch ids to search for reactants to buy to complete the synthesis

    Returns
    -------
    reactants_to_buy: list
        The SMILES of the reactants that need to be bought. Excludes reactants
        made in previous reaction steps
    """
    reactants_to_buy = []
    for batch_id in batch_ids:
        batchobj = Batch.objects.get(id=batch_id)
        targetqueryset = batchobj.targets.all()
        methods = [target.methods.all() for target in targetqueryset]
        method_sublist = [item for sublist in methods for item in sublist]
        reactions = [method.reactions.all() for method in method_sublist]
        reaction_sublist = [item for sublist in reactions for item in sublist]
        reactants = [reaction.reactants.all() for reaction in reaction_sublist]
        reactants_sublist = [item for sublist in reactants for item in sublist]
        reactants_batch_to_buy = list(
            set(
                [
                    reactant.smiles
                    for reactant in reactants_sublist
                    if reactant.previousreactionproduct == False
                ]
            )
        )
        reactants_to_buy = reactants_to_buy + reactants_batch_to_buy
    return list(set(reactants_to_buy))


def wellIndexToWellName(wellindex: int, platesize: int) -> str:
    """Converts a well index to a human readable well name

    Parameters
    ----------
    wellindex: int
        The well index to convert to a well name
    platesize: int
        The plate size to convert the well index to a well name

    Returns
    -------
    wellname: str
        The well name eg. A01
    """
    platewellindexconversions = {
        4: "A" + "%02d" % ((wellindex) + 1,),
        24: "ABCD"[(wellindex) % 4] + "%02d" % ((wellindex) // 4 + 1,),
        96: "ABCDEFGH"[(wellindex) % 8] + "%02d" % ((wellindex) // 8 + 1,),
        384: "ABCDEFGHIJKLMNOP"[(wellindex) % 16] + "%02d" % ((wellindex) // 16 + 1,),
    }
    wellname = platewellindexconversions[platesize]
    return wellname


def getBatchTargetMWs(batch_id: int) -> list[float]:
    """Gets the molecular weights of the final target compounds
    for batches

    Parameters
    ----------
    batch_id: int
        The batch id to get the target molecular weights for

    Returns
    -------
    target_MWs: list
        The molecular weights of the targets for a batch
    """

    batchobj = Batch.objects.get(id=batch_id)
    targetqueryset = batchobj.targets.all()
    smiles = [targetobj.smiles for targetobj in targetqueryset]
    target_MWs = getMWs(smiles=smiles)
    return target_MWs


def getBatchTargetSmiles(batch_id: int) -> list[float]:
    """Gets the SMILES of the final target compounds
    for a batch

    Parameters
    ----------
    batch_id: int
        The batch id to get the target molecular weights for

    Returns
    -------
    target_SMILES: list
        The target SMILES for a batch
    """

    batchobj = Batch.objects.get(id=batch_id)
    targetqueryset = batchobj.targets.all()
    target_SMILES = [targetobj.smiles for targetobj in targetqueryset]
    return target_SMILES


def getBatchReactionIDs(batch_id: int, reaction_number: int) -> list[float]:
    """Gets the reaction ids for a reaction number in
       a batch

    Parameters
    ----------
    batch_id: int
        The batch id to get the target molecular weights for
    reaction_number: int
        The reactions to find product SMILES for

    Returns
    -------
    reaction_IDs: list
        The reaction IDs for a reaction step in a batch
    """
    reaction_IDs = []
    batchobj = Batch.objects.get(id=batch_id)
    targetqueryset = batchobj.targets.all().order_by("id")
    for targetobj in targetqueryset:
        methodqueryset = targetobj.methods.all().order_by("id")
        for methodobj in methodqueryset:
            reactionqueryset = (
                methodobj.reactions.all().filter(number=reaction_number).order_by("id")
            )
            for reactionobj in reactionqueryset:
                reaction_IDs.append(reactionobj.id)
    return reaction_IDs


def updateTargetMols(batch_id: int, concentration: float, volume: float) -> int:
    """Updates the Target mols in a Batch - using concentartion and volume

    Parameters
    ----------
    batch_id: int
        The batch model object id that the targets are linked to
    concentration: float
        The concentration (mM) to be updated for the targets
    volume: float
        The volume (uL) to be updated for the targets

    Returns
    -------
    """
    target_qs = Batch.objects.get(id=batch_id).targets.all()
    for target_obj in target_qs:
        smiles = target_obj.smiles
        mols = calculateMolsFromConc(
            target_concentration=concentration, target_volume=volume
        )
        mass = calculateMassFromMols(mols=mols, SMILES=smiles)
        target_obj.mols = mols
        target_obj.concentration = concentration
        target_obj.volume = volume
        target_obj.mass = mass
        target_obj.save()


def updateReactionSuccessToFail(reaction_ids: list[int]):
    """Updates reactions to be failures

    Parameters
    ----------
    reactions_ids: list[int]
        The reactions to update as synthetic failures
    """
    if Reaction.objects.filter(id__in=reaction_ids).exists():
        Reaction.objects.filter(id__in=reaction_ids).update(success=False)


def updateBatchMethodOTFriendly(batch_id: int):
    """Updates a batch of methods to all be OT friendly

    Parameters
    ----------
    batch_id: int
        The batch id to get the target molecular weights for
    """
    batchobj = Batch.objects.get(id=batch_id)
    targetqueryset = batchobj.targets.all().order_by("id")
    for targetobj in targetqueryset:
        methodqueryset = targetobj.methods.all().order_by("id")
        for methodobj in methodqueryset:
            methodobj.otchem = True
            methodobj.save()


def deleteBatchActionSessions(batch_id: int):
    """Deletes all actions for a batch

    Parameters
    ----------
    batch_id: int
        The batch id to delete all actions for
    """
    reaction_qs = getBatchReactions(batchid=batch_id)
    getActionSessionQuerySet(reaction_ids=reaction_qs).delete()


def updateRecipeType(
    batch_id: int, reaction_class: str, current_recipe: str, recipe_to_use: str
):
    """Updates the recipe type for a Reactions in a Batch

    Parameters
    ----------
    batch_id: int
        The batch id to update the recipe type for
    reaction_class: str
        The reaction class to update the recipe for
    current_recipe: str
        The current recipe to update
    recipe_to_use: str
        The recipe to update to
    """
    reaction_qs = getBatchReactions(batchid=batch_id).filter(
        reactionclass=reaction_class, recipe=current_recipe
    )
    for reaction_obj in reaction_qs:
        reaction_obj.recipe = recipe_to_use
        reaction_obj.save()


# def deleteBatchActionSessions(batch_id: int):
#     """Delete the Action Sessions related to a Batch

#     Parameters
#     ----------
#     batch_id: int
#         The batch id to delete all actions for
#     """
#     reaction_qs = getBatchReactions(batchid=batch_id)
#     getActionSessionQuerySet(reaction_ids=reaction_qs).delete()


def getBatchReactionProductSmiles(batch_id: int, reaction_number: int) -> list[float]:
    """Gets the MWs of the products for a reaction in
       a batch

    Parameters
    ----------
    batch_id: int
        The batch id to get the target molecular weights for
    reaction_number: int
        The reactions to find product SMILES for

    Returns
    -------
    product_MWs: list
        The product MWs for a reaction step in a batch
    """
    product_SMILES = []
    batchobj = Batch.objects.get(id=batch_id)
    targetqueryset = batchobj.targets.all().order_by("id")
    for targetobj in targetqueryset:
        methodqueryset = targetobj.methods.all().order_by("id")
        for methodobj in methodqueryset:
            reactionqueryset = (
                methodobj.reactions.all().filter(number=reaction_number).order_by("id")
            )
            for reactionobj in reactionqueryset:
                product_SMILES = product_SMILES + list(
                    reactionobj.products.all().values_list("smiles", flat=True)
                )
    return product_SMILES


def getPlateMap(plate_ids: list, out_dir: str):
    """Generates a Plate Map for a list of plate ids

    Parameters
    ----------
    plate_ids: list
        The plate ids to generate the platemap info for
    out_dir: str
        The directory to write the csv to

    Returns
    -------
    platemap_csv: File
        The csv files in tmp-files
    """

    plate_info = {
        "plate_id": [],
        "well_index": [],
        "target_ids": [],
        "target_names": [],
        "target_smiles": [],
        "target_MWs": [],
        "reactant_1_smiles": [],
        "reactant_2_smiles": [],
        "reactant_1_MWs": [],
        "reactant_2_MWs": [],
    }

    try:
        for plate_id in plate_ids:
            plate = Plate.objects.get(id=plate_id)
            wells = plate.well_set.all().order_by("id")

            for well in wells:
                well_index = well.index
                target_id = well.method_id.target_id.id
                target_name = well.method_id.target_id.name
                target_smi = well.smiles
                target_mw = getMWs(smiles=[target_smi])[0]
                reactant_smiles = well.reaction_id.reactants.values_list(
                    "smiles", flat=True
                )
                reactant_mws = getMWs(smiles=reactant_smiles)
                if len(reactant_smiles) == 1:
                    reactant_1_smi = reactant_smiles[0]
                    reactant_1_mw = reactant_mws[0]
                    reactant_2_smi = ""
                    reactant_2_mw = ""
                if len(reactant_smiles) == 2:
                    reactant_1_smi = reactant_smiles[0]
                    reactant_1_mw = reactant_mws[0]
                    reactant_2_smi = reactant_smiles[1]
                    reactant_2_mw = reactant_mws[1]
                plate_info["plate_id"].append(plate_id)
                plate_info["well_index"].append(well_index)
                plate_info["target_ids"].append(target_id)
                plate_info["target_names"].append(target_name)
                plate_info["target_smiles"].append(target_smi)
                plate_info["target_MWs"].append(target_mw)
                plate_info["reactant_1_smiles"].append(reactant_1_smi)
                plate_info["reactant_2_smiles"].append(reactant_2_smi)
                plate_info["reactant_1_MWs"].append(reactant_1_mw)
                plate_info["reactant_2_MWs"].append(reactant_2_mw)

        filename = "platemapids-" + "-".join(map(str, plate_ids))

        with open("{}{}.csv".format(out_dir, filename), "w") as f:
            writer = csv.writer(f)
            limit = len(plate_info["well_index"])
            writer.writerow(plate_info.keys())
            for i in range(limit):
                writer.writerow([plate_info[x][i] for x in plate_info.keys()])
        f.close()

    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def getBatchReactionProductMWs(batch_id: int, reaction_number: int) -> list[float]:
    """Gets the SMILES of the products for a list of reactions
       for a batch

    Parameters
    ----------
    batch_id: int
        The batch id to get the target molecular weights for
    reaction_number: int
        The reactions to find product SMILES for

    Returns
    -------
    product_SMILES: list
        The product SMILES for a batch of reactions
    """
    product_SMILES = []
    batchobj = Batch.objects.get(id=batch_id)
    targetqueryset = batchobj.targets.all().order_by("id")
    for targetobj in targetqueryset:
        methodqueryset = targetobj.methods.all().order_by("id")
        for methodobj in methodqueryset:
            reactionqueryset = (
                methodobj.reactions.all().filter(number=reaction_number).order_by("id")
            )
            for reactionobj in reactionqueryset:
                product_SMILES = product_SMILES + list(
                    reactionobj.products.all()
                    .order_by("id")
                    .values_list("smiles", flat=True)
                )
    product_MWs = getMWs(smiles=product_SMILES)
    return product_MWs


def getPlateQuerySet(plate_id: int = None, otsession_id: int = None) -> QuerySet[Plate]:
    if plate_id:
        platequeryset = Plate.objects.filter(id=plate_id)
    if otsession_id:
        platequeryset = Plate.objects.filter(otsession_id=otsession_id)
    return platequeryset


def getProductQuerySet(reaction_ids: list) -> QuerySet[Product]:
    """Get product queryset for reaction ids

    Parameters
    ----------
    reaction_ids: list
        The reaction ids to search for related products

    Returns
    -------
    productqueryset: QuerySet[Product]
        The product queryset related to the reaction ids
    """
    productqueryset = Product.objects.filter(reaction_id__in=reaction_ids)
    return productqueryset


def getProduct(reaction_id: int) -> Product:
    """Get product object

    Parameters
    ----------
    reaction_id: int
        The reaction id to search for a matching product

    Returns
    -------
    productobj: Product
        The product Django model object
    """
    productobj = Product.objects.get(reaction_id=reaction_id)
    return productobj


def getProductSmiles(reaction_ids: list) -> list:
    """Get product smiles of reactions

    Parameters
    ----------
    reaction_ids: list
        The reactions to get product smiles for

    Returns
    -------
    productsmiles: list
        The list of product smiles
    """

    productsmiles = Product.objects.filter(reaction_id__in=reaction_ids).values_list(
        "smiles", flat=True
    )
    if productsmiles:
        return list(productsmiles)
    else:
        return None


def getReaction(reaction_id: int) -> Reaction:
    """Get reaction object

    Parameters
    ----------
    reaction_id: int
        The reaction id to search for a reaction

    Returns
    -------
    reactionobj: Reaction
        The reaction Django model object
    """
    reactionobj = Reaction.objects.get(id=reaction_id)
    return reactionobj


def getReactionQuerySet(
    reaction_ids: list = None, method_id: int = None
) -> QuerySet[Reaction]:
    """Get a  synthesis methods reactions

    Parameters
    ----------
    reaction_id: int or Reaction
        The reaction ids to find reactions for
    method_id: int
        The optional synthesis method's id to get reactions for

    Returns
    -------
    reactionqueryset: QuerySet[Reaction]
        The reactions of a synthesis method
    """
    if reaction_ids:
        reactionqueryset = Reaction.objects.filter(id__in=reaction_ids).order_by("id")
    if method_id:
        reactionqueryset = Reaction.objects.filter(method_id=method_id).order_by("id")
    return reactionqueryset


def checkProceedingReactions(reaction_id: int) -> QuerySet[Reaction]:
    """Checks if there are any reactions that proceed the reaction

    Parameters
    ----------
    reaction_id: int
        The reaction id of the Django model object to search for
        all relative proceeding reactions objects

    Returns
    -------
    proceedingreactionqueryset: QuerySet[Reaction]
        Returns the reactions that proceed the reaction
    """
    reactionobj = getReaction(reaction_id=reaction_id)
    proceedingreactionqueryset = Method.objects.get(
        id=reactionobj.method_id.id
    ).reactions.filter(id__gt=reaction_id)
    return proceedingreactionqueryset


def getReactionYields(reactionclasslist: list, recipelist) -> list[int]:
    """Gets the reaction yields

    Parameters
    ----------
    reactionclasslist: list
        The reaction classes to find yields for
    recipelist: list
        The list of recipes to find yields for

    Returns
    -------
    reactionyields: list[float]
        Returns the reaction yields eg. 0.80
    """
    reactionyields = [
        (encoded_recipes[reactionclass]["recipes"][recipe]["yield"] / 100)
        for reactionclass, recipe in zip(reactionclasslist, recipelist)
    ]
    return reactionyields


def checkPreviousReactionProducts(reaction_id: int, smiles: str) -> bool:
    """Checks if any previous reactions had a product matching the smiles

    Parameters
    ----------
    reaction_id: int
        The reaction id of the Django model object to search for
        all relative previous reactions objects. The previous reactions may
        have products that are this reaction's reactant input
    smiles: str
        The SMILES of the reaction's reactant and previous reaction products

    Returns
    -------
    status: bool
        The status is True if a match is found
    """
    reactionobj = getReaction(reaction_id=reaction_id)
    reactionqueryset = getReactionQuerySet(method_id=reactionobj.method_id.id)
    prevreactionqueryset = getPreviousObjEntries(
        queryset=reactionqueryset, obj=reactionobj
    )
    productmatches = []
    if prevreactionqueryset:
        for reactionobj in prevreactionqueryset:
            productobj = getProduct(reaction_id=reactionobj)
            if productobj.smiles == smiles:
                productmatches.append(productobj)
        if productmatches:
            return True
        else:
            return False
    else:
        return False


def getPreviousReactionQuerySets(reaction_id: int, smiles: str) -> QuerySet[Reaction]:
    """Checks if any previous reactions had a product matching the smiles

    Parameters
    ----------
    reaction_id: int
        The reaction id of the Django model object to search for
        all relative previous reactions objects. The previous reactions may
        have products that are this reaction's reactant input
    smiles: str
        The SMILES of the reaction's reactant and previous reaction products

    Returns
    -------
    previousreactionqueryset: QuerySet[Reaction]
        Returns the reactions that yiled products that match the SMILES searched
    """
    reactionobj = getReaction(reaction_id=reaction_id)
    previousreactionqueryset = Method.objects.get(
        id=reactionobj.method_id.id
    ).reactions.filter(id__lt=reaction_id, products__smiles=smiles)
    return previousreactionqueryset


def getMWs(smiles: list[str]) -> list[float]:
    """Gets the molecular weights of a list of compounds SMILES

    Parameters
    ----------
    smiles: list[str]
        The SMILES to calculate molecular weights for
    Returns
    -------
    MWs: list[float]
        The list of molecular weights
    """
    try:
        MWs = [
            Descriptors.MolWt(Chem.MolFromSmiles(smi)) for smi in smiles if smi != ""
        ]
        return MWs
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def getInchiKey(smiles: str) -> str:
    """Gets the inchikey for a compound SMILES

    Parameters
    ----------
    smiles: str
        The SMILESs to convert to an inchikey

    Returns
    -------
    inchikeys: str
        The inchikeys
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Could not create molecule from SMILES: {smiles}")
            return None

        # Use the correct function from the inchi module
        inchikey = Chem.inchi.MolToInchiKey(mol)
        return inchikey
    except Exception as e:
        logger.info(f"{inspect.stack()[0][3]} yielded error: {str(e)}")
        print(f"Failed to generate InChiKey for {smiles}: {e}")
        return None


def calculateMolsFromConc(target_concentration: float, target_volume: float) -> object:
    """Function to calculate product mols of reaction using a target mass

    Parameters
    ----------
    target_concentration: float
        The target concentration (mM) of the product
    target_volume: float
        The target volume (uL) of the product

    Returns
    -------
    product_moles: rdkit mol object
        The product mols
    """
    try:
        target_mols = (target_volume / 1e6) * (target_concentration / 1e3)
        return target_mols
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def calculateMassFromMols(mols: float, SMILES: str) -> object:
    """Function to calculate mass from mols

    Parameters
    ----------
    mols: float
        The mols of the compound
    SMILES: str
        The SMILES of the compound

    Returns
    -------
    mass: float
        The mass (mg) of the compound
    """
    try:
        MW = Descriptors.MolWt(Chem.MolFromSmiles(SMILES))
        mass = (mols * MW) * 1e3
        return mass
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def canonSmiles(smiles: str) -> str:
    """Function to canonicalise SMILES

    Parameters
    ----------
    smiles: str
        The SMILES to be canonicalised

    Returns
    -------
    canon_smiles: str
        The canonicalised SMILES
    status: bool
        Returns False if the input smiles could not be canonicalised
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            canon_smiles = Chem.MolToSmiles(mol)
            return canon_smiles
        else:
            return False
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def combiChem(
    reactant_1_SMILES: list, reactant_2_SMILES: list, are_product_SMILES: bool = False
) -> list:
    """Gets all possible combinations between two uneven lists of
       reactants

    Parameters
    ----------
    reactant_1_SMILES: list
        The list of reactant one smiles
    reactant_2_SMILES: list
        The second list of reactant two smiles
    are_product_SMILES: boolean
        Set to True if reactant_2_SMILES is list of products from
        previous reaction step

    Returns
    -------
    all_possible_combinations: list
        All possible reactant combinations possible
        between reactat 1 and reactant two lists
        as a list of tuples
    """
    if len(reactant_1_SMILES) == 0:
        reactant_2_SMILES_canon = [canonSmiles(smi) for smi in reactant_2_SMILES]
        if not are_product_SMILES:
            reactant_2_SMILES_canon = list(dict.fromkeys(reactant_2_SMILES_canon))
        all_possible_combinations = list(
            itertools.product([""], reactant_2_SMILES_canon)
        )
    if len(reactant_2_SMILES) == 0:
        reactant_1_SMILES_canon = [canonSmiles(smi) for smi in reactant_1_SMILES]
        if not are_product_SMILES:
            reactant_1_SMILES_canon = list(dict.fromkeys(reactant_1_SMILES_canon))
        all_possible_combinations = list(
            itertools.product([""], reactant_1_SMILES_canon)
        )
    if len(reactant_1_SMILES) != 0 and len(reactant_2_SMILES) != 0:
        reactant_1_SMILES_canon = [canonSmiles(smi) for smi in reactant_1_SMILES]
        reactant_2_SMILES_canon = [canonSmiles(smi) for smi in reactant_2_SMILES]
        if not are_product_SMILES:
            reactant_2_SMILES_canon = list(dict.fromkeys(reactant_2_SMILES_canon))
        all_possible_combinations = list(
            itertools.product(
                list(dict.fromkeys(reactant_1_SMILES_canon)), reactant_2_SMILES_canon
            )
        )
    return all_possible_combinations


def createCombiChemCSV(csv_input_file: str, out_dir: str):
    """Creates a .csv file for all the combinations possible for a given input
        of reactant SMILES pairs

    Parameters
    ----------
    csv_input_file: str
        The path to .csv file to read the reactant SMILES pairs and reactant classes from
    out_dir: str
        The directory to write the csv to
    """
    try:
        output_list = []
        input_df = pd.read_csv(csv_input_file)
        grouped_df = input_df.groupby(["reactant_class", "reaction_recipe"])
        for group in grouped_df:
            reaction_classes = group[1]["reactant_class"].tolist()
            reaction_recipes = group[1]["reaction_recipe"].tolist()
            reactant_1_SMILES = group[1]["reactant_1"].tolist()
            reactant_2_SMILES = group[1]["reactant_2"].tolist()
            reaction_SMARTS = encoded_recipes[reaction_classes[0]][reaction_recipes[0]][
                "reactionSMARTS"
            ]

            all_possible_combinations = combiChem(reactant_1_SMILES, reactant_2_SMILES)
            product_smiles = []
            for reactant_pair in all_possible_combinations:
                product_mols = checkReactantSMARTS(
                    reactant_SMILES=reactant_pair, reaction_SMARTS=reaction_SMARTS
                )
                product_smiles.append(Chem.MolToSmiles(product_mols[0]))

            all_possible_reactant_1_smiles = [x[0] for x in all_possible_combinations]
            all_possible_reactant_2_smiles = [x[1] for x in all_possible_combinations]

            output_list.append(
                [
                    all_possible_reactant_1_smiles,
                    all_possible_reactant_2_smiles,
                    product_smiles,
                    reaction_classes,
                    reaction_recipes,
                ]
            )
        out_df = pd.DataFrame(
            output_list,
            columns=[
                "reactant_1_SMILES",
                "reactant_2_SMILES",
                "product_SMILES",
                "reaction_class",
                "reaction_recipe",
            ],
        )
        out_df.write_csv(
            out_dir
            + "{}-combi-chem.csv".format(datetime.now().strftime("%Y-%m-%d-%H-%M-%S"))
        )

    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def createSVGString(smiles: str) -> str:
    """Function that creates a SVG image string from smiles string

    Parameters
    ----------
    smiles: string
        The SMILES to create an SVG image string from

    Returns
    -------
    svg_string: string
        The SVG image string
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(100, 50)
        drawer.SetFontSize(8)
        drawer.SetLineWidth(1)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_string = drawer.GetDrawingText()
        return svg_string
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def createReactionSVGString(smarts: str) -> str:
    """Function that creates a SVG image string from smarts string

    Parameters
    ----------
    smarts: string
        The SMARTS reaction pattern to create an SVG image
        string from

    Returns
    -------
    svg_string: string
        The SVG image string of the SMARTS pattern
    """
    try:
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(900, 200)
        drawer.DrawReaction(smarts)
        drawer.FinishDrawing()
        svg_string = drawer.GetDrawingText()
        return svg_string
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def getAddtionOrder(
    product_smi: str, reactant_SMILES: tuple, reaction_SMARTS: list
) -> list:
    """Gets reactant pair addition order as SMILES that yields the expected
       product via the reaction SMARTS pattern

    Parameters
    ----------
    product_smi: str
        The product SMILES
    reactant_SMILES_pair: tuple
        The reactant SMILES pair for a reaction
    reaction_SMARTS: lists
        The reaction SMARTS pattern. List for multiple SMARTS patterns/reaction transformations.

    Returns
    -------
    reactant_SMILES_pair: list
        The list of ordered reactant smiles
    status: None
        None if no order can create the input product
    """
    try:
        reaction_SMARTS = reaction_SMARTS[0]
        reaction_SMARTS_reactants = reaction_SMARTS.split(">>")[0].split(".")
        rxn = AllChem.ReactionFromSmarts(reaction_SMARTS)
        reactant_mols = [
            Chem.MolFromSmiles(smi) for smi in reactant_SMILES if smi != ""
        ]

        if reaction_SMARTS_reactants == "":
            ordered_smis = [canonSmiles(smi) for smi in reactant_SMILES if smi != ""]

        if len(reactant_mols) == 1:
            ordered_smis = [canonSmiles(smi) for smi in reactant_SMILES if smi != ""]

        if len(reactant_mols) > 1:

            if Chem.MolFromSmarts(reaction_SMARTS_reactants[0]) == Chem.MolFromSmarts(
                reaction_SMARTS_reactants[1]
            ):
                ordered_smis = [
                    canonSmiles(smi) for smi in reactant_SMILES if smi != ""
                ]
            else:
                for reactant_permutation in list(itertools.permutations(reactant_mols)):
                    try:
                        products = rxn.RunReactants(reactant_permutation)
                        product_mols = [product[0] for product in products]
                        if not product_mols:
                            continue  # reactants were in wrong order so no product
                    except Exception as e:
                        logger.info(
                            inspect.stack()[0][3] + " yielded error: {}".format(e)
                        )
                        print(e)
                        print(reactant_permutation)
                        continue
                    product_smis = [
                        Chem.MolToSmiles(m) for m in product_mols if m is not None
                    ]
                    if product_smi in product_smis:
                        ordered_smis = [
                            Chem.MolToSmiles(m) for m in reactant_permutation
                        ]
        if "ordered_smis" in locals():
            return ordered_smis
        else:
            print("Addition order not found for product SMILES".format(product_smi))
            print("The reaction SMARTS pattern is {}".format(reaction_SMARTS))
            print("The reactant SMILES are {}".format(reactant_SMILES))
            return None
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def checkReactantSMARTS(reactant_SMILES: tuple, reaction_SMARTS: list) -> list:
    """Checks if reactant pair can produce a product

    Parameters
    ----------
    reactant_SMILES_pair: tuple
        The pair of reactant smiles to check
    reaction_SMARTS: str
        The reaction SMARTS pattern used to check the reactant SMILES

    Returns
    -------
    product_mols: list
        The list of product mols formed between reactant SMILES from SMARTS pattern
    status: None
        Returns None if no product mols are formed
    """
    reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactant_SMILES if smi != ""]
    if len(reaction_SMARTS) == 1:
        SMARTS_pattern = reaction_SMARTS[0]
        rxn = AllChem.ReactionFromSmarts(SMARTS_pattern)
        if len(reactant_mols) == 1:
            try:
                products = rxn.RunReactants(reactant_mols)
                product_mols = [product[0] for product in products]
                if product_mols:
                    product_smiles = set(
                        [Chem.MolToSmiles(mol) for mol in product_mols]
                    )
                    product_mols = [Chem.MolFromSmiles(smi) for smi in product_smiles]
            except Exception as e:
                logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
                print(e)

        if len(reactant_mols) > 1:
            print("The number of reactant mols is > 1")
            for reactant_permutation in list(itertools.permutations(reactant_mols)):
                try:
                    products = rxn.RunReactants(reactant_permutation)
                    product_mols = [product[0] for product in products]
                    if product_mols:
                        print("The product mols are {}".format(product_mols))
                        product_smiles = set(
                            [Chem.MolToSmiles(mol) for mol in product_mols]
                        )
                        product_mols = [
                            Chem.MolFromSmiles(smi) for smi in product_smiles
                        ]
                        break
                    if not product_mols:
                        continue  # reactants were in wrong order so no product
                except Exception as e:
                    logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
                    print(e)
                    print(reactant_permutation)
                    continue
        if "product_mols" in locals() and len(product_mols) != 0:
            return [product_mols[0]]
        else:
            print(SMARTS_pattern)
            print(reactant_SMILES)
            return None

    if len(reaction_SMARTS) > 1:
        product_mols = []
        if len(reactant_mols) == 1:
            print("The number of reactant mols is 1")
            try:
                for SMARTS_pattern in reaction_SMARTS:
                    rxn = AllChem.ReactionFromSmarts(SMARTS_pattern)
                    if not "product_mol" in locals():
                        products = rxn.RunReactants(reactant_mols)
                    if "product_mol" in locals():
                        products = rxn.RunReactants(product_mol)
                    product_mol = [products[0][0]]
                    print("The product mol is {}".format(product_mol))
                    print(
                        "The product smiles is {}".format(
                            Chem.MolToSmiles(product_mol[0])
                        )
                    )
                    if product_mol:
                        product_mols.append(product_mol[-1])
            except Exception as e:
                print("The reactant smiles are {}".format(reactant_SMILES))
                print("The SMARTS pattern is {}".format(SMARTS_pattern))
                logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
                print(e)
        if len(reactant_mols) > 1:
            print("The number of reactant mols is > 1")
            try:
                for reactant_permutation in list(itertools.permutations(reactant_mols)):
                    for SMARTS_pattern in reaction_SMARTS:
                        rxn = AllChem.ReactionFromSmarts(SMARTS_pattern)
                        if not "product_mol" in locals:
                            products = rxn.RunReactants(reactant_mols)
                        if "product_mol" in locals():
                            products = rxn.RunReactants(product_mol)
                        product_mol = [products[0][0]]
                        print("The product mol is {}".format(product_mol))
                        print(
                            "The product smiles is {}".format(
                                Chem.MolToSmiles(product_mol)
                            )
                        )
                        if product_mol:
                            product_mols.extend(product_mol[-1])
                        if not product_mol:
                            continue
            except Exception as e:
                logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
                print(e)
                print(reactant_permutation)
        if len(product_mols) != 0:
            return product_mols
        else:
            print(SMARTS_pattern)
            print(reactant_SMILES)
            return None


def getPubChemCompound(inchikey: str) -> object:
    """Searches PubChem for compound using inchi key

    Parameters
    ----------
    inchikey: str
        The inchikey of the compound to search the PubChem DB for

    Returns
    -------
    compound: object
        The PuBChem compound class object
    status: None
        Returns None if no compound is found or an error occurs
    """
    try:
        compound = pcp.get_compounds(inchikey, "inchikey")[0]
        if not compound.cid:
            return None
        else:
            return compound
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(
            "Pubchempy could not retrieve compound entry for input inchikey: {} with error {}".format(
                inchikey, e
            )
        )


def getPubChemCAS(compound: object) -> str:
    """Get CAS identifier for PubChem compound synonyms

    Parameters
    ----------
    compound: object
        A PuBChem compound object

    Returns
    -------
    cas: str
        The CAS id of the compound
    """
    try:
        synonyms = compound.synonyms
        if synonyms:
            for syn in synonyms:
                match = re.match("(\d{1,7}-\d{1,2}-\d)", syn)
                if match:
                    cas = match.group(1)
                    return cas
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def getMolecularFormula(smiles: list) -> list:
    """Gets the molecular formula of a list of compounds SMILES
    Parameters
    ----------
    smiles: list[str]
        The SMILES to calculate molecular formula for
    Returns
    -------
    formula: list[str]
    """
    try:
        formula = [
            rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi)) for smi in smiles
        ]
        return formula
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def getChemicalName(inchikey: str) -> str:
    """Searches PubChem for compound using SMILES

    Parameters
    ----------
    inchikey: str
        The inchiley of the compound to search the PubChem DB for
        it's IUPAC name

    Returns
    -------
    name: str
        The IUPAC name of the compound
    status: None
        Returns None if no compound IUPAC name is found or if an error
        occurs
    """
    try:
        name = pcp.get_compounds(inchikey, "inchikey")[0].iupac_name
        if not name:
            return None
        else:
            return name
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(
            "Pubchempy could not convert SMILES to a IUPAC name with error {}".format(e)
        )


def matchSMARTS(smiles: str, smarts: str) -> bool:
    """Matches a SMILES pattern to a SMARTS pattern

    Parameters
    ----------
    smiles: str
        The SMILES pattern to match
    smarts: str
        The SMARTS pattern to match

    Returns
    -------
    status: bool
        Returns True if the SMILES pattern matches the SMARTS pattern
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return True
        else:
            return False
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def atomRemover(mol, rxn):
    """Remove atoms from a molecule using a reaction pattern.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule to modify
    rxn : rdkit.Chem.rdChemReactions.ChemicalReaction
        The reaction pattern to apply to the molecule

    Returns
    -------
    rdkit.Chem.rdchem.Mol or None
        The modified molecule with atoms removed, or None if the reaction fails

    Examples
    --------
    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> mol = Chem.MolFromSmiles('CC(=O)O')
    >>> rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])[OH]>>[C:1](=[O:2])')
    >>> result = atomRemover(mol, rxn)
    """
    try:
        ps = rxn.RunReactants((mol,))

        logger.debug(f"Attempting to run reaction on molecule: {Chem.MolToSmiles(mol)}")

        if not ps:
            logger.warning("Could not run the reaction, returning original molecule")
            return Chem.Mol(mol)

        for p in ps:
            res = Chem.RemoveHs(p[0])
            logger.info(f"Successfully removed atoms, result: {Chem.MolToSmiles(res)}")
            return res

    except Exception as e:
        logger.error(f"Error in atomRemover: {str(e)}")
        return None


def getFrags(mols: list, smarts: str) -> list:
    """Get the fragments of a list of molecules"
    Parameters
    ----------
    frags: list[rdkit.Chem.rdchem.Mol]
        The molecules to fragment

    Returns
    -------
    frags: list[rdkit.Chem.rdchem.Mol]
        The fragments of the input molecules
    """
    frag_mols = []
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        for mol in mols:
            try:
                ps = rxn.RunReactants((mol,))
                if not ps:
                    frag_mols.append(None)
                    continue
                for p in ps:
                    res = Chem.RemoveHs(p[0])
                    frag_mols.append(res)
            except Exception as e:
                logger.error(f"Error in getFrags: {str(e)}")
                frag_mols.append(None)
                continue
        return frag_mols
    except Exception as e:
        logger.error(f"Error in getFrags: {str(e)}")
        return None


def removeRadicals(mol):
    """Remove radicals from a molecule by adding hydrogens.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule that may contain radicals

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        The molecule with radicals removed, or None if the operation fails

    Examples
    --------
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles('[CH2]CC')
    >>> result = removeRadicals(mol)
    >>> Chem.MolToSmiles(result)
    'CCC'
    """
    try:
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                atom.SetNumRadicalElectrons(0)
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)

        Chem.SanitizeMol(mol)

        return mol

    except Exception as e:
        logger.error(f"Error removing radicals: {str(e)}")
        return None


def stripSalts(smiles: str, return_details: bool = False):
    """Strips salts from a SMILES string by returning the largest molecular fragment.

    Parameters
    ----------
    smiles : str
        The input SMILES string potentially containing salts
    return_details : bool, optional
        If True, returns additional information about salt stripping

    Returns
    -------
    str or tuple
        If return_details is False: The SMILES string with salts removed
        If return_details is True: A tuple (desalted_smiles, salts_removed, salt_fragments)
        Returns the original SMILES if processing fails

    Examples
    --------
    >>> strip_salts("CC(=O)O.Na")  # Sodium acetate
    'CC(=O)O'
    >>> strip_salts("CC(=O)O.[Na+]", True)  # Sodium acetate with details
    ('CC(=O)O', True, ['[Na+]'])
    """
    try:
        # Canonicalize input SMILES first
        canonical_smiles = canonSmiles(smiles)
        if not canonical_smiles:
            logger.warning(f"Could not canonicalize SMILES: {smiles}")
            return (smiles, False, []) if return_details else smiles

        # Convert to RDKit molecule
        mol = Chem.MolFromSmiles(canonical_smiles)
        if mol is None:
            logger.warning(f"Could not parse SMILES: {canonical_smiles}")
            return (canonical_smiles, False, []) if return_details else canonical_smiles

        # Get fragments
        fragments = Chem.GetMolFrags(mol, asMols=True)
        if len(fragments) <= 1:
            # No salts present
            return (canonical_smiles, False, []) if return_details else canonical_smiles

        # Find the largest fragment by molecular weight
        fragment_weights = [Descriptors.MolWt(frag) for frag in fragments]
        largest_idx = fragment_weights.index(max(fragment_weights))
        main_fragment = fragments[largest_idx]

        # Get salt fragments
        salt_fragments = []
        for i, frag in enumerate(fragments):
            if i != largest_idx:
                salt_fragments.append(Chem.MolToSmiles(frag))

        # Convert to canonical SMILES
        desalted_smiles = Chem.MolToSmiles(main_fragment)

        logger.info(f"Removed salts from {canonical_smiles} -> {desalted_smiles}")
        logger.info(f"Salt fragments: {salt_fragments}")

        if return_details:
            return desalted_smiles, True, salt_fragments
        else:
            return desalted_smiles

    except Exception as e:
        logger.error(f"Error stripping salts from {smiles}: {str(e)}")
        return (smiles, False, []) if return_details else smiles
