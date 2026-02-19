"""Database query and update utilities for Django models.

This module contains functions for querying and updating Django ORM models
related to batches, reactions, targets, methods, actions, and plates.
"""

from django.db.models import QuerySet
from django.db.models import Q, Max

import csv
import inspect
import logging

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

logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Action QuerySets
# -----------------------------------------------------------------------------


def get_add_action_query_set(
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


def get_extract_action_query_set(
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


def get_ot_batch_protocol_query_set(batch_id: int) -> QuerySet[OTBatchProtocol]:
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


# -----------------------------------------------------------------------------
# ActionSession utilities
# -----------------------------------------------------------------------------


def get_action_session_sequence_numbers(
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


def get_action_session_types(actionsessionqueryset: QuerySet[ActionSession]) -> QuerySet:
    """Set of action session types

    Returns
    ------
    actionsessiontypes: QuerySet
        The set of action session types in a queryset
        eg. ["reaction", "workup", "stir"]
    """
    actionsessiontypes = set(list(actionsessionqueryset.values_list("type", flat=True)))
    return actionsessiontypes


def get_grouped_action_session_sequences(
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


def get_grouped_action_session_types(
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


def get_action_session_query_set(
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


# -----------------------------------------------------------------------------
# Query helpers
# -----------------------------------------------------------------------------


def get_previous_obj_entries(queryset: list, obj: object) -> QuerySet:
    """Finds all previous objects relative to obj of queryset"""
    previousqueryset = queryset.filter(pk__lt=obj.pk).order_by("-pk")
    return previousqueryset


# -----------------------------------------------------------------------------
# Reaction checks and queries
# -----------------------------------------------------------------------------


def check_previous_reaction_failures(reactionobj: Reaction) -> bool:
    """Check if any previous reaction failures for a method"""
    reactionqueryset = get_reactions(method_ids=[reactionobj.method_id.id])
    previousreactionqueryset = get_previous_obj_entries(
        queryset=reactionqueryset, obj=reactionobj
    )
    failedreactions = previousreactionqueryset.filter(success=False)
    if failedreactions.exists():
        return True
    else:
        return False


def check_no_method_steps(reactionobj: Reaction) -> bool:
    """Check no reaction steps in method is > 1"""
    methodobj = reactionobj.method_id
    noreactionsteps = methodobj.nosteps
    if noreactionsteps > 1:
        return True
    else:
        return False


def get_reactions_to_do(groupreactionqueryset: QuerySet[Reaction]) -> QuerySet[Reaction]:
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
        if check_no_method_steps(reactionobj=reactionobj):
            if not check_previous_reaction_failures(reactionobj=reactionobj):
                reactionstodo.append(reactionobj.id)
    groupreactiontodoqueryset = groupreactionqueryset.filter(id__in=reactionstodo)
    return groupreactiontodoqueryset


# -----------------------------------------------------------------------------
# Basic model queries
# -----------------------------------------------------------------------------


def get_targets(batch_ids: QuerySet[Batch]) -> QuerySet[Target]:
    targetqueryset = Target.objects.filter(batch_id__in=batch_ids).order_by("id")
    return targetqueryset


def get_methods(target_ids: QuerySet[Target]) -> QuerySet[Method]:
    methodqueryset = (
        Method.objects.filter(target_id__in=target_ids)
        .filter(otchem=True)
        .order_by("id")
    )
    return methodqueryset


def get_reactions(method_ids: QuerySet[Method]) -> QuerySet[Reaction]:
    reactionqueryset = Reaction.objects.filter(method_id__in=method_ids).order_by("id")
    return reactionqueryset


def get_batch_tag(batchid):
    batch_obj = Batch.objects.get(id=batchid)
    batchtag = batch_obj.batchtag
    return batchtag


def get_batch_reactions(batchid: int) -> QuerySet[Reaction]:
    targetqueryset = get_targets(batch_ids=[batchid])
    if targetqueryset:
        methodqueryset = get_methods(target_ids=targetqueryset)
        if methodqueryset:
            reactionqueryset = get_reactions(method_ids=methodqueryset)
            if reactionqueryset:
                return reactionqueryset


def get_max_reaction_number(reactionqueryset: QuerySet[Reaction]) -> int:
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


def group_reactions(reactionqueryset: QuerySet[Reaction], maxreactionnumber: int):
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


def get_reactants_to_buy(batch_ids: list[int]) -> list:
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


# -----------------------------------------------------------------------------
# Batch data retrieval
# -----------------------------------------------------------------------------


def get_batch_target_smiles(batch_id: int) -> list[float]:
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


def get_batch_reaction_ids(batch_id: int, reaction_number: int) -> list[float]:
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


def get_batch_reaction_product_smiles(batch_id: int, reaction_number: int) -> list[float]:
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


# -----------------------------------------------------------------------------
# Update operations
# -----------------------------------------------------------------------------


def update_target_mols(batch_id: int, concentration: float, volume: float) -> int:
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
    # Import here to avoid circular imports
    from .conversions import calculate_mols_from_conc, calculate_mass_from_mols

    target_qs = Batch.objects.get(id=batch_id).targets.all()
    for target_obj in target_qs:
        smiles = target_obj.smiles
        mols = calculate_mols_from_conc(
            target_concentration=concentration, target_volume=volume
        )
        mass = calculate_mass_from_mols(mols=mols, SMILES=smiles)
        target_obj.mols = mols
        target_obj.concentration = concentration
        target_obj.volume = volume
        target_obj.mass = mass
        target_obj.save()


def update_reaction_success_to_fail(reaction_ids: list[int]):
    """Updates reactions to be failures

    Parameters
    ----------
    reactions_ids: list[int]
        The reactions to update as synthetic failures
    """
    if Reaction.objects.filter(id__in=reaction_ids).exists():
        Reaction.objects.filter(id__in=reaction_ids).update(success=False)


def update_batch_method_ot_friendly(batch_id: int):
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


def delete_batch_action_sessions(batch_id: int):
    """Deletes all actions for a batch

    Parameters
    ----------
    batch_id: int
        The batch id to delete all actions for
    """
    reaction_qs = get_batch_reactions(batchid=batch_id)
    get_action_session_query_set(reaction_ids=reaction_qs).delete()


def update_recipe_type(
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
    reaction_qs = get_batch_reactions(batchid=batch_id).filter(
        reactionclass=reaction_class, recipe=current_recipe
    )
    for reaction_obj in reaction_qs:
        reaction_obj.recipe = recipe_to_use
        reaction_obj.save()


# -----------------------------------------------------------------------------
# Plate queries
# -----------------------------------------------------------------------------


def get_plate_query_set(plate_id: int = None, otsession_id: int = None) -> QuerySet[Plate]:
    if plate_id:
        platequeryset = Plate.objects.filter(id=plate_id)
    if otsession_id:
        platequeryset = Plate.objects.filter(otsession_id=otsession_id)
    return platequeryset


def get_plate_map(plate_ids: list, out_dir: str):
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
    # Import here to avoid circular imports
    from .chem_utils import get_mws

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
                target_mw = get_mws(smiles=[target_smi])[0]
                reactant_smiles = well.reaction_id.reactants.values_list(
                    "smiles", flat=True
                )
                reactant_mws = get_mws(smiles=reactant_smiles)
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


# -----------------------------------------------------------------------------
# Product queries
# -----------------------------------------------------------------------------


def get_product_query_set(reaction_ids: list) -> QuerySet[Product]:
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


def get_product(reaction_id: int) -> Product:
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


def get_product_smiles(reaction_ids: list) -> list:
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


# -----------------------------------------------------------------------------
# Reaction queries
# -----------------------------------------------------------------------------


def get_reaction(reaction_id: int) -> Reaction:
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


def get_reaction_temperature(reaction_id: int) -> float:
    """Get reaction temperature

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
    return reactionobj.temperature


def get_reaction_class(reaction_id: int) -> str:
    """Get reaction class

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
    return reactionobj.reactionclass


def get_reaction_recipe(reaction_id: int) -> str:
    """Get reaction recipe

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
    return reactionobj.recipe


def get_reaction_query_set(
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


def check_proceeding_reactions(reaction_id: int) -> QuerySet[Reaction]:
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
    reactionobj = get_reaction(reaction_id=reaction_id)
    proceedingreactionqueryset = Method.objects.get(
        id=reactionobj.method_id.id
    ).reactions.filter(id__gt=reaction_id)
    return proceedingreactionqueryset


def get_reaction_yields(reactionclasslist: list, recipelist) -> list[int]:
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
    from .recipe_utils import get_recipe_yield

    reactionyields = [
        get_recipe_yield(reactionclass, recipe)
        for reactionclass, recipe in zip(reactionclasslist, recipelist)
    ]
    return reactionyields


def check_previous_reaction_products(reaction_id: int, smiles: str) -> bool:
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
    reactionobj = get_reaction(reaction_id=reaction_id)
    reactionqueryset = get_reaction_query_set(method_id=reactionobj.method_id.id)
    prevreactionqueryset = get_previous_obj_entries(
        queryset=reactionqueryset, obj=reactionobj
    )
    productmatches = []
    if prevreactionqueryset:
        for reactionobj in prevreactionqueryset:
            productobj = get_product(reaction_id=reactionobj)
            if productobj.smiles == smiles:
                productmatches.append(productobj)
        if productmatches:
            return True
        else:
            return False
    else:
        return False


def get_previous_reaction_query_sets(reaction_id: int, smiles: str) -> QuerySet[Reaction]:
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
    reactionobj = get_reaction(reaction_id=reaction_id)
    previousreactionqueryset = Method.objects.get(
        id=reactionobj.method_id.id
    ).reactions.filter(id__lt=reaction_id, products__smiles=smiles)
    return previousreactionqueryset


# -----------------------------------------------------------------------------
# Batch MW calculations (requires chem_utils)
# -----------------------------------------------------------------------------


def get_batch_target_mws(batch_id: int) -> list[float]:
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
    from .chem_utils import get_mws

    batchobj = Batch.objects.get(id=batch_id)
    targetqueryset = batchobj.targets.all()
    smiles = [targetobj.smiles for targetobj in targetqueryset]
    target_MWs = get_mws(smiles=smiles)
    return target_MWs


def get_batch_reaction_product_mws(batch_id: int, reaction_number: int) -> list[float]:
    """Gets the MWs of the products for a list of reactions
       for a batch

    Parameters
    ----------
    batch_id: int
        The batch id to get the target molecular weights for
    reaction_number: int
        The reactions to find product SMILES for

    Returns
    -------
    product_MWs: list
        The product MWs for a batch of reactions
    """
    from .chem_utils import get_mws

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
    product_MWs = get_mws(smiles=product_SMILES)
    return product_MWs
