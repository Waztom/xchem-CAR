"""CAR's celery tasks"""
from __future__ import annotations
from celery import shared_task, current_task
from django.conf import settings
from django.db import transaction
import logging

logger = logging.getLogger(__name__)


from zipfile import ZipFile
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from itertools import count
import os

from .models import (
    Batch,
    SolventPrep,
    CompoundOrder,
    OTProject,
    OTBatchProtocol,
    OTScript,
    OTSession,
)

from .validate import ValidateFile
from .createmodels import (
    createProjectModel,
    createBatchModel,
    createTargetModel,
    createCatalogEntryModel,
    createMethodModel,
    createReactionModel,
    createProductModel,
    createReactantModel,
    CreateEncodedActionModels,
)

from .manifold.apicalls import (
    getExactSearch,
    getManifoldRetrosynthesisBatch,
)
from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import (
    checkPreviousReactionProducts,
    getActionSessionQuerySet,
    getActionSessionSequenceNumbers,
    getActionSessionTypes,
    getAddtionOrder,
    canonSmiles,
    getBatchReactions,
    getBatchTag,
    getGroupedActionSessionSequences,
    getGroupedActionSessionTypes,
    getMaxReactionNumber,
    getOTBatchProtocolQuerySet,
    getReactionsToDo,
    groupReactions,
)

from .opentrons.otsession import CreateOTSession
from .opentrons.otwrite import OTWrite


def delete_tmp_file(filepath):
    os.remove(filepath)


@shared_task
def validateFileUpload(
    csv_fp, validate_type: str = None, project_info=None, validate_only=True
):
    """Celery task to process validate the uploaded files for retrosynthesis planning.

    Parameters
    ----------
    csv_fp: str
        filepath of the uploaded csv file, which is saved to temporary storage by `viewer.views.UploadCSV`
    validate_type: validate different types of upload files
    project_info: dict
        dictionary of project details (name, email and project_name) that will be used to create the project model
        if the csv file is validated
    validate_only: boolean
        set to True to delete tmp file saved for validation and set to False to save tmp file and keep
        for uploading and creating model objects

    Returns
    -------
    validate_output: tuple
        contains the following:
            - validate dict (dict): dict containing any errors found during the validation step
            - validated (bool): True if the file(s) were validated, False if not
            - filename (str): name of the uploaded csv file
    """

    validation = ValidateFile(csv_to_validate=csv_fp, validate_type=validate_type)

    if validate_only:
        delete_tmp_file(csv_fp)
        csv_fp = None

    uploaded_dict = validation.df.to_dict("list")
    validated = validation.validated
    validate_dict = validation.validate_dict

    return (
        validate_type,
        validate_dict,
        validated,
        project_info,
        csv_fp,
        uploaded_dict,
    )


@shared_task
def uploadManifoldReaction(validate_output):

    _, validate_dict, validated, project_info, csv_fp, uploaded_dict = validate_output
    uploaded_df = pd.DataFrame(uploaded_dict)

    if not validated:
        delete_tmp_file(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        project_id = createProjectModel(project_info)
        project_info["project_id"] = project_id

        grouped_targets = uploaded_df.groupby("batch-tag")

        for batchtag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batchtag=batchtag,
            )
            target_names = list(group["target-names"])
            target_smiles = list(group["target-SMILES"])
            target_concentrations = list(group["concentration-required-mM"])
            target_volumes = list(group["amount-required-uL"])
            for i in range(
                0, len(target_smiles), 10
            ):  # Manifold can do 10 smiles in one batch query
                target_names_10 = target_names[i : i + 10]
                target_smiles_10 = target_smiles[i : i + 10]
                target_cocentrations_10 = target_concentrations[i : i + 10]
                target_volumes_10 = target_volumes[i : i + 10]

                retrosynthesis_results = getManifoldRetrosynthesisBatch(
                    smiles=target_smiles_10
                )
                if "results" in retrosynthesis_results:
                    batchrouteresults = retrosynthesis_results["results"]

                    for (
                        target_name,
                        target_smi,
                        target_concentration,
                        target_volume,
                        routeresult,
                    ) in zip(
                        target_names_10,
                        target_smiles_10,
                        target_cocentrations_10,
                        target_volumes_10,
                        batchrouteresults,
                    ):
                        # Need to capture when route is not found
                        if "routes" in routeresult:
                            routes = routeresult["routes"]
                            if routes:
                                target_id = createTargetModel(
                                    batch_id=batch_id,
                                    name=target_name,
                                    smiles=target_smi,
                                    concentration=target_concentration,
                                    volume=target_volume,
                                )
                                first_route = routes[0]

                                if first_route["molecules"][0]["isBuildingBlock"]:
                                    catalog_entries = first_route["molecules"][0][
                                        "catalogEntries"
                                    ]
                                    for catalog_entry in catalog_entries:
                                        createCatalogEntryModel(
                                            catalog_entry=catalog_entry,
                                            target_id=target_id,
                                        )

                                for route in routes[1:]:
                                    no_steps = len(route["reactions"])
                                    reactions = route["reactions"]
                                    encoded_reactions_found = [
                                        reaction
                                        for reaction in reactions
                                        if reaction["name"] in encoded_recipes
                                    ]

                                    if len(encoded_reactions_found) != no_steps:
                                        # Introduce a combinatorial calc of all recipe combos possible
                                        # between reactions and then create many recipe methods
                                        method_id = createMethodModel(
                                            target_id=target_id,
                                            nosteps=no_steps,
                                            otchem=False,
                                        )

                                        for index, reaction in enumerate(
                                            reversed(reactions)
                                        ):
                                            reaction_name = reaction["name"]
                                            reactant_smiles = reaction["reactantSmiles"]
                                            product_smiles = reaction["productSmiles"]
                                            reaction_smarts = (
                                                AllChem.ReactionFromSmarts(
                                                    "{}>>{}".format(
                                                        ".".join(reactant_smiles),
                                                        product_smiles,
                                                    ),
                                                    useSmiles=True,
                                                )
                                            )
                                            if len(reactant_smiles) == 1:
                                                intramolecular = True
                                            if len(reactant_smiles) == 2:
                                                intramolecular = False

                                            reaction_id = createReactionModel(
                                                method_id=method_id,
                                                reaction_class=reaction_name,
                                                reaction_number=index + 1,
                                                intramolecular=intramolecular,
                                                reaction_smarts=reaction_smarts,
                                            )

                                            createProductModel(
                                                reaction_id=reaction_id,
                                                product_smiles=product_smiles,
                                            )

                                            for reactant_smi in reactant_smiles:
                                                previousreactionqueryset = (
                                                    checkPreviousReactionProducts(
                                                        reaction_id=reaction_id,
                                                        smiles=reactant_smi,
                                                    )
                                                )
                                                if previousreactionqueryset:
                                                    previous_reaction_product = True
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    createCatalogEntryModel(
                                                        reactant_id=reactant_id,
                                                        previous_reaction_product=True,
                                                    )

                                                if not previousreactionqueryset:
                                                    previous_reaction_product = False
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    catalog_entries = [
                                                        molecule["catalogEntries"]
                                                        for molecule in route[
                                                            "molecules"
                                                        ]
                                                        if molecule["smiles"]
                                                        == reactant_smi
                                                    ][0]
                                                    for (
                                                        catalog_entry
                                                    ) in catalog_entries:
                                                        createCatalogEntryModel(
                                                            catalog_entry=catalog_entry,
                                                            reactant_id=reactant_id,
                                                        )

                                    if len(encoded_reactions_found) == no_steps:
                                        method_id = createMethodModel(
                                            target_id=target_id,
                                            nosteps=no_steps,
                                            otchem=True,
                                        )

                                        for index, reaction in enumerate(
                                            reversed(reactions)
                                        ):
                                            reaction_name = reaction["name"]
                                            reactant_smiles = reaction["reactantSmiles"]
                                            product_smiles = reaction["productSmiles"]
                                            reaction_temperature = [
                                                actionsession["actions"][0]["content"][
                                                    "temperature"
                                                ]["value"]
                                                for actionsession in encoded_recipes[
                                                    reaction_name
                                                ]["recipes"]["standard"][
                                                    "actionsessions"
                                                ]
                                                if actionsession["type"] == "stir"
                                            ][0]
                                            intramolecular_possible = encoded_recipes[
                                                reaction_name
                                            ]["intramolecular"]

                                            # recipe_rxn_smarts = encoded_recipes[
                                            #     reaction_name
                                            # ]["recipes"]["standard"]["reactionSMARTS"]

                                            if (
                                                len(reactant_smiles) == 1
                                                and intramolecular_possible
                                            ):
                                                # reactant_smiles_ordered = (
                                                #     reactant_smiles
                                                # )
                                                intramolecular = True
                                            else:
                                                # reactant_smiles_ordered = getAddtionOrder(
                                                #     product_smi=product_smiles,
                                                #     reactant_SMILES=reactant_smiles,
                                                #     reaction_SMARTS=recipe_rxn_smarts,
                                                # )
                                                intramolecular = False
                                                # if not reactant_smiles_ordered:
                                                #     continue

                                            reaction_smarts = (
                                                AllChem.ReactionFromSmarts(
                                                    "{}>>{}".format(
                                                        ".".join(reactant_smiles),
                                                        product_smiles,
                                                    ),
                                                    useSmiles=True,
                                                )
                                            )

                                            reaction_id = createReactionModel(
                                                method_id=method_id,
                                                reaction_class=reaction_name,
                                                reaction_number=index + 1,
                                                intramolecular=intramolecular,
                                                reaction_recipe="standard",  # Need to change when we get multiple recipes runnning!
                                                reaction_temperature=reaction_temperature,
                                                reaction_smarts=reaction_smarts,
                                            )

                                            createProductModel(
                                                reaction_id=reaction_id,
                                                product_smiles=product_smiles,
                                            )

                                            for reactant_smi in reactant_smiles:
                                                previousreactionqueryset = (
                                                    checkPreviousReactionProducts(
                                                        reaction_id=reaction_id,
                                                        smiles=reactant_smi,
                                                    )
                                                )
                                                if previousreactionqueryset:
                                                    previous_reaction_product = True
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    createCatalogEntryModel(
                                                        reactant_id=reactant_id,
                                                        previous_reaction_product=True,
                                                    )

                                                if not previousreactionqueryset:
                                                    previous_reaction_product = False
                                                    reactant_id = createReactantModel(
                                                        reaction_id=reaction_id,
                                                        reactant_smiles=reactant_smi,
                                                        previous_reaction_product=previous_reaction_product,
                                                    )
                                                    catalog_entries = [
                                                        molecule["catalogEntries"]
                                                        for molecule in route[
                                                            "molecules"
                                                        ]
                                                        if molecule["smiles"]
                                                        == reactant_smi
                                                    ][0]
                                                    for (
                                                        catalog_entry
                                                    ) in catalog_entries:
                                                        createCatalogEntryModel(
                                                            catalog_entry=catalog_entry,
                                                            reactant_id=reactant_id,
                                                        )

        delete_tmp_file(csv_fp)

        return validate_dict, validated, project_info


@shared_task
def uploadCustomReaction(validate_output):
    (
        _,
        validate_dict,
        validated,
        project_info,
        csv_fp,
        uploaded_dict,
    ) = validate_output
    uploaded_df = pd.DataFrame(uploaded_dict)

    if not validated:
        delete_tmp_file(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        project_id = createProjectModel(project_info)
        project_info["project_id"] = project_id

        grouped_targets = uploaded_df.groupby("batch-tag")
        for batchtag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batchtag=batchtag,
            )

            for (
                target_name,
                target_smiles,
                target_concentration,
                target_volume,
                no_steps,
                reactant_pair_smiles_tuples,
                reaction_groupby_column_tuples,
                reaction_name_tuples,
                reaction_recipe_tuples,
                reaction_product_smiles_tuples,
            ) in zip(
                group["target-names"],
                group["target-SMILES"],
                group["concentration-required-mM"],
                group["amount-required-uL"],
                group["no-steps"],
                group["reactant-pair-smiles"],
                group["reaction-groupby-column"],
                group["reaction-name"],
                group["reaction-recipe"],
                group["product-smiles"],
            ):
                target_id = createTargetModel(
                    batch_id=batch_id,
                    name=target_name,
                    smiles=target_smiles,
                    concentration=target_concentration,
                    volume=target_volume,
                )
                method_id = createMethodModel(
                    target_id=target_id,
                    nosteps=no_steps,
                    otchem=True,
                )

                for (
                    index,
                    reactant_pair_smiles,
                    reaction_groupby_column,
                    reaction_name,
                    reaction_recipe,
                    reaction_product_smiles,
                ) in zip(
                    count(),
                    reactant_pair_smiles_tuples,
                    reaction_groupby_column_tuples,
                    reaction_name_tuples,
                    reaction_recipe_tuples,
                    reaction_product_smiles_tuples,
                ):
                    print("Reactant pair smiles: ", reactant_pair_smiles)
                    print("Reaction product smiles: ", reaction_product_smiles)
                    reaction_smarts = AllChem.ReactionFromSmarts(
                        "{}>>{}".format(
                            ".".join(reactant_pair_smiles), reaction_product_smiles
                        ),
                        useSmiles=True,
                    )
                    # Need to move this to OT session - this will enable changing reaction temps in recipes
                    # and not have to upload again...
                    reaction_temperature = [
                        actionsession["actions"][0]["content"]["temperature"]["value"]
                        for actionsession in encoded_recipes[reaction_name]["recipes"][
                            reaction_recipe
                        ]["actionsessions"]
                        if actionsession["type"] == "stir"
                    ][0]

                    reaction_id = createReactionModel(
                        method_id=method_id,
                        reaction_class=reaction_name,
                        reaction_number=index + 1,
                        intramolecular=False,
                        reaction_recipe=reaction_recipe,
                        reaction_temperature=reaction_temperature,
                        reaction_smarts=reaction_smarts,
                        groupby_column=reaction_groupby_column,
                    )

                    createProductModel(
                        reaction_id=reaction_id,
                        product_smiles=reaction_product_smiles,
                    )

                    for reactant_smi in reactant_pair_smiles:
                        previousreactionqueryset = checkPreviousReactionProducts(
                            reaction_id=reaction_id,
                            smiles=reactant_smi,
                        )
                        if previousreactionqueryset:
                            reactant_id = createReactantModel(
                                reaction_id=reaction_id,
                                reactant_smiles=reactant_smi,
                                previous_reaction_product=True,
                            )
                            createCatalogEntryModel(
                                reactant_id=reactant_id,
                                previous_reaction_product=True,
                            )

                        if not previousreactionqueryset:
                            reactant_id = createReactantModel(
                                reaction_id=reaction_id,
                                reactant_smiles=reactant_smi,
                                previous_reaction_product=False,
                            )
                            #### Creating catalog entries takes very long!
                            createCatalogEntryModel(
                                reactant_id=reactant_id,
                                previous_reaction_product=False,
                                lab_inventory=True,
                            )
                            catalog_entries = getExactSearch(smiles=reactant_smi)
                            if "results" in catalog_entries:
                                for catalog_entry in catalog_entries["results"]:
                                    createCatalogEntryModel(
                                        catalog_entry=catalog_entry,
                                        reactant_id=reactant_id,
                                        previous_reaction_product=False,
                                    )

    delete_tmp_file(csv_fp)

    return validate_dict, validated, project_info


@shared_task
def uploadCombiCustomReaction(validate_output):
    (
        _,
        validate_dict,
        validated,
        project_info,
        csv_fp,
        uploaded_dict,
    ) = validate_output
    uploaded_df = pd.DataFrame(uploaded_dict)

    if not validated:
        delete_tmp_file(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        project_id = createProjectModel(project_info)
        project_info["project_id"] = project_id

        grouped_targets = uploaded_df.groupby("batch-tag")
        for batchtag, group in grouped_targets:
            batch_id = createBatchModel(
                project_id=project_id,
                batchtag=batchtag,
            )

            for (
                target_name,
                target_smiles,
                target_concentration,
                target_volume,
                no_steps,
                reactant_pair_smiles_tuples,
                reaction_name_tuples,
                reaction_recipe_tuples,
                reaction_product_smiles_tuples,
            ) in zip(
                group["target-names"],
                group["target-SMILES"],
                group["concentration-required-mM"],
                group["amount-required-uL"],
                group["no-steps"],
                group["reactant-pair-smiles"],
                group["reaction-name"],
                group["reaction-recipe"],
                group["product-smiles"],
            ):
                target_id = createTargetModel(
                    batch_id=batch_id,
                    name=target_name,
                    smiles=target_smiles,
                    concentration=target_concentration,
                    volume=target_volume,
                )
                method_id = createMethodModel(
                    target_id=target_id,
                    nosteps=no_steps,
                    otchem=True,
                )

                for (
                    index,
                    reactant_pair_smiles,
                    reaction_name,
                    reaction_recipe,
                    reaction_product_smiles,
                ) in zip(
                    count(),
                    reactant_pair_smiles_tuples,
                    reaction_name_tuples,
                    reaction_recipe_tuples,
                    reaction_product_smiles_tuples,
                ):
                    reaction_smarts = AllChem.ReactionFromSmarts(
                        "{}>>{}".format(
                            ".".join(reactant_pair_smiles), reaction_product_smiles
                        ),
                        useSmiles=True,
                    )
                    # Need to move this to OT session - this will enable changing reaction temps in recipes
                    # and not have to upload again...
                    reaction_temperature = [
                        actionsession["actions"][0]["content"]["temperature"]["value"]
                        for actionsession in encoded_recipes[reaction_name]["recipes"][
                            reaction_recipe
                        ]["actionsessions"]
                        if actionsession["type"] == "stir"
                    ][0]

                    reaction_id = createReactionModel(
                        method_id=method_id,
                        reaction_class=reaction_name,
                        reaction_number=index + 1,
                        intramolecular=False,
                        reaction_recipe=reaction_recipe,
                        reaction_temperature=reaction_temperature,
                        reaction_smarts=reaction_smarts,
                    )

                    createProductModel(
                        reaction_id=reaction_id,
                        product_smiles=reaction_product_smiles,
                    )

                    for reactant_smi in reactant_pair_smiles:
                        previousreactionqueryset = checkPreviousReactionProducts(
                            reaction_id=reaction_id,
                            smiles=reactant_smi,
                        )
                        if previousreactionqueryset:
                            reactant_id = createReactantModel(
                                reaction_id=reaction_id,
                                reactant_smiles=reactant_smi,
                                previous_reaction_product=True,
                            )
                            createCatalogEntryModel(
                                reactant_id=reactant_id,
                                previous_reaction_product=True,
                            )

                        if not previousreactionqueryset:
                            reactant_id = createReactantModel(
                                reaction_id=reaction_id,
                                reactant_smiles=reactant_smi,
                                previous_reaction_product=False,
                            )
                            #### Creating catalog entries takes very long!
                            createCatalogEntryModel(
                                reactant_id=reactant_id,
                                previous_reaction_product=False,
                                lab_inventory=True,
                            )
                            catalog_entries = getExactSearch(smiles=reactant_smi)
                            if "results" in catalog_entries:
                                for catalog_entry in catalog_entries["results"]:
                                    createCatalogEntryModel(
                                        catalog_entry=catalog_entry,
                                        reactant_id=reactant_id,
                                        previous_reaction_product=False,
                                    )

    delete_tmp_file(csv_fp)

    return validate_dict, validated, project_info


@shared_task
def createOTScript(batchids: list, protocol_name: str, custom_SM_files: dict = None):
    """
    Create otscripts and starting plates for a list of batch ids

    Parameters
    ----------
    batchids: list
        list of batch ids to create otscripts for
    protocol_name: str
        name of the protocol to create
    custom_SM_files: dict
        dictionary of optional custom SM files to use for the protocol
        Keys are batch IDs (as strings) and values are paths to CSV files
        Format: {"batch_id": "path/to/csv"}
    """
    task_summary = {}
    otprojectobj = OTProject()
    projectobj = Batch.objects.get(id=batchids[0]).project_id
    otprojectobj.project_id = projectobj
    otprojectobj.name = protocol_name
    otprojectobj.save()

    for batchid in batchids:
        reactionqueryset = getBatchReactions(batchid=batchid)
        actionsessionqueryset = getActionSessionQuerySet(reaction_ids=reactionqueryset)
        if not actionsessionqueryset:
            for reactionobj in reactionqueryset:
                reaction_id = reactionobj.id
                reactionclass = reactionobj.reactionclass
                reactionrecipe = reactionobj.recipe
                target_id = reactionobj.method_id.target_id.id
                reactant_pair_smiles = reactionobj.reactants.all().values_list(
                    "smiles", flat=True
                )
                recipe = encoded_recipes[reactionclass]["recipes"][reactionrecipe]
                intramolecular_possible = encoded_recipes[reactionclass][
                    "intramolecular"
                ]

                if intramolecular_possible and len(reactant_pair_smiles) == 1:
                    intramolecular_possible = True
                else:
                    intramolecular_possible = False

                actionsessions = recipe["actionsessions"]

                CreateEncodedActionModels(
                    intramolecular=intramolecular_possible,
                    actionsessions=actionsessions,
                    target_id=target_id,
                    reaction_id=reaction_id,
                    reactant_pair_smiles=list(reactant_pair_smiles),
                    reaction_name=reactionclass,
                )

        batchtag = getBatchTag(batchid=batchid)
        otbatchprotocolqueryset = getOTBatchProtocolQuerySet(batch_id=batchid)

        if otbatchprotocolqueryset and otbatchprotocolqueryset[0].zipfile:
            otbatchprotocolobj = otbatchprotocolqueryset[0]
            otbatchprotocolobj.otproject_id = otprojectobj
            otbatchprotocolobj.celery_taskid = current_task.request.id
            otbatchprotocolobj.save()
            task_summary[batchid] = True
        else:
            otbatchprotocolobj = OTBatchProtocol()
            otbatchprotocolobj.batch_id = Batch.objects.get(id=batchid)
            otbatchprotocolobj.otproject_id = otprojectobj
            otbatchprotocolobj.celery_taskid = current_task.request.id
            otbatchprotocolobj.save()
            maxreactionnumber = getMaxReactionNumber(reactionqueryset=reactionqueryset)
            groupedreactionquerysets = groupReactions(
                reactionqueryset=reactionqueryset, maxreactionnumber=maxreactionnumber
            )
            for index, groupreactionqueryset in enumerate(groupedreactionquerysets):
                if index == 0:
                    reaction_ids = [reaction.id for reaction in groupreactionqueryset]
                    actionsessionqueryset = getActionSessionQuerySet(
                        reaction_ids=reaction_ids
                    )
                    sessionnumbers = getActionSessionSequenceNumbers(
                        actionsessionqueryset=actionsessionqueryset
                    )
                    groupedactionsessionsequences = getGroupedActionSessionSequences(
                        sessionnumbers=sessionnumbers,
                        actionsessionqueryset=actionsessionqueryset,
                    )
                    for groupactionsession in groupedactionsessionsequences:
                        actionsessiontypes = getActionSessionTypes(
                            actionsessionqueryset=groupactionsession
                        )
                        groupedactionsessiontypes = getGroupedActionSessionTypes(
                            actionsessiontypes=actionsessiontypes,
                            actionsessionqueryset=groupactionsession,
                        )
                        for groupactionsessiontype in groupedactionsessiontypes:
                            human_actionsessionqueryset = groupactionsessiontype.filter(
                                driver="human"
                            )
                            robot_actionsessionqueryset = groupactionsessiontype.filter(
                                driver="robot"
                            )
                            if human_actionsessionqueryset:
                                pass
                            if robot_actionsessionqueryset:
                                actionsession_ids = (
                                    robot_actionsessionqueryset.values_list(
                                        "id", flat=True
                                    )
                                )

                                logger.info(
                                    f"The custom starting material files are: {custom_SM_files}"
                                )
                                # Get custom starting material CSV path if it exists for this batch
                                custom_sm_csv_path = None
                                if custom_SM_files and str(batchid) in custom_SM_files:
                                    custom_sm_csv_path = custom_SM_files[str(batchid)]
                                    logger.info(
                                        f"Using custom SMILES CSV for batch {batchid}: {custom_sm_csv_path}"
                                    )

                                # Create multiple sessions for reactions
                                sessions = createMultipleOTSessions(
                                    reactionstep=index + 1,
                                    otbatchprotocolobj=otbatchprotocolobj,
                                    actionsessionqueryset=robot_actionsessionqueryset,
                                    customSMcsvpath=custom_sm_csv_path,
                                    batchtag=batchtag,  # Pass batchtag
                                )

                if index > 0:
                    groupreactiontodoqueryset = getReactionsToDo(
                        groupreactionqueryset=groupreactionqueryset
                    )
                    if len(groupreactiontodoqueryset) == 0:
                        break
                    else:
                        reaction_ids = [
                            reaction.id for reaction in groupreactiontodoqueryset
                        ]
                        actionsessionqueryset = getActionSessionQuerySet(
                            reaction_ids=reaction_ids
                        )
                        sessionnumbers = getActionSessionSequenceNumbers(
                            actionsessionqueryset=actionsessionqueryset
                        )
                        groupedactionsessionsequences = (
                            getGroupedActionSessionSequences(
                                sessionnumbers=sessionnumbers,
                                actionsessionqueryset=actionsessionqueryset,
                            )
                        )
                        for groupactionsession in groupedactionsessionsequences:
                            actionsessiontypes = getActionSessionTypes(
                                actionsessionqueryset=groupactionsession
                            )
                            groupedactionsessiontypes = getGroupedActionSessionTypes(
                                actionsessiontypes=actionsessiontypes,
                                actionsessionqueryset=groupactionsession,
                            )
                            for groupactionsessiontype in groupedactionsessiontypes:
                                human_actionsessionqueryset = (
                                    groupactionsessiontype.filter(driver="human")
                                )
                                robot_actionsessionqueryset = (
                                    groupactionsessiontype.filter(driver="robot")
                                )
                                if human_actionsessionqueryset:
                                    pass
                                if robot_actionsessionqueryset:
                                    actionsession_ids = (
                                        robot_actionsessionqueryset.values_list(
                                            "id", flat=True
                                        )
                                    )

                                    # Get custom starting material CSV path if it exists for this batch
                                    custom_sm_csv_path = None
                                    if (
                                        custom_SM_files
                                        and str(batchid) in custom_SM_files
                                    ):
                                        custom_sm_csv_path = custom_SM_files[
                                            str(batchid)
                                        ]

                                    # Create multiple sessions for reactions
                                    sessions = createMultipleOTSessions(
                                        reactionstep=index + 1,
                                        otbatchprotocolobj=otbatchprotocolobj,
                                        actionsessionqueryset=robot_actionsessionqueryset,
                                        customSMcsvpath=custom_sm_csv_path,
                                        batchtag=batchtag,  # Pass batchtag
                                    )

            createZipOTBatchProtocol = ZipOTBatchProtocol(
                otbatchprotocolobj=otbatchprotocolobj, batchtag=batchtag
            )

            if createZipOTBatchProtocol.errors:
                task_summary[batchid] = False
            else:
                task_summary[batchid] = True
    else:
        task_summary[batchid] = False

    # Clean up temporary CSV files if they exist
    if custom_SM_files:
        for filepath in custom_SM_files.values():
            try:
                if os.path.exists(filepath):
                    os.remove(filepath)
            except Exception as e:
                print(f"Error cleaning up file {filepath}: {str(e)}")

    return task_summary, otprojectobj.id


def createMultipleOTSessions(
    reactionstep: int,
    otbatchprotocolobj: OTBatchProtocol,
    actionsessionqueryset: QuerySet[ActionSession],
    customSMcsvpath: str = None,
    max_reactions_per_session: int = None,
    batchtag: str = None,  # Add batchtag parameter
):
    """Splits a large reaction set into multiple OT sessions to prevent deck overflow.

    For each session, immediately performs OTWrite to ensure database changes from
    one session are available to subsequent sessions.
    """
    # Get unique reaction IDs from action sessions
    reaction_ids = list(
        set(actionsessionqueryset.values_list("reaction_id", flat=True))
    )
    total_reactions = len(reaction_ids)

    logger.info(
        f"Processing {total_reactions} unique reactions across {actionsessionqueryset.count()} action sessions"
    )

    if max_reactions_per_session is None:
        # Start with a smaller default to avoid issues
        max_reactions_per_session = 600 if total_reactions > 600 else total_reactions

    # Calculate initial number of sessions needed
    num_sessions = (
        total_reactions + max_reactions_per_session - 1
    ) // max_reactions_per_session

    logger.info(
        f"Initial plan: Creating {num_sessions} sessions with max {max_reactions_per_session} reactions per session"
    )

    created_sessions = []
    reaction_groups = []

    # Create initial reaction groups
    for i in range(0, total_reactions, max_reactions_per_session):
        batch = reaction_ids[i : i + max_reactions_per_session]
        reaction_groups.append(batch)

    # Try to create sessions for each reaction group
    for group_index, reaction_group in enumerate(reaction_groups):
        try:
            # Find action sessions related to this reaction group
            group_action_sessions = actionsessionqueryset.filter(
                reaction_id__in=reaction_group
            )
            group_action_session_ids = list(
                group_action_sessions.values_list("id", flat=True)
            )

            logger.info(
                f"Attempting session {group_index + 1} with {len(reaction_group)} reactions"
            )

            # Create a session for this group
            session = CreateOTSession(
                reactionstep=reactionstep,
                otbatchprotocolobj=otbatchprotocolobj,
                actionsessionqueryset=group_action_sessions,
                customSMcsvpath=customSMcsvpath,
            )

            # IMPORTANT: Immediately run OTWrite for this session before proceeding
            # This ensures database changes from this write are available to next sessions
            session_batchtag = (
                f"{batchtag}_session_{session.otsessionobj.id}"
                if batchtag
                else f"session_{session.otsessionobj.id}"
            )

            logger.info(f"Running OTWrite for session {group_index + 1}")
            otwrite_result = OTWrite(
                batchtag=session_batchtag,
                otsessionobj=session.otsessionobj,
                actionsession_ids=group_action_session_ids,
            )

            # If we get here, everything worked - store session info
            created_sessions.append(
                {
                    "session": session,
                    "action_session_ids": group_action_session_ids,
                    "otwrite_result": otwrite_result,
                }
            )

            logger.info(f"Successfully created and wrote session {group_index + 1}")

        except ValueError as ve:
            # Handle deck slot issues by splitting the group
            error_msg = str(ve)
            logger.error(f"ValueError in session {group_index + 1}: {error_msg}")
            if (
                error_msg == "No deck slots available - cannot create more plates"
                and len(reaction_group) > 1
            ):
                logger.warning(f"Deck full error detected: {str(ve)}")
                logger.warning(
                    f"Splitting group of {len(reaction_group)} reactions in half"
                )

                # Recursive call with smaller max_reactions_per_session for just this group
                smaller_max = max(1, len(reaction_group) // 2)
                logger.info(f"Trying with batch size: {smaller_max}")

                sub_action_sessions = actionsessionqueryset.filter(
                    reaction_id__in=reaction_group
                )
                sub_sessions = createMultipleOTSessions(
                    reactionstep=reactionstep,
                    otbatchprotocolobj=otbatchprotocolobj,
                    actionsessionqueryset=sub_action_sessions,
                    customSMcsvpath=customSMcsvpath,
                    max_reactions_per_session=smaller_max,
                    batchtag=batchtag,  # Pass the batchtag
                )
                created_sessions.extend(sub_sessions)
            else:
                # If not a deck slot issue or already at minimum group size, propagate the error
                logger.error(f"ValueError in session {group_index + 1}: {str(ve)}")
                raise

    logger.info(
        f"Created {len(created_sessions)} sessions to handle {total_reactions} reactions"
    )
    return created_sessions


class ZipOTBatchProtocol(object):
    """
    Creates a ZipBatchProtocol object for writing compound order csvs and otscripts
    for a batch
    """

    def __init__(self, otbatchprotocolobj: object, batchtag: str):
        """
        zipOTBatchProtocol constructor
        Args:
            otbatchprotocolobj (Django object): OT batch protocol object created for a batch
            batchtag (str): Batch tag for naming zip file
        """
        self.errors = {"function": [], "errorwarning": []}
        self.otbatchprotocolobj = otbatchprotocolobj
        self.mediaroot = settings.MEDIA_ROOT
        self.otsessionqueryset = self.getOTSessionQuerySet()
        self.zipfn = "batch-{}-protocol.zip".format(batchtag)
        self.ziptmpfp = os.path.join(settings.MEDIA_ROOT, "tmp", "batchprotocoltmp.zip")
        self.ziparchive = ZipFile(self.ziptmpfp, "w")

        for otsession_obj in self.otsessionqueryset:
            solventprepqueryset = self.getSolventPrepQuerySet(
                otsessionobj=otsession_obj
            )
            compoundorderqueryset = self.getCompoundOrderQuerySet(
                otsessionobj=otsession_obj
            )
            otscriptqueryset = self.getOTScriptQuerySet(otsessionobj=otsession_obj)

            if solventprepqueryset:
                for solventprepobj in solventprepqueryset:
                    filepath = self.getSolventPrepFilePath(
                        solventprepobj=solventprepobj
                    )
                    destdir = "solventprep"
                    self.writeZip(destdir=destdir, filepath=filepath)

            if compoundorderqueryset:
                for compoundorderobj in compoundorderqueryset:
                    filepath = self.getCompoundOrderFilePath(
                        compoundorderobj=compoundorderobj
                    )
                    destdir = "compoundorders"
                    self.writeZip(destdir=destdir, filepath=filepath)

            if otscriptqueryset:
                for otscriptobj in otscriptqueryset:
                    filepath = self.getOTScriptFilePath(otscriptobj=otscriptobj)
                    destdir = "otscripts"
                    self.writeZip(destdir=destdir, filepath=filepath)

        self.ziparchive.close()
        self.writeZipToMedia()
        self.deleteTmpZip()

    def addWarning(self, function: str, errorwarning: str):
        self.errors["function"].append(function)
        self.errors["errorwarning"].append(errorwarning)

    def getOTSessionQuerySet(self):
        """Retrieve OTSession model queryset"""
        otsessionqueryset = OTSession.objects.filter(
            otbatchprotocol_id=self.otbatchprotocolobj
        )

        if not otsessionqueryset:
            self.addWarning(
                function=self.getOTSessionQuerySet.__name__,
                errorwarning="No queryset found",
            )
        else:
            return otsessionqueryset

    def getSolventPrepQuerySet(self, otsessionobj: object):
        """Retrieve SolventPrep model queryset
        Args:
            otsessionobj (Django obj): OTSession Django object
        """
        solventprepqueryset = SolventPrep.objects.filter(otsession_id=otsessionobj)

        if not solventprepqueryset:
            self.addWarning(
                function=self.getSolventPrepQuerySet.__name__,
                errorwarning="No queryset found",
            )
            return None
        else:
            return solventprepqueryset

    def getCompoundOrderQuerySet(self, otsessionobj: object):
        """Retrieve CompoundOrder model queryset
        Args:
            otsessionobj (Django obj): OTSession Django object
        """
        compoundorderqueryset = CompoundOrder.objects.filter(otsession_id=otsessionobj)

        if not compoundorderqueryset:
            self.addWarning(
                function=self.getCompoundOrderQuerySet.__name__,
                errorwarning="No queryset found",
            )
        else:
            return compoundorderqueryset

    def getOTScriptQuerySet(self, otsessionobj: object):
        """Retrieve OTScript model queryset
        Args:
            otsessionobj (Django obj): OTScript Django object
        """
        otscriptqueryset = OTScript.objects.filter(otsession_id=otsessionobj)

        if not otscriptqueryset:
            self.addWarning(
                function=self.getOTScriptQuerySet.__name__,
                errorwarning="No queryset found",
            )
        else:
            return otscriptqueryset

    def getSolventPrepFilePath(self, solventprepobj: object):
        """Retrieve SolventPrep csv file path
        Args:
            solventprepobj (Django obj): SolventPrep Django object
        """
        filepath = os.path.join(self.mediaroot, solventprepobj.solventprepcsv.name)
        return filepath

    def getCompoundOrderFilePath(self, compoundorderobj: object):
        """Retrieve CompoundOrder csv file path
        Args:
            compoundorderobj (Django obj): OTSession Django object
        """
        filepath = os.path.join(self.mediaroot, compoundorderobj.ordercsv.name)
        return filepath

    def getOTScriptFilePath(self, otscriptobj: object):
        """Retrieve OTScript Python file path
        Args:
            otscriptobj (Django obj): OTScript Django object
        """
        filepath = os.path.join(self.mediaroot, otscriptobj.otscript.name)
        return filepath

    def writeZip(self, destdir: str, filepath: str):
        """Add the requested file to the zip archive.

        Args:
            destdir (str): directory to write file to in ziparchive
            filepath (str): filepath from record
        """
        arcname = os.path.join(destdir, filepath.split("/")[-1])
        self.ziparchive.write(filename=filepath, arcname=arcname)

    def writeZipToMedia(self):
        """Write the ziparchive to medida for the
        OTBatchProtocol Django object
        """
        zf = open(self.ziptmpfp, "rb")
        self.otbatchprotocolobj.zipfile.save(self.zipfn, zf)
        self.otbatchprotocolobj.save()

    def deleteTmpZip(self):
        """ "Delete the temporary zip archive created by ZipFile"""
        os.remove(self.ziptmpfp)


@shared_task
def canonicalizeSmiles(csvfile: str = None, smiles: list = None):
    """ "
    Canonicalizes smiles from csv file uploaded from frontend
    """
    validated = True

    if csvfile:
        csvdf = pd.read_csv(csvfile, encoding="utf8")
        smiles = [smi for smi in csvdf["SMILES"]]
        delete_tmp_file(csvfile)
    if smiles:
        smiles = smiles

    molcheck = [Chem.MolFromSmiles(smi) for smi in smiles]

    if None not in molcheck:
        canonicalizedsmiles = [canonSmiles(smi) for smi in smiles]
        return validated, canonicalizedsmiles
    else:
        validated = False
        indexerrors = [i for i, v in enumerate(molcheck) if v is None]
        errorsummary = "There was an error with the smiles csv at index: {}".format(
            indexerrors
        )
        return validated, errorsummary
