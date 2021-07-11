from celery import shared_task
from django.core.files.storage import default_storage

import pandas as pd

from .validate import ValidateFile

from .IBM.createibmmodels import (
    createProjectModel,
    createTargetModel,
    createMethodModel,
    createReactionModel,
    createProductModel,
    createActionModel,
)
from .IBM.apicalls import createIBMProject, getIBMRetroSyn, collectIBMReactionInfo
from .IBM.filtermethod import filtermethod

from .manifold.apicalls import getManifoldretrosynthesis
from .recipebuilder.createencodedmodels import createEncodedActionModel
from .recipebuilder.encodedrecipes import encoded_recipes
from rdkit.Chem import AllChem


def delete_tmp_file(filepath):
    default_storage.delete(filepath)


# Celery validate task
@shared_task
def validateFileUpload(csv_fp, validate_type=None, project_info=None, validate_only=True):
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

    # Delete tempory file if only validate selected
    if validate_only:
        default_storage.delete(csv_fp)
        csv_fp = None

    # Convert dataframe to dictionary to make it JSON serializable
    uploaded_dict = validation.df.to_dict("list")
    validated = validation.validated
    validate_dict = validation.validate_dict

    return (validate_dict, validated, csv_fp, project_info, uploaded_dict)


@shared_task
def uploadIBMReaction(validate_output):
    # Validate output is a list - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    validate_dict, validated, csv_fp, project_info, uploaded_dict = validate_output

    if not validated:
        # Delete tempory file if only validate selected
        default_storage.delete(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        # Create project model and return project id
        project_id, project_name = createProjectModel(project_info)
        # Add project name to project info dict for emailing when upload complete
        project_info["project_name"] = project_name

        # Create an IBM project
        IBM_project_id, rxn4chemistry_wrapper = createIBMProject(project_name)

        target_no = 1
        for smiles, target_mass in zip(
            uploaded_dict["Targets"], uploaded_dict["Ammount_required (mg)"]
        ):

            # Number of steps/name of reactants from IBM API???????
            # Create a Target model and return target id
            target_id = createTargetModel(
                project_id=project_id,
                smiles=smiles,
                target_no=target_no,
                target_mass=target_mass,
            )

            # Create IBM
            # Run IBM API cal to get retrosyn info
            max_steps = 3
            results = getIBMRetroSyn(
                rxn4chemistry_wrapper=rxn4chemistry_wrapper,
                smiles=smiles,
                max_steps=max_steps,
            )

            # Check if IBM API has been successful
            if results:
                # Set maximum number of attempts. Sometimes variation in routes is
                # very similar and wastes times
                max_attempts = 3
                attempts_dict = {}

                # Set maximum number of methods/pathways to collect
                no_pathways_found = len(results["retrosynthetic_paths"])

                if no_pathways_found <= 3:
                    max_pathways = no_pathways_found
                if no_pathways_found > 3:
                    max_pathways = 3

                pathway_no = 1

                pathway_filter = []

                for pathway in results["retrosynthetic_paths"]:
                    print(pathway_no)
                    try:
                        attempts_dict[pathway_no] += 1
                    except:
                        attempts_dict[pathway_no] = 1

                    no_attempts = attempts_dict[pathway_no]

                    if no_attempts > max_attempts:
                        break

                    if pathway_no <= max_pathways and pathway["confidence"] > 0.90:

                        # Get reaction info about pathway
                        reaction_info = collectIBMReactionInfo(
                            rxn4chemistry_wrapper=rxn4chemistry_wrapper, pathway=pathway
                        )

                        # Check if reaction info call has been successful
                        if reaction_info:
                            # Need to add check if reaction class and reactants already exist -
                            # in other words looking for diversity and not dulpication
                            # of reactions. This ignores actions of the method
                            method_integer = filtermethod(reaction_info=reaction_info)

                            if method_integer not in pathway_filter:
                                pathway_no += 1

                                pathway_filter.append(method_integer)

                                # Create a Method model
                                method_id = createMethodModel(
                                    target_id=target_id,
                                    nosteps=max_steps,
                                )

                                product_no = 1
                                for product_smiles, reaction_class, actions, reaction_smarts in zip(
                                    reaction_info["product_smiles"],
                                    reaction_info["rclass"],
                                    reaction_info["actions"],
                                    reaction_info["reactions"],
                                ):
                                    # Product_smiles and reaction class is a list of individual elements
                                    # Reactants and actions is a list of lists

                                    # Create a Reaction model
                                    reaction_id = createReactionModel(
                                        method_id=method_id,
                                        reaction_class=reaction_class,
                                        reaction_smarts=reaction_smarts,
                                    )

                                    # Create a Product model
                                    createProductModel(
                                        reaction_id=reaction_id,
                                        project_name=project_name,
                                        target_no=target_no,
                                        pathway_no=pathway_no,
                                        product_no=product_no,
                                        product_smiles=product_smiles,
                                    )

                                    # Create action models
                                    action_no = 1
                                    for action in actions:
                                        created_model = createActionModel(
                                            reaction_id=reaction_id,
                                            action_no=action_no,
                                            action=action,
                                        )
                                        if created_model:
                                            action_no += 1

                                    product_no += 1

                    if pathway_no > max_pathways:
                        break

            target_no += 1

    # Delete tempory file

    default_storage.delete(csv_fp)

    return validate_dict, validated, project_info


@shared_task
def uploadManifoldReaction(validate_output):
    # Validate output is a list - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    validate_dict, validated, csv_fp, project_info, uploaded_dict = validate_output

    if not validated:
        # Delete tempory file if only validate selected
        default_storage.delete(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        # Create project model and return project id
        project_id, project_name = createProjectModel(project_info)
        # Add project name to project info dict for emailing when upload complete
        project_info["project_name"] = project_name

        # Do Postera stuff
        target_no = 1
        for target_smiles, target_mass in zip(
            uploaded_dict["Targets"], uploaded_dict["Ammount_required (mg)"]
        ):

            retrosynthesis_result = getManifoldretrosynthesis(target_smiles)
            routes = retrosynthesis_result["routes"]

            target_id = createTargetModel(
                project_id=project_id,
                smiles=target_smiles,
                target_no=target_no,
                target_mass=target_mass,
            )

            target_no += 1

            pathway_no = 1

            for route in routes:
                no_steps = len(route["reactions"])

                if no_steps > 0:
                    # Check if reactions are OT friendly
                    reactions = route["reactions"]

                    reactions_found = [
                        reaction for reaction in reactions if reaction["name"] in encoded_recipes
                    ]

                    if len(reactions_found) == no_steps:

                        method_id = createMethodModel(
                            target_id=target_id,
                            nosteps=no_steps,
                        )

                        # Then loop over route for synthetic steps
                        product_no = 1

                        for reaction in reactions:
                            reaction_class = reaction["name"]
                            # Check if OT friendly and in encoded recipes
                            print(reaction_class)
                            if reaction_class in encoded_recipes:
                                encoded_recipe = encoded_recipes[reaction_class]["recipe"]
                                reactant_SMILES = reaction["reactantSmiles"]
                                product_smiles = reaction["productSmiles"]

                                # Create a Reaction model
                                reaction_smarts = AllChem.ReactionFromSmarts(
                                    "{}>>{}".format(".".join(reactant_SMILES), product_smiles),
                                    useSmiles=True,
                                )
                                reaction_id = createReactionModel(
                                    method_id=method_id,
                                    reaction_class=reaction_class,
                                    reaction_smarts=reaction_smarts,
                                )

                                # Create a Product model
                                createProductModel(
                                    reaction_id=reaction_id,
                                    project_name=project_name,
                                    target_no=target_no,
                                    pathway_no=pathway_no,
                                    product_no=product_no,
                                    product_smiles=product_smiles,
                                )

                                for action in encoded_recipe:
                                    createEncodedActionModel(
                                        reaction_id=reaction_id,
                                        action=action,
                                        reactants=reactant_SMILES,
                                        target_id=target_id,
                                    )

                                product_no += 1

                            else:
                                pass

                pathway_no += 1

    default_storage.delete(csv_fp)

    return validate_dict, validated, project_info


@shared_task
def uploadCustomReaction(validate_output):
    # Validate output is a list - this is one way to get
    # Celery chaining to work where second function uses list output
    # from first function (validate) called
    validate_dict, validated, csv_fp, project_info, uploaded_dict = validate_output

    if not validated:
        # Delete tempory file if only validate selected
        default_storage.delete(csv_fp)
        return (validate_dict, validated, project_info)

    if validated:
        # Create project model and return project id
        project_id, project_name = createProjectModel(project_info)
        # Add project name to project info dict for emailing when upload complete
        project_info["project_name"] = project_name

        # Do custom chem stuff
        target_no = 1
        pathway_no = 1
        product_no = 1

        for reactant_pair_smiles, reaction_name, target_smiles, target_mass in zip(
            uploaded_dict["reactant_pair_smiles"],
            uploaded_dict["Reaction-name"],
            uploaded_dict["target-smiles"],
            uploaded_dict["Ammount_required (mg)"],
        ):
            target_id = createTargetModel(
                project_id=project_id,
                smiles=target_smiles,
                target_no=target_no,
                target_mass=target_mass,
            )

            target_no += 1

            method_id = createMethodModel(
                target_id=target_id,
                nosteps=1,
            )

            reaction_smarts = AllChem.ReactionFromSmarts(
                "{}>>{}".format(".".join(reactant_pair_smiles), target_smiles),
                useSmiles=True,
            )

            reaction_id = createReactionModel(
                method_id=method_id,
                reaction_class=reaction_name,
                reaction_smarts=reaction_smarts,
            )

            # Create a Product model
            createProductModel(
                reaction_id=reaction_id,
                project_name=project_name,
                target_no=target_no,
                pathway_no=pathway_no,
                product_no=product_no,
                product_smiles=target_smiles,
            )

            encoded_recipe = encoded_recipes[reaction_name]["recipe"]

            for action in encoded_recipe:
                createEncodedActionModel(
                    reaction_id=reaction_id,
                    action=action,
                    reactants=reactant_pair_smiles,
                    target_id=target_id,
                )

    default_storage.delete(csv_fp)

    return validate_dict, validated, project_info