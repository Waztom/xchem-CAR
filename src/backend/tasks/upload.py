"""Upload tasks: Manifold, custom-chem, and combi-custom-chem."""
from __future__ import annotations

from itertools import count

import pandas as pd
from celery import shared_task
from django.db import transaction
from rdkit.Chem import AllChem

from ..createmodels import (
    create_batch_model,
    create_catalog_entry_model,
    create_method_model,
    create_product_model,
    create_project_model,
    create_reactant_model,
    create_reaction_model,
    create_target_model,
)
from ..db_utils import check_previous_reaction_products
from ..manifold.apicalls import (
    get_exact_search,
    get_manifold_retrosynthesis_batch,
)
from ..recipe_utils import (
    get_recipe_intramolecular,
    get_recipe_stir_temperature,
    recipe_exists,
)
from .validation import delete_tmp_file


# ---------------------------------------------------------------------------
# Shared upload helpers
# ---------------------------------------------------------------------------


def _unpack_validate_output(validate_output):
    """Unpack the standard 6-element validation tuple.

    Returns
    -------
    tuple
        (validate_dict, validated, project_info, csv_fp, uploaded_df)
    """
    _, validate_dict, validated, project_info, csv_fp, uploaded_dict = validate_output
    uploaded_df = pd.DataFrame(uploaded_dict)
    return validate_dict, validated, project_info, csv_fp, uploaded_df


def _handle_not_validated(validate_dict, project_info, csv_fp):
    """Short-circuit return when validation failed: delete tmp file and return."""
    delete_tmp_file(csv_fp)
    return (validate_dict, False, project_info)


def _create_reactant_and_catalog(
    reaction_id,
    reactant_smi,
    fetch_pubchem=True,
    fetch_catalogue=False,
    manifold_molecules=None,
):
    """Create a Reactant model and its associated CatalogEntry models.

    Handles both the "previous reaction product" and "new reactant" paths,
    and optionally fetches catalog entries from Manifold molecule data or
    the exact-search API.

    Parameters
    ----------
    reaction_id : int
        The parent Reaction model id.
    reactant_smi : str
        Canonical SMILES for the reactant.
    fetch_pubchem : bool
        Whether to fetch PubChem data for the reactant.
    fetch_catalogue : bool
        Whether to call get_exact_search for vendor catalog entries
        (used by custom-chem and combi uploads).
    manifold_molecules : list[dict] or None
        If provided (Manifold uploads), catalog entries are looked up from
        this list by SMILES match instead of using get_exact_search.
    """
    previousreactionqueryset = check_previous_reaction_products(
        reaction_id=reaction_id,
        smiles=reactant_smi,
    )
    if previousreactionqueryset:
        reactant_id = create_reactant_model(
            reaction_id=reaction_id,
            reactant_smiles=reactant_smi,
            previous_reaction_product=True,
            fetch_pubchem=fetch_pubchem,
        )
        create_catalog_entry_model(
            reactant_id=reactant_id,
            previous_reaction_product=True,
        )
    else:
        reactant_id = create_reactant_model(
            reaction_id=reaction_id,
            reactant_smiles=reactant_smi,
            previous_reaction_product=False,
            fetch_pubchem=fetch_pubchem,
        )
        if manifold_molecules is not None:
            # Manifold path: look up catalog entries from the route's molecule list
            catalog_entries = [
                molecule["catalogEntries"]
                for molecule in manifold_molecules
                if molecule["smiles"] == reactant_smi
            ][0]
            for catalog_entry in catalog_entries:
                create_catalog_entry_model(
                    catalog_entry=catalog_entry,
                    reactant_id=reactant_id,
                )
        else:
            # Custom / combi path: lab inventory entry + optional vendor search
            create_catalog_entry_model(
                reactant_id=reactant_id,
                previous_reaction_product=False,
                lab_inventory=True,
            )
            if fetch_catalogue:
                catalog_entries = get_exact_search(smiles=reactant_smi)
                if "results" in catalog_entries:
                    for catalog_entry in catalog_entries["results"]:
                        create_catalog_entry_model(
                            catalog_entry=catalog_entry,
                            reactant_id=reactant_id,
                            previous_reaction_product=False,
                        )


def _process_manifold_route_reactions(
    route, target_id, fetch_pubchem=True,
):
    """Process all reactions in a single Manifold route.

    Determines whether all reactions are encoded (otchem=True) or not
    (otchem=False), creates the Method, Reaction, Product, Reactant, and
    CatalogEntry models accordingly.

    Parameters
    ----------
    route : dict
        A Manifold route dict with 'reactions' and 'molecules' keys.
    target_id : int
        The parent Target model id.
    fetch_pubchem : bool
        Whether to fetch PubChem data for products/reactants.
    """
    reactions = route["reactions"]
    no_steps = len(reactions)
    molecules = route["molecules"]

    encoded_reactions_found = [
        reaction for reaction in reactions if recipe_exists(reaction["name"])
    ]
    all_encoded = len(encoded_reactions_found) == no_steps

    method_id = create_method_model(
        target_id=target_id,
        nosteps=no_steps,
        otchem=all_encoded,
    )

    for index, reaction in enumerate(reactions):
        reaction_name = reaction["name"]
        reactant_smiles = reaction["reactantSmiles"]
        product_smiles = reaction["productSmiles"]

        reaction_smarts = AllChem.ReactionFromSmarts(
            "{}>>{}".format(".".join(reactant_smiles), product_smiles),
            useSmiles=True,
        )

        # Determine intramolecular flag
        if all_encoded:
            intramolecular_possible = get_recipe_intramolecular(reaction_name)
            intramolecular = (
                len(reactant_smiles) == 1 and intramolecular_possible
            )
        else:
            if len(reactant_smiles) == 1:
                intramolecular = True
            elif len(reactant_smiles) == 2:
                intramolecular = False

        # Build create_reaction_model kwargs
        reaction_kwargs = dict(
            method_id=method_id,
            reaction_class=reaction_name,
            reaction_number=index + 1,
            intramolecular=intramolecular,
            reaction_smarts=reaction_smarts,
        )
        if all_encoded:
            reaction_kwargs["reaction_recipe"] = "standard"

        reaction_id = create_reaction_model(**reaction_kwargs)

        create_product_model(
            reaction_id=reaction_id,
            product_smiles=product_smiles,
            fetch_pubchem=fetch_pubchem,
        )

        for reactant_smi in reactant_smiles:
            _create_reactant_and_catalog(
                reaction_id=reaction_id,
                reactant_smi=reactant_smi,
                fetch_pubchem=fetch_pubchem,
                manifold_molecules=molecules,
            )


def _upload_custom_or_combi(
    validate_output,
    fetch_pubchem=True,
    fetch_catalogue=True,
    has_groupby_column=False,
):
    """Shared implementation for upload_custom_reaction and upload_combi_custom_reaction.

    Parameters
    ----------
    validate_output : tuple
        The standard 6-element validation output tuple.
    fetch_pubchem : bool
        Whether to fetch PubChem data for products/reactants.
    fetch_catalogue : bool
        Whether to call get_exact_search for vendor catalog entries.
    has_groupby_column : bool
        True for custom-chem uploads (which include a groupby column),
        False for combi-custom-chem uploads.
    """
    validate_dict, validated, project_info, csv_fp, uploaded_df = (
        _unpack_validate_output(validate_output)
    )

    if not validated:
        return _handle_not_validated(validate_dict, project_info, csv_fp)

    project_id = create_project_model(project_info)
    project_info["project_id"] = project_id

    grouped_targets = uploaded_df.groupby("batch-tag")
    for batchtag, group in grouped_targets:
        batch_id = create_batch_model(
            project_id=project_id,
            batchtag=batchtag,
        )

        # Build the per-target column lists; groupby column is optional
        target_columns = [
            group["target-names"],
            group["target-SMILES"],
            group["concentration-required-mM"],
            group["amount-required-uL"],
            group["no-steps"],
            group["reactant-pair-smiles"],
            group["reaction-name"],
            group["reaction-recipe"],
            group["product-smiles"],
        ]
        if has_groupby_column:
            target_columns.append(group["reaction-groupby-column"])

        for target_row in zip(*target_columns):
            # Unpack common fields
            (
                target_name,
                target_smiles,
                target_concentration,
                target_volume,
                no_steps,
                reactant_pair_smiles_tuples,
                reaction_name_tuples,
                reaction_recipe_tuples,
                reaction_product_smiles_tuples,
            ) = target_row[:9]

            reaction_groupby_column_tuples = (
                target_row[9] if has_groupby_column else None
            )

            target_id = create_target_model(
                batch_id=batch_id,
                name=target_name,
                smiles=target_smiles,
                concentration=target_concentration,
                volume=target_volume,
            )
            method_id = create_method_model(
                target_id=target_id,
                nosteps=no_steps,
                otchem=True,
            )

            # Build per-step iterables
            step_iterables = [
                reactant_pair_smiles_tuples,
                reaction_name_tuples,
                reaction_recipe_tuples,
                reaction_product_smiles_tuples,
            ]
            if has_groupby_column:
                step_iterables.append(reaction_groupby_column_tuples)

            for step_row in zip(count(), *step_iterables):
                index = step_row[0]
                reactant_pair_smiles = step_row[1]
                reaction_name = step_row[2]
                reaction_recipe = step_row[3]
                reaction_product_smiles = step_row[4]
                reaction_groupby_column = (
                    step_row[5] if has_groupby_column else None
                )

                reaction_smarts = AllChem.ReactionFromSmarts(
                    "{}>>{}".format(
                        ".".join(reactant_pair_smiles), reaction_product_smiles
                    ),
                    useSmiles=True,
                )
                reaction_temperature = get_recipe_stir_temperature(
                    reaction_name, reaction_recipe
                )

                reaction_kwargs = dict(
                    method_id=method_id,
                    reaction_class=reaction_name,
                    reaction_number=index + 1,
                    intramolecular=False,
                    reaction_recipe=reaction_recipe,
                    reaction_temperature=reaction_temperature,
                    reaction_smarts=reaction_smarts,
                )
                if has_groupby_column:
                    reaction_kwargs["groupby_column"] = reaction_groupby_column

                reaction_id = create_reaction_model(**reaction_kwargs)

                create_product_model(
                    reaction_id=reaction_id,
                    product_smiles=reaction_product_smiles,
                    fetch_pubchem=fetch_pubchem,
                )

                for reactant_smi in reactant_pair_smiles:
                    _create_reactant_and_catalog(
                        reaction_id=reaction_id,
                        reactant_smi=reactant_smi,
                        fetch_pubchem=fetch_pubchem,
                        fetch_catalogue=fetch_catalogue,
                    )

    delete_tmp_file(csv_fp)
    return validate_dict, validated, project_info


# ---------------------------------------------------------------------------
# Upload tasks
# ---------------------------------------------------------------------------


@shared_task
def upload_manifold_reaction(validate_output, fetch_pubchem=True):
    validate_dict, validated, project_info, csv_fp, uploaded_df = (
        _unpack_validate_output(validate_output)
    )

    if not validated:
        return _handle_not_validated(validate_dict, project_info, csv_fp)

    with transaction.atomic():
        project_id = create_project_model(project_info)
        project_info["project_id"] = project_id

        grouped_targets = uploaded_df.groupby("batch-tag")
        for batchtag, group in grouped_targets:
            batch_id = create_batch_model(
                project_id=project_id,
                batchtag=batchtag,
            )
            target_names = list(group["target-names"])
            target_smiles = list(group["target-SMILES"])
            target_concentrations = list(group["concentration-required-mM"])
            target_volumes = list(group["amount-required-uL"])

            # Manifold can do 10 smiles in one batch query
            for i in range(0, len(target_smiles), 10):
                target_names_10 = target_names[i : i + 10]
                target_smiles_10 = target_smiles[i : i + 10]
                target_concentrations_10 = target_concentrations[i : i + 10]
                target_volumes_10 = target_volumes[i : i + 10]

                retrosynthesis_results = get_manifold_retrosynthesis_batch(
                    smiles=target_smiles_10
                )
                if "results" not in retrosynthesis_results:
                    continue

                for (
                    target_name,
                    target_smi,
                    target_concentration,
                    target_volume,
                    routeresult,
                ) in zip(
                    target_names_10,
                    target_smiles_10,
                    target_concentrations_10,
                    target_volumes_10,
                    retrosynthesis_results["results"],
                ):
                    if "routes" not in routeresult:
                        continue
                    routes = routeresult["routes"]
                    if not routes:
                        continue

                    target_id = create_target_model(
                        batch_id=batch_id,
                        name=target_name,
                        smiles=target_smi,
                        concentration=target_concentration,
                        volume=target_volume,
                    )

                    # First route: building-block catalog entries
                    first_route = routes[0]
                    if first_route["molecules"][0]["isBuildingBlock"]:
                        for catalog_entry in first_route["molecules"][0]["catalogEntries"]:
                            create_catalog_entry_model(
                                catalog_entry=catalog_entry,
                                target_id=target_id,
                            )

                    # Remaining routes: reaction methods
                    for route in routes[1:]:
                        _process_manifold_route_reactions(
                            route=route,
                            target_id=target_id,
                            fetch_pubchem=fetch_pubchem,
                        )

    delete_tmp_file(csv_fp)
    return validate_dict, validated, project_info


@shared_task
def upload_custom_reaction(validate_output, fetch_pubchem=True, fetch_catalogue=True):
    return _upload_custom_or_combi(
        validate_output,
        fetch_pubchem=fetch_pubchem,
        fetch_catalogue=fetch_catalogue,
        has_groupby_column=True,
    )


@shared_task
def upload_combi_custom_reaction(validate_output, fetch_pubchem=True, fetch_catalogue=True):
    return _upload_custom_or_combi(
        validate_output,
        fetch_pubchem=fetch_pubchem,
        fetch_catalogue=fetch_catalogue,
        has_groupby_column=False,
    )
