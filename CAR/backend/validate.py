from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .recipebuilder.encodedrecipes import encoded_recipes


def add_warning(field, warning_string, validate_dict):
    validate_dict["field"].append(field)
    validate_dict["warning_string"].append(warning_string)

    return validate_dict


def checkColumnNames(columns, validated, validate_dict, expected_no_columns):
    if expected_no_columns == 2:
        column_names = ["Targets", "Ammount_required (mg)"]
    if expected_no_columns == 4:
        column_names = ["Reactant-1", "Reactant-2", "Reaction-name", "Ammount_required (mg)"]
    if not all(columns == column_names):
        validate_dict = add_warning(
            field="name_columns",
            warning_string="Column names should be set to: {}".format(column_names),
            validate_dict=validate_dict,
        )
        validated = False
    return validated, validate_dict


def checkNumberColumns(columns, validated, validate_dict, expected_no_columns):
    no_columns = len(columns)

    if no_columns != expected_no_columns:
        validate_dict = add_warning(
            field="number_columns",
            warning_string="Found {} column names. Set and name columns to 'Targets' only".format(
                no_columns
            ),
            validate_dict=validate_dict,
        )
        validated = False

    if no_columns == expected_no_columns:
        validated, validate_dict = checkColumnNames(
            columns, validated, validate_dict, expected_no_columns
        )

    return validated, validate_dict


def checkSMILES(target_smiles, column_name, index, validated, validate_dict):

    mol = Chem.MolFromSmiles(target_smiles)

    if mol is None:
        validate_dict = add_warning(
            field="check_smiles",
            warning_string="Input target smiles: '{}' at index {} in column {} is not a valid smiles".format(
                target_smiles, index, column_name
            ),
            validate_dict=validate_dict,
        )
        validated = False

    return validated, validate_dict


def checkIsNumber(amount, index, validated, validate_dict):
    if not isinstance(amount, (int, float)):
        validate_dict = add_warning(
            field="check_number",
            warning_string="Target mass {} at index {} is not a valid number".format(amount, index),
            validate_dict=validate_dict,
        )
        validated = False

    return validated, validate_dict


def checkReaction(reactant_pair, reaction_name, index, validated, validate_dict):
    reaction_smarts = encoded_recipes[reaction_name]["reactionSMARTS"]
    reacts = [Chem.MolFromSmiles(smi) for smi in reactant_pair]

    for smarts in reaction_smarts:
        print(smarts)
        reaction = rdChemReactions.ReactionFromSmarts(smarts)
        # Perfom reaction
        products = reaction.RunReactants(reacts)
        if len(products) != 0:
            print(Chem.MolToSmiles(products[0][0]))

    if len(products) == 0:
        validate_dict = add_warning(
            field="check_reaction",
            warning_string="Reaction for reactants at index {} is not a valid reaction".format(
                index
            ),
            validate_dict=validate_dict,
        )
        validated = False
        product_smiles = "None"

    if len(products) != 0:
        product_mol = products[0][0]
        product_smiles = Chem.MolToSmiles(product_mol)

    return validated, validate_dict, product_smiles
