"""Checks validation of file for uploading to CAR"""
import math
from typing import BinaryIO
import inspect
import pandas as pd
from rdkit import Chem

from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import (
    getAddtionOrder,
    checkReactantSMARTS,
    combiChem,
)

import logging

logger = logging.getLogger(__name__)


class ValidateFile(object):
    """
    Creates a validate object for checking file validation for upload
    """

    def __init__(self, csv_to_validate: BinaryIO, validate_type: str):
        """ValidateFile constructor

        Parameters
        ----------

        csv_to_validate: IO
            The uploaded csv file for validating
        """
        self.df = pd.read_csv(csv_to_validate, encoding="utf8")
        self.df_columns = self.df.columns
        self.no_df_columns = len(self.df_columns)
        self.index_df_rows = range(0, len(self.df), 1)
        self.upload_type = validate_type
        self.validate_dict = {"field": [], "warning_string": []}
        self.validated = True

        if self.upload_type == "custom-chem":
            self.validateCustomChem()
        if self.upload_type == "combi-custom-chem":
            self.validateCustomCombiChem()

        if self.upload_type == "retro-API":
            expected_column_names = [
                "target-SMILES",
                "target-names",
                "concentration-required-mM",
                "amount-required-uL",
                "batch-tag",
            ]
            expected_no_columns = len(expected_column_names)
            self.checkNumberColumns(
                expected_no_columns=expected_no_columns,
                expected_column_names=expected_column_names,
            )
            if self.validated:
                self.checkColumnNames(expected_column_names=expected_column_names)
            if self.validated:
                self.target_smiles = [smi.strip() for smi in self.df["target-SMILES"]]
                self.df["target-SMILES"] = self.target_smiles
                self.checkSMILES(
                    df_rows_index=self.index_df_rows,
                    smiles=self.target_smiles,
                    smiles_type="target",
                )

    def validateCustomChem(self):
        max_no_steps = max(self.df["no-steps"])
        reaction_numbers = list(range(1, max_no_steps + 1))
        expected_groupby_column_names = [
            "reaction-groupby-column-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reactant_1_column_names = [
            "reactant-1-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reactant_2_column_names = [
            "reactant-2-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reaction_name_column_names = [
            "reaction-name-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reaction_recipe_column_names = [
            "reaction-recipe-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_product_column_names = [
            "reaction-product-smiles-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_column_names = (
            [
                "target-names",
                "no-steps",
                "concentration-required-mM",
                "amount-required-uL",
                "batch-tag",
            ]
            + expected_groupby_column_names
            + expected_reactant_1_column_names
            + expected_reactant_2_column_names
            + expected_product_column_names
            + expected_reaction_name_column_names
            + expected_reaction_recipe_column_names
        )
        expected_no_columns = len(expected_column_names)
        self.checkNumberColumns(
            expected_no_columns=expected_no_columns,
            expected_column_names=expected_column_names,
        )
        if self.validated:
            self.checkColumnNames(expected_column_names=expected_column_names)
        if self.validated:
            self.target_names = self.df["target-names"].tolist()
            self.batchtags = self.df["batch-tag"].tolist()
            self.concentrations = self.df["concentration-required-mM"].tolist()
            self.amounts = self.df["amount-required-uL"].tolist()
            self.nosteps = self.df["no-steps"].tolist()
            self.target_smiles = []
            self.product_smiles = []
            self.reactant_pair_smiles = []
            self.reaction_groupby_column = []
            self.reaction_names = []
            self.reaction_recipes = []
            all_reaction_info = {}
            for reaction_number in reaction_numbers:
                reaction_info = {}
                reaction_groupby_column = (
                    self.df["reaction-groupby-column-{}".format(reaction_number)]
                    .apply(bool)
                    .tolist()
                )
                reaction_names = self.df[
                    "reaction-name-{}".format(reaction_number)
                ].tolist()
                reaction_recipes = self.df[
                    "reaction-recipe-{}".format(reaction_number)
                ].tolist()
                product_smiles = self.df[
                    "reaction-product-smiles-{}".format(reaction_number)
                ].tolist()

                reactant_1_SMILES = [
                    "" if str(reactant) == "nan" else reactant.strip()
                    for reactant in self.df["reactant-1-{}".format(reaction_number)]
                ]

                reactant_2_SMILES = [
                    "" if str(reactant) == "nan" else reactant.strip()
                    for reactant in self.df["reactant-2-{}".format(reaction_number)]
                ]
                reactant_pair_smiles = list(zip(reactant_1_SMILES, reactant_2_SMILES))
                if self.validated:
                    reactant_pair_smiles_ordered, product_smiles = self.checkReaction(
                        reactant_pair_smiles=reactant_pair_smiles,
                        reaction_names=reaction_names,
                        reaction_recipes=reaction_recipes,
                        product_smiles=product_smiles,
                    )
                    print(
                        "The reactant pair smiles ordered are: ",
                        reactant_pair_smiles_ordered,
                    )
                    if reaction_number == max_no_steps:
                        self.target_smiles = self.target_smiles + product_smiles
                    reaction_info[
                        "reaction-groupby-column-{}".format(reaction_number)
                    ] = reaction_groupby_column
                    reaction_info[
                        "reaction-name-{}".format(reaction_number)
                    ] = reaction_names
                    reaction_info[
                        "reaction-recipe-{}".format(reaction_number)
                    ] = reaction_recipes
                    reaction_info[
                        "reaction-reactant-pair-smiles-{}".format(reaction_number)
                    ] = reactant_pair_smiles_ordered
                    reaction_info[
                        "reaction-product-smiles-{}".format(reaction_number)
                    ] = product_smiles
                    all_reaction_info.update(reaction_info)
        if self.validated:
            products = list(
                zip(
                    *[
                        all_reaction_info[
                            "reaction-product-smiles-{}".format(reactionnumber)
                        ]
                        for reactionnumber in reaction_numbers
                    ]
                )
            )
            reactant_pair_smiles = list(
                zip(
                    *[
                        all_reaction_info[
                            "reaction-reactant-pair-smiles-{}".format(reactionnumber)
                        ]
                        for reactionnumber in reaction_numbers
                    ]
                )
            )
            reaction_groupby_column = list(
                zip(
                    *[
                        all_reaction_info[
                            "reaction-groupby-column-{}".format(reactionnumber)
                        ]
                        for reactionnumber in reaction_numbers
                    ]
                )
            )
            reaction_names = list(
                zip(
                    *[
                        all_reaction_info["reaction-name-{}".format(reactionnumber)]
                        for reactionnumber in reaction_numbers
                    ]
                )
            )
            reaction_recipes = list(
                zip(
                    *[
                        all_reaction_info["reaction-recipe-{}".format(reactionnumber)]
                        for reactionnumber in reaction_numbers
                    ]
                )
            )

            self.product_smiles = self.product_smiles + products
            self.reactant_pair_smiles = self.reactant_pair_smiles + reactant_pair_smiles
            self.reaction_groupby_column = (
                self.reaction_groupby_column + reaction_groupby_column
            )
            self.reaction_names = self.reaction_names + reaction_names
            self.reaction_recipes = self.reaction_recipes + reaction_recipes

            self.df = pd.DataFrame()
            self.df["batch-tag"] = self.batchtags
            self.df["target-names"] = self.target_names
            self.df["target-SMILES"] = self.target_smiles
            self.df["concentration-required-mM"] = self.concentrations
            self.df["amount-required-uL"] = self.amounts
            self.df["no-steps"] = self.nosteps
            self.df["reactant-pair-smiles"] = self.reactant_pair_smiles
            self.df["reaction-groupby-column"] = self.reaction_groupby_column
            self.df["reaction-name"] = self.reaction_names
            self.df["reaction-recipe"] = self.reaction_recipes
            self.df["product-smiles"] = self.product_smiles
            self.checkIsNumber(values=self.concentrations)
            self.checkIsNumber(values=self.amounts)

    def validateCustomCombiChem(self):
        max_no_steps = int(max(self.df["no-steps"]))
        reaction_numbers = list(range(1, max_no_steps + 1))
        expected_reactant_1_column_names = [
            "reactant-1-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reaction_name_column_names = [
            "reaction-name-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_reaction_recipe_column_names = [
            "reaction-recipe-{}".format(reaction_number)
            for reaction_number in reaction_numbers
        ]
        expected_column_names = (
            [
                "combi-group",
                "no-steps",
                "concentration-required-mM",
                "amount-required-uL",
                "batch-tag",
                "reactant-2-1",
            ]
            + expected_reactant_1_column_names
            + expected_reaction_name_column_names
            + expected_reaction_recipe_column_names
        )
        expected_no_columns = len(expected_column_names)
        self.checkNumberColumns(
            expected_no_columns=expected_no_columns,
            expected_column_names=expected_column_names,
        )
        if self.validated:
            self.checkColumnNames(expected_column_names=expected_column_names)

        if self.validated:
            self.target_names = []
            self.target_smiles = []
            self.nosteps = []
            self.product_smiles = []
            self.reactant_pair_smiles = []
            self.reaction_names = []
            self.reaction_recipes = []
            self.batchtags = []
            self.concentrations = []
            self.amounts = []
            combi_grouped = self.df.groupby(["combi-group"])
            for combi_group_name, combi_group in combi_grouped:
                combi_group_info = {}
                combi_group = combi_group.reset_index()
                max_no_steps_combi_group = int(max(combi_group["no-steps"]))
                reaction_numbers_group = list(range(1, max_no_steps_combi_group + 1))
                columns_count = combi_group.nunique(
                    axis="rows", dropna=True
                )  # NB "nan" values (empty row values) not counted
                number_reactant_1s = [
                    columns_count["reactant-1-{}".format(reaction_number)]
                    for reaction_number in reaction_numbers_group
                    if "reactant-1-{}".format(reaction_number) in columns_count
                    and columns_count["reactant-1-{}".format(reaction_number)] != 0
                ]
                number_reactant_2s = [
                    columns_count["reactant-2-{}".format(reaction_number)]
                    for reaction_number in reaction_numbers_group
                    if "reactant-2-{}".format(reaction_number) in columns_count
                    and columns_count["reactant-1-{}".format(reaction_number)] != 0
                ]
                no_targets = int(math.prod(number_reactant_1s + number_reactant_2s))
                target_names = [
                    "{}-{}".format(combi_group_name, i) for i in range(no_targets)
                ]
                batch_tags = [combi_group.at[0, "batch-tag"]] * no_targets
                concentrations = [
                    float(combi_group.at[0, "concentration-required-mM"])
                ] * no_targets
                amounts = [float(combi_group.at[0, "amount-required-uL"])] * no_targets
                no_steps = [combi_group.at[0, "no-steps"]] * no_targets
                self.target_names = self.target_names + target_names
                self.batchtags = self.batchtags + batch_tags
                self.concentrations = self.concentrations + concentrations
                self.amounts = self.amounts + amounts
                self.nosteps = self.nosteps + no_steps
                for reaction_number in reaction_numbers_group:
                    reaction_combi_group_info = {}
                    if reaction_number == 1:
                        reactant_1_SMILES = [
                            reactant.strip()
                            for reactant in combi_group[
                                "reactant-1-{}".format(reaction_number)
                            ]
                            if str(reactant) != "nan"
                        ]

                        reactant_2_SMILES = [
                            reactant.strip()
                            for reactant in combi_group[
                                "reactant-2-{}".format(reaction_number)
                            ]
                            if str(reactant) != "nan"
                        ]
                        are_product_SMILES = False

                    if reaction_number > 1:
                        reactant_1_SMILES = [
                            reactant.strip()
                            for reactant in combi_group[
                                "reactant-1-{}".format(reaction_number)
                            ]
                            if str(reactant) != "nan"
                        ]
                        reactant_2_SMILES = product_smiles[:number_reactant_pair_smiles]
                        are_product_SMILES = True

                    reactant_pair_smiles = combiChem(
                        reactant_1_SMILES=reactant_1_SMILES,
                        reactant_2_SMILES=reactant_2_SMILES,
                        are_product_SMILES=are_product_SMILES,
                    )

                    number_reactant_pair_smiles = len(reactant_pair_smiles)
                    if number_reactant_pair_smiles != no_targets:
                        reactant_pair_smiles = reactant_pair_smiles * (
                            no_targets // len(reactant_pair_smiles)
                        )
                    reaction_names = [
                        combi_group.at[0, "reaction-name-{}".format(reaction_number)]
                    ] * no_targets
                    reaction_recipes = [
                        combi_group.at[0, "reaction-recipe-{}".format(reaction_number)]
                    ] * no_targets
                    reactant_pair_smiles_ordered, product_smiles = self.checkReaction(
                        reactant_pair_smiles=reactant_pair_smiles,
                        reaction_names=reaction_names,
                        reaction_recipes=reaction_recipes,
                    )
                    if reaction_number == max_no_steps_combi_group:
                        self.target_smiles = self.target_smiles + product_smiles
                    reaction_combi_group_info[
                        "reaction-name-{}".format(reaction_number)
                    ] = reaction_names
                    reaction_combi_group_info[
                        "reaction-recipe-{}".format(reaction_number)
                    ] = reaction_recipes
                    reaction_combi_group_info[
                        "reaction-reactant-pair-smiles-{}".format(reaction_number)
                    ] = reactant_pair_smiles_ordered
                    reaction_combi_group_info[
                        "reaction-product-smiles-{}".format(reaction_number)
                    ] = product_smiles
                    combi_group_info.update(reaction_combi_group_info)

                products = list(
                    zip(
                        *[
                            combi_group_info[
                                "reaction-product-smiles-{}".format(reactionnumber)
                            ]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                reactant_pair_smiles = list(
                    zip(
                        *[
                            combi_group_info[
                                "reaction-reactant-pair-smiles-{}".format(
                                    reactionnumber
                                )
                            ]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                reaction_names = list(
                    zip(
                        *[
                            combi_group_info["reaction-name-{}".format(reactionnumber)]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                reaction_recipes = list(
                    zip(
                        *[
                            combi_group_info[
                                "reaction-recipe-{}".format(reactionnumber)
                            ]
                            for reactionnumber in reaction_numbers_group
                        ]
                    )
                )
                self.product_smiles = self.product_smiles + products
                self.reactant_pair_smiles = (
                    self.reactant_pair_smiles + reactant_pair_smiles
                )
                self.reaction_names = self.reaction_names + reaction_names
                self.reaction_recipes = self.reaction_recipes + reaction_recipes

            if self.validated:
                self.df = pd.DataFrame()
                self.df["batch-tag"] = self.batchtags
                self.df["target-names"] = self.target_names
                self.df["target-SMILES"] = self.target_smiles
                self.df["concentration-required-mM"] = self.concentrations
                self.df["amount-required-uL"] = self.amounts
                self.df["no-steps"] = self.nosteps
                self.df["reactant-pair-smiles"] = self.reactant_pair_smiles
                self.df["reaction-name"] = self.reaction_names
                self.df["reaction-recipe"] = self.reaction_recipes
                self.df["product-smiles"] = self.product_smiles
                self.checkIsNumber(values=self.concentrations)
                self.checkIsNumber(values=self.amounts)

    def add_warning(self, field, warning_string):
        self.validate_dict["field"].append(field)
        self.validate_dict["warning_string"].append(warning_string)

    def checkColumnNames(self, expected_column_names):
        try:
            if not set(self.df_columns) == set(expected_column_names):
                self.add_warning(
                    field="name_columns",
                    warning_string="Column names should be set to: {}".format(
                        expected_column_names
                    ),
                )
                self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="name_columns",
                warning_string="Column names check failed with error: {}".format(e),
            )
            self.validated = False

    def checkNumberColumns(self, expected_no_columns, expected_column_names):
        try:
            if self.no_df_columns != expected_no_columns:
                self.add_warning(
                    field="number_columns",
                    warning_string="Found {} columns. Expected {} columns. Set and name columns to {} only".format(
                        self.no_df_columns,
                        expected_no_columns,
                        expected_column_names,
                    ),
                )
                self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="number_columns",
                warning_string="Number columns check failed with error: {}".format(e),
            )
            self.validated = False

    def checkSMILES(
        self, df_rows_index: list[int], smiles: list[str], smiles_type: str
    ):
        """Checks if input SMILES from df is valid

        Parameters
        ----------
        df_rows_index: list[int]
            The index of the df rows being tested - for error reporting
        smiles: list[str]
            The SMILES being tested eg. Target or reactant pair SMILES
        smiles_type: str
            The type of SMILES being tested eg. target or reactant_pair

        """
        try:
            for index, smi in zip(df_rows_index, smiles):
                if all(isinstance(item, tuple) for item in smi):
                    mol_test = [Chem.MolFromSmiles(smi) for smi in smi]
                else:
                    mol_test = [Chem.MolFromSmiles(smi)]
                if None in mol_test:
                    self.add_warning(
                        field="check_smiles",
                        warning_string="Input {} smiles: '{}' at index {} is not a valid smiles".format(
                            smiles_type,
                            smi,
                            index,
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_smiles",
                warning_string="Input {} smiles check failed for smiles: {} and with error: {}".format(
                    smiles_type, smiles, e
                ),
            )
            self.validated = False

    def checkIsNumber(self, values):
        try:
            # self.target_amounts = [amount for amount in self.df["volume-required-uL"]]
            for index, value in zip(self.index_df_rows, values):
                if not isinstance(value, (int, float)):
                    self.add_warning(
                        field="check_number",
                        warning_string="Value {} at index {} is not a valid number".format(
                            value, index
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_number",
                warning_string="Input amount check failed with error: {}".format(e),
            )
            self.validated = False

    def checkIsString(self):
        try:
            self.batchtags = [tag.strip() for tag in self.df["batch-tag"]]
            for index, tag in zip(self.index_df_rows, self.batchtags):
                if not type(tag) == str:
                    self.add_warning(
                        field="check_string",
                        warning_string="Batch tag {} at index {} is not a valid string format".format(
                            tag, index
                        ),
                    )
                    self.validated = False
        except Exception as e:
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_string",
                warning_string="Input batch tag check failed with error: {}".format(e),
            )
            self.validated = False

    def checkReaction(
        self,
        reactant_pair_smiles: list,
        reaction_names: list[str],
        reaction_recipes: list[str],
        product_smiles: list[str] = None,
    ):
        try:
            product_created_smiles = []
            reactant_pair_smiles_ordered = []
            for index, (reactant_pair, reaction_name, reaction_recipe) in enumerate(
                zip(reactant_pair_smiles, reaction_names, reaction_recipes)
            ):
                smarts = encoded_recipes[reaction_name]["recipes"][reaction_recipe][
                    "reactionSMARTS"
                ]
                if not all(smarts):
                    print(
                        "Warning ignoring smarts pattern for reaction: ", reaction_name
                    )
                    print("This can only be used for custom chem uploads")
                    product_created_smiles.append(product_smiles[index])
                    reactant_pair_smiles_ordered.append(reactant_pair)
                    continue

                product_mols = checkReactantSMARTS(
                    reactant_SMILES=reactant_pair, reaction_SMARTS=smarts
                )
                if not product_mols:
                    self.add_warning(
                        field="check_reaction",
                        warning_string="Reaction for reactants: {} and reaction: {} is not a valid reaction".format(
                            reactant_pair, reaction_name
                        ),
                    )
                    self.validated = False

                if product_mols:
                    if product_smiles:
                        product_mol = Chem.MolFromSmiles(product_smiles[index])
                    else:
                        product_mol = product_mols[-1]
                    product_smi = Chem.MolToSmiles(product_mol)
                    reactant_smis = getAddtionOrder(
                        product_smi=product_smi,
                        reactant_SMILES=reactant_pair,
                        reaction_SMARTS=smarts,
                    )
                    product_created_smiles.append(product_smi)
                    reactant_pair_smiles_ordered.append(reactant_smis)
            return reactant_pair_smiles_ordered, product_created_smiles

        except Exception as e:
            logger.info(
                inspect.stack()[0][3]
                + "Reaction check failed with error: {} for reaction: {}, recipe: {}, reactant pair: {} ".format(
                    e, reaction_name, reaction_recipe, reactant_pair
                ),
            )
            logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
            self.add_warning(
                field="check_reaction",
                warning_string="Reaction check failed with error: {} for reaction: {}, recipe: {}, reactant pair: {} ".format(
                    e, reaction_name, reaction_recipe, reactant_pair
                ),
            )
            self.validated = False
