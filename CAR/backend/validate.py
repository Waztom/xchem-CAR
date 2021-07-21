"""Checks validation of file for uploading to CAR"""
from __future__ import annotations
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdChemReactions

from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import canonSmiles


class ValidateFile(object):
    """
    Creates a validate object for checking file validation for upload
    """

    def __init__(self, csv_to_validate: csvFile, validate_type: str):
        """
        ValidateFile constructor
        Args:
            csv_to_validate (.csv): Uploaded .csv file for testing validation
        """
        self.df = pd.read_csv(csv_to_validate, encoding="utf8")
        self.df_columns = self.df.columns
        self.no_df_columns = len(self.df_columns)
        self.index_df_rows = range(0, len(self.df), 1)
        self.target_amounts = [amount for amount in self.df["Ammount_required (mg)"]]
        self.upload_type = validate_type
        self.validate_dict = {"field": [], "warning_string": []}
        self.validated = True

        if self.upload_type == "custom-chem":
            self.expected_no_columns = 4
            self.expected_column_names = [
                "Reactant-1",
                "Reactant-2",
                "Reaction-name",
                "Ammount_required (mg)",
            ]
            self.checkNumberColumns()
            if self.validated:
                self.checkColumnNames()
            if self.validated:
                self.reactant_pair_smiles = [
                    reactants for reactants in zip(self.df["Reactant-1"], self.df["Reactant-2"])
                ]
                self.df["reactant_pair_smiles"] = self.reactant_pair_smiles
                self.reaction_names = [reaction_name for reaction_name in self.df["Reaction-name"]]
                self.checkReaction()
                if self.validated:
                    self.df["target-smiles"] = self.product_smiles
                    self.checkIsNumber()

        if self.upload_type == "retro-API":
            self.expected_no_columns = 2
            self.expected_column_names = ["Targets", "Ammount_required (mg)"]
            self.checkNumberColumns()
            if self.validated:
                self.checkColumnNames()
            if self.validated:
                self.target_smiles = [canonSmiles(smi.strip()) for smi in self.df["Targets"]]
                self.df["Targets"] = self.target_smiles
                self.checkTargetSMILES()
                if self.validated:
                    self.checkIsNumber()

    def add_warning(self, field, warning_string):
        self.validate_dict["field"].append(field)
        self.validate_dict["warning_string"].append(warning_string)

    def checkColumnNames(self):
        if not all(self.df_columns == self.expected_column_names):
            self.add_warning(
                field="name_columns",
                warning_string="Column names should be set to: {}".format(
                    self.expected_column_names
                ),
            )
            self.validated = False

    def checkNumberColumns(self):
        if self.no_df_columns != self.expected_no_columns:
            self.add_warning(
                field="number_columns",
                warning_string="Found {} columns. Expected {} columns. Set and name columns to {} only".format(
                    self.no_df_columns,
                    self.expected_no_columns,
                    self.expected_column_names,
                ),
            )
            self.validated = False

    def checkTargetSMILES(self):
        for index, smi in zip(self.index_df_rows, self.target_smiles):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                self.add_warning(
                    field="check_smiles",
                    warning_string="Input target smiles: '{}' at index {} is not a valid smiles".format(
                        smi,
                        index,
                    ),
                )
                self.validated = False

    def checkIsNumber(self):
        for index, amount in zip(self.index_df_rows, self.target_amounts):
            if not isinstance(amount, (int, float)):
                self.add_warning(
                    field="check_number",
                    warning_string="Target mass {} at index {} is not a valid number".format(
                        amount, index
                    ),
                )
                self.validated = False

    def checkReaction(self):
        self.product_smiles = []

        for index, reactant_pair, reaction_name in zip(
            self.index_df_rows, self.reactant_pair_smiles, self.reaction_names
        ):
            reaction_smarts = encoded_recipes[reaction_name]["reactionSMARTS"]
            reacts = [Chem.MolFromSmiles(smi) for smi in reactant_pair]

            for smarts in reaction_smarts:
                reaction = rdChemReactions.ReactionFromSmarts(smarts)
                products = reaction.RunReactants(reacts)
                if len(products) != 0:
                    print(Chem.MolToSmiles(products[0][0]))

            if len(products) == 0:
                self.add_warning(
                    field="check_reaction",
                    warning_string="Reaction for reactants at index {} is not a valid reaction".format(
                        index
                    ),
                )
                self.validated = False

            if len(products) != 0:
                product_mol = products[0][0]
                self.product_smiles.append(Chem.MolToSmiles(product_mol))
