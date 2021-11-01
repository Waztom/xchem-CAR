"""Checks validation of file for uploading to CAR"""
from __future__ import annotations
from itertools import product
import pandas as pd
from rdkit import Chem

from .recipebuilder.encodedrecipes import encoded_recipes
from .utils import canonSmiles, getAddtionOrder, checkReactantSMARTS, combichem


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
        self.upload_type = validate_type
        self.validate_dict = {"field": [], "warning_string": []}
        self.validated = True

        if self.upload_type == "custom-chem":
            self.validatecustomchem()
        if self.upload_type == "combi-custom-chem":
            self.validatecustomcombichem()

        if self.upload_type == "retro-API":
            self.expected_no_columns = 2
            self.expected_column_names = ["targets", "amount-required-mg"]
            self.checkNumberColumns()
            if self.validated:
                self.checkColumnNames()
            if self.validated:
                self.target_smiles = [canonSmiles(smi.strip()) for smi in self.df["targets"]]
                self.df["targets"] = self.target_smiles
                self.checkTargetSMILES()
                if self.validated:
                    self.checkIsNumber()

    def validatecustomchem(self):
        self.expected_no_columns = 4
        self.expected_column_names = [
            "reactant-1",
            "reactant-2",
            "reaction-name",
            "amount-required-mg",
        ]
        self.checkNumberColumns()
        if self.validated:
            self.checkColumnNames()
        if self.validated:
            self.reactant_pair_smiles = [
                reactants for reactants in zip(self.df["reactant-1"], self.df["reactant-2"])
            ]
            self.df["reactant-pair-smiles"] = self.reactant_pair_smiles
            self.checkReactantSMILES()
            if self.validated:
                self.reactant_pair_smiles = [
                    (canonSmiles(smi[0]), canonSmiles(smi[1])) for smi in self.reactant_pair_smiles
                ]
                self.reaction_names = self.df["reaction-name"]
                self.checkReaction()
                if self.validated:
                    self.df["reactant-pair-smiles"] = self.reactant_pair_smiles_ordered
                    self.df["target-smiles"] = self.product_smiles
                    self.checkIsNumber()

    def validatecustomcombichem(self):
        self.expected_no_columns = 4
        self.expected_column_names = [
            "reactant-1",
            "reactant-2",
            "reaction-name",
            "amount-required-mg",
        ]
        self.checkNumberColumns()
        if self.validated:
            self.checkColumnNames()
        if self.validated:
            self.reactant_pair_smiles = []
            self.reaction_names = []
            grouped = self.df.groupby("reaction-name")
            for name, group in grouped:
                reactant_1_SMILES = set(
                    [reactant for reactant in group["reactant-1"] if str(reactant) != "nan"]
                )
                reactant_2_SMILES = set(
                    [reactant for reactant in group["reactant-2"] if str(reactant) != "nan"]
                )
                reactant_pair_smiles = combichem(
                    reactant_1_SMILES=reactant_1_SMILES, reactant_2_SMILES=reactant_2_SMILES
                )
                reaction_names = [name] * len(reactant_pair_smiles)
                self.reactant_pair_smiles = self.reactant_pair_smiles + reactant_pair_smiles
                self.reaction_names = self.reaction_names + reaction_names

            self.checkReactantSMILES()
            if self.validated:
                self.reactant_pair_smiles = [
                    (canonSmiles(smi[0]), canonSmiles(smi[1])) for smi in self.reactant_pair_smiles
                ]
                self.checkReaction()
                if self.validated:
                    amount_required_mg = self.df.at[0, "amount-required-mg"]
                    self.df = pd.DataFrame()
                    self.df["reactant-pair-smiles"] = self.reactant_pair_smiles_ordered
                    self.df["target-smiles"] = self.product_smiles
                    self.df["reaction-name"] = self.reaction_names
                    self.df["amount-required-mg"] = [amount_required_mg] * len(
                        self.reactant_pair_smiles
                    )
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

    def checkReactantSMILES(self):
        for index, smi_pair in zip(self.index_df_rows, self.reactant_pair_smiles):
            mols = [Chem.MolFromSmiles(smi) for smi in smi_pair]
            if None in mols:
                none_test_indices = [index for index, mol in enumerate(mols) if mol == None]
                invalid_smiles = [smi_pair[index] for index in none_test_indices]
                self.add_warning(
                    field="check_smiles",
                    warning_string="Input reactant smiles: ".join(
                        "{} ".format(*smi) for smi in invalid_smiles
                    )
                    + "at index {} is not a valid smiles".format(
                        index,
                    ),
                )
                self.validated = False

    def checkIsNumber(self):
        self.target_amounts = [amount for amount in self.df["amount-required-mg"]]
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
        self.reactant_pair_smiles_ordered = []
        no_reaction_tests = len(self.reaction_names)

        for index, reactant_pair, reaction_name in zip(
            range(no_reaction_tests), self.reactant_pair_smiles, self.reaction_names
        ):
            smarts = encoded_recipes[reaction_name]["reactionSMARTS"]
            product_mols = checkReactantSMARTS(
                reactant_SMILES=reactant_pair, reaction_SMARTS=smarts
            )

            if not product_mols:
                print(smarts)
                print(reactant_pair)
                self.add_warning(
                    field="check_reaction",
                    warning_string="Reaction for reactants at index {} is not a valid reaction".format(
                        index
                    ),
                )
                self.validated = False

            if product_mols:
                product_mol = product_mols[
                    0
                ]  # Need to build in something to show muttiple products and then let user choose product!
                product_smi = Chem.MolToSmiles(product_mol)
                reactant_smis = getAddtionOrder(
                    product_smi=product_smi,
                    reactant_SMILES=reactant_pair,
                    reaction_SMARTS=smarts,
                )
                self.product_smiles.append(product_smi)
                self.reactant_pair_smiles_ordered.append(reactant_smis)
