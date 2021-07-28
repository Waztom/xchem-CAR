"""Create IBMAPI object for calling IBM API"""
from __future__ import annotations
import functools
from rxn4chemistry import RXN4ChemistryWrapper
import os
import time
from rdkit.Chem import AllChem
import requests
import json

import sys
from .commonsolvents import common_solvents

sys.path.append("..")


from ..utils import canonSmiles


class IBMAPI(object):
    """
    Creates an IBMAPI object for calling IBM API
    """

    def __init__(self, project_name: str):
        """
        ValidateFile constructor
        Args:
            project_name (str): Project name used in project model
        """
        self.project_name = project_name
        self.api_key = os.environ["IBM_API_KEY"]
        self.createIBMProject()

    def createIBMProject(self):
        rxn4chemistry_wrapper = RXN4ChemistryWrapper(api_key=self.api_key)
        rxn4chemistry_wrapper.create_project(self.project_name)
        self.rxn4chemistry_wrapper = rxn4chemistry_wrapper

    def convertIBMNameToSmiles(self, chemical_name):
        try:
            data = [chemical_name]
            headers = {
                "Authorization": self.api_key,
                "Content-Type": "application/json",
                "Accept": "application/json",
            }
            url = "https://rxn.res.ibm.com/rxn/api/api/v1/actions/convert-material-to-smiles"
            r = requests.post(url=url, data=json.dumps(data), headers=headers, cookies={})
            response_dict = r.json()
            smiles = response_dict["payload"][chemical_name]
            return smiles
        except:
            return False

    def getIBMRetroSyn(self, smiles, max_steps):
        """
        Use the IBM API to get some possible retrosynthesis routes
        """
        # Create dummy dictionary to create while loop to catch when status is a SUCCESS
        results = {}
        results["status"] = None
        results["results"] = None

        while results["results"] is None:
            try:
                time.sleep(10)
                response = self.rxn4chemistry_wrapper.predict_automatic_retrosynthesis(
                    product=smiles, max_steps=max_steps
                )
                while results["status"] != "SUCCESS":
                    time.sleep(130)
                    results = (
                        self.rxn4chemistry_wrapper.get_predict_automatic_retrosynthesis_results(
                            response["prediction_id"]
                        )
                    )
                    results["results"] = results
            except Exception as e:
                return None
        return results["results"]

    def collect_actions(self, tree):
        actions = []

        if "children" in tree and len(tree["children"]):
            actions.append(tree["actions"])

        for node in tree["children"]:
            actions.extend(self.collect_actions(node))

        return actions

    def collect_rclass(self, tree):
        rclass = []
        if "children" in tree and len(tree["children"]):
            rclass.append(tree["rclass"])

        for node in tree["children"]:
            rclass.extend(self.collect_rclass(node))
        return rclass

    def collect_products(self, tree):
        products = []
        if "children" in tree and len(tree["children"]):
            products.append(canonSmiles(tree["smiles"]))

        for node in tree["children"]:
            products.extend(self.collect_products(node))
        return products

    def collect_reactants(self, tree):
        reactants = []
        if "children" in tree and len(tree["children"]):
            reactants.append([canonSmiles(node["smiles"]) for node in tree["children"]])

        for node in tree["children"]:
            reactants.extend(self.collect_reactants(node))
        return reactants

    def collect_reactions(self, tree):
        reactions = []
        if "children" in tree and len(tree["children"]):
            reactions.append(
                AllChem.ReactionFromSmarts(
                    "{}>>{}".format(
                        ".".join(
                            [
                                canonSmiles(node["smiles"])
                                for node in tree["children"]
                                if canonSmiles(node["smiles"]) not in common_solvents
                            ]
                        ),
                        tree["smiles"],
                    ),
                    useSmiles=True,
                )
            )

        for node in tree["children"]:
            reactions.extend(self.collect_reactions(node))
        return reactions

    def collectIBMReactionInfo(self, pathway):
        reaction_info = {}

        try:
            time.sleep(10)
            pathway_synthesis_response = self.rxn4chemistry_wrapper.create_synthesis_from_sequence(
                sequence_id=pathway["sequenceId"]
            )
            pathway_synthesis_id = pathway_synthesis_response["synthesis_id"]
            time.sleep(10)
            synthesis_tree, reactions, actions = self.rxn4chemistry_wrapper.get_synthesis_plan(
                synthesis_id=pathway_synthesis_id
            )

            reaction_info["actions"] = self.collect_actions(synthesis_tree)
            reaction_info["rclass"] = self.collect_rclass(pathway)
            reaction_info["product_smiles"] = self.collect_products(pathway)
            reaction_info["reactants"] = self.collect_reactants(pathway)
            reaction_info["reactions"] = self.collect_reactions(pathway)

            return reaction_info
        except:
            return None

    def convert(self, bytes_list):
        res = functools.reduce(lambda total, d: 10 * total + d, bytes_list, 0)
        return res

    def convertbytes(self, reaction_class):
        list_bytes = list(reaction_class.encode("utf8"))
        return list_bytes

    def filtermethod(self, reaction_info):
        reaction_classes = reaction_info["rclass"]
        reactants = reaction_info["reactants"]
        flat_list = [item for sublist in reactants for item in sublist]
        reactants = [reactant for reactant in flat_list if reactant not in common_solvents]
        reactants = list(dict.fromkeys(reactants))
        rxn_classes_integers = [
            self.convert(self.convertbytes(rxn_class)) for rxn_class in reaction_classes
        ]
        reactant_integers = [self.convert(self.convertbytes(reactant)) for reactant in reactants]
        method_integer = sum(reactant_integers) + sum(rxn_classes_integers)
        return method_integer
