import functools
from .common_solvents import common_solvents


def convert(bytes_list):
    # multiply each integer element with its
    # corresponding power and perform summation
    res = functools.reduce(lambda total, d: 10 * total + d, bytes_list, 0)

    return res


def convertbytes(reaction_class):
    list_bytes = list(reaction_class.encode("utf8"))

    return list_bytes


def filtermethod(reaction_info):
    reaction_classes = reaction_info["rclass"]
    reactants = reaction_info["reactants"]

    # Flatten reactants list and exclude solvent reactants
    flat_list = [item for sublist in reactants for item in sublist]
    reactants = [reactant for reactant in flat_list if reactant not in common_solvents]

    # IBM API sometimes yields duplicate reactants
    reactants = list(dict.fromkeys(reactants))

    # Try encode to bytes to convert strings to numbers
    rxn_classes_integers = [convert(convertbytes(rxn_class)) for rxn_class in reaction_classes]
    reactant_integers = [convert(convertbytes(reactant)) for reactant in reactants]

    method_integer = sum(reactant_integers) + sum(rxn_classes_integers)

    return method_integer
