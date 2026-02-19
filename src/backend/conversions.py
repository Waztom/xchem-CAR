"""Conversion and calculation utilities.

This module contains simple conversion functions for wells, molar calculations,
mass calculations, and string sanitization.
"""

from rdkit.Chem import Descriptors
from rdkit import Chem

import re
import inspect
import logging

logger = logging.getLogger(__name__)


def well_index_to_well_name(wellindex: int, platesize: int) -> str:
    """Converts a well index to a human readable well name

    Parameters
    ----------
    wellindex: int
        The well index to convert to a well name
    platesize: int
        The plate size to convert the well index to a well name

    Returns
    -------
    wellname: str
        The well name eg. A01
    """
    platewellindexconversions = {
        4: "A" + "%02d" % ((wellindex) + 1,),
        24: "ABCD"[(wellindex) % 4] + "%02d" % ((wellindex) // 4 + 1,),
        96: "ABCDEFGH"[(wellindex) % 8] + "%02d" % ((wellindex) // 8 + 1,),
        384: "ABCDEFGHIJKLMNOP"[(wellindex) % 16] + "%02d" % ((wellindex) // 16 + 1,),
    }
    wellname = platewellindexconversions[platesize]
    return wellname


def calculate_mols_from_conc(target_concentration: float, target_volume: float) -> object:
    """Function to calculate product mols of reaction using a target mass

    Parameters
    ----------
    target_concentration: float
        The target concentration (mM) of the product
    target_volume: float
        The target volume (uL) of the product

    Returns
    -------
    product_moles: rdkit mol object
        The product mols
    """
    try:
        target_mols = (target_volume / 1e6) * (target_concentration / 1e3)
        return target_mols
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def calculate_mass_from_mols(mols: float, SMILES: str) -> object:
    """Function to calculate mass from mols

    Parameters
    ----------
    mols: float
        The mols of the compound
    SMILES: str
        The SMILES of the compound

    Returns
    -------
    mass: float
        The mass (mg) of the compound
    """
    try:
        MW = Descriptors.MolWt(Chem.MolFromSmiles(SMILES))
        mass = (mols * MW) * 1e3
        return mass
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def sanitize_for_python_var(name):
    """
    Convert a string to a valid Python variable name.

    Parameters
    ----------
    name: str
        The string to be converted to a valid Python variable name

    Returns
    -------
    str
        A string that conforms to Python variable naming rules
    """
    # Replace any non-alphanumeric characters (except underscore) with underscore
    name = re.sub(r"[^\w]", "_", name)

    # Ensure it starts with a letter or underscore
    if name and not (name[0].isalpha() or name[0] == "_"):
        name = "plate_" + name

    return name
