"""PubChem API utilities for compound lookups.

This module contains functions for querying the PubChem database
to retrieve compound information, CAS numbers, and chemical names.
"""

import pubchempy as pcp
import re
import inspect
import logging

logger = logging.getLogger(__name__)


def get_pubchem_compound(inchikey: str) -> object:
    """Searches PubChem for compound using inchi key

    Parameters
    ----------
    inchikey: str
        The inchikey of the compound to search the PubChem DB for

    Returns
    -------
    compound: object
        The PuBChem compound class object
    status: None
        Returns None if no compound is found or an error occurs
    """
    try:
        compound = pcp.get_compounds(inchikey, "inchikey")[0]
        if not compound.cid:
            return None
        else:
            return compound
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(
            "Pubchempy could not retrieve compound entry for input inchikey: {} with error {}".format(
                inchikey, e
            )
        )


def get_pubchem_cas(compound: object) -> str:
    """Get CAS identifier for PubChem compound synonyms

    Parameters
    ----------
    compound: object
        A PuBChem compound object

    Returns
    -------
    cas: str
        The CAS id of the compound
    """
    try:
        synonyms = compound.synonyms
        if synonyms:
            for syn in synonyms:
                match = re.match("(\d{1,7}-\d{1,2}-\d)", syn)
                if match:
                    cas = match.group(1)
                    return cas
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(e)


def get_chemical_name(inchikey: str) -> str:
    """Searches PubChem for compound using SMILES

    Parameters
    ----------
    inchikey: str
        The inchiley of the compound to search the PubChem DB for
        it's IUPAC name

    Returns
    -------
    name: str
        The IUPAC name of the compound
    status: None
        Returns None if no compound IUPAC name is found or if an error
        occurs
    """
    try:
        name = pcp.get_compounds(inchikey, "inchikey")[0].iupac_name
        if not name:
            return None
        else:
            return name
    except Exception as e:
        logger.info(inspect.stack()[0][3] + " yielded error: {}".format(e))
        print(
            "Pubchempy could not convert SMILES to a IUPAC name with error {}".format(e)
        )
