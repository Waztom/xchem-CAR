import requests
from ratelimit import limits, sleep_and_retry
import os

api_key = os.environ["MANIFOLD_API_KEY"]


@sleep_and_retry
@limits(calls=1000, period=60)
def getManifoldRetrosynthesis(smiles: str):
    """Call Manifold API to search for a retrosynthesis for a given smiles

    Parameters
    ----------
    smiles: str
        SMILES for Manifold retrosynthesis search

    Returns
    -------
    response.json(): json
        The route matches
    """

    data = {
        "smiles": smiles,
        "maxLeadTimeWeeks": 12,
        "maxSearchDepth": 3,
        "maxNumRoutesToReturn": 3,
        "catalogs": [
            "enamine_bb",
            "molport",
            "mcule",
            "mcule_ultimate",
        ],
    }

    response = requests.post(
        url="https://api.postera.ai/api/v1/retrosynthesis/",
        headers={
            "X-API-KEY": api_key,
        },
        json=data,
    )
    return response.json()


@sleep_and_retry
@limits(calls=1000, period=60)
def getManifoldRetrosynthesisBatch(smiles: list):
    """Call Manifold API to search for a retrosynthesis for a given list
    of target compound SMILES

    Parameters
    ----------
    smiles: list
        The list of SMILES for the Manifold retrosynthesis search

    Returns
    -------
    response.json(): json
        The route matches to the list of SMILES
    """

    data = {
        "smilesList": smiles,
        "maxLeadTimeWeeks": 12,
        "maxSearchDepth": 3,
        "maxNumRoutesToReturn": 3,
        "catalogs": [
            "enamine_bb",
            "molport",
            "mcule",
            "mcule_ultimate",
        ],
    }

    response = requests.post(
        url="https://api.postera.ai/api/v1/retrosynthesis/batch/",
        headers={
            "X-API-KEY": api_key,
        },
        json=data,
    )
    return response.json()


@sleep_and_retry
@limits(calls=1000, period=60)
def getExactSearch(smiles: str):
    """Searches for exact compound match for catalogue info

    Parameters
    ----------
    smiles: str
        The SMILES of the compound to search the manifold API
        for catalogue info

    Returns
    -------
    response.json(): json
        The catalogue matches
    """
    data = {
        "smiles": smiles,
        "maxLeadTimeWeeks": 12,
        "patentDatabases": [],
    }

    response = requests.post(
        "https://api.postera.ai/api/v1/exact/",
        headers={
            "X-API-KEY": api_key,
        },
        json=data,
    )

    return response.json()


@sleep_and_retry
@limits(calls=1000, period=60)
def getExactSearchBatch(smilesList: list):
    """Searches for exact compound match for catalogue info

    Parameters
    ----------
    smiles: lits
        The list of SMILES of compounds to search the manifold API
        for catalogue info

    Returns
    -------
    response.json(): json
        The catalogue matches for the list of SMILES
    """
    data = {
        "smilesList": smilesList,
        "maxLeadTimeWeeks": 12,
        "vendors": ["enamine_bb", "enamine_bb_EU-US", "enamine_real"],
        "patentDatabases": [],
    }

    response = requests.post(
        "https://api.postera.ai/api/v1/exact/batch/",
        headers={
            "X-API-KEY": api_key,
        },
        json=data,
    )

    return response.json()
