import json
import requests
import os

api_key = os.environ["MANIFOLD_API_KEY"]


def getManifoldretrosynthesis(target_smiles):
    """
    Function to call the Manifold API to search for a retrosynthesis for a given smiles
    """
    data = {
        "smiles": target_smiles,
        "catalogs": [
            "mcule",
            "mcule_ultimate",
        ],
        "maxLeadTimeWeeks": 12,
        "maxSearchDepth": 3,
        "maxNumRoutesToReturn": 30,
        "reactionTag": "diamond_robotic_synthesis",
    }

    response = requests.post(
        url="https://api.postera.ai/api/v1/retrosynthesis/",
        headers={
            "X-API-KEY": api_key,
        },
        json=data,
    )

    return response.json()


def getexactsearch(target_smiles):
    response = requests.post(
        "https://api.postera.ai/api/v1/exact/",
        headers={
            "X-API-KEY": api_key,
        },
        json={"smiles": target_smiles},
    )

    return response.json()