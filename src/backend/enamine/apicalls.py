import requests
from ratelimit import limits, sleep_and_retry
import os
from typing import List, Dict, Any, Optional
from enum import Enum

ENAMINE_API_KEY = os.environ["ENAMINE_REAL_TOOLS_API_KEY"]
ENAMINE_BASE_URL = "https://real.enamine.net"

# Human-readable reaction type mappings
class ReactionType(Enum):
    VERY_TRACTABLE = 0          # Most straightforward, often 2-component reactions, typically cheapest
    MODERATELY_TRACTABLE = 1050 # More complex reactions, products typically more expensive
    LEGACY = 1051               # Legacy reactions, not part of current REAL version
    CUSTOM_SYNTHESIS = 1055     # Candidates that didn't pass, custom synthesis available

REACTION_TYPE_DESCRIPTIONS = {
    0: "Very tractable",
    1050: "Moderately tractable", 
    1051: "Legacy",
    1055: "Custom synthesis"
}

# Default headers for API requests
def _get_headers() -> Dict[str, str]:
    """Get standard headers for REAL Tools API requests."""
    return {
        "Content-Type": "application/json",
        "X-API-KEY": ENAMINE_API_KEY
    }

@sleep_and_retry
@limits(calls=100, period=60)  # Rate limiting
def search_smiles_single(
    smiles: str, 
    reaction_types: Optional[List[int]] = None,
    any_stereo: bool = True
) -> List[Dict[str, Any]]:
    """
    Search for a single SMILES in REAL space.
    
    Args:
        smiles: Target molecule SMILES string
        reaction_types: List of reaction type integers (default: [0, 1050] for very/moderately tractable)
        any_stereo: If True, search for any stereoisomer; if False, search exact stereo match
        
    Returns:
        List of found structures with synthesis information
        
    Raises:
        requests.RequestException: If API request fails
    """
    if reaction_types is None:
        reaction_types = [ReactionType.VERY_TRACTABLE.value, ReactionType.MODERATELY_TRACTABLE.value]
    
    endpoint = "any-stereo" if any_stereo else "exact-stereo"
    url = f"{ENAMINE_BASE_URL}/api/v1/space/real/search-structure/single/{endpoint}"
    
    payload = {
        "smiles": smiles,
        "reaction_types": reaction_types
    }
    
    response = requests.post(url, json=payload, headers=_get_headers())
    response.raise_for_status()
    
    return response.json()

@sleep_and_retry
@limits(calls=50, period=60)  # Lower rate limit for batch requests
def search_smiles_batch(
    smiles_list: List[str],
    reaction_types: Optional[List[int]] = None, 
    any_stereo: bool = True
) -> List[Dict[str, Any]]:
    """
    Search for multiple SMILES in REAL space.
    
    Args:
        smiles_list: List of SMILES strings to search
        reaction_types: List of reaction type integers (default: [0, 1050] for very/moderately tractable)
        any_stereo: If True, search for any stereoisomer; if False, search exact stereo match
        
    Returns:
        List of found structures with synthesis information
        
    Raises:
        requests.RequestException: If API request fails
    """
    if reaction_types is None:
        reaction_types = [ReactionType.VERY_TRACTABLE.value, ReactionType.MODERATELY_TRACTABLE.value]
        
    endpoint = "any-stereo" if any_stereo else "exact-stereo"
    url = f"{ENAMINE_BASE_URL}/api/v1/space/real/search-structure/batch/{endpoint}"
    
    payload = {
        "smiles_list": smiles_list,
        "reaction_types": reaction_types
    }
    
    response = requests.post(url, json=payload, headers=_get_headers())
    response.raise_for_status()
    
    return response.json()

@sleep_and_retry
@limits(calls=100, period=60)
def get_synthons_by_id(synthon_ids: List[int]) -> List[Dict[str, Any]]:
    """
    Get synthon SMILES by their numeric IDs.
    
    Args:
        synthon_ids: List of numeric synthon IDs
        
    Returns:
        List of dictionaries containing:
        - id: synthon ID
        - roleInReaction: role number in the reaction
        - sSmiles: synthon SMILES string
        
    Raises:
        requests.RequestException: If API request fails
    """
    url = f"{ENAMINE_BASE_URL}/api/v1/space/real/utils/get-synthons-by-id"
    
    response = requests.post(url, json=synthon_ids, headers=_get_headers())
    response.raise_for_status()
    
    return response.json()


def extract_building_blocks(search_result: Dict[str, Any]) -> List[Dict[str, str]]:
    """
    Extract building block information from a REAL search result.
    
    Args:
        search_result: Single result from REAL search
        
    Returns:
        List of building blocks with code and SMILES
    """
    building_blocks = []
    
    if 'vSynt' in search_result:
        for synthesis in search_result['vSynt']:
            if 'rgn' in synthesis:
                for reagent in synthesis['rgn']:
                    building_blocks.append({
                        'code': reagent.get('code', ''),
                        'smiles': reagent.get('smiles', '')
                    })
    
    return building_blocks

def extract_reaction_codes(search_result: Dict[str, Any]) -> List[str]:
    """
    Extract RSN reaction codes from a REAL search result.
    
    Args:
        search_result: Single result from REAL search
        
    Returns:
        List of RSN codes
    """
    reaction_codes = []
    
    if 'vSynt' in search_result:
        for synthesis in search_result['vSynt']:
            if 'rsn' in synthesis:
                reaction_codes.extend(synthesis['rsn'])
                
    return reaction_codes

# Example usage functions (for documentation/testing)
def example_search_single():
    """
    Example of how to search for a single target molecule.
    """
    # Example target molecule SMILES
    target_smiles = "CNC(=O)C1Cc2ccccc2N1C(=O)c1cccc(OCC(F)(F)F)n1"
    
    print("Searching for target molecule...")
    print(f"SMILES: {target_smiles}")
    
    # Search with default settings (very tractable + moderately tractable, any stereo)
    results = search_smiles_single(target_smiles)
    
    print(f"Found {len(results)} structures")
    
    # Extract building blocks for first result if available
    if results:
        first_result = results[0]
        building_blocks = extract_building_blocks(first_result)
        reaction_codes = extract_reaction_codes(first_result)
        print(f"Building blocks needed: {len(building_blocks)}")
        for bb in building_blocks:
            print(f"  - {bb['code']}: {bb['smiles']}")
        print(f"Reaction codes: {reaction_codes[:2]}...")  # Show first 2
    
    return results

def example_search_batch():
    """
    Example of how to search for multiple target molecules.
    """
    # Example target molecules
    target_smiles_list = [
        "CNC(=O)C1Cc2ccccc2N1C(=O)c1cccc(OCC(F)(F)F)n1",
        "O=C(Cn1cc(-c2ccccc2)nn1)N1CCC(N2CCCC2)C1"
    ]
    
    print("Searching for multiple target molecules...")
    
    # Search batch
    results = search_smiles_batch(target_smiles_list)
    
    structures_found = len([r for r in results if 'error' not in r])
    errors = len([r for r in results if 'error' in r])
    
    print(f"Total structures searched: {len(results)}")
    print(f"Structures found: {structures_found}")
    print(f"Errors: {errors}")
    
    return results

