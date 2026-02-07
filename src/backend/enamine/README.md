# Enamine REAL Tools API Integration for CAR

This module provides integration with the Enamine REAL Tools API to search for synthetic routes and building blocks for target molecules in the CAR system.

## Features

- **Single and batch SMILES searches** in the REAL space
- **Human-readable reaction type mappings** for easy selection
- **Building block extraction** from search results  
- **Synthon ID correlation** to get SMILES from synthon IDs
- **Any stereo search by default** as requested
- **Rate limiting** to respect API quotas
- **Comprehensive error handling** for invalid SMILES

## Dependencies

- `requests`: HTTP client for API calls
- `ratelimit`: Rate limiting decorator
- `typing`: Type hints for better code quality

## Notes

- **Default search mode**: Any stereo (as requested)
- **Default reaction types**: Very tractable (0) + Moderately tractable (1050)
- **Building blocks**: Extracted from `rgn` field in API responses
- **Reaction codes**: Extracted from `rsn` field for synthesis tracking

## Setup

1. Obtain an API key from [Enamine REAL Tools](https://real.enamine.net/api/docs)
2. Set the environment variable:
   ```bash
   export ENAMINE_REAL_TOOLS_API_KEY="your_api_key_here"
   ```

## Reaction Types

The API supports four reaction types with human-readable mappings:

| Code | Description | Usage in CAR |
|------|-------------|---------------|
| 0 | Very tractable | Most straightforward reactions, typically cheapest |
| 1050 | Moderately tractable | More complex reactions, typically more expensive |
| 1051 | Legacy | Legacy reactions, not in current REAL version |
| 1055 | Custom synthesis | Candidates for custom synthesis by Enamine |

**Default for CAR**: Types 0 and 1050 (very and moderately tractable)

## Basic Usage

### Single Molecule Search

```python
from enamine.apicalls import search_smiles_single, extract_building_blocks

# Search for a target molecule  
target_smiles = "CNC(=O)C1Cc2ccccc2N1C(=O)c1cccc(OCC(F)(F)F)n1"
results = search_smiles_single(target_smiles)

print(f"Found {len(results)} structures")
if results:
    building_blocks = extract_building_blocks(results[0])
    print(f"Building blocks needed: {len(building_blocks)}")
```

### Batch Search

```python
from enamine.apicalls import search_smiles_batch

# Search multiple molecules at once
smiles_list = [
    "CNC(=O)C1Cc2ccccc2N1C(=O)c1cccc(OCC(F)(F)F)n1",
    "O=C(Cn1cc(-c2ccccc2)nn1)N1CCC(N2CCCC2)C1"
]

results = search_smiles_batch(smiles_list)
structures_found = len([r for r in results if 'error' not in r])
print(f"Total structures found: {structures_found}")
```

### Extract Building Blocks

```python
from enamine.apicalls import extract_building_blocks, extract_reaction_codes

# Get building blocks from a search result
for structure in results:
    if 'error' not in structure:  # Skip error results
        building_blocks = extract_building_blocks(structure)
        reaction_codes = extract_reaction_codes(structure)
        
        print(f"Structure: {structure['smiles']}")
        print("Building blocks:")
        for bb in building_blocks:
            print(f"  - {bb['code']}: {bb['smiles']}")
        print(f"Reaction codes: {reaction_codes}")
```

### Synthon ID Correlation

```python
from enamine.apicalls import get_synthons_by_id

# Get SMILES for specific synthon IDs
synthon_ids = [23844188, 18715594]
synthons = get_synthons_by_id(synthon_ids)

for synthon in synthons:
    print(f"Synthon {synthon['id']}: {synthon['sSmiles']}")
    print(f"  Role in reaction: {synthon['roleInReaction']}")
```

## API Response Structure

Each search result contains:

- `smiles`: Found structure in ChemAxon format
- `sntDt`: Estimated availability date
- `vSynt`: List of virtual syntheses
  - `rgn`: Reagents/building blocks with code and SMILES
  - `rsn`: RSN codes encoding the molecule
- `query_smiles`: Original query SMILES
- `space`: Category estimation (e.g., "sREAL")
- `exact`: (any-stereo only) Whether structure is exact match

## Error Handling

The API handles various error conditions:

- **Invalid SMILES**: Returns error field with details
- **Charged molecules**: Returns "No charged structures in REAL"
- **Rate limiting**: Automatic retry with exponential backoff
- **Network errors**: Raises `requests.RequestException`

## Integration with CAR Workflow

1. **Target Selection**: Use `search_smiles_single()` to find synthesis routes
2. **Route Analysis**: Extract building blocks with `extract_building_blocks()`
3. **Cost Estimation**: Use reaction type info for cost planning
4. **Batch Processing**: Use `search_smiles_batch()` for multiple targets
5. **Synthesis Planning**: Use RSN codes and building blocks for protocol design

## Rate Limiting

The integration includes built-in rate limiting:
- Single searches: 100 calls/minute
- Batch searches: 50 calls/minute  
- Synthon queries: 100 calls/minute

