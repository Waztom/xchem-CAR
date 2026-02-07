#!/usr/bin/env python3
"""
Enamine REAL Tools Target Processing Prototype

This prototype script processes a CSV file of target molecules using the Enamine REAL Tools API
to find synthesis routes, building blocks, and reaction information.

Usage:
    python process_targets_prototype.py input_file.csv [options]

Requirements:
    - Input CSV with 'target-SMILES' column
    - ENAMINE_REAL_TOOLS_API_KEY environment variable
    - enamine-real-tools-rxns.csv file for reaction mapping
"""

import pandas as pd
import argparse
import sys
import os
import re
from typing import List, Dict, Any, Optional
from pathlib import Path

# Import our Enamine API functions
try:
    from apicalls import search_smiles_batch, extract_building_blocks, extract_reaction_codes, ReactionType, get_synthons_by_id
except ImportError:
    print("Error: Could not import apicalls module. Make sure apicalls.py is in the same directory.")
    sys.exit(1)


def parse_rsn_code(rsn_code: str) -> Dict[str, Any]:
    """
    Parse RSN code to extract reaction ID and synthon IDs.
    
    Example: "s487____14925160____28909612" -> 
    {
        'reaction_code': 's487',
        'reaction_id': 487,
        'synthon_ids': [14925160, 28909612]
    }
    
    Args:
        rsn_code: RSN code string from REAL Tools API
        
    Returns:
        Dictionary with parsed information
    """
    try:
        parts = rsn_code.split('____')
        if len(parts) < 2:
            return {'reaction_code': rsn_code, 'reaction_id': None, 'synthon_ids': []}
        
        reaction_code = parts[0]
        # Strip letters to get numeric reaction ID
        reaction_id_match = re.search(r'\d+', reaction_code)
        reaction_id = int(reaction_id_match.group()) if reaction_id_match else None
        
        # Extract synthon IDs (all parts after reaction code)
        synthon_ids = []
        for part in parts[1:]:
            try:
                synthon_ids.append(int(part))
            except ValueError:
                continue
                
        return {
            'reaction_code': reaction_code,
            'reaction_id': reaction_id,
            'synthon_ids': synthon_ids
        }
    except Exception as e:
        print(f"Warning: Failed to parse RSN code '{rsn_code}': {e}")
        return {'reaction_code': rsn_code, 'reaction_id': None, 'synthon_ids': []}


def fetch_synthon_smarts_batch(all_synthon_ids: List[int]) -> Dict[int, str]:
    """
    Fetch synthon SMARTS for a list of synthon IDs with efficient batching.
    
    Args:
        all_synthon_ids: List of unique synthon IDs to fetch
        
    Returns:
        Dictionary mapping synthon ID to SMARTS string
    """
    synthon_smarts_map = {}
    batch_size = 4  # API maximum limit for synthon queries
    
    print(f"Fetching SMARTS for {len(all_synthon_ids)} unique synthon IDs...")
    
    # Process in batches to avoid overwhelming the API
    for i in range(0, len(all_synthon_ids), batch_size):
        batch = all_synthon_ids[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(all_synthon_ids) + batch_size - 1) // batch_size
        
        print(f"Fetching synthon batch {batch_num}/{total_batches} ({len(batch)} synthons)")
        print(f"  Batch IDs: {batch[:5]}...")  # Show first 5 IDs
        
        try:
            # Ensure all IDs are proper integers
            clean_batch = [int(synthon_id) for synthon_id in batch if isinstance(synthon_id, (int, float)) and synthon_id > 0]
            
            if not clean_batch:
                print(f"  Warning: No valid synthon IDs in batch {batch_num}")
                continue
                
            synthon_data = get_synthons_by_id(clean_batch)
            print(f"  API returned {len(synthon_data)} synthon records")
            
            # Map synthon ID to SMARTS
            for synthon in synthon_data:
                synthon_id = synthon.get('id')
                smarts = synthon.get('sSmiles', '')  # API returns 'sSmiles' field
                if synthon_id is not None:
                    synthon_smarts_map[synthon_id] = smarts
                    
        except Exception as e:
            print(f"Warning: Failed to fetch synthon batch {batch_num}: {e}")
            print(f"  Batch content: {batch}")
            # Add empty SMARTS for failed batch
            for synthon_id in batch:
                if isinstance(synthon_id, (int, float)) and synthon_id > 0:
                    synthon_smarts_map[int(synthon_id)] = ''
    
    return synthon_smarts_map


def collect_all_synthon_ids(processed_results: List[Dict[str, Any]]) -> List[int]:
    """
    Collect all unique synthon IDs from processed results.
    
    Args:
        processed_results: List of processed result dictionaries
        
    Returns:
        List of unique synthon IDs
    """
    all_synthon_ids = set()
    
    for result in processed_results:
        if result.get('error'):
            continue
            
        # Collect from the synthon_ids string (semicolon-separated)
        synthon_ids_str = result.get('synthon_ids', '')
        if synthon_ids_str and isinstance(synthon_ids_str, str):
            # Parse semicolon-separated synthon IDs
            for synthon_id_str in synthon_ids_str.split(';'):
                if synthon_id_str.strip():
                    try:
                        # Try int first, fall back to float->int for legacy CSVs with .0 values
                        try:
                            synthon_id = int(synthon_id_str.strip())
                        except ValueError:
                            synthon_id = int(float(synthon_id_str.strip()))
                            
                        if synthon_id > 0:  # Only add positive valid IDs
                            all_synthon_ids.add(synthon_id)
                    except (ValueError, TypeError):
                        print(f"Warning: Invalid synthon ID '{synthon_id_str}' - skipping")
                        continue
        elif isinstance(synthon_ids_str, list):
            # Handle case where it's already a list
            for synthon_id in synthon_ids_str:
                try:
                    synthon_id_int = int(synthon_id) if synthon_id else 0
                    if synthon_id_int > 0:
                        all_synthon_ids.add(synthon_id_int)
                except (ValueError, TypeError):
                    continue
    
    return sorted(list(all_synthon_ids))


def load_reaction_mapping(reaction_csv_path: str) -> pd.DataFrame:
    """
    Load the reaction mapping CSV file.
    
    Args:
        reaction_csv_path: Path to enamine-real-tools-rxns.csv
        
    Returns:
        DataFrame with reaction information
    """
    try:
        reaction_df = pd.read_csv(reaction_csv_path)
        print(f"Loaded reaction mapping with {len(reaction_df)} reactions")
        return reaction_df
    except FileNotFoundError:
        print(f"Warning: Reaction mapping file '{reaction_csv_path}' not found.")
        print("Will proceed without reaction mapping.")
        return pd.DataFrame()
    except Exception as e:
        print(f"Error loading reaction mapping: {e}")
        return pd.DataFrame()


def process_search_results(results: List[Dict[str, Any]], reaction_df: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Process REAL Tools search results to extract relevant information.
    
    Args:
        results: List of search results from REAL Tools API
        reaction_df: DataFrame with reaction mapping information
        
    Returns:
        List of processed result dictionaries
    """
    processed_results = []
    
    for result in results:
        if 'error' in result:
            processed_results.append({
                'target_smiles': result.get('query_smiles', ''),
                'error': result['error'],
                'building_blocks': [],
                'reaction_info': [],
                'synthon_ids': []
            })
            continue
            
        target_smiles = result.get('query_smiles', result.get('smiles', ''))
        building_blocks = extract_building_blocks(result)
        reaction_codes = extract_reaction_codes(result)
        
        # Process reaction codes and synthon IDs
        reaction_info = []
        all_synthon_ids = []
        
        for rsn_code in reaction_codes:
            parsed_rsn = parse_rsn_code(rsn_code)
            
            # Look up reaction information if mapping available
            reaction_details = {'rsn_code': rsn_code}
            if not reaction_df.empty and parsed_rsn['reaction_id'] is not None:
                matching_reactions = reaction_df[reaction_df['reaction-ID'] == parsed_rsn['reaction_id']]
                if not matching_reactions.empty:
                    reaction_details.update(matching_reactions.iloc[0].to_dict())
            
            reaction_details.update(parsed_rsn)
            reaction_info.append(reaction_details)
            all_synthon_ids.extend(parsed_rsn['synthon_ids'])
        
        processed_results.append({
            'target_smiles': target_smiles,
            'found_smiles': result.get('smiles', ''),
            'space': result.get('space', ''),
            'synthesis_date': result.get('sntDt', ''),
            'exact_match': result.get('exact', None),
            'building_blocks': building_blocks,
            'reaction_info': reaction_info,
            'synthon_ids': all_synthon_ids,
            'num_building_blocks': len(building_blocks),
            'vSynt': result.get('vSynt', []),  # Keep original vSynt for detailed processing
            'error': None
        })
    
    return processed_results


def create_output_dataframe(processed_results: List[Dict[str, Any]], reaction_df: pd.DataFrame, synthon_smarts_map: Optional[Dict[int, str]] = None) -> pd.DataFrame:
    """
    Create output DataFrame from processed results.
    Each reaction gets its own row with all building blocks for that reaction.
    
    Args:
        processed_results: List of processed result dictionaries
        reaction_df: DataFrame with reaction mapping information
        
    Returns:
        DataFrame ready for CSV export
    """
    output_rows = []
    
    for result in processed_results:
        if result['error']:
            # Add error row
            output_rows.append({
                'target_smiles': result['target_smiles'],
                'status': 'ERROR',
                'error': result['error'],
                'found_smiles': '',
                'space': '',
                'synthesis_date': '',
                'exact_match': '',
                'reaction_code': '',
                'reaction_id': '',
                'reaction_name': '',
                'reaction_components': '',
                'reaction_space': '',
                'reaction_description': '',
                'reaction_tags': '',
                'synthon_ids': '',
                'num_building_blocks': 0,
                'reactant_1_code': '',
                'reactant_1_smiles': '',
                'synthon_id_1': '',
                'reactant_2_code': '',
                'reactant_2_smiles': '',
                'synthon_id_2': '',
                'reactant_3_code': '',
                'reactant_3_smiles': '',
                'synthon_id_3': '',
                'reactant_4_code': '',
                'reactant_4_smiles': '',
                'synthon_id_4': '',
            })
            continue
            
        if not result['building_blocks']:
            # No building blocks found
            output_rows.append({
                'target_smiles': result['target_smiles'],
                'status': 'NOT_FOUND',
                'error': '',
                'found_smiles': result['found_smiles'],
                'space': result['space'],
                'synthesis_date': result['synthesis_date'],
                'exact_match': result['exact_match'],
                'method_no': '',
                'reaction_code': '',
                'reaction_id': '',
                'reaction_name': '',
                'reaction_components': '',
                'reaction_space': '',
                'reaction_description': '',
                'reaction_tags': '',
                'synthon_ids': '',
                'num_building_blocks': 0,
                'reactant_1_code': '',
                'reactant_1_smiles': '',
                'synthon_id_1': '',
                'reactant_2_code': '',
                'reactant_2_smiles': '',
                'synthon_id_2': '',
                'reactant_3_code': '',
                'reactant_3_smiles': '',
                'synthon_id_3': '',
                'reactant_4_code': '',
                'reactant_4_smiles': '',
                'synthon_id_4': '',
            })
            continue
        
        # Group building blocks by reaction (using vSynt index)
        if 'vSynt' in result and len(result['vSynt']) > 0:
            method_counter = 1
            # Process each synthesis route separately
            for synt_idx, synthesis in enumerate(result['vSynt']):
                # Get building blocks for this synthesis route
                route_building_blocks = []
                if 'rgn' in synthesis:
                    for reagent in synthesis['rgn']:
                        route_building_blocks.append({
                            'code': reagent.get('code', ''),
                            'smiles': reagent.get('smiles', '')
                        })
                
                # Check if too many building blocks
                if len(route_building_blocks) > 4:
                    output_rows.append({
                        'target_smiles': result['target_smiles'],
                        'status': '>FOUR_BBs',
                        'error': f'Found {len(route_building_blocks)} building blocks (max 4)',
                        'found_smiles': result['found_smiles'],
                        'space': result['space'],
                        'synthesis_date': result['synthesis_date'],
                        'exact_match': result['exact_match'],
                        'method_no': method_counter,
                        'reaction_code': '',
                        'reaction_id': '',
                        'reaction_name': '',
                        'reaction_components': '',
                        'reaction_space': '',
                        'reaction_description': '',
                        'reaction_tags': '',
                        'synthon_ids': '',
                        'num_building_blocks': len(route_building_blocks),
                        'reactant_1_code': '',
                        'reactant_1_smiles': '',
                        'reactant_synthon_id_1': 0,
                        'reactant_synthon_smarts_1': '',
                        'reactant_2_code': '',
                        'reactant_2_smiles': '',
                        'reactant_synthon_id_2': 0,
                        'reactant_synthon_smarts_2': '',
                        'reactant_3_code': '',
                        'reactant_3_smiles': '',
                        'reactant_synthon_id_3': 0,
                        'reactant_synthon_smarts_3': '',
                        'reactant_4_code': '',
                        'reactant_4_smiles': '',
                        'reactant_synthon_id_4': 0,
                        'reactant_synthon_smarts_4': '',
                    })
                    method_counter += 1
                    continue
                
                # Get reaction codes for this synthesis route
                route_reaction_codes = synthesis.get('rsn', [])
                
                # Process each reaction code separately (each gets its own row)
                for reaction_code in route_reaction_codes:
                    route_synthon_ids = []
                    reaction_id = ''
                    reaction_name = ''
                    reaction_components = ''
                    reaction_space = ''
                    reaction_description = ''
                    reaction_tags = ''
                    
                    if reaction_code:
                        parsed_rsn = parse_rsn_code(reaction_code)
                        reaction_id = parsed_rsn.get('reaction_id', '')
                        route_synthon_ids = parsed_rsn.get('synthon_ids', [])
                        
                        # Look up reaction details from mapping DataFrame
                        if not reaction_df.empty and reaction_id:
                            matching_reactions = reaction_df[reaction_df['reaction-ID'] == reaction_id]
                            if not matching_reactions.empty:
                                reaction_details = matching_reactions.iloc[0]
                                reaction_name = reaction_details.get('reaction-name', '')
                                reaction_components = reaction_details.get('no-components', '')
                                reaction_space = reaction_details.get('reaction-space', '')
                                reaction_description = reaction_details.get('description', '')
                                reaction_tags = reaction_details.get('tags', '')
                    
                    # Create row for this reaction route
                    row = {
                        'target_smiles': result['target_smiles'],
                        'status': 'FOUND',
                        'error': '',
                        'found_smiles': result['found_smiles'],
                        'space': result['space'],
                        'synthesis_date': result['synthesis_date'],
                        'exact_match': result['exact_match'],
                        'method_no': method_counter,
                        'reaction_code': reaction_code,
                        'reaction_id': reaction_id,
                        'reaction_name': reaction_name,
                        'reaction_components': reaction_components,
                        'reaction_space': reaction_space,
                        'reaction_description': reaction_description,
                        'reaction_tags': reaction_tags,
                        'synthon_ids': ';'.join(map(str, route_synthon_ids)),
                        'num_building_blocks': len(route_building_blocks),
                    }
                    
                    # Add reactants and synthon IDs in pairs (up to 4 per reaction)
                    for i in range(4):
                        if i < len(route_building_blocks):
                            row[f'reactant_{i+1}_code'] = route_building_blocks[i]['code']
                            row[f'reactant_{i+1}_smiles'] = route_building_blocks[i]['smiles']
                        else:
                            row[f'reactant_{i+1}_code'] = ''
                            row[f'reactant_{i+1}_smiles'] = ''
                            
                        if i < len(route_synthon_ids):
                            synthon_id = int(route_synthon_ids[i])  # Ensure integer
                            row[f'reactant_synthon_id_{i+1}'] = synthon_id
                            # Add SMARTS if available
                            if synthon_smarts_map and synthon_id in synthon_smarts_map:
                                row[f'reactant_synthon_smarts_{i+1}'] = synthon_smarts_map[synthon_id]
                            else:
                                row[f'reactant_synthon_smarts_{i+1}'] = ''
                        else:
                            row[f'reactant_synthon_id_{i+1}'] = 0  # Use 0 instead of empty string
                            row[f'reactant_synthon_smarts_{i+1}'] = ''
                    
                    output_rows.append(row)
                    method_counter += 1
        else:
            # Fallback: create single row with all building blocks (if <= 4)
            if len(result['building_blocks']) > 4:
                # Too many building blocks
                output_rows.append({
                    'target_smiles': result['target_smiles'],
                    'status': '>FOUR_BBs',
                    'error': f'Found {len(result["building_blocks"])} building blocks (max 4)',
                    'found_smiles': result['found_smiles'],
                    'space': result['space'],
                    'synthesis_date': result['synthesis_date'],
                    'exact_match': result['exact_match'],
                    'method_no': 1,
                    'reaction_code': '',
                    'reaction_id': '',
                    'reaction_name': '',
                    'reaction_components': '',
                    'reaction_space': '',
                    'reaction_description': '',
                    'reaction_tags': '',
                    'synthon_ids': ';'.join(map(str, result['synthon_ids'])),
                    'num_building_blocks': len(result['building_blocks']),
                    'reactant_1_code': '',
                    'reactant_1_smiles': '',
                    'reactant_synthon_id_1': 0,
                    'reactant_synthon_smarts_1': '',
                    'reactant_2_code': '',
                    'reactant_2_smiles': '',
                    'reactant_synthon_id_2': 0,
                    'reactant_synthon_smarts_2': '',
                    'reactant_3_code': '',
                    'reactant_3_smiles': '',
                    'reactant_synthon_id_3': 0,
                    'reactant_synthon_smarts_3': '',
                    'reactant_4_code': '',
                    'reactant_4_smiles': '',
                    'reactant_synthon_id_4': 0,
                    'reactant_synthon_smarts_4': '',
                })
            else:
                row = {
                    'target_smiles': result['target_smiles'],
                    'status': 'FOUND',
                    'error': '',
                    'found_smiles': result['found_smiles'],
                    'space': result['space'],
                    'synthesis_date': result['synthesis_date'],
                    'exact_match': result['exact_match'],
                    'method_no': 1,
                    'reaction_code': '',
                    'reaction_id': '',
                    'reaction_name': '',
                    'reaction_components': '',
                    'reaction_space': '',
                    'reaction_description': '',
                    'reaction_tags': '',
                    'synthon_ids': ';'.join(map(str, result['synthon_ids'])),
                    'num_building_blocks': len(result['building_blocks']),
                }
                
                # Add reactants and synthon IDs in pairs (up to 4)
                for i in range(4):
                    if i < len(result['building_blocks']):
                        row[f'reactant_{i+1}_code'] = result['building_blocks'][i]['code']
                        row[f'reactant_{i+1}_smiles'] = result['building_blocks'][i]['smiles']
                    else:
                        row[f'reactant_{i+1}_code'] = ''
                        row[f'reactant_{i+1}_smiles'] = ''
                        
                    if i < len(result['synthon_ids']):
                        synthon_id = int(result['synthon_ids'][i])  # Ensure integer
                        row[f'reactant_synthon_id_{i+1}'] = synthon_id
                        # Add SMARTS if available
                        if synthon_smarts_map and synthon_id in synthon_smarts_map:
                            row[f'reactant_synthon_smarts_{i+1}'] = synthon_smarts_map[synthon_id]
                        else:
                            row[f'reactant_synthon_smarts_{i+1}'] = ''
                    else:
                        row[f'reactant_synthon_id_{i+1}'] = 0  # Use 0 instead of empty string
                        row[f'reactant_synthon_smarts_{i+1}'] = ''
                
                output_rows.append(row)
    
    return pd.DataFrame(output_rows)


def process_targets_batch(
    smiles_list: List[str], 
    reaction_df: pd.DataFrame,
    reaction_types: Optional[List[int]] = None
) -> List[Dict[str, Any]]:
    """
    Process a batch of target SMILES using REAL Tools API.
    
    Args:
        smiles_list: List of target SMILES to process
        reaction_df: DataFrame with reaction mapping
        reaction_types: List of reaction types to search (default: [0, 1050])
        
    Returns:
        List of processed results
    """
    if reaction_types is None:
        reaction_types = [ReactionType.VERY_TRACTABLE.value, ReactionType.MODERATELY_TRACTABLE.value]
    
    print(f"Processing batch of {len(smiles_list)} targets...")
    print(f"Reaction types: {reaction_types}")
    
    try:
        # Search using REAL Tools API
        search_results = search_smiles_batch(smiles_list, reaction_types=reaction_types, any_stereo=True)
        
        # Process results
        processed_results = process_search_results(search_results, reaction_df)
        
        return processed_results
        
    except Exception as e:
        print(f"Error processing batch: {e}")
        # Return error results for all SMILES in batch
        return [{
            'target_smiles': smiles,
            'error': f"Batch processing failed: {str(e)}",
            'building_blocks': [],
            'reaction_info': [],
            'synthon_ids': []
        } for smiles in smiles_list]


def main():
    """Main function to process target molecules."""
    parser = argparse.ArgumentParser(description='Process target molecules using Enamine REAL Tools API')
    parser.add_argument('input_file', help='Input CSV file with target-SMILES column')
    parser.add_argument('-o', '--output', help='Output CSV file (default: input_file_processed.csv)')
    parser.add_argument('-b', '--batch-size', type=int, default=1000, 
                       help='Batch size for API requests (default: 1000)')
    parser.add_argument('-r', '--reactions-file', default='enamine-real-tools-rxns.csv',
                       help='Reaction mapping CSV file (default: enamine-real-tools-rxns.csv)')
    parser.add_argument('--reaction-types', nargs='+', type=int, 
                       default=[0, 1050], help='Reaction types to search (default: 0 1050)')
    parser.add_argument('--include-legacy', action='store_true',
                       help='Include legacy reactions (type 1051)')
    parser.add_argument('--include-custom', action='store_true', 
                       help='Include custom synthesis (type 1055)')
    
    args = parser.parse_args()
    
    # Check API key
    if not os.environ.get('ENAMINE_REAL_TOOLS_API_KEY'):
        print("Error: ENAMINE_REAL_TOOLS_API_KEY environment variable not set")
        sys.exit(1)
    
    # Prepare reaction types
    reaction_types = list(args.reaction_types)
    if args.include_legacy:
        reaction_types.append(ReactionType.LEGACY.value)
    if args.include_custom:
        reaction_types.append(ReactionType.CUSTOM_SYNTHESIS.value)
    
    print(f"Using reaction types: {reaction_types}")
    
    # Load input file
    try:
        input_df = pd.read_csv(args.input_file)
        print(f"Loaded {len(input_df)} targets from {args.input_file}")
    except Exception as e:
        print(f"Error loading input file: {e}")
        sys.exit(1)
    
    # Check for required column
    if 'target-SMILES' not in input_df.columns:
        print("Error: Input CSV must contain 'target-SMILES' column")
        print(f"Available columns: {list(input_df.columns)}")
        sys.exit(1)
    
    # Load reaction mapping
    reaction_df = load_reaction_mapping(args.reactions_file)
    
    # Get target SMILES list
    target_smiles = input_df['target-SMILES'].dropna().tolist()
    print(f"Processing {len(target_smiles)} target molecules...")
    
    # Process in batches
    all_results = []
    batch_size = args.batch_size
    
    for i in range(0, len(target_smiles), batch_size):
        batch = target_smiles[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(target_smiles) + batch_size - 1) // batch_size
        
        print(f"\nProcessing batch {batch_num}/{total_batches} ({len(batch)} targets)")
        
        batch_results = process_targets_batch(batch, reaction_df, reaction_types)
        all_results.extend(batch_results)
        
        print(f"Completed batch {batch_num}/{total_batches}")
    
    # Collect all unique synthon IDs and fetch their SMARTS
    print("\nCollecting unique synthon IDs...")
    all_synthon_ids = collect_all_synthon_ids(all_results)
    print(f"Found {len(all_synthon_ids)} unique synthon IDs: {all_synthon_ids[:10]}...")  # Show first 10
    
    synthon_smarts_map = {}
    if all_synthon_ids:
        synthon_smarts_map = fetch_synthon_smarts_batch(all_synthon_ids)
        print(f"Fetched SMARTS for {len(synthon_smarts_map)} synthons")
        print(f"Sample SMARTS mapping: {dict(list(synthon_smarts_map.items())[:3])}")  # Show first 3
    else:
        print("No synthon IDs found to fetch SMARTS")

    # Create output DataFrame
    print("\nCreating output DataFrame...")
    output_df = create_output_dataframe(all_results, reaction_df, synthon_smarts_map)
    # Set output filename
    if args.output:
        output_file = args.output
    else:
        input_path = Path(args.input_file)
        output_file = input_path.parent / f"{input_path.stem}_processed.csv"
    
    # Save results
    output_df.to_csv(output_file, index=False)
    
    # Print summary
    total_targets = len(all_results)
    found_targets = len([r for r in all_results if not r.get('error')])
    error_targets = total_targets - found_targets
    
    print(f"\n{'='*50}")
    print("PROCESSING SUMMARY")
    print(f"{'='*50}")
    print(f"Total targets processed: {total_targets}")
    print(f"Targets with synthesis routes found: {found_targets}")
    print(f"Targets with errors: {error_targets}")
    print(f"Success rate: {found_targets/total_targets*100:.1f}%")
    print(f"\nOutput saved to: {output_file}")
    print(f"Total output rows: {len(output_df)}")


if __name__ == "__main__":
    main()