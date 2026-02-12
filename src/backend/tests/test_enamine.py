"""
Unit tests for the Enamine REAL Tools integration.

Tests cover:
- API wrapper functions (with mocked responses)
- Helper functions for parsing and extracting data
- Data processing utilities
"""

import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
import os
import sys

# Add backend directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from backend.enamine.apicalls import (
    extract_building_blocks,
    extract_reaction_codes,
    ReactionType,
    REACTION_TYPE_DESCRIPTIONS,
)


class TestReactionTypes(unittest.TestCase):
    """Test reaction type enums and mappings."""
    
    def test_reaction_type_values(self):
        """Test that reaction type enum values are correct."""
        self.assertEqual(ReactionType.VERY_TRACTABLE.value, 0)
        self.assertEqual(ReactionType.MODERATELY_TRACTABLE.value, 1050)
        self.assertEqual(ReactionType.LEGACY.value, 1051)
        self.assertEqual(ReactionType.CUSTOM_SYNTHESIS.value, 1055)
    
    def test_reaction_type_descriptions(self):
        """Test that all reaction types have descriptions."""
        self.assertIn(0, REACTION_TYPE_DESCRIPTIONS)
        self.assertIn(1050, REACTION_TYPE_DESCRIPTIONS)
        self.assertIn(1051, REACTION_TYPE_DESCRIPTIONS)
        self.assertIn(1055, REACTION_TYPE_DESCRIPTIONS)
        
        self.assertEqual(REACTION_TYPE_DESCRIPTIONS[0], "Very tractable")
        self.assertEqual(REACTION_TYPE_DESCRIPTIONS[1050], "Moderately tractable")


class TestExtractBuildingBlocks(unittest.TestCase):
    """Test the extract_building_blocks function."""
    
    def test_extract_building_blocks_with_valid_data(self):
        """Test extracting building blocks from a valid search result."""
        search_result = {
            'smiles': 'CNC(=O)C1Cc2ccccc2N1C(=O)c1cccc(OCC(F)(F)F)n1',
            'vSynt': [
                {
                    'rgn': [
                        {'code': 'EN300-123456', 'smiles': 'CC(C)N'},
                        {'code': 'EN300-789012', 'smiles': 'c1ccccc1Br'}
                    ],
                    'rsn': ['s487____14925160____28909612']
                }
            ]
        }
        
        building_blocks = extract_building_blocks(search_result)
        
        self.assertEqual(len(building_blocks), 2)
        self.assertEqual(building_blocks[0]['code'], 'EN300-123456')
        self.assertEqual(building_blocks[0]['smiles'], 'CC(C)N')
        self.assertEqual(building_blocks[1]['code'], 'EN300-789012')
        self.assertEqual(building_blocks[1]['smiles'], 'c1ccccc1Br')
    
    def test_extract_building_blocks_empty_vsynt(self):
        """Test extracting building blocks when vSynt is empty."""
        search_result = {'smiles': 'CC', 'vSynt': []}
        
        building_blocks = extract_building_blocks(search_result)
        
        self.assertEqual(len(building_blocks), 0)
    
    def test_extract_building_blocks_no_vsynt(self):
        """Test extracting building blocks when vSynt is missing."""
        search_result = {'smiles': 'CC'}
        
        building_blocks = extract_building_blocks(search_result)
        
        self.assertEqual(len(building_blocks), 0)
    
    def test_extract_building_blocks_multiple_syntheses(self):
        """Test extracting building blocks from multiple synthesis routes."""
        search_result = {
            'smiles': 'CNC(=O)C1Cc2ccccc2N1',
            'vSynt': [
                {
                    'rgn': [{'code': 'EN1', 'smiles': 'CC'}],
                    'rsn': ['s100____111']
                },
                {
                    'rgn': [{'code': 'EN2', 'smiles': 'CCC'}, {'code': 'EN3', 'smiles': 'CCCC'}],
                    'rsn': ['s200____222']
                }
            ]
        }
        
        building_blocks = extract_building_blocks(search_result)
        
        # Should get all building blocks from all synthesis routes
        self.assertEqual(len(building_blocks), 3)
        codes = [bb['code'] for bb in building_blocks]
        self.assertIn('EN1', codes)
        self.assertIn('EN2', codes)
        self.assertIn('EN3', codes)


class TestExtractReactionCodes(unittest.TestCase):
    """Test the extract_reaction_codes function."""
    
    def test_extract_reaction_codes_single(self):
        """Test extracting a single reaction code."""
        search_result = {
            'smiles': 'CC',
            'vSynt': [
                {
                    'rgn': [],
                    'rsn': ['s487____14925160____28909612']
                }
            ]
        }
        
        codes = extract_reaction_codes(search_result)
        
        self.assertEqual(len(codes), 1)
        self.assertEqual(codes[0], 's487____14925160____28909612')
    
    def test_extract_reaction_codes_multiple(self):
        """Test extracting multiple reaction codes."""
        search_result = {
            'smiles': 'CC',
            'vSynt': [
                {'rgn': [], 'rsn': ['s100____111', 's200____222']},
                {'rgn': [], 'rsn': ['s300____333']}
            ]
        }
        
        codes = extract_reaction_codes(search_result)
        
        self.assertEqual(len(codes), 3)
        self.assertIn('s100____111', codes)
        self.assertIn('s200____222', codes)
        self.assertIn('s300____333', codes)
    
    def test_extract_reaction_codes_empty(self):
        """Test extracting reaction codes when none exist."""
        search_result = {'smiles': 'CC', 'vSynt': []}
        
        codes = extract_reaction_codes(search_result)
        
        self.assertEqual(len(codes), 0)


class TestParseRsnCode(unittest.TestCase):
    """Test the parse_rsn_code function from enamine_real_tools."""
    
    def setUp(self):
        """Import parse_rsn_code function."""
        from backend.enamine.enamine_real_tools import parse_rsn_code
        self.parse_rsn_code = parse_rsn_code
    
    def test_parse_rsn_code_standard(self):
        """Test parsing a standard RSN code."""
        rsn_code = "s487____14925160____28909612"
        
        result = self.parse_rsn_code(rsn_code)
        
        self.assertEqual(result['reaction_code'], 's487')
        self.assertEqual(result['reaction_id'], 487)
        self.assertEqual(result['synthon_ids'], [14925160, 28909612])
    
    def test_parse_rsn_code_single_synthon(self):
        """Test parsing RSN code with single synthon."""
        rsn_code = "s100____12345"
        
        result = self.parse_rsn_code(rsn_code)
        
        self.assertEqual(result['reaction_code'], 's100')
        self.assertEqual(result['reaction_id'], 100)
        self.assertEqual(result['synthon_ids'], [12345])
    
    def test_parse_rsn_code_three_synthons(self):
        """Test parsing RSN code with three synthons."""
        rsn_code = "s250____111____222____333"
        
        result = self.parse_rsn_code(rsn_code)
        
        self.assertEqual(result['reaction_code'], 's250')
        self.assertEqual(result['reaction_id'], 250)
        self.assertEqual(result['synthon_ids'], [111, 222, 333])
    
    def test_parse_rsn_code_invalid(self):
        """Test parsing an invalid RSN code."""
        rsn_code = "invalid"
        
        result = self.parse_rsn_code(rsn_code)
        
        self.assertEqual(result['reaction_code'], 'invalid')
        self.assertIsNone(result['reaction_id'])
        self.assertEqual(result['synthon_ids'], [])


class TestCollectSynthonIds(unittest.TestCase):
    """Test the collect_all_synthon_ids function."""
    
    def setUp(self):
        """Import collect_all_synthon_ids function."""
        from backend.enamine.enamine_real_tools import collect_all_synthon_ids
        self.collect_all_synthon_ids = collect_all_synthon_ids
    
    def test_collect_synthon_ids_from_string(self):
        """Test collecting synthon IDs from semicolon-separated string."""
        processed_results = [
            {'synthon_ids': '111;222;333', 'error': None}
        ]
        
        ids = self.collect_all_synthon_ids(processed_results)
        
        self.assertEqual(len(ids), 3)
        self.assertIn(111, ids)
        self.assertIn(222, ids)
        self.assertIn(333, ids)
    
    def test_collect_synthon_ids_multiple_results(self):
        """Test collecting synthon IDs from multiple results."""
        processed_results = [
            {'synthon_ids': '111;222', 'error': None},
            {'synthon_ids': '333;444', 'error': None}
        ]
        
        ids = self.collect_all_synthon_ids(processed_results)
        
        self.assertEqual(len(ids), 4)
    
    def test_collect_synthon_ids_with_duplicates(self):
        """Test that duplicate synthon IDs are deduplicated."""
        processed_results = [
            {'synthon_ids': '111;222', 'error': None},
            {'synthon_ids': '222;333', 'error': None}
        ]
        
        ids = self.collect_all_synthon_ids(processed_results)
        
        self.assertEqual(len(ids), 3)  # 111, 222, 333 (no duplicate 222)
    
    def test_collect_synthon_ids_skip_errors(self):
        """Test that error results are skipped."""
        processed_results = [
            {'synthon_ids': '111;222', 'error': None},
            {'synthon_ids': '999', 'error': 'Some error'}
        ]
        
        ids = self.collect_all_synthon_ids(processed_results)
        
        self.assertEqual(len(ids), 2)
        self.assertNotIn(999, ids)
    
    def test_collect_synthon_ids_empty_string(self):
        """Test with empty synthon_ids string."""
        processed_results = [
            {'synthon_ids': '', 'error': None}
        ]
        
        ids = self.collect_all_synthon_ids(processed_results)
        
        self.assertEqual(len(ids), 0)


class TestSearchSmilesAPI(unittest.TestCase):
    """Test API functions with mocked responses."""
    
    @patch('backend.enamine.apicalls.requests.post')
    @patch.dict(os.environ, {'ENAMINE_REAL_TOOLS_API_KEY': 'test_api_key'})
    def test_search_smiles_single_success(self, mock_post):
        """Test successful single SMILES search."""
        # Need to reimport after patching environment
        from backend.enamine.apicalls import search_smiles_single
        
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {
                'smiles': 'CNC(=O)C1Cc2ccccc2N1',
                'vSynt': [
                    {
                        'rgn': [{'code': 'EN300-123', 'smiles': 'CC'}],
                        'rsn': ['s100____111']
                    }
                ]
            }
        ]
        mock_response.raise_for_status = MagicMock()
        mock_post.return_value = mock_response
        
        results = search_smiles_single('CNC(=O)C1Cc2ccccc2N1')
        
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]['smiles'], 'CNC(=O)C1Cc2ccccc2N1')
        mock_post.assert_called_once()
    
    @patch('backend.enamine.apicalls.requests.post')
    @patch.dict(os.environ, {'ENAMINE_REAL_TOOLS_API_KEY': 'test_api_key'})
    def test_search_smiles_batch_success(self, mock_post):
        """Test successful batch SMILES search."""
        from backend.enamine.apicalls import search_smiles_batch
        
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {'smiles': 'CC', 'query_smiles': 'CC'},
            {'smiles': 'CCC', 'query_smiles': 'CCC'}
        ]
        mock_response.raise_for_status = MagicMock()
        mock_post.return_value = mock_response
        
        results = search_smiles_batch(['CC', 'CCC'])
        
        self.assertEqual(len(results), 2)
        mock_post.assert_called_once()
    
    @patch('backend.enamine.apicalls.requests.post')
    @patch.dict(os.environ, {'ENAMINE_REAL_TOOLS_API_KEY': 'test_api_key'})
    def test_get_synthons_by_id_success(self, mock_post):
        """Test successful synthon ID lookup."""
        from backend.enamine.apicalls import get_synthons_by_id
        
        mock_response = MagicMock()
        mock_response.json.return_value = [
            {'id': 12345, 'sSmiles': 'CC[N:1]', 'roleInReaction': 1},
            {'id': 67890, 'sSmiles': 'c1ccc([Br:1])cc1', 'roleInReaction': 2}
        ]
        mock_response.raise_for_status = MagicMock()
        mock_post.return_value = mock_response
        
        results = get_synthons_by_id([12345, 67890])
        
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0]['id'], 12345)
        self.assertEqual(results[0]['sSmiles'], 'CC[N:1]')


class TestLoadReactionMapping(unittest.TestCase):
    """Test the load_reaction_mapping function."""
    
    def setUp(self):
        """Import load_reaction_mapping function."""
        from backend.enamine.enamine_real_tools import load_reaction_mapping
        self.load_reaction_mapping = load_reaction_mapping
    
    def test_load_reaction_mapping_file_not_found(self):
        """Test loading when file doesn't exist."""
        result = self.load_reaction_mapping('/nonexistent/path/reactions.csv')
        
        self.assertTrue(result.empty)
    
    @patch('backend.enamine.enamine_real_tools.pd.read_csv')
    def test_load_reaction_mapping_success(self, mock_read_csv):
        """Test successful loading of reaction mapping."""
        mock_df = pd.DataFrame({
            'reaction-ID': [100, 200, 300],
            'reaction-name': ['Amide coupling', 'SNAr', 'Suzuki']
        })
        mock_read_csv.return_value = mock_df
        
        result = self.load_reaction_mapping('reactions.csv')
        
        self.assertEqual(len(result), 3)
        self.assertIn('reaction-ID', result.columns)


class TestProcessSearchResults(unittest.TestCase):
    """Test the process_search_results function."""
    
    def setUp(self):
        """Import process_search_results function."""
        from backend.enamine.enamine_real_tools import process_search_results
        self.process_search_results = process_search_results
        self.empty_reaction_df = pd.DataFrame()
    
    def test_process_results_with_error(self):
        """Test processing results that contain errors."""
        results = [
            {'query_smiles': 'invalid_smiles', 'error': 'Invalid SMILES'}
        ]
        
        processed = self.process_search_results(results, self.empty_reaction_df)
        
        self.assertEqual(len(processed), 1)
        self.assertEqual(processed[0]['error'], 'Invalid SMILES')
        self.assertEqual(processed[0]['target_smiles'], 'invalid_smiles')
    
    def test_process_results_success(self):
        """Test processing successful search results."""
        results = [
            {
                'query_smiles': 'CC',
                'smiles': 'CC',
                'space': 'sREAL',
                'sntDt': '2024-01-01',
                'exact': True,
                'vSynt': [
                    {
                        'rgn': [{'code': 'EN1', 'smiles': 'C'}],
                        'rsn': ['s100____111']
                    }
                ]
            }
        ]
        
        processed = self.process_search_results(results, self.empty_reaction_df)
        
        self.assertEqual(len(processed), 1)
        self.assertIsNone(processed[0]['error'])
        self.assertEqual(processed[0]['target_smiles'], 'CC')
        self.assertEqual(processed[0]['space'], 'sREAL')
        self.assertEqual(len(processed[0]['building_blocks']), 1)


if __name__ == '__main__':
    unittest.main()
