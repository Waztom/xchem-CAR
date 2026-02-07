#!/usr/bin/env python3
"""
Quick test script for the Enamine REAL Tools processor.

This script demonstrates how to run the prototype with a small subset of data
to test the functionality before processing the full dataset.
"""

import os
import sys
import pandas as pd
from pathlib import Path

def test_small_batch():
    """Test with just the first 5 molecules."""
    
    # Check if we're in the right directory
    current_dir = Path.cwd()
    if not (current_dir / "postera-ver1.5.csv").exists():
        print("Please run this script from the enamine directory containing postera-ver1.5.csv")
        return
    
    # Check API key
    if not os.environ.get('ENAMINE_REAL_TOOLS_API_KEY'):
        print("⚠️  Warning: ENAMINE_REAL_TOOLS_API_KEY not set")
        print("   Set it with: export ENAMINE_REAL_TOOLS_API_KEY='your_key_here'")
        return
    
    print("🧪 Testing Enamine REAL Tools processor with small batch...")
    
    # Load the input file and take just first 5 rows
    input_df = pd.read_csv("postera-ver1.5.csv")
    test_df = input_df.head(5)
    
    # Save small test file
    test_file = "test_small.csv"
    test_df.to_csv(test_file, index=False)
    
    print(f"📁 Created test file with {len(test_df)} molecules:")
    for i, row in test_df.iterrows():
        print(f"   {row['target-names']}: {row['target-SMILES']}")
    
    # Run the processor
    print("\n🚀 Running processor...")
    cmd = f"python process_targets_prototype.py {test_file} --batch-size 5 --output test_results.csv"
    print(f"📝 Command: {cmd}")
    
    # Execute the command
    exit_code = os.system(cmd)
    
    if exit_code == 0:
        print("\n✅ Processing completed successfully!")
        
        # Show results summary
        if Path("test_results.csv").exists():
            results_df = pd.read_csv("test_results.csv")
            print(f"\n📊 Results Summary:")
            print(f"   Total output rows: {len(results_df)}")
            print(f"   Statuses: {results_df['status'].value_counts().to_dict()}")
            
            # Show first few results
            print(f"\n🔍 First few results:")
            for i, row in results_df.head(3).iterrows():
                if row['status'] == 'FOUND':
                    print(f"   ✅ {row['target_smiles'][:50]}...")
                    print(f"      Building block: {row['building_block_code']} - {row['building_block_smiles'][:40]}...")
                    print(f"      Reaction: {row['reaction_name']} (ID: {row['reaction_id']})")
                elif row['status'] == 'ERROR':
                    print(f"   ❌ {row['target_smiles'][:50]}...")
                    print(f"      Error: {row['error']}")
                else:
                    print(f"   ⚪ {row['target_smiles'][:50]}... (NOT_FOUND)")
        
        print(f"\n📄 Full results saved to: test_results.csv")
        
    else:
        print(f"\n❌ Processing failed with exit code {exit_code}")
    
    # Cleanup small test file
    if Path(test_file).exists():
        os.remove(test_file)

def run_full_batch():
    """Run the full dataset processing."""
    
    print("🚀 Processing full Postera dataset...")
    
    cmd = "python process_targets_prototype.py postera-ver1.5.csv --batch-size 1000 --output postera_results.csv"
    print(f"📝 Command: {cmd}")
    
    print("⏳ This may take several minutes depending on API rate limits...")
    exit_code = os.system(cmd) 
    
    if exit_code == 0:
        print("✅ Full processing completed!")
        print("📄 Results saved to: postera_results.csv")
    else:
        print(f"❌ Processing failed with exit code {exit_code}")

def main():
    """Main function with user choice."""
    
    print("🧬 Enamine REAL Tools Processor Test")
    print("=" * 40)
    
    while True:
        print("\nOptions:")
        print("1. Test with 5 molecules (recommended first)")
        print("2. Process full dataset (103 molecules)")
        print("3. Exit")
        
        choice = input("\nSelect option (1/2/3): ").strip()
        
        if choice == "1":
            test_small_batch()
            break
        elif choice == "2":
            confirm = input("Process all 103 molecules? This may take time (y/n): ").strip().lower()
            if confirm in ['y', 'yes']:
                run_full_batch()
            break
        elif choice == "3":
            print("👋 Goodbye!")
            break
        else:
            print("❌ Invalid option. Please enter 1, 2, or 3.")

if __name__ == "__main__":
    main()