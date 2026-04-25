#!/usr/bin/env python3
"""
SMIRKS Regression Test
======================
For each FOUND row in the Enamine REAL Tools processed CSVs, this script:

1. Looks up the SMIRKS for the reaction_name via the mappings in
   prepare_syndirella_input.py.
2. Desalts the reactant SMILES (same as the main pipeline).
3. Runs _validate_smirks() which applies the RDKit reaction on the
   reactants and checks against target, found_smiles, and deprotected
   forms.
4. Reports per-file and per-reaction-type pass/fail statistics.
"""

import os
import sys
import pandas as pd
from collections import defaultdict

from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.ERROR)

# Import from the module under test
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from prepare_syndirella_input import (
    ENAMINE_REACTION_TO_CLASS,
    REACTION_CLASS_SMIRKS,
    _desalt_smiles,
    _refine_reaction_class,
    _validate_smirks,
    get_smirks_for_reaction,
)


def test_row(row):
    """
    Test a single row from a processed CSV using _validate_smirks().
    Returns (status, detail).
    """
    rxn_name = str(row.get('reaction_name', '')).strip()
    if not rxn_name or rxn_name == 'nan':
        return 'no_reaction', 'empty reaction_name'

    # Get reactants
    r1_raw = str(row.get('reactant_1_smiles', '')).strip()
    r2_raw = str(row.get('reactant_2_smiles', '')).strip()
    r3_raw = str(row.get('reactant_3_smiles', '')).strip()
    r4_raw = str(row.get('reactant_4_smiles', '')).strip()

    # Target and found SMILES
    target_smi = str(row.get('target_smiles', '')).strip()
    found_smi = str(row.get('found_smiles', '')).strip()

    # Desalt reactants (same as pipeline)
    r1 = _desalt_smiles(r1_raw)
    r2 = _desalt_smiles(r2_raw)

    if not r1 or r1 == 'nan' or not r2 or r2 == 'nan':
        return 'bad_reactants', f'missing reactants: r1={r1_raw}, r2={r2_raw}'

    # Look up SMIRKS
    rxn_class, smirks = get_smirks_for_reaction(rxn_name)
    if not rxn_class:
        return 'no_reaction', f'unmapped reaction: {rxn_name}'

    # Refine class based on actual reactants
    refined = _refine_reaction_class(rxn_class, r1, r2)
    if refined != rxn_class:
        smirks = REACTION_CLASS_SMIRKS.get(refined, smirks)

    has_r3 = r3_raw and r3_raw != 'nan'
    has_r4 = r4_raw and r4_raw != 'nan'
    r3 = _desalt_smiles(r3_raw) if has_r3 else None
    r4 = _desalt_smiles(r4_raw) if has_r4 else None

    status, product = _validate_smirks(
        smirks, r1, r2, target_smi,
        r3=r3, r4=r4,
        found_smi=found_smi if found_smi and found_smi != 'nan' else None,
    )

    # Map _validate_smirks statuses to regression test statuses
    status_map = {
        'pass': 'pass',
        'pass_found': 'pass',
        'pass_deprotected': 'pass',
        'mismatch': 'mismatch',
        'no_products': 'no_products',
        'no_smirks': 'no_smirks',
        'multi_step1_pass': 'multi_fires',
        'multi_step1_fail': 'multi_no_product',
        'error': 'bad_reactants',
    }
    mapped = status_map.get(status, status)

    detail = f'rxn={rxn_name} class={refined} status={status}'
    if mapped in ('no_products', 'mismatch', 'multi_no_product'):
        detail += (
            f'\n  r1={r1}\n  r2={r2}\n  target={target_smi}'
            + (f'\n  product={product}' if product else '')
        )

    return mapped, detail


def test_file(filepath):
    """Test all FOUND rows in a processed CSV. Returns stats dict."""
    df = pd.read_csv(filepath)
    found = df[df['status'] == 'FOUND'].copy()

    stats = defaultdict(int)
    rxn_stats = defaultdict(lambda: defaultdict(int))
    failures = []

    for idx, row in found.iterrows():
        status, detail = test_row(row)
        stats[status] += 1
        rxn_name = str(row.get('reaction_name', '')).strip()
        rxn_stats[rxn_name][status] += 1

        if status in ('no_products', 'mismatch', 'multi_no_product'):
            failures.append((idx, rxn_name, status, detail))

    return dict(stats), dict(rxn_stats), failures, len(found)


def main():
    processed_files = [
        '/container/src/output-files/CpKRS-resynthesis-sanity-check_processed.csv',
        '/container/src/output-files/EV2A-crudes-activity-cliff_processed.csv',
        '/container/src/output-files/Zika-crudes-activity-cliff_processed.csv',
        '/container/src/output-files/postera-ver1.5_processed.csv',
    ]

    grand_stats = defaultdict(int)
    grand_total = 0
    all_failures = []

    for filepath in processed_files:
        basename = os.path.basename(filepath)
        if not os.path.exists(filepath):
            print(f"SKIP: {basename} (not found)")
            continue

        print(f"\n{'='*70}")
        print(f"Testing: {basename}")
        print(f"{'='*70}")

        stats, rxn_stats, failures, total = test_file(filepath)
        grand_total += total

        # Per-reaction summary
        print(f"\n  {'Reaction':<40s} {'pass':>5s} {'mis':>5s} "
              f"{'noprod':>6s} {'nosmi':>5s} {'multi':>5s}")
        print(f"  {'-'*40} {'-'*5} {'-'*5} {'-'*6} {'-'*5} {'-'*5}")
        for rxn_name in sorted(rxn_stats.keys()):
            rs = rxn_stats[rxn_name]
            p = rs.get('pass', 0)
            m = rs.get('mismatch', 0)
            np_ = rs.get('no_products', 0)
            ns = rs.get('no_smirks', 0)
            mf = rs.get('multi_fires', 0)
            mnp = rs.get('multi_no_product', 0)
            multi_str = f"{mf}" + (f"/{mf+mnp}" if mnp else "")
            print(f"  {rxn_name:<40s} {p:>5d} {m:>5d} "
                  f"{np_:>6d} {ns:>5d} {multi_str:>5s}")

        # File summary
        n_pass = stats.get('pass', 0)
        n_mismatch = stats.get('mismatch', 0)
        n_no_prod = stats.get('no_products', 0)
        n_no_smirks = stats.get('no_smirks', 0)
        n_no_rxn = stats.get('no_reaction', 0)
        n_multi_fires = stats.get('multi_fires', 0)
        n_multi_noprod = stats.get('multi_no_product', 0)
        n_bad = stats.get('bad_reactants', 0)

        two_comp_total = n_pass + n_mismatch + n_no_prod
        two_comp_rate = (n_pass / two_comp_total * 100) if two_comp_total > 0 else 0

        print(f"\n  SUMMARY ({basename}):")
        print(f"    Total FOUND rows     : {total}")
        print(f"    2-comp pass          : {n_pass}")
        print(f"    2-comp mismatch      : {n_mismatch}")
        print(f"    2-comp no products   : {n_no_prod}")
        print(f"    2-comp pass rate     : {two_comp_rate:.1f}%")
        print(f"    Multi-comp fires     : {n_multi_fires}")
        print(f"    Multi-comp no product: {n_multi_noprod}")
        print(f"    No SMIRKS defined    : {n_no_smirks}")
        print(f"    No reaction mapping  : {n_no_rxn}")
        print(f"    Bad reactants        : {n_bad}")

        for k, v in stats.items():
            grand_stats[k] += v
        all_failures.extend(failures)

    # === Grand summary ===
    print(f"\n{'='*70}")
    print(f"GRAND SUMMARY")
    print(f"{'='*70}")
    print(f"  Total FOUND rows tested: {grand_total}")
    for k in ['pass', 'mismatch', 'no_products', 'multi_fires',
              'multi_no_product', 'no_smirks', 'no_reaction', 'bad_reactants']:
        print(f"  {k:<25s}: {grand_stats.get(k, 0)}")

    two_comp = grand_stats.get('pass', 0) + grand_stats.get('mismatch', 0) + grand_stats.get('no_products', 0)
    if two_comp > 0:
        rate = grand_stats.get('pass', 0) / two_comp * 100
        print(f"\n  2-component pass rate: {grand_stats.get('pass',0)}/{two_comp} = {rate:.1f}%")

    multi = grand_stats.get('multi_fires', 0) + grand_stats.get('multi_no_product', 0)
    if multi > 0:
        mrate = grand_stats.get('multi_fires', 0) / multi * 100
        print(f"  Multi-comp fire rate : {grand_stats.get('multi_fires',0)}/{multi} = {mrate:.1f}%")

    # Show some failure examples
    if all_failures:
        print(f"\n{'='*70}")
        print(f"FAILURE EXAMPLES (first 20)")
        print(f"{'='*70}")
        for idx, rxn_name, status, detail in all_failures[:20]:
            print(f"\n  [{status}] row {idx}, {rxn_name}:")
            for line in detail.split('\n'):
                print(f"    {line}")

    # Exit code: 0 if pass rate >= 80%, else 1
    if two_comp > 0 and rate < 50:
        print(f"\nFAIL: pass rate {rate:.1f}% < 50%")
        sys.exit(1)
    else:
        print(f"\nOK")
        sys.exit(0)


if __name__ == '__main__':
    main()
