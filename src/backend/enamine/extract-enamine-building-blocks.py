"""
Extract Building Blocks and Reaction Centre SMARTS from Enamine REAL Output
===========================================================================
Parses the Enamine REAL Tools synthesis CSV output and extracts:
  - Unique building blocks (Enamine catalogue codes + SMILES)
  - The reaction type each building block participates in
  - Reaction centre SMARTS for each building block (the functional group
    in the original BB that reacts), suitable for finding analogues

The reaction centre SMARTS are derived by:
  1. Mapping each Enamine reaction type (all 190 in enamine-real-tools-rxns.csv)
     to a reaction class with known functional group roles
  2. Analysing the synthon SMARTS pattern ([U] attachment point) to determine
     the role of each building block in the reaction
  3. Validating via RDKit substructure matching against the BB SMILES

Output is formatted for downstream analogue-searching workflows (e.g. Enamine
REAL synthon searching, catalogue substructure search, etc.).

Usage:
    python extract-enamine-building-blocks.py <input_csv> [output_csv]

    If output_csv is not given, it defaults to <input_basename>_building-blocks.csv
"""

import sys
import os
import csv
import re
import argparse
from collections import OrderedDict

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.logger().setLevel(RDLogger.ERROR)


# ═══════════════════════════════════════════════════════════════════════════
# Reaction Centre SMARTS — Functional Group Definitions
# ═══════════════════════════════════════════════════════════════════════════
# Each "role" a building block can play is defined here with:
#   - rxn_centre_smarts: SMARTS for the reactive group in the ORIGINAL BB
#   - label: human-readable description
#   - synthon_indicators: substrings in the synthon SMARTS that identify
#     this role (checked in order; first match wins)

ROLE_DEFINITIONS = {
    # ── Nitrogen nucleophiles ─────────────────────────────────────────────
    'amine': {
        'rxn_centre_smarts': '[NX3;H1,H2,H3;!$(NC=O);!$(NS(=O)=O)]',
        'label': 'Primary or secondary amine (incl. ammonia)',
        'synthon_indicators': ['N[U]', 'N([U])', '[U]N'],
    },
    'primary_amine': {
        'rxn_centre_smarts': '[NX3;H2;!$(NC=O);!$(NS(=O)=O)]',
        'label': 'Primary amine',
        'synthon_indicators': ['N[U]'],
    },
    'secondary_amine': {
        'rxn_centre_smarts': '[NX3;H1;!$(NC=O);!$(NS(=O)=O)]',
        'label': 'Secondary amine',
        'synthon_indicators': ['N([U])', 'N[U]'],
    },
    'aromatic_amine': {
        'rxn_centre_smarts': '[NX3;H2,H1;!$(NC=O);!$(NS(=O)=O)]',
        'label': 'Aromatic amine',
        'synthon_indicators': ['N[U]', 'N([U])'],
    },
    'amino_acid': {
        'rxn_centre_smarts': '[NX3;H2,H1]',
        'label': 'Amino acid (amine group reacts)',
        'synthon_indicators': ['N[U]', 'N([U])'],
    },
    'hydrazine': {
        'rxn_centre_smarts': '[NX3;H2,H1][NX3]',
        'label': 'Hydrazine or hydrazide',
        'synthon_indicators': ['N[U]'],
    },
    'nh_sulfamide': {
        'rxn_centre_smarts': '[NX3;H1]S(=O)(=O)',
        'label': 'NH-Sulfamide (NH reacts)',
        'synthon_indicators': ['N[U]', 'N([U])'],
    },
    # ── Oxygen nucleophiles ───────────────────────────────────────────────
    'alcohol': {
        'rxn_centre_smarts': '[OX2H1]',
        'label': 'Alcohol or phenol (-OH)',
        'synthon_indicators': ['O[U]', '[U]O'],
    },
    'phenol': {
        'rxn_centre_smarts': '[OX2H1][c]',
        'label': 'Phenol (ArOH)',
        'synthon_indicators': ['O[U]', '[U]O'],
    },
    # ── Sulfur nucleophiles ───────────────────────────────────────────────
    'thiol': {
        'rxn_centre_smarts': '[SX2H1]',
        'label': 'Thiol (-SH)',
        'synthon_indicators': ['S[U]', '[U]S'],
    },
    # ── Carboxylic acid & acyl electrophiles ──────────────────────────────
    'carboxylic_acid': {
        'rxn_centre_smarts': '[$([CX3](=O)[OX2H1]),$([CX3](=O)[O-])]',
        'label': 'Carboxylic acid or carboxylate (-COOH / -COO-)',
        'synthon_indicators': ['C(=O)[U]', 'C(=O)O[U]', '[U]C(=O)'],
    },
    'acid_chloride': {
        'rxn_centre_smarts': '[$([CX3](=O)[Cl]),$([CX3](=O)[OX2][CX3]=O)]',
        'label': 'Acyl chloride or anhydride (-COCl / -CO-O-CO-)',
        'synthon_indicators': ['C(=O)[U]'],
    },
    'acryloyl': {
        'rxn_centre_smarts': '[CH2]=[CH][CX3](=O)[Cl]',
        'label': 'Acryloyl chloride (covalent warhead)',
        'synthon_indicators': ['C(=O)[U]'],
    },
    'cyclic_anhydride': {
        'rxn_centre_smarts': '[CX3](=O)[OX2][CX3](=O)',
        'label': 'Cyclic anhydride',
        'synthon_indicators': ['C(=O)[U]'],
    },
    'oxalic_acid_derivative': {
        'rxn_centre_smarts': '[CX3](=O)[CX3](=O)',
        'label': 'Oxalic acid derivative (oxamide precursor)',
        'synthon_indicators': ['C(=O)[U]', 'C(=O)C(=O)[U]'],
    },
    # ── Sulfonyl electrophiles ────────────────────────────────────────────
    'sulfonyl_halide': {
        'rxn_centre_smarts': '[SX4](=O)(=O)[Cl,F,$([OX2][c,C])]',
        'label': 'Sulfonyl halide or activated ester (-SO2Cl/-SO2F/-SO2OAr)',
        'synthon_indicators': ['S(=O)(=O)[U]', '[U]S(=O)(=O)'],
    },
    'vinyl_sulfonyl': {
        'rxn_centre_smarts': '[CH2]=[CH][SX4](=O)(=O)',
        'label': 'Vinyl sulfonyl (covalent warhead)',
        'synthon_indicators': ['S(=O)(=O)[U]'],
    },
    'sulfonyl_fluoride_acyl': {
        'rxn_centre_smarts': '[CX3](=O)[OX2H1]',
        'label': 'Carboxylic acid with sulfonyl fluoride',
        'synthon_indicators': ['C(=O)[U]'],
    },
    # ── Isocyanate / Carbamate / Urea formers ────────────────────────────
    'isocyanate': {
        'rxn_centre_smarts': '[NX2]=[CX2]=[OX1]',
        'label': 'Isocyanate (-N=C=O)',
        'synthon_indicators': ['N[U]', 'N([U])'],
    },
    'carbamoyl_chloride': {
        'rxn_centre_smarts': '[Cl][CX3](=O)[NX3]',
        'label': 'Carbamoyl chloride (ClCONR)',
        'synthon_indicators': ['N[U]', 'N([U])'],
    },
    'reactive_carbamate': {
        'rxn_centre_smarts': '[$([OX2][CX3](=O)[NX3]),$([OX2;H0][CX3](=O)),$([n][CX3](=O))]',
        'label': 'Reactive carbamate/carbonate/activated ester/CDI-adduct',
        'synthon_indicators': ['N[U]', 'O[U]', 'N([U])', '[U]N', '[U]O'],
    },
    'chloroformate': {
        'rxn_centre_smarts': '[Cl][CX3](=O)[OX2]',
        'label': 'Chloroformate (ClCOOR)',
        'synthon_indicators': ['O[U]', 'N[U]'],
    },
    # ── Isothiocyanate ────────────────────────────────────────────────────
    'isothiocyanate': {
        'rxn_centre_smarts': '[NX2]=[CX2]=[SX1]',
        'label': 'Isothiocyanate (-N=C=S)',
        'synthon_indicators': ['N[U]'],
    },
    # ── Carbonyl (for reductive amination / condensation) ─────────────────
    'aldehyde': {
        'rxn_centre_smarts': '[CX3H1](=O)',
        'label': 'Aldehyde (-CHO)',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    'ketone': {
        'rxn_centre_smarts': '[#6][CX3](=O)[#6]',
        'label': 'Ketone (R-CO-R)',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    'aldehyde_or_ketone': {
        'rxn_centre_smarts': '[CX3](=O)',
        'label': 'Aldehyde or ketone',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    'methylene_active': {
        'rxn_centre_smarts': '[CX4H2]([CX3]=O)',
        'label': 'Active methylene compound',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    # ── Alkyl electrophiles ───────────────────────────────────────────────
    'alkyl_halide': {
        'rxn_centre_smarts': '[CX4][Cl,Br,I,$([OX2]S(=O)(=O))]',
        'label': 'Alkyl halide/sulfonate (-CH2X / -CH2OMs/OTs)',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    'epoxide': {
        'rxn_centre_smarts': '[CX4]1[OX2][CX4]1',
        'label': 'Epoxide (oxirane)',
        'synthon_indicators': ['C[U]', '[U]C', 'O[U]'],
    },
    # ── Aryl halides ──────────────────────────────────────────────────────
    'aryl_halide': {
        'rxn_centre_smarts': '[c][F,Cl,Br,I]',
        'label': 'Aryl halide (Ar-X)',
        'synthon_indicators': [],  # classified by BB SMILES
    },
    'aryl_dihalide': {
        'rxn_centre_smarts': '[c][F,Cl,Br,I]',
        'label': 'Aryl dihalide (for double arylation)',
        'synthon_indicators': [],
    },
    # ── Boronic acids / esters (Suzuki) ───────────────────────────────────
    'boronic_acid': {
        'rxn_centre_smarts': '[$([#6]B([OX2])[OX2]),$([#6][B-](F)(F)F),$([#6]B(O)O)]',
        'label': 'Boronic acid, pinacol ester, or trifluoroborate',
        'synthon_indicators': [],  # classified by BB SMILES (has B)
    },
    # ── Azide / Alkyne (Click) ────────────────────────────────────────────
    'azide': {
        'rxn_centre_smarts': '[NX1]~[NX2]~[NX2]',
        'label': 'Azide (-N3)',
        'synthon_indicators': ['N[U]'],
    },
    'terminal_alkyne': {
        'rxn_centre_smarts': '[CX2]#[CX2H1]',
        'label': 'Terminal alkyne (-C≡CH)',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    # ── Amidoxime / Nitrile (oxadiazole) ──────────────────────────────────
    'amidoxime': {
        'rxn_centre_smarts': '[NX3]C(=[NX2])[OX2H1]',
        'label': 'Amidoxime (for oxadiazole formation)',
        'synthon_indicators': ['N[U]', 'C[U]'],
    },
    'nitrile': {
        'rxn_centre_smarts': '[CX2]#[NX1]',
        'label': 'Nitrile (-C≡N, for oxadiazole/tetrazole)',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    # ── Thioamide / Alpha-haloketone (thiazole) ──────────────────────────
    'thioamide': {
        'rxn_centre_smarts': '[CX3](=S)[NX3,SX2]',
        'label': 'Thioamide or thiourea (for thiazole)',
        'synthon_indicators': ['N[U]', 'S[U]'],
    },
    'alpha_haloketone': {
        'rxn_centre_smarts': '[CX4][Cl,Br][CX3]=O',
        'label': 'Alpha-haloketone (for thiazole)',
        'synthon_indicators': ['C[U]'],
    },
    'thiosemicarbazide': {
        'rxn_centre_smarts': '[NX3]NC(=S)[NX3]',
        'label': 'Thiosemicarbazide',
        'synthon_indicators': ['N[U]'],
    },
    # ── Isonitrile (for Ugi-type MCRs) ────────────────────────────────────
    'isonitrile': {
        'rxn_centre_smarts': '[CX1-]#[NX2+]',
        'label': 'Isonitrile/isocyanide (for Ugi)',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
    # ── Alpha-amino ester (hydantoin) ─────────────────────────────────────
    'alpha_amino_ester': {
        'rxn_centre_smarts': '[NX3][CX4][CX3](=O)[OX2]',
        'label': 'Alpha-amino ester (for hydantoin)',
        'synthon_indicators': ['N[U]', 'C[U]'],
    },
    # ── Lactim ester (triazole) ───────────────────────────────────────────
    'lactim_ester': {
        'rxn_centre_smarts': '[OX2][CX3]=[NX2]',
        'label': 'Lactim ester (for triazole)',
        'synthon_indicators': ['N[U]', 'C[U]'],
    },
    'hydrazide': {
        'rxn_centre_smarts': '[NX3;H2,H1]NC(=O)',
        'label': 'Hydrazide (for triazole/heterocycle)',
        'synthon_indicators': ['N[U]'],
    },
    # ── Generic / Fallback ────────────────────────────────────────────────
    'generic_nucleophile': {
        'rxn_centre_smarts': '[NX3;H2,H1,OX2H1,SX2H1]',
        'label': 'Nucleophile (N/O/S)',
        'synthon_indicators': ['N[U]', 'O[U]', 'S[U]', '[U]S', '[U]N', '[U]O', 'N([U])'],
    },
    'generic_electrophile': {
        'rxn_centre_smarts': '[CX4][Cl,Br,I,$([OX2]S(=O)(=O))]',
        'label': 'Electrophile (alkyl halide/sulfonate)',
        'synthon_indicators': ['C[U]', '[U]C'],
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# Reaction Class Definitions
# ═══════════════════════════════════════════════════════════════════════════
# Each class maps role names from ROLE_DEFINITIONS.  The classifier tries
# role_a indicators first, then role_b, etc.  For Suzuki, BB-SMILES analysis
# is used instead.

REACTION_CLASSES = {
    # ── AMIDATION (coupling reagent: acid + amine → amide) ────────────────
    'amidation': {
        'description': 'Amide bond formation (carboxylic acid + amine)',
        'role_a': 'amine',
        'role_b': 'carboxylic_acid',
    },
    # ── ACETYLATION (acid chloride/anhydride + amine → amide) ─────────────
    'acetylation': {
        'description': 'Acylation with acid chloride or anhydride',
        'role_a': 'amine',
        'role_b': 'acid_chloride',
        'classify_by_bb_smiles_acyl': True,
    },
    # ── ACYLATION of O,C nucleophiles ─────────────────────────────────────
    'acylation_oc': {
        'description': 'Acylation of O,C-nucleophiles (ester bond)',
        'role_a': 'alcohol',
        'role_b': 'carboxylic_acid',
    },
    # ── ACRYLAMIDE ────────────────────────────────────────────────────────
    'acrylamide': {
        'description': 'Acylation with acryloyl chloride (covalent binder)',
        'role_a': 'amine',
        'role_b': 'acryloyl',
    },
    # ── ANHYDRIDE ACYLATION ───────────────────────────────────────────────
    'anhydride': {
        'description': 'Acylation of amine with cyclic anhydride',
        'role_a': 'amine',
        'role_b': 'cyclic_anhydride',
    },
    # ── ACYLSULFAMIDE (acid + NH-sulfamide) ───────────────────────────────
    'acylsulfamide': {
        'description': 'Acylation of NH-sulfamides',
        'role_a': 'nh_sulfamide',
        'role_b': 'carboxylic_acid',
    },
    # ── SULFONAMIDE (sulfonyl halide + amine) ─────────────────────────────
    'sulfonamide': {
        'description': 'Sulfonamide formation (sulfonyl halide + amine)',
        'role_a': 'amine',
        'role_b': 'sulfonyl_halide',
    },
    # ── SULFONAMIDE of PHENOL ─────────────────────────────────────────────
    'sulfonamide_phenol': {
        'description': 'Sulfacylation of phenol (sulfonyl halide + ArOH)',
        'role_a': 'phenol',
        'role_b': 'sulfonyl_halide',
    },
    # ── VINYL SULFONAMIDE ─────────────────────────────────────────────────
    'vinyl_sulfonamide': {
        'description': 'Vinylsulfonamide formation (covalent warhead)',
        'role_a': 'amine',
        'role_b': 'vinyl_sulfonyl',
    },
    # ── REDUCTIVE AMINATION ───────────────────────────────────────────────
    'reductive_amination': {
        'description': 'Reductive amination (amine + aldehyde/ketone)',
        'role_a': 'amine',
        'role_b': 'aldehyde_or_ketone',
    },
    'reductive_amination_ketone': {
        'description': 'Reductive amination of sterically hindered ketones',
        'role_a': 'primary_amine',
        'role_b': 'ketone',
    },
    # ── SUZUKI COUPLING ───────────────────────────────────────────────────
    'suzuki': {
        'description': 'Suzuki coupling (aryl halide + boronic acid/ester)',
        'role_a': 'boronic_acid',
        'role_b': 'aryl_halide',
        'classify_by_bb_smiles': True,
    },
    # ── N,O,S-ARYLATION ──────────────────────────────────────────────────
    'arylation': {
        'description': 'N,O,S-arylation (nucleophile + aryl halide)',
        'role_a': 'amine',
        'role_b': 'aryl_halide',
    },
    # ── ALKYLATION of amines ──────────────────────────────────────────────
    'alkylation_amine': {
        'description': 'N-alkylation of amines',
        'role_a': 'amine',
        'role_b': 'alkyl_halide',
    },
    # ── ALKYLATION (general N/O/S nucleophile + alkyl halide) ─────────────
    'alkylation_nos': {
        'description': 'Alkylation of N/O/S nucleophiles',
        'role_a': 'generic_nucleophile',
        'role_b': 'alkyl_halide',
    },
    # ── ALKYLATION of acid / ester formation ──────────────────────────────
    'alkylation_acid': {
        'description': 'Alkylation of carboxylic acid (ester formation)',
        'role_a': 'carboxylic_acid',
        'role_b': 'alkyl_halide',
    },
    # ── ALKYLATION of thiols / sulfides ───────────────────────────────────
    'alkylation_thiol': {
        'description': 'S-alkylation (thiol + alkyl halide)',
        'role_a': 'thiol',
        'role_b': 'alkyl_halide',
    },
    # ── EPOXIDE ring-opening ──────────────────────────────────────────────
    'epoxide_amine': {
        'description': 'Epoxide ring-opening with amine',
        'role_a': 'amine',
        'role_b': 'epoxide',
    },
    'epoxide_phenol': {
        'description': 'Epoxide ring-opening with phenol',
        'role_a': 'phenol',
        'role_b': 'epoxide',
    },
    # ── UREA FORMATION ────────────────────────────────────────────────────
    'urea': {
        'description': 'Urea formation (amine + amine via CDI/carbamate)',
        'role_a': 'amine',
        'role_b': 'amine',
        'classify_by_bb_smiles_urea_carbamate': True,
    },
    'urea_isocyanate': {
        'description': 'Urea formation (amine + isocyanate)',
        'role_a': 'amine',
        'role_b': 'isocyanate',
        'classify_by_bb_smiles_urea_iso': True,
    },
    'urea_carbamate': {
        'description': 'Urea formation from reactive carbamate',
        'role_a': 'amine',
        'role_b': 'reactive_carbamate',
        'classify_by_bb_smiles_urea_carbamate': True,
    },
    # ── CARBAMATE / URETHANE FORMATION ────────────────────────────────────
    'carbamate_cdi': {
        'description': 'Carbamate formation (amine + alcohol via CDI)',
        'role_a': 'amine',
        'role_b': 'alcohol',
    },
    'carbamate_chloroformate': {
        'description': 'Carbamate from amine + chloroformate',
        'role_a': 'amine',
        'role_b': 'chloroformate',
        'classify_by_bb_smiles_urea_carbamate': True,
    },
    # ── THIOUREA ──────────────────────────────────────────────────────────
    'thiourea': {
        'description': 'Thiourea formation (amine + isothiocyanate)',
        'role_a': 'amine',
        'role_b': 'isothiocyanate',
    },
    # ── OXAMIDE ───────────────────────────────────────────────────────────
    'oxamide': {
        'description': 'Disubstituted oxamide formation',
        'role_a': 'amine',
        'role_b': 'oxalic_acid_derivative',
    },
    # ── HYDRAZONE CONDENSATION ────────────────────────────────────────────
    'hydrazone': {
        'description': 'Hydrazone formation (hydrazine/hydrazide + carbonyl)',
        'role_a': 'hydrazine',
        'role_b': 'aldehyde_or_ketone',
    },
    # ── KNOEVENAGEL CONDENSATION ──────────────────────────────────────────
    'knoevenagel': {
        'description': 'Knoevenagel condensation (active methylene + aldehyde)',
        'role_a': 'methylene_active',
        'role_b': 'aldehyde',
    },
    # ── MANNICH REACTION ──────────────────────────────────────────────────
    'mannich': {
        'description': 'Mannich reaction',
        'role_a': 'amine',
        'role_b': 'aldehyde_or_ketone',
    },
    # ── CLICK / TRIAZOLE ──────────────────────────────────────────────────
    'click_triazole': {
        'description': '1,2,3-Triazole (azide + terminal alkyne, CuAAC)',
        'role_a': 'azide',
        'role_b': 'terminal_alkyne',
    },
    'click_from_amine': {
        'description': 'Triazole: in-situ azide from amine + terminal alkyne',
        'role_a': 'primary_amine',
        'role_b': 'terminal_alkyne',
    },
    # ── TRIAZOLE from lactim ester + hydrazide ────────────────────────────
    'triazole_lactim': {
        'description': 'Triazole from lactim ester + hydrazide',
        'role_a': 'lactim_ester',
        'role_b': 'hydrazide',
    },
    # ── OXADIAZOLE ────────────────────────────────────────────────────────
    'oxadiazole_amidoxime': {
        'description': '1,2,4-Oxadiazole from amidoxime + carboxylic acid',
        'role_a': 'amidoxime',
        'role_b': 'carboxylic_acid',
    },
    'oxadiazole_nitrile': {
        'description': '1,2,4-Oxadiazole from nitrile + carboxylic acid',
        'role_a': 'nitrile',
        'role_b': 'carboxylic_acid',
    },
    'oxadiazole_cyclization': {
        'description': 'Oxadiazole cyclization',
        'role_a': 'generic_nucleophile',
        'role_b': 'carboxylic_acid',
    },
    'aminoxadiazole': {
        'description': '2-Amino-1,3,4-oxadiazole from hydrazide + amine (Thio-CDI)',
        'role_a': 'hydrazide',
        'role_b': 'amine',
    },
    # ── THIAZOLE FORMATION ────────────────────────────────────────────────
    'thiazole_2comp': {
        'description': 'Thiazole from thioamide/thiourea + alpha-haloketone',
        'role_a': 'thioamide',
        'role_b': 'alpha_haloketone',
    },
    'thiazole_3comp': {
        'description': 'Thiazole from thiosemicarbazide + aldehyde + haloketone',
        'role_a': 'thiosemicarbazide',
        'role_b': 'aldehyde',
        'role_c': 'alpha_haloketone',
    },
    'thiazole_thiourea_3comp': {
        'description': 'Thiourea + thiazole cyclization (3-component)',
        'role_a': 'isothiocyanate',
        'role_b': 'amine',
        'role_c': 'alpha_haloketone',
    },
    # ── HYDANTOIN ─────────────────────────────────────────────────────────
    'hydantoin': {
        'description': 'Hydantoin from amine + alpha-amino ester',
        'role_a': 'primary_amine',
        'role_b': 'alpha_amino_ester',
    },
    # ── BIGINELLI ─────────────────────────────────────────────────────────
    'biginelli': {
        'description': 'Biginelli reaction (urea + aldehyde + beta-ketoester)',
        'role_a': 'amine',
        'role_b': 'aldehyde',
        'role_c': 'methylene_active',
    },
    # ── CYCLIC GUANIDINE ──────────────────────────────────────────────────
    'cyclic_guanidine': {
        'description': 'Cyclic guanidine (thiouronium) formation',
        'role_a': 'amine',
        'role_b': 'thioamide',
    },
    # ── TETRAZOLE ─────────────────────────────────────────────────────────
    'tetrazole_acylation': {
        'description': 'Acylation + tetrazole cyclization of nitrile',
        'role_a': 'amine',
        'role_b': 'nitrile',
    },
    'tetrazole_alkylation': {
        'description': 'Alkylation + tetrazole cyclization',
        'role_a': 'amine',
        'role_b': 'alkyl_halide',
    },
    'tetrazole_sulfonamide': {
        'description': 'Sulfonamide + tetrazole cyclization',
        'role_a': 'amine',
        'role_b': 'sulfonyl_halide',
    },
    # ── UGI MCR ───────────────────────────────────────────────────────────
    'ugi_4cr': {
        'description': 'Ugi 4-component reaction',
        'role_a': 'amine',
        'role_b': 'aldehyde',
        'role_c': 'carboxylic_acid',
        'role_d': 'isonitrile',
    },
    'ugi_3cr': {
        'description': 'Ugi 3-component reaction',
        'role_a': 'amine',
        'role_b': 'aldehyde',
        'role_c': 'isonitrile',
    },
    # ── GBB MCR ───────────────────────────────────────────────────────────
    'gbb': {
        'description': 'Groebke-Blackburn-Bienayme 3CR',
        'role_a': 'amine',
        'role_b': 'aldehyde',
        'role_c': 'isonitrile',
    },
    # ── MULTI-COMPONENT sequential reactions ──────────────────────────────
    'multi_amide_arylation': {
        'description': 'Amidation + arylation',
        'role_a': 'amine', 'role_b': 'carboxylic_acid', 'role_c': 'aryl_halide',
    },
    'multi_amide_amide': {
        'description': 'Double amidation',
        'role_a': 'amine', 'role_b': 'carboxylic_acid', 'role_c': 'carboxylic_acid',
    },
    'multi_sulfonamide_arylation': {
        'description': 'Sulfonamide + arylation',
        'role_a': 'amine', 'role_b': 'sulfonyl_halide', 'role_c': 'aryl_halide',
    },
    'multi_sulfonamide_suzuki': {
        'description': 'Sulfonamide + Suzuki',
        'role_a': 'amine', 'role_b': 'sulfonyl_halide', 'role_c': 'boronic_acid',
    },
    'multi_amide_suzuki': {
        'description': 'Amidation + Suzuki',
        'role_a': 'amine', 'role_b': 'carboxylic_acid', 'role_c': 'boronic_acid',
    },
    'multi_amide_click': {
        'description': 'Amidation + click triazole',
        'role_a': 'amine', 'role_b': 'carboxylic_acid', 'role_c': 'terminal_alkyne',
    },
    'multi_alkylation_click': {
        'description': 'Alkylation + click triazole',
        'role_a': 'amine', 'role_b': 'alkyl_halide', 'role_c': 'terminal_alkyne',
    },
    'multi_amide_alkylation': {
        'description': 'Amidation + alkylation',
        'role_a': 'amine', 'role_b': 'carboxylic_acid', 'role_c': 'alkyl_halide',
    },
    'multi_reductive_amination_amide': {
        'description': 'Reductive amination + amidation',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'carboxylic_acid',
    },
    'multi_reductive_amination_arylation': {
        'description': 'Reductive amination + arylation',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'aryl_halide',
    },
    'multi_reductive_amination_sulfonamide': {
        'description': 'Reductive amination + sulfonamide',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'sulfonyl_halide',
    },
    'multi_reductive_amination_urea': {
        'description': 'Reductive amination + urea',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'amine',
    },
    'multi_reductive_amination_carbamate': {
        'description': 'Reductive amination + carbamate',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'chloroformate',
    },
    'multi_reductive_amination_reductive_amination': {
        'description': 'Double reductive amination',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'aldehyde',
    },
    'multi_reductive_amination_acyl_halide': {
        'description': 'Reductive amination + acylation with acid halide',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'acid_chloride',
    },
    'multi_arylation_arylation': {
        'description': 'Double arylation',
        'role_a': 'amine', 'role_b': 'aryl_dihalide', 'role_c': 'amine',
    },
    'multi_arylation_suzuki': {
        'description': 'Arylation + Suzuki',
        'role_a': 'amine', 'role_b': 'aryl_halide', 'role_c': 'boronic_acid',
    },
    'multi_triazole_amide': {
        'description': 'Triazole + amidation',
        'role_a': 'lactim_ester', 'role_b': 'hydrazide', 'role_c': 'carboxylic_acid',
    },
    'multi_triazole_urea': {
        'description': 'Triazole + urea',
        'role_a': 'lactim_ester', 'role_b': 'hydrazide', 'role_c': 'amine',
    },
    'multi_triazole_alkylation': {
        'description': 'Triazole + alkylation',
        'role_a': 'lactim_ester', 'role_b': 'hydrazide', 'role_c': 'alkyl_halide',
    },
    'multi_triazole_sulfonamide': {
        'description': 'Triazole + sulfonamide',
        'role_a': 'lactim_ester', 'role_b': 'hydrazide', 'role_c': 'sulfonyl_halide',
    },
    'multi_thiazole_amide': {
        'description': 'Thiazole + amidation',
        'role_a': 'thioamide', 'role_b': 'alpha_haloketone', 'role_c': 'carboxylic_acid',
    },
    'multi_thiazole_sulfonamide': {
        'description': 'Thiazole + sulfonamide',
        'role_a': 'thioamide', 'role_b': 'alpha_haloketone', 'role_c': 'sulfonyl_halide',
    },
    'multi_anhydride_amide_arylation': {
        'description': 'Anhydride acylation + arylation',
        'role_a': 'amine', 'role_b': 'cyclic_anhydride', 'role_c': 'aryl_halide',
    },
    'multi_sulfonamide_amide': {
        'description': 'Sulfonamide + amidation',
        'role_a': 'amine', 'role_b': 'sulfonyl_halide', 'role_c': 'carboxylic_acid',
    },
    'multi_acylsulfamide_alkylation': {
        'description': 'Acyl-sulfamide + amine alkylation',
        'role_a': 'nh_sulfamide', 'role_b': 'acid_chloride', 'role_c': 'amine',
    },
    'multi_click_amide': {
        'description': 'Click triazole + amidation',
        'role_a': 'azide', 'role_b': 'terminal_alkyne', 'role_c': 'carboxylic_acid',
    },
    'multi_oxoisoindole_amide': {
        'description': 'Oxoisoindole formation + amidation',
        'role_a': 'amine', 'role_b': 'aldehyde', 'role_c': 'carboxylic_acid',
    },
    'dihydrouracil': {
        'description': 'Dihydrouracyl synthesis from two amines',
        'role_a': 'amine', 'role_b': 'amine', 'role_c': 'amine',
    },
    'purine_substitution': {
        'description': '6,9-disubstituted purine formation',
        'role_a': 'amine', 'role_b': 'aryl_halide', 'role_c': 'alkyl_halide',
    },
    # Generic fallbacks
    'generic_2comp': {
        'description': 'Unknown 2-component reaction',
        'role_a': 'generic_nucleophile', 'role_b': 'generic_electrophile',
    },
    'generic_3comp': {
        'description': 'Unknown 3-component reaction',
        'role_a': 'generic_nucleophile', 'role_b': 'generic_electrophile',
        'role_c': 'generic_electrophile',
    },
    'generic_4comp': {
        'description': 'Unknown 4-component reaction',
        'role_a': 'generic_nucleophile', 'role_b': 'generic_electrophile',
        'role_c': 'generic_electrophile', 'role_d': 'generic_electrophile',
    },
}

# ═══════════════════════════════════════════════════════════════════════════
# Product Bond SMARTS — Expected substructure in the target molecule
# ═══════════════════════════════════════════════════════════════════════════
# For each reaction class, the SMARTS pattern that MUST be present in the
# target/product molecule if the reaction happened correctly.

PRODUCT_BOND_SMARTS = {
    # ── Amide bond formers ─────────────────────────────────────────────────
    'amidation':                        '[#7][CX3](=O)',
    'acetylation':                      '[#7][CX3](=O)',
    'acylation_oc':                     '[OX2][CX3](=O)',
    'acrylamide':                       '[#7][CX3](=O)',
    'anhydride':                        '[#7][CX3](=O)',
    'acylsulfamide':                    '[#7]S(=O)(=O)[#7]',
    # ── Sulfonamide ───────────────────────────────────────────────────────
    'sulfonamide':                      '[#7]S(=O)(=O)',
    'sulfonamide_phenol':               '[OX2]S(=O)(=O)',
    'vinyl_sulfonamide':                '[#7]S(=O)(=O)',
    # ── Reductive amination ───────────────────────────────────────────────
    'reductive_amination':              '[#7][CX4]',
    'reductive_amination_ketone':       '[#7][CX4]',
    # ── Suzuki ────────────────────────────────────────────────────────────
    'suzuki':                           '[c]~[c]',
    # ── Arylation ─────────────────────────────────────────────────────────
    'arylation':                        '[#7,#8,#16][c]',
    # ── Alkylation ────────────────────────────────────────────────────────
    'alkylation_amine':                 '[#7][CX4]',
    'alkylation_nos':                   '[#7,#8,#16][CX4]',
    'alkylation_acid':                  '[OX2][CX4]',
    'alkylation_thiol':                 '[#16][CX4]',
    # ── Epoxide ───────────────────────────────────────────────────────────
    'epoxide_amine':                    '[#7][CX4][CX4][OX2]',
    'epoxide_phenol':                   '[OX2][CX4][CX4][OX2]',
    # ── Urea / Carbamate / Thiourea ───────────────────────────────────────
    'urea':                             '[#7]C(=O)[#7]',
    'urea_isocyanate':                  '[#7]C(=O)[#7]',
    'urea_carbamate':                   '[#7]C(=O)[#7,#8]',
    'carbamate_cdi':                    '[#7]C(=O)[OX2]',
    'carbamate_chloroformate':          '[#7]C(=O)[OX2]',
    'thiourea':                         '[#7]C(=S)[#7]',
    # ── Oxamide ───────────────────────────────────────────────────────────
    'oxamide':                          '[#7]C(=O)C(=O)[#7]',
    # ── Condensation ──────────────────────────────────────────────────────
    'hydrazone':                        '[#7][#7]=[CX3]',
    'knoevenagel':                      '[CX3]=[CX3]',
    'mannich':                          '[#7][CX4]',
    # ── Heterocycle formation ─────────────────────────────────────────────
    'click_triazole':                   '[n]:[n]:[n]',
    'click_from_amine':                 '[n]:[n]:[n]',
    'triazole_lactim':                  '[n]:[n]:[n]',
    'oxadiazole_amidoxime':             '[nX2]:[c]:[nX2]',
    'oxadiazole_nitrile':               '[nX2]:[c]:[nX2]',
    'oxadiazole_cyclization':           '[nX2]:[c]:[nX2]',
    'aminoxadiazole':                   '[nX2]:[c]:[nX2]',
    'thiazole_2comp':                   '[sX2]:[c]:[nX2]',
    'thiazole_3comp':                   '[sX2]:[c]:[nX2]',
    'thiazole_thiourea_3comp':          '[sX2]:[c]:[nX2]',
    'hydantoin':                        '[#7]C(=O)[#7]C(=O)',
    'biginelli':                        '[#7]C(=O)[#7]',
    'cyclic_guanidine':                 '[#7]C(=[#7])[#7]',
    'tetrazole_acylation':              '[n]:[n]:[n]:[n]',
    'tetrazole_alkylation':             '[n]:[n]:[n]:[n]',
    'tetrazole_sulfonamide':            '[n]:[n]:[n]:[n]',
    # ── MCR ───────────────────────────────────────────────────────────────
    'ugi_4cr':                          '[#7][CX3](=O)',
    'ugi_3cr':                          '[#7][CX3](=O)',
    'gbb':                              '[nX2]:[c]:[nX2]',
    # ── Multi-component (check the primary bond) ─────────────────────────
    'multi_amide_arylation':            '[#7][CX3](=O)',
    'multi_amide_amide':                '[#7][CX3](=O)',
    'multi_sulfonamide_arylation':      '[#7]S(=O)(=O)',
    'multi_sulfonamide_suzuki':         '[#7]S(=O)(=O)',
    'multi_amide_suzuki':               '[#7][CX3](=O)',
    'multi_amide_click':                '[#7][CX3](=O)',
    'multi_alkylation_click':           '[#7][CX4]',
    'multi_amide_alkylation':           '[#7][CX3](=O)',
    'multi_reductive_amination_amide':  '[#7][CX4]',
    'multi_reductive_amination_arylation': '[#7][CX4]',
    'multi_reductive_amination_sulfonamide': '[#7][CX4]',
    'multi_reductive_amination_urea':   '[#7][CX4]',
    'multi_reductive_amination_carbamate': '[#7][CX4]',
    'multi_reductive_amination_reductive_amination': '[#7][CX4]',
    'multi_reductive_amination_acyl_halide': '[#7][CX4]',
    'multi_arylation_arylation':        '[#7,#8,#16][c]',
    'multi_arylation_suzuki':           '[#7,#8,#16][c]',
    'multi_triazole_amide':             '[n]:[n]:[n]',
    'multi_triazole_urea':              '[n]:[n]:[n]',
    'multi_triazole_alkylation':        '[n]:[n]:[n]',
    'multi_triazole_sulfonamide':       '[n]:[n]:[n]',
    'multi_thiazole_amide':             '[sX2]:[c]:[nX2]',
    'multi_thiazole_sulfonamide':       '[sX2]:[c]:[nX2]',
    'multi_anhydride_amide_arylation':  '[#7][CX3](=O)',
    'multi_sulfonamide_amide':          '[#7]S(=O)(=O)',
    'multi_acylsulfamide_alkylation':   '[#7]S(=O)(=O)',
    'multi_click_amide':                '[n]:[n]:[n]',
    'multi_oxoisoindole_amide':         '[#7][CX3](=O)',
    'dihydrouracil':                    '[#7]C(=O)[#7]',
    'purine_substitution':              '[#7][c]',
}


# ═══════════════════════════════════════════════════════════════════════════
# Protecting Group Definitions
# ═══════════════════════════════════════════════════════════════════════════
# Enamine REAL reactions with prefixes Boc-, B-, K-, BBoc-, R039-Me3Si-
# involve a protecting group removal step.  These definitions let us
# detect the PG in the original BB SMILES and report a SMARTS pattern
# that enumeration workflows can use to exclude the PG moiety.
#
# pg_smarts  : SMARTS that matches the protecting group in the BB
# label      : human-readable description of the protecting group

PROTECTING_GROUP_INFO = {
    'boc': {
        'label': 'Boc (tert-butyloxycarbonyl)',
        'pg_smarts': '[#7]C(=O)OC([CH3])([CH3])[CH3]',
    },
    'tbu_ester': {
        'label': 'tert-butyl ester',
        'pg_smarts': '[CX3](=O)OC([CH3])([CH3])[CH3]',
    },
    'saponification': {
        'label': 'Alkyl ester (saponification)',
        'pg_smarts': '[CX3](=O)O[CX4;H1,H2,H3]',
    },
    'tms': {
        'label': 'TMS (trimethylsilyl)',
        'pg_smarts': '[Si]([CH3])([CH3])[CH3]',
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# Mapping: Enamine Reaction Name → Reaction Class
# ═══════════════════════════════════════════════════════════════════════════
# Every reaction from enamine-real-tools-rxns.csv is mapped here.

ENAMINE_REACTION_TO_CLASS = {
    # ── AMIDATION variants (acid + amine → amide) ─────────────────────────
    'AMIDE':              'amidation',       # 11  Mukaiyama acylation
    'AMIDE1':             'amidation',       # 22  EDC acylation
    'AMIDE2':             'amidation',       # 527  CDI acylation
    'AMIDEMSCl':          'amidation',       # 1626  Acylation with MsCl
    'K-AMIDE':            'amidation',       # 270112  Mukaiyama + saponification
    'K-AMIDE1':           'amidation',       # 188690  EDC + saponification
    'B-AMIDE':            'amidation',       # 272270  Mukaiyama + t-Bu cleavage
    'B-AMIDE1':           'amidation',       # 269946  EDC + t-Bu cleavage
    'B-AMIDE2':           'amidation',       # 273652  CDI + t-Bu cleavage
    'Boc-AMIDE':          'amidation',       # 272430  Mukaiyama + Boc deprotection
    'Boc-AMIDE1':         'amidation',       # 240690  EDC + Boc deprotection
    'Boc-AMIDE2':         'amidation',       # 270006  CDI + Boc deprotection
    'Boc-AMIDEMSCl':      'amidation',       # 282770  MsCl + Boc deprotection
    'BBoc-AMIDE1':        'amidation',       # 273450  EDC + Boc + t-Bu
    'AMINOACID1':         'amidation',       # 272610  Acylation of amino acids
    'R119-ACYLATION':     'amidation',       # 276750  Amides with sulfonyl fluoride
    'R039-Me3Si-Ac':      'acetylation',     # 274079  Acyl chloride + amino acid

    # ── ACETYLATION (acid chloride/anhydride + amine) ─────────────────────
    'ACETYLATION':        'acetylation',     # 270062

    # ── ACYLATION of O,C nucleophiles ─────────────────────────────────────
    'R111-ACYLATION':     'acylation_oc',    # 276436

    # ── ACRYLAMIDE ────────────────────────────────────────────────────────
    'R017-ACRYLAMIDE':    'acrylamide',      # 273270

    # ── ANHYDRIDE ─────────────────────────────────────────────────────────
    'ANHYDRIDE1':         'anhydride',       # 58

    # ── ACYLSULFAMIDE ─────────────────────────────────────────────────────
    'ACYLSULFAMIDE':      'acylsulfamide',   # 2740
    'Boc-ACYLSULFAMIDE':  'acylsulfamide',   # 273692
    'B-ACYLSULFAMIDE':    'acylsulfamide',   # 274614

    # ── SULFONAMIDE (sulfonyl halide + amine) ─────────────────────────────
    'SULFAMIDE':          'sulfonamide',     # 20
    'SULFAMIDE1':         'sulfonamide',     # 40
    'SULFAMIDE2':         'sulfonamide',     # 232682
    'SULFAMIDE4-ShB':     'sulfonamide',     # 271082
    'K-SULFAMIDE1':       'sulfonamide',     # 196680
    'B-SULFAMIDE1':       'sulfonamide',     # 270188
    'Boc-SULFAMIDE1':     'sulfonamide',     # 270084
    'Boc-SULFAMIDE2':     'sulfonamide',     # 273578
    'R039-Me3Si-Sac':     'sulfonamide',     # 274078
    'R165-F-SULFAMIDE':   'sulfonamide',     # 278578

    # ── SULFACYLATION of PHENOL ───────────────────────────────────────────
    'SULFACYLATION-O':    'sulfonamide_phenol',  # 87

    # ── VINYL SULFONAMIDE ─────────────────────────────────────────────────
    'VINYLSULFAMIDE':     'vinyl_sulfonamide',   # 282270

    # ── REDUCTIVE AMINATION ───────────────────────────────────────────────
    'REDUCTION':          'reductive_amination',         # 207
    'REDUCTION3':         'reductive_amination',         # 270004
    'REDUCTION-HET':      'reductive_amination',         # 270302
    'REDUCTION-AMINOACID': 'reductive_amination',        # 271304
    'R070-R3-AMINOACID':  'reductive_amination',         # 275532
    'Boc-REDUCTION3':     'reductive_amination',         # 271302
    'K-REDUCTION3':       'reductive_amination',         # 274370
    'B-REDUCTION3':       'reductive_amination',         # 273456
    'REDUCTION-TI':       'reductive_amination_ketone',  # 269862
    'Boc-REDUCTION-TI':   'reductive_amination_ketone',  # 272690

    # ── SUZUKI COUPLING ───────────────────────────────────────────────────
    'SUZUKI':             'suzuki',          # 271570
    'Boc-SUZUKI':         'suzuki',          # 273584
    'AC-SUZUKI':          'suzuki',          # 273712
    'K-SUZUKI':           'suzuki',          # 275652

    # ── N,O,S-ARYLATION ──────────────────────────────────────────────────
    'ARYLATION':          'arylation',       # 27
    'ARYLATION60':        'arylation',       # 88
    'Boc-ARYLATION':      'arylation',       # 270166
    'K-ARYLATION':        'arylation',       # 270344
    'B-ARYLATION':        'arylation',       # 272710

    # ── N-ALKYLATION of amines ────────────────────────────────────────────
    'MA1':                'alkylation_amine', # 38
    'MA3':                'alkylation_amine', # 44
    'MA2-POLY':           'alkylation_amine', # 2230
    'Boc-MA2':            'alkylation_amine', # 269982
    'K-MA2':              'alkylation_amine', # 270122
    'B-MA2':              'alkylation_amine', # 269956

    # ── ALKYLATION of N/O/S nucleophiles ──────────────────────────────────
    'ESTER-NS':           'alkylation_nos',  # 7
    'OXYME':              'alkylation_nos',   # 34
    'Boc-OXYME':          'alkylation_nos',   # 272692
    'B-OXYME':            'alkylation_nos',   # 273580
    'ESTER-PH':           'alkylation_nos',   # 2624
    'Boc-ESTER-NS':       'alkylation_nos',   # 273654
    'B-ESTER-NS':         'alkylation_nos',   # 270034

    # ── ALKYLATION of acid (ester formation) ──────────────────────────────
    'ESTER':              'alkylation_acid',  # 1458
    'B-ESTER':            'alkylation_acid',  # 279392

    # ── S-ALKYLATION / SULFIDE / SULFOXIDE / SULFONE ──────────────────────
    'SULFIDE':            'alkylation_thiol', # 62
    'SULFONE':            'alkylation_thiol', # 2714
    'SULFOXIDE':          'alkylation_thiol', # 265282
    'SULFIDE-OX1':        'alkylation_thiol', # 1498
    'SULFIDE-OX2':        'alkylation_thiol', # 1500

    # ── EPOXIDE RING-OPENING ──────────────────────────────────────────────
    'EPOXYDE':            'epoxide_amine',    # 63
    'EPOXYDE2':           'epoxide_amine',    # 273610
    'EPOXYDE3':           'epoxide_phenol',   # 273790

    # ── UREA FORMATION ────────────────────────────────────────────────────
    'UREA':               'urea_isocyanate',  # 68
    'UREA-CDI':           'urea',             # 512
    'URETANE':            'urea_carbamate',   # 487
    'URETANE1':           'urea_isocyanate',  # 2430
    'URETANE2':           'urea_carbamate',   # 2554
    'URETANE4':           'urea_isocyanate',  # 2708
    'URETANE-T':          'urea',             # 272164
    'Boc-URETANE':        'urea_carbamate',   # 273574
    'Boc-URETANE1':       'urea_isocyanate',  # 273458
    'Boc-URETANE4':       'urea_carbamate',   # 273464
    'B-URETANE1':         'urea_isocyanate',  # 273390
    'B-URETANE4':         'urea_isocyanate',  # 273392
    'R039-Me3Si-Urea':    'urea_isocyanate',  # 274080

    # ── CARBAMATE / URETHANE ──────────────────────────────────────────────
    'R157-CARBAMATE':     'carbamate_cdi',           # 277990
    'R192-CARBAMATE':     'carbamate_chloroformate', # 282070
    'R039-Me3Si-Uretane': 'carbamate_chloroformate', # 274081

    # ── THIOUREA ──────────────────────────────────────────────────────────
    'THYOUREA':           'thiourea',        # 8

    # ── OXAMIDE ───────────────────────────────────────────────────────────
    'OXAMIDE':            'oxamide',         # 2718
    'OXAMIDE1Pri':        'oxamide',         # 271948
    'OXAMIDE1Sec':        'oxamide',         # 271949
    'B-OXAMIDE':          'oxamide',         # 273460
    'Boc-OXAMIDE':        'oxamide',         # 273462
    'B-OXAMIDE1Pri':      'oxamide',         # 273492
    'B-OXAMIDE1Sec':      'oxamide',         # 273494
    'Boc-OXAMIDE1Pri':    'oxamide',         # 273496
    'Boc-OXAMIDE1Sec':    'oxamide',         # 273498

    # ── HYDRAZONE ─────────────────────────────────────────────────────────
    'HYDRAZON':           'hydrazone',       # 4

    # ── KNOEVENAGEL ───────────────────────────────────────────────────────
    'METHYLEN':           'knoevenagel',     # 6
    'METHYLEN1':          'knoevenagel',     # 28

    # ── MANNICH ───────────────────────────────────────────────────────────
    'MANNICH':            'mannich',         # 26

    # ── CLICK / TRIAZOLE ──────────────────────────────────────────────────
    'R030-TRIAZOLE':      'click_triazole',          # 273910
    'R040-TRIAZOLE':      'click_from_amine',        # 274090
    'Boc-R030-TRIAZOLE':  'click_triazole',          # 275090
    'K-R040-TRIAZOLE':    'click_from_amine',        # 275530
    'B-R040-TRIAZOLE':    'click_from_amine',        # 275534
    'Boc-R040-TRIAZOLE':  'click_from_amine',        # 275536
    'Radd-178-AT-Click':  'click_from_amine',        # 282512

    # ── TRIAZOLE from lactim ──────────────────────────────────────────────
    'R036-Lactim':        'triazole_lactim', # 274052

    # ── OXADIAZOLE ────────────────────────────────────────────────────────
    'AMIDOXIME1':         'oxadiazole_amidoxime',    # 265764
    'Boc-AMIDOXIME1':     'oxadiazole_amidoxime',    # 270036
    'AMIDOXIME1-CN':      'oxadiazole_nitrile',      # 270196
    'Boc-AMIDOXIME1-CN':  'oxadiazole_nitrile',      # 271362
    'ANONYMOUS3':         'oxadiazole_cyclization',   # 271722
    'R037a-AMINOXADIAZOLE': 'aminoxadiazole',         # 274134
    'Radd-117-Ac-CN2Oxdzln': 'oxadiazole_nitrile',   # 280212

    # ── THIAZOLE ──────────────────────────────────────────────────────────
    'THIAZOLE':           'thiazole_3comp',           # 1
    'THIAZOLE2':          'thiazole_thiourea_3comp',  # 2
    'THIAZOLE3':          'thiazole_thiourea_3comp',  # 3
    'THIAZOLE4':          'thiazole_2comp',           # 43

    # ── HYDANTOIN ─────────────────────────────────────────────────────────
    'HYDANTOIN4':         'hydantoin',       # 272104
    'HYDANTOIN6':         'hydantoin',       # 272212

    # ── BIGINELLI ─────────────────────────────────────────────────────────
    'BIGINELLI':          'biginelli',       # 12

    # ── CYCLIC GUANIDINE / THIOURONIUM ────────────────────────────────────
    'R113-CyclicTHIOURONIUM': 'cyclic_guanidine',    # 276492
    'THIOURONIUM1':       'cyclic_guanidine',         # 264724

    # ── TETRAZOLE ─────────────────────────────────────────────────────────
    'T-ACYLATION':        'tetrazole_acylation',     # 280130
    'T-MA2':              'tetrazole_alkylation',     # 280178
    'Radd-139-T-SO2N':    'tetrazole_sulfonamide',   # 280750

    # ── UGI MCR ───────────────────────────────────────────────────────────
    'Radd-152-Ugi':       'ugi_4cr',                 # 281210
    'Radd-173-Split-Ugi': 'ugi_4cr',                 # 282410
    'Radd-176-Ugi-MF':    'ugi_3cr',                 # 282490
    'Radd-160-Ugi-T':     'ugi_3cr',                 # 281630
    'Radd-203-Ugi-3CR-CA': 'ugi_3cr',                # 283090
    'Radd-204-Ugi-3CR-AA': 'ugi_3cr',                # 283110

    # ── GBB MCR ───────────────────────────────────────────────────────────
    'Radd-190-GBB':       'gbb',                      # 282730

    # ── MULTI-COMP: Double amidation ──────────────────────────────────────
    'Radd-022-AA-Ac-Boc-Ac':  'multi_amide_amide',   # 274552
    'Radd-049-DA-Ac-Boc-Ac':  'multi_amide_amide',   # 275592
    'Radd-138-Ac-tBu-Ac':     'multi_amide_amide',   # 280710

    # ── MULTI-COMP: Amide + arylation ─────────────────────────────────────
    'Radd-062-DA-Ar-Boc-Ac':  'multi_amide_arylation',   # 276450
    'Radd-010-DA-Ac-Boc-Ar':  'multi_amide_arylation',   # 274310

    # ── MULTI-COMP: Sulfonamide + amide ───────────────────────────────────
    'Radd-063-DA-SO2N-Boc-Ac': 'multi_sulfonamide_amide', # 276452

    # ── MULTI-COMP: Amide + alkylation ────────────────────────────────────
    'Radd-039-DA-Ac-Boc-Alk':  'multi_amide_alkylation',  # 275030
    'Radd-065-DA-Alk-Boc-Ac':  'multi_amide_alkylation',  # 276456

    # ── MULTI-COMP: Amide + click ─────────────────────────────────────────
    'Radd-036-NAc-Click':      'multi_amide_click',    # 274860
    'Radd-093-Ac-Click#':      'multi_amide_click',    # 279130
    'Radd-058-Click-Boc-Ac':   'multi_click_amide',    # 276010
    'Radd-115-NAc-Click-OH':   'multi_amide_click',    # 280190

    # ── MULTI-COMP: Alkylation + click ────────────────────────────────────
    'R051-Alk-Click':          'multi_alkylation_click',  # 274952

    # ── MULTI-COMP: Amide + Suzuki ────────────────────────────────────────
    'Radd-053-Ac-Suzuki':      'multi_amide_suzuki',   # 275730

    # ── MULTI-COMP: Sulfonamide + Suzuki ──────────────────────────────────
    'Radd-051-SO2N-Suzuki':    'multi_sulfonamide_suzuki',  # 275650

    # ── MULTI-COMP: Sulfonamide + arylation ───────────────────────────────
    'Radd-130-SO2N-Ar':        'multi_sulfonamide_arylation',  # 280530

    # ── MULTI-COMP: Red. amin. + various ──────────────────────────────────
    'R071-Rd-Acl':             'multi_reductive_amination_acyl_halide',  # 275550
    'R071b-Rd-Ac':             'multi_reductive_amination_amide',        # 282151
    'Radd-078-Ac-Boc-RAT':    'multi_reductive_amination_amide',        # 277970
    'R071a-Rd-Ar':             'multi_reductive_amination_arylation',    # 282150
    'R071c-Rd-Sac':            'multi_reductive_amination_sulfonamide',  # 282152
    'R071d-Rd-Ur':             'multi_reductive_amination_urea',         # 282153
    'R071g-Rd-Ur1':            'multi_reductive_amination_urea',         # 282156
    'R071h-Rd-Ur4':            'multi_reductive_amination_urea',         # 282157
    'R071e-Rd-Carbamate':      'multi_reductive_amination_carbamate',    # 282154
    'R071f-Rd-Rd':             'multi_reductive_amination_reductive_amination',  # 282155

    # ── MULTI-COMP: Arylation + Suzuki / Red. amin. + Suzuki ─────────────
    'Radd-158-RedAmin-Suzuki': 'multi_arylation_suzuki',    # 281510
    'Radd-159-Arylation-Suzuki': 'multi_arylation_suzuki',  # 281570
    'Radd-054-Alk-Suzuki':     'multi_arylation_suzuki',    # 275732

    # ── MULTI-COMP: Double arylation ──────────────────────────────────────
    'Radd-168-Ar-Ar':          'multi_arylation_arylation',  # 282030
    'ARYLATION2':              'multi_arylation_arylation',  # 270942
    'ARYLATION3':              'multi_arylation_arylation',  # 273170

    # ── MULTI-COMP: Anhydride + arylation ─────────────────────────────────
    'Radd-111-Anh-Ac-Ar':     'multi_anhydride_amide_arylation',  # 280150

    # ── MULTI-COMP: Ac-SAc ────────────────────────────────────────────────
    'Radd-113-Ac-SAc':         'multi_sulfonamide_amide',    # 280176

    # ── MULTI-COMP: Acylsulfamide + alkylation ────────────────────────────
    'Radd-179-SAc-Alk':        'multi_acylsulfamide_alkylation',  # 282530

    # ── MULTI-COMP: Triazole (lactim) + X ─────────────────────────────────
    'R194a-Lactim-Boc-Ac':     'multi_triazole_amide',       # 282890
    'R194b-Lactim-Boc-Ur1':    'multi_triazole_urea',        # 282970
    'R194c-Lactim-Boc-Alk':    'multi_triazole_alkylation',  # 282971
    'R194d-Lactim-Boc-Sac':    'multi_triazole_sulfonamide', # 282972

    # ── MULTI-COMP: Triazole (click) + heterocycle ────────────────────────
    'Radd-068-135Subst124Triazole': 'multi_click_amide',     # 276630

    # ── MULTI-COMP: Thiazole + amide / sulfonamide ────────────────────────
    'Radd-196-Thzl-Boc-Ac':   'multi_thiazole_amide',        # 282910
    'Radd-202-Thzl-Boc-SO2N':  'multi_thiazole_sulfonamide', # 283070

    # ── MULTI-COMP: Oxoisoindole ──────────────────────────────────────────
    'R198a-OXOISOINDOLE-Boc-Ac': 'multi_oxoisoindole_amide', # 283351

    # ── MULTI-COMP: Dihydrouracyl ─────────────────────────────────────────
    'R142-DIHYDROURACYL':      'dihydrouracil',               # 277130

    # ── MULTI-COMP: Purine ────────────────────────────────────────────────
    'R055':                    'purine_substitution',          # 275054

    # ── MISC: No-description reactions (generic fallback) ─────────────────
    'R003':                    'generic_3comp',     # 271242
    'R014':                    'generic_3comp',     # 273230
    'R015':                    'generic_3comp',     # 273232
}


# ═══════════════════════════════════════════════════════════════════════════
# Utility Functions
# ═══════════════════════════════════════════════════════════════════════════

def get_u_neighbour_info(synthon_smarts):
    """
    Parse a synthon SMARTS to find the atom directly bonded to [U].
    Returns (symbol, is_aromatic) or (None, None).
    """
    mol = Chem.MolFromSmiles(synthon_smarts)
    if mol is None:
        mol = Chem.MolFromSmarts(synthon_smarts)
    if mol is None:
        return None, None

    u_idx = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 92:  # Uranium = [U]
            u_idx = atom.GetIdx()
            break
    if u_idx is None:
        return None, None

    neighbours = mol.GetAtomWithIdx(u_idx).GetNeighbors()
    if not neighbours:
        return None, None

    nb = neighbours[0]
    return nb.GetSymbol(), nb.GetIsAromatic()


def classify_suzuki_role(bb_smiles):
    """For Suzuki, determine boronic acid vs aryl halide from BB SMILES."""
    mol = Chem.MolFromSmiles(bb_smiles)
    if mol is None:
        return 'unknown'
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 5:  # Boron
            return 'boronic_acid'
    halide_pat = Chem.MolFromSmarts('[c][Br,I,Cl]')
    if halide_pat and mol.HasSubstructMatch(halide_pat):
        return 'aryl_halide'
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in (17, 35, 53):  # Cl, Br, I
            return 'aryl_halide'
    return 'unknown'


# ── Protecting Group Detection ────────────────────────────────────────────

def get_reaction_pg_types(rxn_name):
    """
    Determine protecting group type(s) from the Enamine reaction name.

    Enamine uses naming prefixes to indicate a deprotection/cleavage step:
      Boc-   → tert-butyloxycarbonyl (Boc) removal from amine
      B-     → tert-butyl ester cleavage (acid deprotection)
      K-     → saponification (alkyl ester hydrolysis on acid)
      BBoc-  → both Boc + tert-butyl ester
      Me3Si  → trimethylsilyl (TMS) removal
      -Boc-  → Boc deprotection in multi-component reactions
      -tBu-  → tert-butyl ester in multi-component reactions

    Returns a list of PG type keys (from PROTECTING_GROUP_INFO).
    Order is significant for BBoc- (check Boc before tBu so that the
    amine BB matches Boc first, leaving tBu for the acid BB).
    """
    if not rxn_name:
        return []

    # BBoc- means both Boc (on amine) + tBu ester (on acid)
    if rxn_name.startswith('BBoc-'):
        return ['boc', 'tbu_ester']

    # Boc- prefix = Boc deprotection of amine
    if rxn_name.startswith('Boc-'):
        return ['boc']

    # B- prefix = tert-butyl ester cleavage
    if rxn_name.startswith('B-'):
        return ['tbu_ester']

    # K- prefix = saponification (base hydrolysis of alkyl ester)
    if rxn_name.startswith('K-'):
        return ['saponification']

    # R039-Me3Si- = TMS deprotection
    if 'Me3Si' in rxn_name:
        return ['tms']

    # Multi-component reactions with Boc in the name
    # e.g., Radd-022-AA-Ac-Boc-Ac, R194a-Lactim-Boc-Ac, etc.
    if '-Boc-' in rxn_name:
        return ['boc']

    # Multi-component with tBu in the name
    # e.g., Radd-138-Ac-tBu-Ac
    if '-tBu-' in rxn_name:
        return ['tbu_ester']

    return []


def detect_protecting_group(bb_smiles, pg_types):
    """
    Check if a building block contains any of the specified protecting groups.

    For each PG type (in order), check whether the BB SMILES contains the
    corresponding PG pattern via RDKit substructure matching.  Returns the
    PG SMARTS string on the first match, or '' if none detected.

    The order matters for BBoc- reactions: 'boc' is checked before 'tbu_ester'
    so the amine BB (with N-Boc) captures 'boc' and the acid BB (with
    C(=O)OC(C)(C)C but no N-attached version) captures 'tbu_ester'.
    """
    if not pg_types or not bb_smiles:
        return ''

    clean = str(bb_smiles).strip()
    if not clean or clean == 'nan':
        return ''

    # Strip salts — keep largest fragment
    if '.' in clean:
        clean = max(clean.split('.'), key=len)

    mol = Chem.MolFromSmiles(clean)
    if mol is None:
        return ''

    for pg_type in pg_types:
        pg_info = PROTECTING_GROUP_INFO.get(pg_type)
        if pg_info is None:
            continue
        pg_smarts = pg_info['pg_smarts']
        pat = Chem.MolFromSmarts(pg_smarts)
        if pat is None:
            continue
        if mol.HasSubstructMatch(pat):
            return pg_smarts

    return ''


def _classify_acyl_bb(synthon_smarts, bb_smiles):
    """
    For ACETYLATION reactions, classify a BB as either the acyl electrophile
    (acid_chloride / anhydride) or the amine nucleophile based on BB SMILES.
    """
    mol = Chem.MolFromSmiles(bb_smiles) if bb_smiles else None

    # ── acid chloride  R-C(=O)Cl ─────────────────────────────────────
    acyl_cl = Chem.MolFromSmarts('[CX3](=O)[Cl]')
    if mol and acyl_cl and mol.HasSubstructMatch(acyl_cl):
        rd = ROLE_DEFINITIONS['acid_chloride']
        return ('acid_chloride', rd['rxn_centre_smarts'], rd['label'])

    # ── anhydride  R-C(=O)O-C(=O)R ───────────────────────────────────
    anhydride = Chem.MolFromSmarts('[CX3](=O)[OX2][CX3](=O)')
    if mol and anhydride and mol.HasSubstructMatch(anhydride):
        rd = ROLE_DEFINITIONS['acid_chloride']          # same role def
        return ('acid_chloride', rd['rxn_centre_smarts'], rd['label'])

    # ── sulfonyl chloride  R-S(=O)(=O)Cl ─────────────────────────────
    sulfonyl_cl = Chem.MolFromSmarts('[SX4](=O)(=O)[Cl]')
    if mol and sulfonyl_cl and mol.HasSubstructMatch(sulfonyl_cl):
        rd = ROLE_DEFINITIONS['sulfonyl_chloride']
        return ('sulfonyl_chloride', rd['rxn_centre_smarts'], rd['label'])

    # ── amine (default nucleophile) ───────────────────────────────────
    # Also check synthon for `N([U])` indicator as a safety net
    rd = ROLE_DEFINITIONS['amine']
    return ('amine', rd['rxn_centre_smarts'], rd['label'])


def _classify_urea_isocyanate_bb(synthon_smarts, bb_smiles):
    """
    For urea_isocyanate reactions, classify a BB as either the isocyanate
    electrophile or the amine nucleophile based on BB SMILES.
    """
    mol = Chem.MolFromSmiles(bb_smiles) if bb_smiles else None

    # ── isocyanate  R-N=C=O ───────────────────────────────────────────
    iso_pat = Chem.MolFromSmarts('[NX2]=[CX2]=[OX1]')
    if mol and iso_pat and mol.HasSubstructMatch(iso_pat):
        rd = ROLE_DEFINITIONS['isocyanate']
        return ('isocyanate', rd['rxn_centre_smarts'], rd['label'])

    # ── carbamoyl chloride  Cl-C(=O)-N ───────────────────────────────
    carbamoyl_pat = Chem.MolFromSmarts('[Cl][CX3](=O)[NX3]')
    if mol and carbamoyl_pat and mol.HasSubstructMatch(carbamoyl_pat):
        rd = ROLE_DEFINITIONS['carbamoyl_chloride']
        return ('carbamoyl_chloride', rd['rxn_centre_smarts'], rd['label'])

    # ── default: amine ────────────────────────────────────────────────
    rd = ROLE_DEFINITIONS['amine']
    return ('amine', rd['rxn_centre_smarts'], rd['label'])


def _classify_urea_carbamate_bb(synthon_smarts, bb_smiles):
    """
    For urea_carbamate reactions, classify a BB as either the reactive
    carbamate/carbonate electrophile or the amine nucleophile.
    Checks for free amine first; if absent, checks for ester/carbamate.
    """
    mol = Chem.MolFromSmiles(bb_smiles) if bb_smiles else None
    if mol is None:
        rd = ROLE_DEFINITIONS['amine']
        return ('amine', rd['rxn_centre_smarts'], rd['label'])

    # ── Check for free amine first — if present, it's the nucleophile ─
    free_amine = Chem.MolFromSmarts('[NX3;H1,H2,H3;!$(NC=O);!$(NS(=O)=O)]')
    if free_amine and mol.HasSubstructMatch(free_amine):
        rd = ROLE_DEFINITIONS['amine']
        return ('amine', rd['rxn_centre_smarts'], rd['label'])

    # ── Check for aromatic NH (pyrazole, imidazole, etc.) ─────────────
    arom_nh = Chem.MolFromSmarts('[nH]')
    if arom_nh and mol.HasSubstructMatch(arom_nh):
        return ('amine', '[nH]',
                'Aromatic NH nucleophile (pyrazole, imidazole, etc.)')

    # ── reactive carbamate  R-O-C(=O)-NR or R-O-C(=O)-OR ─────────────
    carbamate_pat = Chem.MolFromSmarts('[OX2][CX3](=O)[NX3,OX2]')
    if carbamate_pat and mol.HasSubstructMatch(carbamate_pat):
        rd = ROLE_DEFINITIONS['reactive_carbamate']
        return ('reactive_carbamate', rd['rxn_centre_smarts'], rd['label'])

    # ── chloroformate  Cl-C(=O)-O ────────────────────────────────────
    chloroformate_pat = Chem.MolFromSmarts('[Cl][CX3](=O)[OX2]')
    if chloroformate_pat and mol.HasSubstructMatch(chloroformate_pat):
        rd = ROLE_DEFINITIONS['chloroformate']
        return ('chloroformate', rd['rxn_centre_smarts'], rd['label'])

    # ── isocyanate (sometimes in carbamate reactions) ────────────────
    iso_pat = Chem.MolFromSmarts('[NX2]=[CX2]=[OX1]')
    if iso_pat and mol.HasSubstructMatch(iso_pat):
        rd = ROLE_DEFINITIONS['isocyanate']
        return ('isocyanate', rd['rxn_centre_smarts'], rd['label'])

    # ── activated ester (no free amine → must be electrophile) ────────
    ester_pat = Chem.MolFromSmarts('[OX2;H0][CX3](=O)')
    if ester_pat and mol.HasSubstructMatch(ester_pat):
        rd = ROLE_DEFINITIONS['reactive_carbamate']
        return ('reactive_carbamate', rd['rxn_centre_smarts'], rd['label'])

    # ── CDI-activated amide (C=O bonded to aromatic N, e.g. imidazole) ─
    cdi_pat = Chem.MolFromSmarts('[n][CX3](=O)')
    if cdi_pat and mol.HasSubstructMatch(cdi_pat):
        rd = ROLE_DEFINITIONS['reactive_carbamate']
        return ('reactive_carbamate', rd['rxn_centre_smarts'], rd['label'])

    # ── default: amine ────────────────────────────────────────────────
    rd = ROLE_DEFINITIONS['amine']
    return ('amine', rd['rxn_centre_smarts'], rd['label'])


def match_synthon_indicators(synthon_smarts, indicators):
    """Check if any indicator substring is present in the synthon SMARTS."""
    for ind in indicators:
        if ind in synthon_smarts:
            return True
    return False


def classify_bb_role(synthon_smarts, bb_smiles, reaction_name):
    """
    Determine the role and reaction centre SMARTS for a building block.

    Returns (role_name, rxn_centre_smarts, rxn_centre_label)
    """
    rxn_class_name = ENAMINE_REACTION_TO_CLASS.get(reaction_name)

    # ── Fallback for unmapped reactions ────────────────────────────────
    if rxn_class_name is None:
        return _classify_by_synthon_fallback(synthon_smarts, bb_smiles)

    rxn_class = REACTION_CLASSES.get(rxn_class_name)
    if rxn_class is None:
        return _classify_by_synthon_fallback(synthon_smarts, bb_smiles)

    # ── Special: Suzuki — classify by BB SMILES ───────────────────────
    if rxn_class.get('classify_by_bb_smiles'):
        suzuki_role = classify_suzuki_role(bb_smiles)
        if suzuki_role == 'boronic_acid':
            rd = ROLE_DEFINITIONS['boronic_acid']
            return ('boronic_acid', rd['rxn_centre_smarts'], rd['label'])
        elif suzuki_role == 'aryl_halide':
            rd = ROLE_DEFINITIONS['aryl_halide']
            return ('aryl_halide', rd['rxn_centre_smarts'], rd['label'])
        else:
            return ('unknown', None, 'Unknown Suzuki role')

    # ── Special: ACETYLATION — classify acyl electrophile by BB SMILES ──
    if rxn_class.get('classify_by_bb_smiles_acyl'):
        return _classify_acyl_bb(synthon_smarts, bb_smiles)

    # ── Special: UREA from isocyanate — classify by BB SMILES ────────
    if rxn_class.get('classify_by_bb_smiles_urea_iso'):
        return _classify_urea_isocyanate_bb(synthon_smarts, bb_smiles)

    # ── Special: UREA from reactive carbamate — classify by BB SMILES ─
    if rxn_class.get('classify_by_bb_smiles_urea_carbamate'):
        return _classify_urea_carbamate_bb(synthon_smarts, bb_smiles)

    # ── General: try each role's synthon indicators ───────────────────
    for role_key in ['role_a', 'role_b', 'role_c', 'role_d']:
        role_name = rxn_class.get(role_key)
        if role_name is None:
            continue
        role_def = ROLE_DEFINITIONS.get(role_name)
        if role_def is None:
            continue
        if match_synthon_indicators(synthon_smarts, role_def['synthon_indicators']):
            # For generic roles, refine to a specific nucleophile/electrophile
            if role_name == 'generic_nucleophile':
                return _refine_nucleophile(synthon_smarts, bb_smiles)
            if role_name == 'generic_electrophile':
                return _refine_electrophile(synthon_smarts, bb_smiles)
            return (role_name, role_def['rxn_centre_smarts'], role_def['label'])

    # ── Still no match: use [U]-neighbour analysis as fallback ────────
    return _classify_by_synthon_fallback(synthon_smarts, bb_smiles)


def _refine_nucleophile(synthon_smarts, bb_smiles):
    """
    When a BB is identified as a generic nucleophile, refine to specific type
    based on the [U] neighbor atom in the synthon and BB structure.
    """
    symbol, _ = get_u_neighbour_info(synthon_smarts)
    if symbol is None:
        return ('generic_nucleophile',
                '[$([NX3;H2,H1]),$([NX3;H1]C=O),$([OX2H1]),$([SX2H1])]',
                'Nucleophile (N/O/S)')
    sym = symbol.upper()
    if sym == 'N':
        # Check if the BB only has amide/imide NH (not free amine)
        mol = Chem.MolFromSmiles(bb_smiles)
        if mol:
            free_amine = Chem.MolFromSmarts('[NX3;H1,H2,H3;!$(NC=O);!$(NS(=O)=O)]')
            arom_nh = Chem.MolFromSmarts('[nH]')     # aromatic NH (pyrazole etc.)
            het_nh = Chem.MolFromSmarts('[NX3;H1]')  # broader: any NH
            if free_amine and mol.HasSubstructMatch(free_amine):
                return ('amine', '[NX3;H1,H2,H3;!$(NC=O);!$(NS(=O)=O)]',
                        'Primary or secondary amine (incl. ammonia)')
            elif arom_nh and mol.HasSubstructMatch(arom_nh):
                return ('amine', '[nH]',
                        'Aromatic NH nucleophile (pyrazole, imidazole, etc.)')
            elif het_nh and mol.HasSubstructMatch(het_nh):
                return ('amine', '[NX3;H1]',
                        'NH nucleophile (heterocyclic/amide NH)')
        return ('amine', '[NX3;H2,H1]',
                'Primary or secondary amine')
    elif sym == 'O':
        return ('alcohol', '[OX2H1]', 'Alcohol or phenol (-OH)')
    elif sym == 'S':
        if '[U]S(=O)' in synthon_smarts or 'S(=O)' in synthon_smarts:
            return ('thiol', '[SX2,SX3]',
                    'Thiol or sulfinate (S-nucleophile)')
        return ('thiol', '[SX2H1]', 'Thiol (-SH)')
    return ('generic_nucleophile',
            '[$([NX3;H2,H1]),$([NX3;H1]C=O),$([OX2H1]),$([SX2H1])]',
            'Nucleophile (N/O/S)')


def _refine_electrophile(synthon_smarts, bb_smiles):
    """
    When a BB is identified as a generic electrophile, refine to specific type.
    """
    symbol, is_arom = get_u_neighbour_info(synthon_smarts)
    if symbol is None:
        return ('generic_electrophile',
                '[CX4][Cl,Br,I,$([OX2]S(=O)(=O))]',
                'Electrophile (alkyl halide/sulfonate)')
    sym = symbol.upper()
    if sym == 'C':
        if is_arom:
            return ('aryl_halide', '[c][F,Cl,Br,I]', 'Aryl halide (Ar-X)')
        if '(=O)' in synthon_smarts:
            return ('carboxylic_acid', '[CX3](=O)[OX2H1]',
                    'Carboxylic acid (-COOH)')
        return ('alkyl_halide',
                '[CX4][Cl,Br,I,$([OX2]S(=O)(=O))]',
                'Alkyl halide/sulfonate (-CH2X / -CH2OMs/OTs)')
    elif sym == 'S':
        return ('sulfonyl_halide',
                '[SX4](=O)(=O)[Cl,F,$([OX2][c,C])]',
                'Sulfonyl halide or activated ester')
    return ('generic_electrophile',
            '[CX4][Cl,Br,I,$([OX2]S(=O)(=O))]',
            'Electrophile (alkyl halide/sulfonate)')


def _classify_by_synthon_fallback(synthon_smarts, bb_smiles):
    """
    Fallback classifier based on the atom bonded to [U] in the synthon.
    """
    symbol, is_arom = get_u_neighbour_info(synthon_smarts)
    if symbol is None:
        return ('unknown', None, 'Unclassified')

    sym_upper = symbol.upper()

    if sym_upper == 'N':
        # Check for aromatic NH in BB
        mol = Chem.MolFromSmiles(bb_smiles) if bb_smiles else None
        if mol:
            free_amine = Chem.MolFromSmarts('[NX3;H1,H2,H3;!$(NC=O);!$(NS(=O)=O)]')
            arom_nh = Chem.MolFromSmarts('[nH]')
            if free_amine and mol.HasSubstructMatch(free_amine):
                return ('amine', '[NX3;H1,H2,H3;!$(NC=O);!$(NS(=O)=O)]',
                        'Primary or secondary amine (incl. ammonia)')
            elif arom_nh and mol.HasSubstructMatch(arom_nh):
                return ('amine', '[nH]',
                        'Aromatic NH nucleophile (pyrazole, imidazole, etc.)')
        return ('amine', '[NX3;H1,H2,H3]', 'Primary or secondary amine (incl. ammonia)')
    elif sym_upper == 'O':
        return ('alcohol', '[OX2H1]', 'Alcohol or phenol (-OH)')
    elif sym_upper == 'S':
        # Distinguish sulfonyl electrophile from sulfinate nucleophile
        # Check both indicator directions
        if 'S(=O)(=O)[U]' in synthon_smarts or '[U]S(=O)(=O)' in synthon_smarts:
            return ('sulfonyl_halide',
                    '[SX4](=O)(=O)[Cl,F,$([OX2][c,C])]',
                    'Sulfonyl halide or activated ester')
        if '[U]S(=O)' in synthon_smarts or 'S(=O)[U]' in synthon_smarts:
            return ('thiol', '[SX2,SX3]',
                    'Thiol or sulfinate (S-nucleophile)')
        # Check BB structure for sulfonyl
        mol = Chem.MolFromSmiles(bb_smiles) if bb_smiles else None
        sulfonyl_pat = Chem.MolFromSmarts('[SX4](=O)(=O)[Cl,F]')
        if mol and sulfonyl_pat and mol.HasSubstructMatch(sulfonyl_pat):
            return ('sulfonyl_halide',
                    '[SX4](=O)(=O)[Cl,F,$([OX2][c,C])]',
                    'Sulfonyl halide or activated ester')
        return ('thiol', '[SX2H1]', 'Thiol (-SH)')
    elif sym_upper == 'C':
        if is_arom:
            role = classify_suzuki_role(bb_smiles)
            if role == 'boronic_acid':
                return ('boronic_acid', '[#6]B([OX2])[OX2]',
                        'Boronic acid or pinacol ester')
            return ('aryl_halide', '[c][F,Cl,Br,I]', 'Aryl halide (Ar-X)')
        if '(=O)' in synthon_smarts:
            # Distinguish carboxylic acid from acid chloride / isocyanate
            mol = Chem.MolFromSmiles(bb_smiles) if bb_smiles else None
            if mol:
                acyl_cl = Chem.MolFromSmarts('[CX3](=O)[Cl]')
                iso_pat = Chem.MolFromSmarts('[NX2]=[CX2]=[OX1]')
                cooh = Chem.MolFromSmarts('[$([CX3](=O)[OX2H1]),$([CX3](=O)[O-])]')
                if acyl_cl and mol.HasSubstructMatch(acyl_cl):
                    rd = ROLE_DEFINITIONS['acid_chloride']
                    return ('acid_chloride', rd['rxn_centre_smarts'], rd['label'])
                elif iso_pat and mol.HasSubstructMatch(iso_pat):
                    rd = ROLE_DEFINITIONS['isocyanate']
                    return ('isocyanate', rd['rxn_centre_smarts'], rd['label'])
                elif cooh and mol.HasSubstructMatch(cooh):
                    rd = ROLE_DEFINITIONS['carboxylic_acid']
                    return ('carboxylic_acid', rd['rxn_centre_smarts'], rd['label'])
            return ('carboxylic_acid',
                    '[$([CX3](=O)[OX2H1]),$([CX3](=O)[O-])]',
                    'Carboxylic acid or carboxylate')
        return ('alkyl_halide', '[CX4][Cl,Br,I,$([OX2]S(=O)(=O))]',
                'Alkyl halide/sulfonate (-CH2X / -CH2OMs/OTs)')
    elif sym_upper == 'B':
        return ('boronic_acid', '[#6]B([OX2])[OX2]',
                'Boronic acid or pinacol ester')

    return ('unknown', None, f'Unclassified (atom: {symbol})')


def validate_smarts_match(bb_smiles, rxn_centre_smarts):
    """Check if reaction centre SMARTS matches the BB. Returns True/False/None."""
    if rxn_centre_smarts is None:
        return None
    clean = bb_smiles
    if '.' in clean:
        clean = max(clean.split('.'), key=len)
    mol = Chem.MolFromSmiles(clean)
    if mol is None:
        return None
    pat = Chem.MolFromSmarts(rxn_centre_smarts)
    if pat is None:
        return None
    return mol.HasSubstructMatch(pat)


def validate_synthon_in_target(synthon_smarts, target_smiles):
    """
    Check that the BB scaffold (synthon with [U] → wildcard) is present
    as a substructure of the target product molecule.

    This validates that the building block actually contributed its
    non-reactive scaffold to the final product.

    Returns True/False/None.
    """
    if not synthon_smarts or not target_smiles:
        return None
    synthon_str = str(synthon_smarts).strip()
    target_str = str(target_smiles).strip()
    if synthon_str == '' or target_str == '' or synthon_str == 'nan':
        return None

    # Replace [U] and [Np] attachment points with wildcard [*]
    query_str = synthon_str.replace('[U]', '[*]').replace('[Np]', '[*]')

    target_mol = Chem.MolFromSmiles(target_str)
    if target_mol is None:
        return None

    # Try parsing as SMARTS first (more permissive), then as SMILES
    query_mol = Chem.MolFromSmarts(query_str)
    if query_mol is None:
        query_mol = Chem.MolFromSmiles(query_str)
        if query_mol is None:
            return None

    return target_mol.HasSubstructMatch(query_mol)


def validate_product_bond(reaction_class, target_smiles):
    """
    Check that the expected product bond pattern (from PRODUCT_BOND_SMARTS)
    is present in the target molecule.

    This validates that the reaction type classification is consistent
    with the actual structure of the target — e.g. an amidation target
    must contain an amide bond [#7]C(=O).

    Returns True/False/None.
    """
    if not reaction_class or reaction_class == 'unknown':
        return None
    smarts_str = PRODUCT_BOND_SMARTS.get(reaction_class)
    if smarts_str is None:
        return None

    target_str = str(target_smiles).strip()
    if not target_str or target_str == 'nan':
        return None

    target_mol = Chem.MolFromSmiles(target_str)
    if target_mol is None:
        return None
    pat = Chem.MolFromSmarts(smarts_str)
    if pat is None:
        return None
    return target_mol.HasSubstructMatch(pat)


# ═══════════════════════════════════════════════════════════════════════════
# Main Extraction
# ═══════════════════════════════════════════════════════════════════════════

def extract_building_blocks(input_csv):
    """
    Parse the Enamine REAL output CSV and extract building block info
    with reaction centre SMARTS.  Returns a list of dicts.
    """
    df = pd.read_csv(input_csv)
    found = df[df['status'] == 'FOUND'].copy()

    if found.empty:
        print("WARNING: No FOUND rows in input file.")
        return []

    print(f"Input: {len(df)} total rows, {len(found)} FOUND rows")
    print(f"Unique target compounds: {found['target_smiles'].nunique()}")
    print(f"Reaction types: {sorted(found['reaction_name'].unique())}")

    bb_records = []
    unmapped_rxns = set()

    for _, row in found.iterrows():
        target = row['target_smiles']
        rxn_name = row['reaction_name']
        rxn_desc = row.get('reaction_description', '')
        rxn_id = row.get('reaction_id', '')
        exact_match = row.get('exact_match', '')
        found_smiles = row.get('found_smiles', '')
        method_no = row.get('method_no', '')

        if rxn_name not in ENAMINE_REACTION_TO_CLASS:
            unmapped_rxns.add(rxn_name)

        for i in range(1, 5):
            bb_code = row.get(f'reactant_{i}_code', '')
            bb_smiles = row.get(f'reactant_{i}_smiles', '')
            synthon_smarts = row.get(f'reactant_synthon_smarts_{i}', '')

            if not bb_code or not bb_smiles or pd.isna(bb_code) or pd.isna(bb_smiles):
                continue
            if str(bb_code).strip() == '' or str(bb_smiles).strip() == '':
                continue

            bb_code = str(bb_code).strip()
            bb_smiles = str(bb_smiles).strip()
            synthon_smarts = str(synthon_smarts).strip() if pd.notna(synthon_smarts) else ''

            role, rxn_centre, rxn_label = classify_bb_role(
                synthon_smarts, bb_smiles, rxn_name
            )

            # Canonical SMILES (strip salts)
            clean_smiles = bb_smiles
            if '.' in clean_smiles:
                clean_smiles = max(clean_smiles.split('.'), key=len)
            mol = Chem.MolFromSmiles(clean_smiles)
            if mol:
                clean_smiles = Chem.MolToSmiles(mol)

            match_ok = validate_smarts_match(bb_smiles, rxn_centre)

            rxn_class_name = ENAMINE_REACTION_TO_CLASS.get(rxn_name, 'unknown')
            rxn_class_desc = ''
            if rxn_class_name in REACTION_CLASSES:
                rxn_class_desc = REACTION_CLASSES[rxn_class_name]['description']

            # Validate synthon scaffold is present in the target product
            # Use found_smiles (the actual product) if available, else target
            product_smi = found_smiles if (found_smiles and str(found_smiles).strip()
                                           and str(found_smiles) != 'nan') else target
            synthon_in_target = validate_synthon_in_target(
                synthon_smarts, product_smi
            )

            # Validate expected product bond is in the target
            product_bond_ok = validate_product_bond(rxn_class_name, product_smi)

            # Detect protecting group (if this reaction involves PG removal)
            pg_types = get_reaction_pg_types(rxn_name)
            pg_smarts = detect_protecting_group(bb_smiles, pg_types)

            record = {
                'target_smiles': target,
                'found_smiles': found_smiles,
                'exact_match': exact_match,
                'reaction_name': rxn_name,
                'reaction_id': rxn_id,
                'reaction_description': rxn_desc,
                'reaction_class': rxn_class_name,
                'reaction_class_description': rxn_class_desc,
                'method_no': method_no,
                'reactant_number': i,
                'bb_code': bb_code,
                'bb_smiles_raw': bb_smiles,
                'bb_smiles_clean': clean_smiles,
                'synthon_smarts': synthon_smarts,
                'bb_role': role,
                'rxn_centre_smarts': rxn_centre if rxn_centre else '',
                'rxn_centre_label': rxn_label,
                'smarts_validated': match_ok if match_ok is not None else '',
                'synthon_in_target': synthon_in_target if synthon_in_target is not None else '',
                'product_bond_in_target': product_bond_ok if product_bond_ok is not None else '',
                'protecting_group_SMARTS': pg_smarts,
            }
            bb_records.append(record)

    if unmapped_rxns:
        print(f"\nWARNING: {len(unmapped_rxns)} unmapped reaction name(s) "
              f"(fallback classification used):")
        for r in sorted(unmapped_rxns):
            print(f"  - {r}")

    return bb_records


def deduplicate_building_blocks(records):
    """
    Create deduplicated summary of unique building blocks.
    Returns (full_df, bb_summary_df).
    """
    full_df = pd.DataFrame(records)
    if full_df.empty:
        return full_df, pd.DataFrame()

    dedup = full_df.drop_duplicates(
        subset=['bb_code', 'reaction_name', 'bb_role']
    ).copy()

    bb_summary = dedup[[
        'bb_code', 'bb_smiles_clean', 'reaction_name',
        'reaction_description', 'reaction_class', 'reaction_class_description',
        'bb_role', 'rxn_centre_smarts', 'rxn_centre_label',
        'synthon_smarts', 'smarts_validated',
        'synthon_in_target', 'product_bond_in_target',
        'protecting_group_SMARTS',
    ]].copy()

    bb_summary = bb_summary.sort_values(
        ['reaction_name', 'bb_role', 'bb_code']
    ).reset_index(drop=True)

    return full_df, bb_summary


def print_summary(full_df, bb_summary):
    """Print summary statistics."""
    print("\n" + "=" * 70)
    print("EXTRACTION SUMMARY")
    print("=" * 70)

    if full_df.empty:
        print("No building blocks extracted.")
        return

    print(f"Total BB-target pairs    : {len(full_df)}")
    print(f"Unique building blocks   : {bb_summary['bb_code'].nunique()}")
    print(f"Unique BB+reaction combos: {len(bb_summary)}")

    print(f"\n--- By Reaction Class ---")
    class_counts = bb_summary.groupby('reaction_class')['bb_code'].nunique()
    for cls, count in sorted(class_counts.items()):
        desc = ''
        if cls in REACTION_CLASSES:
            desc = f"  ({REACTION_CLASSES[cls]['description']})"
        print(f"  {cls:40s}: {count:4d} unique BBs{desc}")

    print(f"\n--- By Reaction Name ---")
    rxn_counts = bb_summary.groupby('reaction_name')['bb_code'].nunique()
    for rxn, count in sorted(rxn_counts.items()):
        print(f"  {rxn:30s}: {count} unique BBs")

    print(f"\n--- By BB Role ---")
    role_counts = bb_summary['bb_role'].value_counts()
    for role, count in role_counts.items():
        print(f"  {role:30s}: {count}")

    print(f"\n--- Validation 1: Reaction Centre SMARTS in BB ---")
    print(f"    (Does the BB contain the functional group we assigned?)")
    validated = bb_summary['smarts_validated']
    n_pass = (validated == True).sum()
    n_fail = (validated == False).sum()
    n_na = len(validated) - n_pass - n_fail
    print(f"  Passed  : {n_pass}")
    print(f"  Failed  : {n_fail}")
    print(f"  N/A     : {n_na}")

    if n_fail > 0:
        print(f"\n  Failed:")
        failures = bb_summary[validated == False]
        for _, row in failures.head(20).iterrows():
            print(f"    {row['bb_code']} | {row['bb_smiles_clean'][:50]}")
            print(f"      Role: {row['bb_role']} | SMARTS: {row['rxn_centre_smarts']} | "
                  f"Rxn: {row['reaction_name']}")

    print(f"\n--- Validation 2: BB Synthon Scaffold in Target Product ---")
    print(f"    (Does the target contain the BB scaffold after bond formation?)")
    synthon_val = full_df['synthon_in_target']
    s_pass = (synthon_val == True).sum()
    s_fail = (synthon_val == False).sum()
    s_na = len(synthon_val) - s_pass - s_fail
    print(f"  Passed  : {s_pass}")
    print(f"  Failed  : {s_fail}")
    print(f"  N/A     : {s_na}")

    if s_fail > 0:
        print(f"\n  Failed (showing from full records):")
        failures = full_df[synthon_val == False]
        for _, row in failures.head(20).iterrows():
            prod = row.get('found_smiles', row.get('target_smiles', ''))
            print(f"    {row['bb_code']} | synthon: {str(row['synthon_smarts'])[:60]}")
            print(f"      Target : {str(prod)[:70]}")
            print(f"      Rxn: {row['reaction_name']} | Role: {row['bb_role']}")

    print(f"\n--- Validation 3: Product Bond in Target ---")
    print(f"    (Does the target contain the expected reaction bond, e.g. amide?)")
    # Product bond validation is per-target-row rather than per-BB-deduped,
    # so use the full_df and deduplicate by target+reaction
    bond_rows = full_df.drop_duplicates(
        subset=['found_smiles', 'reaction_class']
    )
    bond_val = bond_rows['product_bond_in_target']
    b_pass = (bond_val == True).sum()
    b_fail = (bond_val == False).sum()
    b_na = len(bond_val) - b_pass - b_fail
    print(f"  Passed  : {b_pass}")
    print(f"  Failed  : {b_fail}")
    print(f"  N/A     : {b_na}")

    if b_fail > 0:
        print(f"\n  Failed:")
        failures = bond_rows[bond_val == False]
        for _, row in failures.head(20).iterrows():
            prod = row.get('found_smiles', row.get('target_smiles', ''))
            rxn_cls = row['reaction_class']
            expected = PRODUCT_BOND_SMARTS.get(rxn_cls, 'N/A')
            print(f"    Class: {rxn_cls} | Expected: {expected}")
            print(f"      Target: {str(prod)[:70]}")

    # Overall summary
    all_pass = (n_fail == 0 and s_fail == 0 and b_fail == 0)
    print(f"\n{'='*70}")
    if all_pass:
        print(f"ALL VALIDATIONS PASSED")
    else:
        total_fail = n_fail + s_fail + b_fail
        print(f"VALIDATION FAILURES: {total_fail} total "
              f"(BB-SMARTS: {n_fail}, synthon-in-target: {s_fail}, "
              f"product-bond: {b_fail})")
    print(f"{'='*70}")

    # ── Protecting Group Summary ──────────────────────────────────────────
    print(f"\n--- Protecting Group Detection ---")
    pg_col = 'protecting-group-SMARTS'
    if pg_col in bb_summary.columns:
        has_pg = bb_summary[pg_col].astype(str).apply(lambda x: x.strip() != '' and x != 'nan')
        n_with_pg = has_pg.sum()
        n_rxns_with_pg = bb_summary.loc[has_pg, 'reaction_name'].nunique() if n_with_pg > 0 else 0

        print(f"  BBs with protecting group   : {n_with_pg} / {len(bb_summary)}")
        print(f"  Reactions involving PG step  : {n_rxns_with_pg}")

        if n_with_pg > 0:
            # Group by PG SMARTS pattern
            pg_counts = bb_summary.loc[has_pg, pg_col].value_counts()
            print(f"\n  By protecting group type:")
            for pg_smarts, count in pg_counts.items():
                # Look up the label
                pg_label = pg_smarts
                for pg_info in PROTECTING_GROUP_INFO.values():
                    if pg_info['pg_smarts'] == pg_smarts:
                        pg_label = pg_info['label']
                        break
                print(f"    {pg_label:40s}: {count:4d} BBs")

            print(f"\n  Sample BBs with protecting groups:")
            sample_pg = bb_summary[has_pg].head(10)
            for _, row in sample_pg.iterrows():
                pg_label = row[pg_col]
                for pg_info in PROTECTING_GROUP_INFO.values():
                    if pg_info['pg_smarts'] == row[pg_col]:
                        pg_label = pg_info['label']
                        break
                print(f"    {row['bb_code']:20s} | {pg_label:30s} | {row['reaction_name']}")
    else:
        print(f"  (column not found in output)")

    print(f"\n--- Sample Output (first 20) ---")
    cols = ['bb_code', 'bb_smiles_clean', 'reaction_name', 'bb_role',
            'rxn_centre_smarts', 'rxn_centre_label']
    pd.set_option('display.max_colwidth', 50)
    print(bb_summary[cols].head(20).to_string(index=False))


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Extract building blocks and reaction centre SMARTS '
                    'from Enamine REAL Tools output CSV.'
    )
    parser.add_argument(
        'input_csv',
        help='Path to the Enamine REAL synthesis output CSV'
    )
    parser.add_argument(
        'output_csv', nargs='?', default=None,
        help='Output CSV path (default: <input>_building-blocks.csv)'
    )
    parser.add_argument(
        '--full', action='store_true',
        help='Also output the full (non-deduplicated) records'
    )
    args = parser.parse_args()

    input_csv = args.input_csv
    if not os.path.exists(input_csv):
        print(f"ERROR: Input file not found: {input_csv}")
        sys.exit(1)

    if args.output_csv:
        output_csv = args.output_csv
    else:
        base, ext = os.path.splitext(input_csv)
        output_csv = f"{base}_building-blocks{ext}"

    print(f"Input  : {os.path.abspath(input_csv)}")
    print(f"Output : {os.path.abspath(output_csv)}")
    print(f"Mapped reactions: {len(ENAMINE_REACTION_TO_CLASS)} "
          f"(of 190 in Enamine catalogue)")
    print()

    records = extract_building_blocks(input_csv)
    full_df, bb_summary = deduplicate_building_blocks(records)
    print_summary(full_df, bb_summary)

    bb_summary.to_csv(output_csv, index=False)
    print(f"\nDeduplicated BB list: {os.path.abspath(output_csv)}")
    print(f"  {len(bb_summary)} rows")

    if args.full:
        full_path = output_csv.replace('_building-blocks', '_building-blocks-full')
        if full_path == output_csv:
            full_path = output_csv.replace('.csv', '_full.csv')
        full_df.to_csv(full_path, index=False)
        print(f"Full records: {os.path.abspath(full_path)}")
        print(f"  {len(full_df)} rows")


if __name__ == '__main__':
    main()
