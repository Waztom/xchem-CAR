"""
Extract Reaction SMIRKS from Enamine REAL Tools Output for Syndirella
=====================================================================
Reads an Enamine REAL Tools synthesis CSV and adds columns containing:
  - reaction_smirks          : SMIRKS for the core transformation
  - syndirella_reaction_name : corresponding Syndirella reaction name
  - reaction_class           : internal classification label

The SMIRKS are written so they are directly usable as Syndirella
(https://github.com/xchem/syndirella) reaction inputs, matching the
format and atom-map conventions of Syndirella's RXN_SMIRKS_CONSTANTS.json.

Multi-step Enamine reactions (Boc-, B-, K- prefixed) are mapped to the
core bond-forming step because Syndirella handles deprotection
internally via its own constants.

Usage:
    python extract-enamine-smirks.py <input_csv> [output_csv]

    If output_csv is not given, it defaults to <input_basename>_smirks.csv
"""

import sys
import os
import argparse

import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.logger().setLevel(RDLogger.ERROR)


# ═══════════════════════════════════════════════════════════════════════════
# Syndirella-Compatible Reaction SMIRKS
# ═══════════════════════════════════════════════════════════════════════════
# Each key is a reaction class (shared with extract-enamine-building-blocks.py).
# The SMIRKS are derived from Syndirella's RXN_SMIRKS_CONSTANTS.json,
# using the parent reaction patterns where possible.

REACTION_CLASS_SMIRKS = {
    # ── Amide bond formers ─────────────────────────────────────────────
    'amidation': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'acetylation': {
        'smirks': '[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3].[#6X3;!$(C-N):1](=[OX1])[#17,#9,#35]>>[#6:1](=[#8:2])-[#7:3]',
        'syndirella_name': 'Amide_Schotten-Baumann',
    },
    'acylation_oc': {
        'smirks': '[Cl,F,Br,I]-[C;H0;D3;+0:1](=[O,S;D1;H0:2])-[*:3].[OH;D1;+0:4]-[*:5]>>[O,S;D1;H0:2]=[C;H0;D3;+0:1](-[*:3])-[O;H0;D2;+0:4]-[*:5]',
        'syndirella_name': 'Schotten-Baumann_to_ester',
    },
    'acrylamide': {
        'smirks': '[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3].[#6X3;!$(C-N):1](=[OX1])[#17,#9,#35]>>[#6:1](=[#8:2])-[#7:3]',
        'syndirella_name': 'Amide_Schotten-Baumann',
    },
    'anhydride': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'acylsulfamide': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[#7;H1:3]S(=O)(=O)>>[#6:1](=[#8:2])-[#7:3]',
        'syndirella_name': 'Amidation',
    },

    # ── Sulfonamide ───────────────────────────────────────────────────────
    'sulfonamide': {
        'smirks': '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
        'syndirella_name': 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)',
    },
    'sulfonamide_phenol': {
        'smirks': '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[OX2H1:5][c:6]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[O:5][c:6]',
        'syndirella_name': 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)',
    },
    'vinyl_sulfonamide': {
        'smirks': '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
        'syndirella_name': 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)',
    },

    # ── Reductive amination ───────────────────────────────────────────────
    'reductive_amination': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'reductive_amination_ketone': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },

    # ── Suzuki coupling ──────────────────────────────────────────────────
    'suzuki': {
        'smirks': '[#6:1]-[Cl,Br,I].[#6:2]-[B](-[O])(-[O])>>[#6:1]-[#6:2]',
        'syndirella_name': 'Suzuki_coupling',
    },

    # ── N,O,S-Arylation ─────────────────────────────────────────────────
    'arylation': {
        'smirks': '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',
        'syndirella_name': 'N-nucleophilic_aromatic_substitution',
    },

    # ── Alkylation of amines ──────────────────────────────────────────────
    'alkylation_amine': {
        'smirks': '[#7&X3;H2,H1;+0:5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[#7&X3;+0:5]-[C&X4;+0:1]',
        'syndirella_name': 'Nucleophilic_substitution_with_amine',
    },

    # ── Alkylation of N/O/S nucleophiles ──────────────────────────────────
    'alkylation_nos': {
        'smirks': '[#7&X3;H2,H1;+0:5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[#7&X3;+0:5]-[C&X4;+0:1]',
        'syndirella_name': 'Nucleophilic_substitution_with_amine',
    },

    # ── Alkylation of acid (ester) ────────────────────────────────────────
    'alkylation_acid': {
        'smirks': '[Cl,Br]-[#6;H2,H1,H0;+0:1].[OH;+0:2]>>[O&X2;H0;D2;+0:2]-[#6&X4;D2;+0:1]',
        'syndirella_name': 'Nucleophilic_substitution',
    },

    # ── S-Alkylation (thiols) ─────────────────────────────────────────────
    'alkylation_thiol': {
        'smirks': '[Cl,Br,I]-[C;H2,H1;+0:1].[SH;+0:2]>>[S&X2;H0;D2;+0:2]-[C&X4;D2;+0:1]',
        'syndirella_name': 'Nucleophilic_substitution_with_thiol',
    },

    # ── Epoxide ring-opening ──────────────────────────────────────────────
    'epoxide_amine': {
        'smirks': '[#6:1]-[CH;D3;+0:2]1-[CH2;D2;+0:3]-[O;H0;D2;+0:4]-1.[N;H2,H1;+0:6]-[#6:7]>>[#6:1]-[C;+0:2](-[OH;D1;+0:4])-[C&X4;+0:3]-[N&X3;+0:6]-[#6:7]',
        'syndirella_name': 'Epoxide_+_amine_coupling',
    },
    'epoxide_phenol': {
        'smirks': '[#6:1]-[CH;D3;+0:2]1-[CH2;D2;+0:3]-[O;H0;D2;+0:4]-1.[OH;D1;+0:6]-[c:7]>>[#6:1]-[C;+0:2](-[OH;D1;+0:4])-[C&X4;+0:3]-[O;+0:6]-[c:7]',
        'syndirella_name': 'Epoxide_+_amine_coupling',
    },

    # ── Urea formation ───────────────────────────────────────────────────
    'urea': {
        'smirks': '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[N&X3;H2,H1;!$(NC=O);!$(NS):4]>>[#7&X3:3]-[#6](=[#8])-[#7&X3:4]',
        'syndirella_name': 'Formation_of_urea_from_two_amines',
    },
    'urea_isocyanate': {
        'smirks': '[*:1]-[N;H0;D2;+0:2]=[C;H0;D2;+0:3]=[O;H0;D1;+0:4].[#6;+0:5]-[N;H2;D1;+0:6]>>[*:1]-[N;H1;D2;+0:2]-[C;H0;D3;+0:3](=[O;H0;D1;+0:4])-[N;H1;D2;+0:6]-[#6;+0:5]',
        'syndirella_name': 'Urea_synthesis_via_isocyanate_and_primary_amine',
    },
    'urea_carbamate': {
        'smirks': '[#7:1]-[C](=[#8X1])[nX3]1[c&H1X3][c&H1X3][nX2][c&H1X3]1.[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):2]>>[#7:1]-[C](=[#8X1])-[#7:2]',
        'syndirella_name': 'Formation_of_urea_by_displacement_of_imidazole_by_amine',
    },

    # ── Carbamate ─────────────────────────────────────────────────────────
    'carbamate_cdi': {
        'smirks': '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[#6:1]-[OH;D1;+0]>>[#7&X3:3]-[#6](=[#8])-[O]-[#6:1]',
        'syndirella_name': 'Formation_of_urea_from_two_amines',
    },
    'carbamate_chloroformate': {
        'smirks': '[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3].[Cl][CX3:1](=[OX1:2])[OX2:4]>>[#7:3]-[CX3:1](=[OX1:2])[OX2:4]',
        'syndirella_name': 'Amide_Schotten-Baumann',
    },

    # ── Thiourea ──────────────────────────────────────────────────────────
    'thiourea': {
        'smirks': '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[NX2:1]=[CX2:2]=[SX1:4]>>[#7&X3:3]-[C:2](=[S:4])-[N:1]',
        'syndirella_name': 'Thiourea_formation',
    },

    # ── Oxamide ───────────────────────────────────────────────────────────
    'oxamide': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },

    # ── Hydrazone condensation ────────────────────────────────────────────
    'hydrazone': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]-[#7]>>[#6:2](=[#7:3])(-[#6:1])',
        'syndirella_name': 'Hydrazone_condensation',
    },

    # ── Knoevenagel condensation ──────────────────────────────────────────
    'knoevenagel': {
        'smirks': '[CX4H2:1]([CX3]=O).[CX3H1:2](=O)>>[C:1]=[C:2]',
        'syndirella_name': 'Knoevenagel_condensation',
    },

    # ── Mannich reaction ──────────────────────────────────────────────────
    'mannich': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },

    # ── Click / Triazole ──────────────────────────────────────────────────
    'click_triazole': {
        'smirks': '[NX1:1]~[NX2:2]~[NX2:3]-[#6:4].[CX2:5]#[CX2:6]-[#6:7]>>[#6:4]-[n:3]1:[n:2]:[n:1]:[c:5](-[#6:7]):[c:6]:1',
        'syndirella_name': 'CuAAC_1,2,3-triazole',
    },
    'click_from_amine': {
        'smirks': '[#7;H2;D1:1]-[#6:4].[CX2:5]#[CX2:6]-[#6:7]>>[#6:4]-[n:1]1:[n]:[n]:[c:5](-[#6:7]):[c:6]:1',
        'syndirella_name': 'CuAAC_1,2,3-triazole',
    },

    # ── Triazole from lactim ──────────────────────────────────────────────
    'triazole_lactim': {
        'smirks': '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
        'syndirella_name': 'Triazole_from_lactim',
    },

    # ── Oxadiazole ────────────────────────────────────────────────────────
    'oxadiazole_amidoxime': {
        'smirks': '[NX3:1]C(=[NX2:2])[OX2H1:3].[#6:4](=[#8])-[#8;H1]>>[c:4]1:[n:2]:[o:3]:[c]:[n:1]:1',
        'syndirella_name': '1,2,4-Oxadiazole_from_amidoxime',
    },
    'oxadiazole_nitrile': {
        'smirks': '[CX2:1]#[NX1:2].[#6:3](=[#8])-[#8;H1]>>[c:1]1:[n:2]:[o]:[c:3]:[n]:1',
        'syndirella_name': '1,2,4-Oxadiazole_from_nitrile',
    },
    'oxadiazole_cyclization': {
        'smirks': '[CX2:1]#[NX1:2].[#6:3](=[#8])-[#8;H1]>>[c:1]1:[n:2]:[o]:[c:3]:[n]:1',
        'syndirella_name': '1,2,4-Oxadiazole_cyclization',
    },
    'aminoxadiazole': {
        'smirks': '[NX3;H2,H1:1]NC(=O).[#7;H2,H1:2]>>[c]1:[n:1]:[o]:[c]:[n:2]:1',
        'syndirella_name': '2-Amino-1,3,4-oxadiazole',
    },

    # ── Thiazole formation ────────────────────────────────────────────────
    'thiazole_2comp': {
        'smirks': '[a,A:1]-[C](=[S:4])[N&H2].[a,A:3]-[C:5](=[O])-[C;H1:6](-[$([Cl,Br,I]),$(OS(=O)(=O)C),$(OS(=O)(=O)c1ccc(C)cc1)])-[a,A:2]>>[c]-1([a,A:1])-[s:4]-[c:6](-[a,A:2])-[c:5](-[a,A:3])-[n]-1',
        'syndirella_name': 'Hantzsch_thiazole_synthesis',
    },
    'thiazole_3comp': {
        'smirks': '[#6:6]-[C;R0:1](=[OD1])-[CH1;R0:5](-[#6:7])-[*;#17,#35,#53].[NH2:2]-[C:3]=[SD1:4]>>[c:1]2(-[#6:6]):[n:2]:[c:3]:[s:4][c:5]([#6:7]):2',
        'syndirella_name': 'Hantzsch_thiazole_synthesis',
    },
    'thiazole_thiourea_3comp': {
        'smirks': '[#6:6]-[C;R0:1](=[OD1])-[CH1;R0:5](-[#6:7])-[*;#17,#35,#53].[NH2:2]-[C:3]=[SD1:4]>>[c:1]2(-[#6:6]):[n:2]:[c:3]:[s:4][c:5]([#6:7]):2',
        'syndirella_name': 'Hantzsch_thiazole_synthesis',
    },

    # ── Hydantoin ─────────────────────────────────────────────────────────
    'hydantoin': {
        'smirks': '[NX3;H2:1].[NX3:2][CX4:3][CX3:4](=O)[OX2]>>[N:1]C(=O)[N:2][C:3][C:4](=O)',
        'syndirella_name': 'Hydantoin_formation',
    },

    # ── Biginelli ─────────────────────────────────────────────────────────
    'biginelli': {
        'smirks': '[N&X3;H2:1]C(=[O,S])[N&X3;H2:2].[CX3H1:3](=O).[CX4H2:4]([CX3]=O)>>[N:1]C(=O)[N:2][C:3][C:4]',
        'syndirella_name': 'Biginelli_reaction',
    },

    # ── Cyclic guanidine ──────────────────────────────────────────────────
    'cyclic_guanidine': {
        'smirks': '[NX3;H2,H1:1].[CX3:2](=S)[NX3:3]>>[N:1][C:2](=[N])[N:3]',
        'syndirella_name': 'Cyclic_guanidine_formation',
    },

    # ── Tetrazole ─────────────────────────────────────────────────────────
    'tetrazole_acylation': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'tetrazole_alkylation': {
        'smirks': '[#7&X3;H2,H1;+0:5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[#7&X3;+0:5]-[C&X4;+0:1]',
        'syndirella_name': 'Nucleophilic_substitution_with_amine',
    },
    'tetrazole_sulfonamide': {
        'smirks': '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
        'syndirella_name': 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)',
    },

    # ── Ugi MCR ───────────────────────────────────────────────────────────
    'ugi_4cr': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[#7;H2,H1:3].[CX3H1:4](=O).[CX1-:5]#[NX2+:6]>>[#6:1](=[#8:2])-[#7:3]-[C:4]-[C:5](=[O])-[N:6]',
        'syndirella_name': 'Ugi_4CR',
    },
    'ugi_3cr': {
        'smirks': '[#7;H2,H1:3].[CX3H1:4](=O).[CX1-:5]#[NX2+:6]>>[#7:3]-[C:4]-[C:5](=[O])-[N:6]',
        'syndirella_name': 'Ugi_3CR',
    },

    # ── GBB MCR ───────────────────────────────────────────────────────────
    'gbb': {
        'smirks': '[#7;H2,H1:3].[CX3H1:4](=O).[CX1-:5]#[NX2+:6]>>[#7:3]-[C:4]-[C:5](=[O])-[N:6]',
        'syndirella_name': 'Groebke-Blackburn-Bienayme_3CR',
    },

    # ── Multi-component reactions ─────────────────────────────────────────
    # For multi-step Enamine reactions, we provide the primary (first)
    # bond-forming SMIRKS.  Syndirella handles multi-step routes natively.
    'multi_amide_arylation': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'multi_amide_amide': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'multi_sulfonamide_arylation': {
        'smirks': '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
        'syndirella_name': 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)',
    },
    'multi_sulfonamide_suzuki': {
        'smirks': '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
        'syndirella_name': 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)',
    },
    'multi_amide_suzuki': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'multi_amide_click': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'multi_alkylation_click': {
        'smirks': '[#7&X3;H2,H1;+0:5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[#7&X3;+0:5]-[C&X4;+0:1]',
        'syndirella_name': 'Nucleophilic_substitution_with_amine',
    },
    'multi_amide_alkylation': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'multi_reductive_amination_amide': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'multi_reductive_amination_arylation': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'multi_reductive_amination_sulfonamide': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'multi_reductive_amination_urea': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'multi_reductive_amination_carbamate': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'multi_reductive_amination_reductive_amination': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'multi_reductive_amination_acyl_halide': {
        'smirks': '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
        'syndirella_name': 'Reductive_amination',
    },
    'multi_arylation_arylation': {
        'smirks': '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',
        'syndirella_name': 'N-nucleophilic_aromatic_substitution',
    },
    'multi_arylation_suzuki': {
        'smirks': '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',
        'syndirella_name': 'N-nucleophilic_aromatic_substitution',
    },
    'multi_triazole_amide': {
        'smirks': '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
        'syndirella_name': 'Triazole_from_lactim',
    },
    'multi_triazole_urea': {
        'smirks': '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
        'syndirella_name': 'Triazole_from_lactim',
    },
    'multi_triazole_alkylation': {
        'smirks': '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
        'syndirella_name': 'Triazole_from_lactim',
    },
    'multi_triazole_sulfonamide': {
        'smirks': '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
        'syndirella_name': 'Triazole_from_lactim',
    },
    'multi_thiazole_amide': {
        'smirks': '[a,A:1]-[C](=[S:4])[N&H2].[a,A:3]-[C:5](=[O])-[C;H1:6](-[$([Cl,Br,I]),$(OS(=O)(=O)C),$(OS(=O)(=O)c1ccc(C)cc1)])-[a,A:2]>>[c]-1([a,A:1])-[s:4]-[c:6](-[a,A:2])-[c:5](-[a,A:3])-[n]-1',
        'syndirella_name': 'Hantzsch_thiazole_synthesis',
    },
    'multi_thiazole_sulfonamide': {
        'smirks': '[a,A:1]-[C](=[S:4])[N&H2].[a,A:3]-[C:5](=[O])-[C;H1:6](-[$([Cl,Br,I]),$(OS(=O)(=O)C),$(OS(=O)(=O)c1ccc(C)cc1)])-[a,A:2]>>[c]-1([a,A:1])-[s:4]-[c:6](-[a,A:2])-[c:5](-[a,A:3])-[n]-1',
        'syndirella_name': 'Hantzsch_thiazole_synthesis',
    },
    'multi_anhydride_amide_arylation': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'multi_sulfonamide_amide': {
        'smirks': '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
        'syndirella_name': 'Sulfonamide_Schotten-Baumann_with_amine_(intermolecular)',
    },
    'multi_acylsulfamide_alkylation': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[#7;H1:3]S(=O)(=O)>>[#6:1](=[#8:2])-[#7:3]',
        'syndirella_name': 'Amidation',
    },
    'multi_click_amide': {
        'smirks': '[NX1:1]~[NX2:2]~[NX2:3]-[#6:4].[CX2:5]#[CX2:6]-[#6:7]>>[#6:4]-[n:3]1:[n:2]:[n:1]:[c:5](-[#6:7]):[c:6]:1',
        'syndirella_name': 'CuAAC_1,2,3-triazole',
    },
    'multi_oxoisoindole_amide': {
        'smirks': '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
        'syndirella_name': 'Amidation',
    },
    'dihydrouracil': {
        'smirks': '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[N&X3;H2,H1;!$(NC=O);!$(NS):4]>>[#7&X3:3]-[#6](=[#8])-[#7&X3:4]',
        'syndirella_name': 'Formation_of_urea_from_two_amines',
    },
    'purine_substitution': {
        'smirks': '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',
        'syndirella_name': 'N-nucleophilic_aromatic_substitution',
    },

    # ── Generic fallbacks ─────────────────────────────────────────────────
    'generic_2comp': {
        'smirks': '',
        'syndirella_name': '',
    },
    'generic_3comp': {
        'smirks': '',
        'syndirella_name': '',
    },
    'generic_4comp': {
        'smirks': '',
        'syndirella_name': '',
    },
}


# ═══════════════════════════════════════════════════════════════════════════
# Post-reaction Deprotection SMIRKS
# ═══════════════════════════════════════════════════════════════════════════
# Enamine reactions with prefixes like Boc-, B-, K- involve a protecting
# group that is removed after the core bond-forming step.  The SMIRKS above
# model only the core reaction, so we need to apply deprotection afterward
# to produce the final product reported by Enamine.

DEPROTECTION_SMIRKS = {
    # Boc (tert-butyloxycarbonyl) on nitrogen  →  free amine
    'boc':        '[#7:1]C(=O)OC(C)(C)C>>[#7:1]',
    # Boc on aromatic nitrogen  →  free [nH]  (explicit H for kekulisation)
    # Used as a fallback when the general 'boc' SMIRKS fails to kekulise.
    'boc_arom':   '[n:1]C(=O)OC(C)(C)C>>[nH:1]',
    # tert-Butyl ester  →  carboxylic acid
    'tbu_ester':  '[C:1](=O)OC(C)(C)C>>[C:1](=O)O',
    # tert-Butyl ether  →  alcohol / phenol  (e.g. C-OtBu → C-OH)
    'tbu_ether':  '[#6:1]OC(C)(C)C>>[#6:1]O',
    # Methyl ester  →  carboxylic acid
    'me_ester':   '[C:1](=O)O[CH3]>>[C:1](=O)O',
    # Ethyl ester  →  carboxylic acid
    'et_ester':   '[C:1](=O)O[CH2][CH3]>>[C:1](=O)O',
    # MOM (methoxymethyl) on nitrogen  →  free N-H
    'mom':        '[#7:1][CH2]OC>>[#7:1]',
    # MOM on aromatic nitrogen  →  free [nH]  (kekulisation-safe)
    'mom_arom':   '[n:1][CH2]OC>>[nH:1]',
}

# Map Enamine reaction-name prefixes to the deprotection(s) needed.
# Order matters for prefixes that overlap (e.g. BBoc- before B-).
PREFIX_TO_DEPROTECTION = [
    ('BBoc-',  ['boc', 'boc_arom', 'tbu_ester', 'tbu_ether']),
    ('Boc-',   ['boc', 'boc_arom', 'tbu_ester', 'tbu_ether', 'mom', 'mom_arom']),
    ('B-',     ['tbu_ester', 'tbu_ether']),
    ('K-',     ['me_ester', 'et_ester']),
]


# ═══════════════════════════════════════════════════════════════════════════
# Mapping: Enamine Reaction Name → Reaction Class
# ═══════════════════════════════════════════════════════════════════════════
# Reused from extract-enamine-building-blocks.py.  Every reaction from
# enamine-real-tools-rxns.csv is mapped here.

ENAMINE_REACTION_TO_CLASS = {
    # ── AMIDATION variants ────────────────────────────────────────────────
    'AMIDE':              'amidation',
    'AMIDE1':             'amidation',
    'AMIDE2':             'amidation',
    'AMIDEMSCl':          'amidation',
    'K-AMIDE':            'amidation',
    'K-AMIDE1':           'amidation',
    'B-AMIDE':            'amidation',
    'B-AMIDE1':           'amidation',
    'B-AMIDE2':           'amidation',
    'Boc-AMIDE':          'amidation',
    'Boc-AMIDE1':         'amidation',
    'Boc-AMIDE2':         'amidation',
    'Boc-AMIDEMSCl':      'amidation',
    'BBoc-AMIDE1':        'amidation',
    'AMINOACID1':         'amidation',
    'R119-ACYLATION':     'amidation',
    'R039-Me3Si-Ac':      'acetylation',

    # ── ACETYLATION ───────────────────────────────────────────────────────
    'ACETYLATION':        'acetylation',

    # ── ACYLATION of O,C nucleophiles ─────────────────────────────────────
    'R111-ACYLATION':     'acylation_oc',

    # ── ACRYLAMIDE ────────────────────────────────────────────────────────
    'R017-ACRYLAMIDE':    'acrylamide',

    # ── ANHYDRIDE ─────────────────────────────────────────────────────────
    'ANHYDRIDE1':         'anhydride',

    # ── ACYLSULFAMIDE ─────────────────────────────────────────────────────
    'ACYLSULFAMIDE':      'acylsulfamide',
    'Boc-ACYLSULFAMIDE':  'acylsulfamide',
    'B-ACYLSULFAMIDE':    'acylsulfamide',

    # ── SULFONAMIDE ───────────────────────────────────────────────────────
    'SULFAMIDE':          'sulfonamide',
    'SULFAMIDE1':         'sulfonamide',
    'SULFAMIDE2':         'sulfonamide',
    'SULFAMIDE4-ShB':     'sulfonamide',
    'K-SULFAMIDE1':       'sulfonamide',
    'B-SULFAMIDE1':       'sulfonamide',
    'Boc-SULFAMIDE1':     'sulfonamide',
    'Boc-SULFAMIDE2':     'sulfonamide',
    'R039-Me3Si-Sac':     'sulfonamide',
    'R165-F-SULFAMIDE':   'sulfonamide',

    # ── SULFACYLATION of PHENOL ───────────────────────────────────────────
    'SULFACYLATION-O':    'sulfonamide_phenol',

    # ── VINYL SULFONAMIDE ─────────────────────────────────────────────────
    'VINYLSULFAMIDE':     'vinyl_sulfonamide',

    # ── REDUCTIVE AMINATION ───────────────────────────────────────────────
    'REDUCTION':          'reductive_amination',
    'REDUCTION3':         'reductive_amination',
    'REDUCTION-HET':      'reductive_amination',
    'REDUCTION-AMINOACID': 'reductive_amination',
    'R070-R3-AMINOACID':  'reductive_amination',
    'Boc-REDUCTION3':     'reductive_amination',
    'K-REDUCTION3':       'reductive_amination',
    'B-REDUCTION3':       'reductive_amination',
    'REDUCTION-TI':       'reductive_amination_ketone',
    'Boc-REDUCTION-TI':   'reductive_amination_ketone',

    # ── SUZUKI COUPLING ───────────────────────────────────────────────────
    'SUZUKI':             'suzuki',
    'Boc-SUZUKI':         'suzuki',
    'AC-SUZUKI':          'suzuki',
    'K-SUZUKI':           'suzuki',

    # ── N,O,S-ARYLATION ──────────────────────────────────────────────────
    'ARYLATION':          'arylation',
    'ARYLATION60':        'arylation',
    'Boc-ARYLATION':      'arylation',
    'K-ARYLATION':        'arylation',
    'B-ARYLATION':        'arylation',

    # ── N-ALKYLATION of amines ────────────────────────────────────────────
    'MA1':                'alkylation_amine',
    'MA3':                'alkylation_amine',
    'MA2-POLY':           'alkylation_amine',
    'Boc-MA2':            'alkylation_amine',
    'K-MA2':              'alkylation_amine',
    'B-MA2':              'alkylation_amine',

    # ── ALKYLATION of N/O/S nucleophiles ──────────────────────────────────
    'ESTER-NS':           'alkylation_nos',
    'OXYME':              'alkylation_nos',
    'Boc-OXYME':          'alkylation_nos',
    'B-OXYME':            'alkylation_nos',
    'ESTER-PH':           'alkylation_nos',
    'Boc-ESTER-NS':       'alkylation_nos',
    'B-ESTER-NS':         'alkylation_nos',

    # ── ALKYLATION of acid (ester) ────────────────────────────────────────
    'ESTER':              'alkylation_acid',
    'B-ESTER':            'alkylation_acid',

    # ── S-ALKYLATION / SULFIDE / SULFOXIDE / SULFONE ──────────────────────
    'SULFIDE':            'alkylation_thiol',
    'SULFONE':            'alkylation_thiol',
    'SULFOXIDE':          'alkylation_thiol',
    'SULFIDE-OX1':        'alkylation_thiol',
    'SULFIDE-OX2':        'alkylation_thiol',

    # ── EPOXIDE RING-OPENING ──────────────────────────────────────────────
    'EPOXYDE':            'epoxide_amine',
    'EPOXYDE2':           'epoxide_amine',
    'EPOXYDE3':           'epoxide_phenol',

    # ── UREA FORMATION ────────────────────────────────────────────────────
    'UREA':               'urea_isocyanate',
    'UREA-CDI':           'urea',
    'URETANE':            'urea_carbamate',
    'URETANE1':           'urea_isocyanate',
    'URETANE2':           'urea_carbamate',
    'URETANE4':           'urea_isocyanate',
    'URETANE-T':          'urea',
    'Boc-URETANE':        'urea_carbamate',
    'Boc-URETANE1':       'urea_isocyanate',
    'Boc-URETANE4':       'urea_carbamate',
    'B-URETANE1':         'urea_isocyanate',
    'B-URETANE4':         'urea_isocyanate',
    'R039-Me3Si-Urea':    'urea_isocyanate',

    # ── CARBAMATE ─────────────────────────────────────────────────────────
    'R157-CARBAMATE':     'carbamate_cdi',
    'R192-CARBAMATE':     'carbamate_chloroformate',
    'R039-Me3Si-Uretane': 'carbamate_chloroformate',

    # ── THIOUREA ──────────────────────────────────────────────────────────
    'THYOUREA':           'thiourea',

    # ── OXAMIDE ───────────────────────────────────────────────────────────
    'OXAMIDE':            'oxamide',
    'OXAMIDE1Pri':        'oxamide',
    'OXAMIDE1Sec':        'oxamide',
    'B-OXAMIDE':          'oxamide',
    'Boc-OXAMIDE':        'oxamide',
    'B-OXAMIDE1Pri':      'oxamide',
    'B-OXAMIDE1Sec':      'oxamide',
    'Boc-OXAMIDE1Pri':    'oxamide',
    'Boc-OXAMIDE1Sec':    'oxamide',

    # ── HYDRAZONE ─────────────────────────────────────────────────────────
    'HYDRAZON':           'hydrazone',

    # ── KNOEVENAGEL ───────────────────────────────────────────────────────
    'METHYLEN':           'knoevenagel',
    'METHYLEN1':          'knoevenagel',

    # ── MANNICH ───────────────────────────────────────────────────────────
    'MANNICH':            'mannich',

    # ── CLICK / TRIAZOLE ──────────────────────────────────────────────────
    'R030-TRIAZOLE':      'click_triazole',
    'R040-TRIAZOLE':      'click_from_amine',
    'Boc-R030-TRIAZOLE':  'click_triazole',
    'K-R040-TRIAZOLE':    'click_from_amine',
    'B-R040-TRIAZOLE':    'click_from_amine',
    'Boc-R040-TRIAZOLE':  'click_from_amine',
    'Radd-178-AT-Click':  'click_from_amine',

    # ── TRIAZOLE from lactim ──────────────────────────────────────────────
    'R036-Lactim':        'triazole_lactim',

    # ── OXADIAZOLE ────────────────────────────────────────────────────────
    'AMIDOXIME1':         'oxadiazole_amidoxime',
    'Boc-AMIDOXIME1':     'oxadiazole_amidoxime',
    'AMIDOXIME1-CN':      'oxadiazole_nitrile',
    'Boc-AMIDOXIME1-CN':  'oxadiazole_nitrile',
    'ANONYMOUS3':         'oxadiazole_cyclization',
    'R037a-AMINOXADIAZOLE': 'aminoxadiazole',
    'Radd-117-Ac-CN2Oxdzln': 'oxadiazole_nitrile',

    # ── THIAZOLE ──────────────────────────────────────────────────────────
    'THIAZOLE':           'thiazole_3comp',
    'THIAZOLE2':          'thiazole_thiourea_3comp',
    'THIAZOLE3':          'thiazole_thiourea_3comp',
    'THIAZOLE4':          'thiazole_2comp',

    # ── HYDANTOIN ─────────────────────────────────────────────────────────
    'HYDANTOIN4':         'hydantoin',
    'HYDANTOIN6':         'hydantoin',

    # ── BIGINELLI ─────────────────────────────────────────────────────────
    'BIGINELLI':          'biginelli',

    # ── CYCLIC GUANIDINE / THIOURONIUM ────────────────────────────────────
    'R113-CyclicTHIOURONIUM': 'cyclic_guanidine',
    'THIOURONIUM1':       'cyclic_guanidine',

    # ── TETRAZOLE ─────────────────────────────────────────────────────────
    'T-ACYLATION':        'tetrazole_acylation',
    'T-MA2':              'tetrazole_alkylation',
    'Radd-139-T-SO2N':    'tetrazole_sulfonamide',

    # ── UGI MCR ───────────────────────────────────────────────────────────
    'Radd-152-Ugi':       'ugi_4cr',
    'Radd-173-Split-Ugi': 'ugi_4cr',
    'Radd-176-Ugi-MF':    'ugi_3cr',
    'Radd-160-Ugi-T':     'ugi_3cr',
    'Radd-203-Ugi-3CR-CA': 'ugi_3cr',
    'Radd-204-Ugi-3CR-AA': 'ugi_3cr',

    # ── GBB MCR ───────────────────────────────────────────────────────────
    'Radd-190-GBB':       'gbb',

    # ── MULTI-COMPONENT ───────────────────────────────────────────────────
    'Radd-022-AA-Ac-Boc-Ac':  'multi_amide_amide',
    'Radd-049-DA-Ac-Boc-Ac':  'multi_amide_amide',
    'Radd-138-Ac-tBu-Ac':     'multi_amide_amide',
    'Radd-062-DA-Ar-Boc-Ac':  'multi_amide_arylation',
    'Radd-010-DA-Ac-Boc-Ar':  'multi_amide_arylation',
    'Radd-063-DA-SO2N-Boc-Ac': 'multi_sulfonamide_amide',
    'Radd-039-DA-Ac-Boc-Alk':  'multi_amide_alkylation',
    'Radd-065-DA-Alk-Boc-Ac':  'multi_amide_alkylation',
    'Radd-036-NAc-Click':      'multi_amide_click',
    'Radd-093-Ac-Click#':      'multi_amide_click',
    'Radd-058-Click-Boc-Ac':   'multi_click_amide',
    'Radd-115-NAc-Click-OH':   'multi_amide_click',
    'R051-Alk-Click':          'multi_alkylation_click',
    'Radd-053-Ac-Suzuki':      'multi_amide_suzuki',
    'Radd-051-SO2N-Suzuki':    'multi_sulfonamide_suzuki',
    'Radd-130-SO2N-Ar':        'multi_sulfonamide_arylation',
    'R071-Rd-Acl':             'multi_reductive_amination_acyl_halide',
    'R071b-Rd-Ac':             'multi_reductive_amination_amide',
    'Radd-078-Ac-Boc-RAT':    'multi_reductive_amination_amide',
    'R071a-Rd-Ar':             'multi_reductive_amination_arylation',
    'R071c-Rd-Sac':            'multi_reductive_amination_sulfonamide',
    'R071d-Rd-Ur':             'multi_reductive_amination_urea',
    'R071g-Rd-Ur1':            'multi_reductive_amination_urea',
    'R071h-Rd-Ur4':            'multi_reductive_amination_urea',
    'R071e-Rd-Carbamate':      'multi_reductive_amination_carbamate',
    'R071f-Rd-Rd':             'multi_reductive_amination_reductive_amination',
    'Radd-158-RedAmin-Suzuki': 'multi_arylation_suzuki',
    'Radd-159-Arylation-Suzuki': 'multi_arylation_suzuki',
    'Radd-054-Alk-Suzuki':     'multi_arylation_suzuki',
    'Radd-168-Ar-Ar':          'multi_arylation_arylation',
    'ARYLATION2':              'multi_arylation_arylation',
    'ARYLATION3':              'multi_arylation_arylation',
    'Radd-111-Anh-Ac-Ar':     'multi_anhydride_amide_arylation',
    'Radd-113-Ac-SAc':         'multi_sulfonamide_amide',
    'Radd-179-SAc-Alk':        'multi_acylsulfamide_alkylation',
    'R194a-Lactim-Boc-Ac':     'multi_triazole_amide',
    'R194b-Lactim-Boc-Ur1':    'multi_triazole_urea',
    'R194c-Lactim-Boc-Alk':    'multi_triazole_alkylation',
    'R194d-Lactim-Boc-Sac':    'multi_triazole_sulfonamide',
    'Radd-068-135Subst124Triazole': 'multi_click_amide',
    'Radd-196-Thzl-Boc-Ac':   'multi_thiazole_amide',
    'Radd-202-Thzl-Boc-SO2N':  'multi_thiazole_sulfonamide',
    'R198a-OXOISOINDOLE-Boc-Ac': 'multi_oxoisoindole_amide',
    'R142-DIHYDROURACYL':      'dihydrouracil',
    'R055':                    'purine_substitution',
    'R003':                    'generic_3comp',
    'R014':                    'generic_3comp',
    'R015':                    'generic_3comp',
}


# ═══════════════════════════════════════════════════════════════════════════
# Core Logic
# ═══════════════════════════════════════════════════════════════════════════

def get_smirks_for_reaction(reaction_name):
    """
    Look up the reaction SMIRKS and Syndirella reaction name for an
    Enamine reaction name.

    Returns (reaction_class, smirks, syndirella_name).
    """
    rxn_class = ENAMINE_REACTION_TO_CLASS.get(reaction_name, '')
    if not rxn_class:
        return ('', '', '')

    info = REACTION_CLASS_SMIRKS.get(rxn_class, {})
    smirks = info.get('smirks', '')
    syndirella_name = info.get('syndirella_name', '')
    return (rxn_class, smirks, syndirella_name)


def validate_smirks_parse(smirks):
    """
    Validate a SMIRKS string can be parsed by RDKit.
    Returns True if valid, False otherwise.
    """
    if not smirks:
        return False
    from rdkit.Chem import AllChem
    rxn = AllChem.ReactionFromSmarts(smirks)
    return rxn is not None


def _canonical_smiles(smiles):
    """Return canonical SMILES or None if unparseable."""
    if not smiles or str(smiles).strip() in ('', 'nan'):
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def _detect_deprotections(reaction_name):
    """
    Given an Enamine reaction name, return a list of deprotection types
    indicated by the prefix (e.g. 'Boc-AMIDE1' → ['boc']).
    """
    for prefix, deprot_types in PREFIX_TO_DEPROTECTION:
        if reaction_name.startswith(prefix):
            return list(deprot_types)
    return []


def _apply_deprotections(mol, deprot_types):
    """
    Apply one or more deprotection SMIRKS to a molecule (in-place iteration).
    Returns a new mol with all specified protecting groups removed.
    """
    from rdkit.Chem import AllChem

    for dt in deprot_types:
        smirks = DEPROTECTION_SMIRKS.get(dt, '')
        if not smirks:
            continue
        rxn = AllChem.ReactionFromSmarts(smirks)
        if rxn is None:
            continue
        # Apply repeatedly until no more matches (handles multiple PGs)
        for _ in range(10):
            products = rxn.RunReactants((mol,))
            if not products:
                break
            # Take the first product
            new_mol = products[0][0]
            try:
                Chem.SanitizeMol(new_mol)
                mol = new_mol
            except Exception:
                break
    return mol


def run_smirks_on_reactants(smirks, reactant_smiles_list):
    """
    Apply a SMIRKS to a list of reactant SMILES and return the set of
    product SMILES (canonical) that RDKit generates.

    Returns a list of canonical SMILES strings, or [] if the reaction
    fails or produces no products.
    """
    from rdkit.Chem import AllChem

    if not smirks or not reactant_smiles_list:
        return []

    rxn = AllChem.ReactionFromSmarts(smirks)
    if rxn is None:
        return []

    # Convert SMILES to mol objects
    reactant_mols = []
    for smi in reactant_smiles_list:
        mol = Chem.MolFromSmiles(str(smi).strip())
        if mol is None:
            return []
        reactant_mols.append(mol)

    # Check reactant count matches SMIRKS expectation
    n_expected = rxn.GetNumReactantTemplates()
    if len(reactant_mols) != n_expected:
        # Try all permutations of the reactants (order matters in SMIRKS)
        from itertools import permutations, combinations
        # If we have more reactants than templates, try all combinations
        if len(reactant_mols) > n_expected:
            combos = combinations(range(len(reactant_mols)), n_expected)
        else:
            return []

        products = set()
        for combo in combos:
            subset = [reactant_mols[i] for i in combo]
            for perm in permutations(subset):
                try:
                    ps = rxn.RunReactants(tuple(perm))
                    for product_set in ps:
                        for p in product_set:
                            try:
                                Chem.SanitizeMol(p)
                                products.add(Chem.MolToSmiles(p))
                            except Exception:
                                pass
                except Exception:
                    pass
        return list(products)

    # Try all permutations of reactants (SMIRKS order may differ from CSV)
    from itertools import permutations
    products = set()
    for perm in permutations(reactant_mols):
        try:
            ps = rxn.RunReactants(tuple(perm))
            for product_set in ps:
                for p in product_set:
                    try:
                        Chem.SanitizeMol(p)
                        products.add(Chem.MolToSmiles(p))
                    except Exception:
                        pass
        except Exception:
            pass
    return list(products)


def validate_smirks_produces_product(smirks, reactant_smiles_list,
                                      expected_product_smiles,
                                      deprot_types=None):
    """
    Check whether applying the SMIRKS to the reactants produces a molecule
    matching the expected product (by canonical SMILES comparison).

    If deprot_types is given, deprotection SMIRKS are applied to each
    product before comparison.  This handles Enamine Boc-, B-, K- etc.
    reactions where the protecting group is removed after the core step.

    Returns:
        (validated: bool or None, actual_product: str or '')
        - True   if at least one product canonical SMILES matches expected
        - False  if products were generated but none match
        - None   if reaction could not be run or product missing
    """
    if not expected_product_smiles or str(expected_product_smiles).strip() in ('', 'nan'):
        return (None, '')
    if not smirks:
        return (None, '')

    expected_can = _canonical_smiles(expected_product_smiles)
    if expected_can is None:
        return (None, '')

    products = run_smirks_on_reactants(smirks, reactant_smiles_list)
    if not products:
        return (None, '')

    if deprot_types is None:
        deprot_types = []

    # Check each product (with optional deprotection)
    for prod_smi in products:
        prod_can = _canonical_smiles(prod_smi)
        if prod_can and prod_can == expected_can:
            return (True, prod_smi)

    # If deprotections specified, try applying them
    if deprot_types:
        for prod_smi in products:
            mol = Chem.MolFromSmiles(prod_smi)
            if mol is None:
                continue
            deprot_mol = _apply_deprotections(mol, deprot_types)
            deprot_smi = Chem.MolToSmiles(deprot_mol)
            deprot_can = _canonical_smiles(deprot_smi)
            if deprot_can and deprot_can == expected_can:
                return (True, deprot_smi)

    # If prefix-based deprotections didn't work, try auto-detecting
    # common protecting groups in the product as a fallback
    if not deprot_types:
        all_types = list(DEPROTECTION_SMIRKS.keys())
        for prod_smi in products:
            mol = Chem.MolFromSmiles(prod_smi)
            if mol is None:
                continue
            deprot_mol = _apply_deprotections(mol, all_types)
            deprot_smi = Chem.MolToSmiles(deprot_mol)
            deprot_can = _canonical_smiles(deprot_smi)
            if deprot_can and deprot_can == expected_can:
                return (True, deprot_smi)

    # Return the first product for inspection even if no match
    return (False, products[0] if products else '')


def add_smirks_to_csv(input_csv):
    """
    Read an Enamine REAL Tools output CSV and add SMIRKS columns.
    Validates each SMIRKS by running the reaction on the reactants.
    Returns the augmented DataFrame.
    """
    df = pd.read_csv(input_csv)

    rxn_classes = []
    smirks_list = []
    syndirella_names = []
    validated_list = []
    actual_products = []

    for _, row in df.iterrows():
        rxn_name = row.get('reaction_name', '')
        if pd.isna(rxn_name):
            rxn_name = ''
        rxn_name = str(rxn_name).strip()

        rxn_class, smirks, syn_name = get_smirks_for_reaction(rxn_name)
        rxn_classes.append(rxn_class)
        smirks_list.append(smirks)
        syndirella_names.append(syn_name)

        # Detect deprotection needed from reaction name prefix
        deprot_types = _detect_deprotections(rxn_name)

        # Collect reactant SMILES from the row
        reactants = []
        for i in range(1, 5):
            bb_smi = row.get(f'reactant_{i}_smiles', '')
            if bb_smi and not pd.isna(bb_smi) and str(bb_smi).strip():
                reactants.append(str(bb_smi).strip())

        # Expected product
        expected = row.get('found_smiles', '')
        if not expected or pd.isna(expected) or str(expected).strip() == '':
            expected = row.get('target_smiles', '')

        validated, actual_prod = validate_smirks_produces_product(
            smirks, reactants, expected, deprot_types=deprot_types
        )
        validated_list.append(validated if validated is not None else '')
        actual_products.append(actual_prod)

    df['reaction_class'] = rxn_classes
    df['reaction_smirks'] = smirks_list
    df['syndirella_reaction_name'] = syndirella_names
    df['smirks_validated'] = validated_list
    df['actual_product_smiles'] = actual_products

    return df


def print_summary(df):
    """Print summary statistics about the SMIRKS extraction."""
    total = len(df)
    found = df[df['status'] == 'FOUND'] if 'status' in df.columns else df
    n_found = len(found)

    has_smirks = (found['reaction_smirks'].astype(str).str.strip() != '')
    n_with_smirks = has_smirks.sum()
    n_without = n_found - n_with_smirks

    print(f"\n{'='*70}")
    print(f"SMIRKS EXTRACTION SUMMARY")
    print(f"{'='*70}")
    print(f"Total rows              : {total}")
    if 'status' in df.columns:
        print(f"FOUND rows              : {n_found}")
    print(f"With reaction SMIRKS    : {n_with_smirks}")
    print(f"Without SMIRKS          : {n_without}")

    if n_found > 0:
        print(f"Coverage                : {n_with_smirks/n_found*100:.1f}%")

    # By reaction class
    if n_with_smirks > 0:
        print(f"\n--- Reactions by Syndirella Name ---")
        syn_counts = found[has_smirks].groupby('syndirella_reaction_name').size()
        for name, count in sorted(syn_counts.items(), key=lambda x: -x[1]):
            print(f"  {name:55s}: {count:4d}")

    # Unmapped reactions
    if n_without > 0:
        unmapped = found[~has_smirks]
        if 'reaction_name' in unmapped.columns:
            missing_rxns = unmapped['reaction_name'].dropna().unique()
            if len(missing_rxns) > 0:
                print(f"\n--- Unmapped Enamine Reaction Names ---")
                for rxn in sorted(missing_rxns):
                    count = (unmapped['reaction_name'] == rxn).sum()
                    print(f"  {rxn:30s}: {count} rows")

    # Validate SMIRKS parse
    unique_smirks = found[has_smirks]['reaction_smirks'].unique()
    n_valid = 0
    n_invalid = 0
    for s in unique_smirks:
        if validate_smirks_parse(s):
            n_valid += 1
        else:
            n_invalid += 1
            print(f"\n  WARNING: Unparseable SMIRKS: {s[:80]}...")

    print(f"\n--- SMIRKS Parse Check ---")
    print(f"  Unique SMIRKS patterns: {len(unique_smirks)}")
    print(f"  RDKit parseable       : {n_valid}")
    print(f"  Failed parse          : {n_invalid}")

    # Product validation (the real test)
    print(f"\n--- SMIRKS Product Validation ---")
    print(f"    (Does applying the SMIRKS to the building blocks produce")
    print(f"     the expected product, compared by canonical SMILES?)")

    val_col = found['smirks_validated']
    n_pass = (val_col == True).sum()
    n_fail = (val_col == False).sum()
    n_na = len(val_col) - n_pass - n_fail
    print(f"  Passed  : {n_pass}")
    print(f"  Failed  : {n_fail}")
    print(f"  N/A     : {n_na}  (reaction not run / no product to compare)")

    if n_pass + n_fail > 0:
        pct = n_pass / (n_pass + n_fail) * 100
        print(f"  Success rate (excl. N/A): {pct:.1f}%")

    if n_fail > 0:
        print(f"\n  Failed rows (first 20):")
        failures = found[val_col == False].head(20)
        for _, row in failures.iterrows():
            expected = row.get('found_smiles', row.get('target_smiles', ''))
            actual = row.get('actual_product_smiles', '')
            rxn = row.get('reaction_name', '')
            print(f"    Rxn: {rxn}")
            print(f"      Expected: {str(expected)[:70]}")
            print(f"      Got     : {str(actual)[:70]}")

            # Show reactants
            reactants = []
            for i in range(1, 5):
                bb = row.get(f'reactant_{i}_smiles', '')
                if bb and not pd.isna(bb) and str(bb).strip():
                    reactants.append(str(bb).strip()[:50])
            print(f"      Reactants: {' + '.join(reactants)}")

    # By reaction class breakdown
    if n_pass > 0 or n_fail > 0:
        print(f"\n  By reaction class:")
        for rxn_cls in sorted(found[has_smirks]['reaction_class'].unique()):
            subset = found[found['reaction_class'] == rxn_cls]
            sub_pass = (subset['smirks_validated'] == True).sum()
            sub_fail = (subset['smirks_validated'] == False).sum()
            sub_na = len(subset) - sub_pass - sub_fail
            total_tested = sub_pass + sub_fail
            rate = f"{sub_pass/total_tested*100:.0f}%" if total_tested > 0 else "N/A"
            print(f"    {rxn_cls:40s}: {sub_pass:3d} pass, "
                  f"{sub_fail:3d} fail, {sub_na:3d} n/a  ({rate})")

    print(f"{'='*70}")


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Add Syndirella-compatible reaction SMIRKS to Enamine '
                    'REAL Tools output CSV.'
    )
    parser.add_argument(
        'input_csv',
        help='Path to the Enamine REAL Tools synthesis output CSV'
    )
    parser.add_argument(
        'output_csv', nargs='?', default=None,
        help='Output CSV path (default: <input>_smirks.csv)'
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
        output_csv = f"{base}_smirks{ext}"

    print(f"Input  : {os.path.abspath(input_csv)}")
    print(f"Output : {os.path.abspath(output_csv)}")
    print(f"Mapped Enamine reactions: {len(ENAMINE_REACTION_TO_CLASS)}")

    df = add_smirks_to_csv(input_csv)
    print_summary(df)

    df.to_csv(output_csv, index=False)
    print(f"\nOutput written: {os.path.abspath(output_csv)}")
    print(f"  {len(df)} rows, added columns: reaction_class, "
          f"reaction_smirks, syndirella_reaction_name, "
          f"smirks_validated, actual_product_smiles")


if __name__ == '__main__':
    main()
