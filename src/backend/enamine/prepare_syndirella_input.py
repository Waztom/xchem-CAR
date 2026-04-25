"""
Prepare Syndirella Manual Input via Enamine REAL Tools
=====================================================
Takes a Syndirella automatic input CSV, calls the Enamine REAL Tools API
for retrosynthesis (or uses a pre-existing output CSV), and produces:

  1. An **Enamine REAL Tools CSV** — the retrosynthesis results from
     the API (target_smiles, status, reaction_name, reactants, …).
  2. A **Syndirella manual input CSV** — with reaction names, reactants,
     products AND hit/template metadata from the auto input.
  3. A **SMIRKS JSON** file — keyed by Enamine reaction name, ready to be
     merged into Syndirella's RXN_SMIRKS_CONSTANTS.json.

The reaction names in the output are the Enamine REAL Tools reaction
types (e.g. AMIDE1, REDUCTION3, Boc-SUZUKI).

Usage:
    # API mode (default) – only the Syndirella auto CSV is required:
    python prepare_syndirella_input.py <syndirella_auto.csv> \\
        [-o manual.csv] [-j smirks.json] [-e enamine_output.csv]

    # Offline mode – skip the API call with a pre-existing Enamine CSV:
    python prepare_syndirella_input.py <syndirella_auto.csv> \\
        --enamine-csv <existing_enamine.csv> [-o manual.csv] [-j smirks.json]
"""

import sys
import os
import json
import argparse

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

# Suppress RDKit warnings
RDLogger.logger().setLevel(RDLogger.ERROR)

# Canonical-tautomer helper (module-level singleton)
_TAUT_ENUMERATOR = rdMolStandardize.TautomerEnumerator()


# ═══════════════════════════════════════════════════════════════════════════
# Reaction SMIRKS by Reaction Class
# ═══════════════════════════════════════════════════════════════════════════
# Each key is a reaction class.  The SMIRKS are derived from Syndirella's
# RXN_SMIRKS_CONSTANTS.json, using parent reaction patterns.

REACTION_CLASS_SMIRKS = {
    # ── Amide bond formers ─────────────────────────────────────────────
    'amidation':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'acetylation':
        '[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3].[#6X3;!$(C-N):1](=[OX1])[#17,#9,#35]>>[#6:1](=[#8:2])-[#7:3]',
    'acylation_oc':
        '[Cl,F,Br,I]-[C;H0;D3;+0:1](=[O,S;D1;H0:2])-[*:3].[OH;D1;+0:4]-[*:5]>>[O,S;D1;H0:2]=[C;H0;D3;+0:1](-[*:3])-[O;H0;D2;+0:4]-[*:5]',
    'acrylamide':
        '[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3].[#6X3;!$(C-N):1](=[OX1])[#17,#9,#35]>>[#6:1](=[#8:2])-[#7:3]',
    'anhydride':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'acylsulfamide':
        '[#6:1](=[#8:2])-[#8;H1].[#7;H2,H1:3][#16:4](=[#8:5])(=[#8:6])>>[#6:1](=[#8:2])-[#7:3][#16:4](=[#8:5])(=[#8:6])',

    # ── Sulfonamide ───────────────────────────────────────────────────────
    'sulfonamide':
        '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
    'sulfonamide_phenol':
        '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[OX2H1:5][c:6]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[O:5][c:6]',
    'vinyl_sulfonamide':
        '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',

    # ── Reductive amination ───────────────────────────────────────────────
    'reductive_amination':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'reductive_amination_ketone':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',

    # ── Suzuki coupling ──────────────────────────────────────────────────
    'suzuki':
        '[#6:1]-[Cl,Br,I].[#6:2]-[B](-[O])(-[O])>>[#6:1]-[#6:2]',

    # ── Stille coupling ──────────────────────────────────────────────────
    'stille':
        '[#6:1]-[#50].[#6:2]-[Cl,Br,I]>>[#6:1]-[#6:2]',

    # ── N,O,S-Arylation ─────────────────────────────────────────────────
    'arylation':
        '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',

    # ── Alkylation of amines ──────────────────────────────────────────────
    'alkylation_amine':
        '[#7&X3;H2,H1;+0:5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[#7&X3;+0:5]-[C&X4;+0:1]',

    # ── Alkylation of N/O/S nucleophiles ──────────────────────────────────
    'alkylation_nos':
        '[$([#7;H2,H1;+0]),$([#8;H1;+0]),$([#16;H1;+0]):5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[*:5]-[C&X4;+0:1]',

    # ── Alkylation of acid (ester) ────────────────────────────────────────
    'alkylation_acid':
        '[Cl,Br]-[#6;H2,H1,H0;+0:1].[OH;+0:2]>>[O&X2;H0;D2;+0:2]-[#6&X4;D2;+0:1]',

    # ── S-Alkylation (thiols) ─────────────────────────────────────────────
    'alkylation_thiol':
        '[Cl,Br,I]-[C;H2,H1;+0:1].[SH;+0:2]>>[S&X2;H0;D2;+0:2]-[C&X4;D2;+0:1]',

    # ── Epoxide ring-opening ──────────────────────────────────────────────
    'epoxide_amine':
        '[#6:1]-[CH;D3;+0:2]1-[CH2;D2;+0:3]-[O;H0;D2;+0:4]-1.[N;H2,H1;+0:6]-[#6:7]>>[#6:1]-[C;+0:2](-[OH;D1;+0:4])-[C&X4;+0:3]-[N&X3;+0:6]-[#6:7]',
    'epoxide_phenol':
        '[#6:1]-[CH;D3;+0:2]1-[CH2;D2;+0:3]-[O;H0;D2;+0:4]-1.[OH;D1;+0:6]-[c:7]>>[#6:1]-[C;+0:2](-[OH;D1;+0:4])-[C&X4;+0:3]-[O;+0:6]-[c:7]',

    # ── Urea formation ───────────────────────────────────────────────────
    'urea':
        '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[N&X3;H2,H1;!$(NC=O);!$(NS):4]>>[#7&X3:3]-[#6](=[#8])-[#7&X3:4]',
    'urea_isocyanate':
        '[*:1]-[N;H0;D2;+0:2]=[C;H0;D2;+0:3]=[O;H0;D1;+0:4].[#6;+0:5]-[N;H2,H1;+0:6]>>[*:1]-[N;+0:2]-[C;H0;D3;+0:3](=[O;H0;D1;+0:4])-[N;+0:6]-[#6;+0:5]',
    'urea_carbamoyl_chloride':
        '[#7;H2,H1;+0:4].[#7:1]-[CX3:2](=[OX1:3])-[Cl]>>[#7:4]-[C:2](=[O:3])-[#7:1]',
    'urea_carbamate':
        '[#7:1]-[C](=[#8X1])[nX3]1[c&H1X3][c&H1X3][nX2][c&H1X3]1.[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):2]>>[#7:1]-[C](=[#8X1])-[#7:2]',
    'urea_phenyl_carbamate':
        '[NX3;H2,H1;+0:1].[#7:2][CX3:3](=[OX1:4])[OX2H0]>>[N:1][C:3](=[O:4])[N:2]',

    # ── Carbamate ─────────────────────────────────────────────────────────
    'carbamate_cdi':
        '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[#6:1]-[OH;D1;+0]>>[#7&X3:3]-[#6](=[#8])-[O]-[#6:1]',
    'carbamate_chloroformate':
        '[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3].[Cl][CX3:1](=[OX1:2])[OX2:4]>>[#7:3]-[CX3:1](=[OX1:2])[OX2:4]',

    # ── Thiourea ──────────────────────────────────────────────────────────
    'thiourea':
        '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[NX2:1]=[CX2:2]=[SX1:4]>>[#7&X3:3]-[C:2](=[S:4])-[N:1]',

    # ── Oxamide ───────────────────────────────────────────────────────────
    'oxamide':
        '[NX3;H2,H1;+0:1].[NX3;H2,H1;+0:2]>>[N:1]C(=O)C(=O)[N:2]',

    # ── Hydrazone condensation ────────────────────────────────────────────
    'hydrazone':
        '[CX3H1:1](=[OX1]).[NX3;H2,H1:2][NX3:3]>>[C:1]=[N:2][N:3]',

    # ── Knoevenagel condensation ──────────────────────────────────────────
    'knoevenagel':
        '[$([CX4H2]([CX3]=O)),$([CX4H2](C#N)):1].[CX3H1:2](=O)>>[C:1]=[C:2]',

    # ── Esterification ────────────────────────────────────────────────────
    'esterification':
        '[OX2H1:1][#6:2].[#6:3][CX3:4](=[OX1:5])[OX2H1]>>[#6:2][O:1][C:4](=[O:5])[#6:3]',

    # ── Mannich reaction ──────────────────────────────────────────────────
    'mannich':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',

    # ── Click / Triazole ──────────────────────────────────────────────────
    'click_triazole':
        '[NX1:1]~[NX2:2]~[NX2:3]-[#6:4].[CX2:5]#[CX2:6]-[#6:7]>>[#6:4]-[n:3]1:[n:2]:[n:1]:[c:5](-[#6:7]):[c:6]:1',
    'click_from_amine':
        '[#7;H2;D1:1]-[#6:4].[CX2:5]#[CX2:6]-[#6:7]>>[#6:4]-[n:1]1:[n]:[n]:[c:5](-[#6:7]):[c:6]:1',

    # ── Triazole from lactim ──────────────────────────────────────────────
    'triazole_lactim':
        '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',

    # ── Oxadiazole ────────────────────────────────────────────────────────
    'oxadiazole_amidoxime':
        '[NX3:1]C(=[NX2:2])[OX2H1:3].[#6:4](=[#8])-[#8;H1]>>[c:4]1:[n:2]:[o:3]:[c]:[n:1]:1',
    'oxadiazole_nitrile':
        '[CX2:1]#[NX1:2].[#6:3](=[#8])-[#8;H1]>>[c:1]1:[n:2]:[o]:[c:3]:[n]:1',
    'oxadiazole_cyclization':
        '[CX2:1]#[NX1:2].[#6:3](=[#8])-[#8;H1]>>[c:1]1:[n:2]:[o]:[c:3]:[n]:1',
    'aminoxadiazole':
        '[NX3;H2,H1:1]NC(=O).[#7;H2,H1:2]>>[c]1:[n:1]:[o]:[c]:[n:2]:1',

    # ── Thiazole formation ────────────────────────────────────────────────
    'thiazole_2comp':
        '[a,A:1]-[C](=[S:4])[N&H2].[a,A:3]-[C:5](=[O])-[C;H1:6](-[$([Cl,Br,I]),$(OS(=O)(=O)C),$(OS(=O)(=O)c1ccc(C)cc1)])-[a,A:2]>>[c]-1([a,A:1])-[s:4]-[c:6](-[a,A:2])-[c:5](-[a,A:3])-[n]-1',
    'thiazole_3comp':
        '[#6:6]-[C;R0:1](=[OD1])-[CH1;R0:5](-[#6:7])-[*;#17,#35,#53].[NH2:2]-[C:3]=[SD1:4]>>[c:1]2(-[#6:6]):[n:2]:[c:3]:[s:4][c:5]([#6:7]):2',
    'thiazole_thiourea_3comp':
        '[#6:6]-[C;R0:1](=[OD1])-[CH1;R0:5](-[#6:7])-[*;#17,#35,#53].[NH2:2]-[C:3]=[SD1:4]>>[c:1]2(-[#6:6]):[n:2]:[c:3]:[s:4][c:5]([#6:7]):2',

    # ── Hydantoin ─────────────────────────────────────────────────────────
    'hydantoin':
        '[NX3;H2:1].[NX3:2][CX4:3][CX3:4](=O)[OX2]>>[N:1]C(=O)[N:2][C:3][C:4](=O)',

    # ── Biginelli ─────────────────────────────────────────────────────────
    'biginelli':
        '[N&X3;H2:1]C(=[O,S])[N&X3;H2:2].[CX3H1:3](=O).[CX4H2:4]([CX3]=O)>>[N:1]C(=O)[N:2][C:3][C:4]',

    # ── Cyclic guanidine ──────────────────────────────────────────────────
    'cyclic_guanidine':
        '[NX3;H2,H1:1].[CX3:2](=S)[NX3:3]>>[N:1][C:2](=[N])[N:3]',

    # ── Tetrazole ─────────────────────────────────────────────────────────
    'tetrazole_acylation':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'tetrazole_alkylation':
        '[#7&X3;H2,H1;+0:5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[#7&X3;+0:5]-[C&X4;+0:1]',
    'tetrazole_sulfonamide':
        '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',

    # ── Ugi MCR ───────────────────────────────────────────────────────────
    'ugi_4cr':
        '[#6:1](=[#8:2])-[#8;H1].[#7;H2,H1:3].[CX3H1:4](=O).[CX1-:5]#[NX2+:6]>>[#6:1](=[#8:2])-[#7:3]-[C:4]-[C:5](=[O])-[N:6]',
    'ugi_3cr':
        '[#7;H2,H1:3].[CX3H1:4](=O).[CX1-:5]#[NX2+:6]>>[#7:3]-[C:4]-[C:5](=[O])-[N:6]',

    # ── GBB MCR ───────────────────────────────────────────────────────────
    'gbb':
        '[#7;H2,H1:3].[CX3H1:4](=O).[CX1-:5]#[NX2+:6]>>[#7:3]-[C:4]-[C:5](=[O])-[N:6]',

    # ── Multi-component (first bond-forming step) ─────────────────────────
    'multi_amide_arylation':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'multi_amide_amide':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'multi_sulfonamide_arylation':
        '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
    'multi_sulfonamide_suzuki':
        '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
    'multi_amide_suzuki':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'multi_amide_click':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'multi_alkylation_click':
        '[#7&X3;H2,H1;+0:5].[Cl,Br,I]-[C;H2,H1;+0:1]>>[#7&X3;+0:5]-[C&X4;+0:1]',
    'multi_amide_alkylation':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'multi_reductive_amination_amide':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'multi_reductive_amination_arylation':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'multi_reductive_amination_sulfonamide':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'multi_reductive_amination_urea':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'multi_reductive_amination_carbamate':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'multi_reductive_amination_reductive_amination':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'multi_reductive_amination_acyl_halide':
        '[#6:2](=[#8])(-[#6:1]).[#7;H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]',
    'multi_arylation_arylation':
        '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',
    'multi_arylation_suzuki':
        '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',
    'multi_triazole_amide':
        '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
    'multi_triazole_urea':
        '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
    'multi_triazole_alkylation':
        '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
    'multi_triazole_sulfonamide':
        '[OX2:1][CX3:2]=[NX2:3].[NX3;H2,H1:4]NC(=O)>>[n:3]1:[n:4]:[n]:[c:2]:[c]:1',
    'multi_thiazole_amide':
        '[a,A:1]-[C](=[S:4])[N&H2].[a,A:3]-[C:5](=[O])-[C;H1:6](-[$([Cl,Br,I]),$(OS(=O)(=O)C),$(OS(=O)(=O)c1ccc(C)cc1)])-[a,A:2]>>[c]-1([a,A:1])-[s:4]-[c:6](-[a,A:2])-[c:5](-[a,A:3])-[n]-1',
    'multi_thiazole_sulfonamide':
        '[a,A:1]-[C](=[S:4])[N&H2].[a,A:3]-[C:5](=[O])-[C;H1:6](-[$([Cl,Br,I]),$(OS(=O)(=O)C),$(OS(=O)(=O)c1ccc(C)cc1)])-[a,A:2]>>[c]-1([a,A:1])-[s:4]-[c:6](-[a,A:2])-[c:5](-[a,A:3])-[n]-1',
    'multi_anhydride_amide_arylation':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'multi_sulfonamide_amide':
        '[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])[#17,#9,#35].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):5]>>[#16X4:1](=[OX1:2])(=[OX1:3])([#6,#7:4])-[#7:5]',
    'multi_acylsulfamide_alkylation':
        '[#6:1](=[#8:2])-[#8;H1].[#7;H1:3]S(=O)(=O)>>[#6:1](=[#8:2])-[#7:3]',
    'multi_click_amide':
        '[NX1:1]~[NX2:2]~[NX2:3]-[#6:4].[CX2:5]#[CX2:6]-[#6:7]>>[#6:4]-[n:3]1:[n:2]:[n:1]:[c:5](-[#6:7]):[c:6]:1',
    'multi_oxoisoindole_amide':
        '[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]',
    'dihydrouracil':
        '[N&X3;H2,H1;!$(NC=O);!$(NS):3].[N&X3;H2,H1;!$(NC=O);!$(NS):4]>>[#7&X3:3]-[#6](=[#8])-[#7&X3:4]',
    'purine_substitution':
        '[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]',

    # ── Generic fallbacks (no SMIRKS) ─────────────────────────────────────
    'generic_2comp': '',
    'generic_3comp': '',
    'generic_4comp': '',
}


# ═══════════════════════════════════════════════════════════════════════════
# Mapping: Enamine Reaction Name → Reaction Class
# ═══════════════════════════════════════════════════════════════════════════

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
    'ACETYLATION':        'acetylation',
    'R111-ACYLATION':     'esterification',
    'R017-ACRYLAMIDE':    'acrylamide',
    'ANHYDRIDE1':         'anhydride',
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
    'SULFACYLATION-O':    'sulfonamide_phenol',
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
    'URETANE':            'urea_phenyl_carbamate',
    'URETANE1':           'urea',
    'URETANE2':           'urea_phenyl_carbamate',
    'URETANE4':           'urea',
    'URETANE-T':          'urea',
    'Boc-URETANE':        'urea_phenyl_carbamate',
    'Boc-URETANE1':       'urea',
    'Boc-URETANE4':       'urea',
    'B-URETANE1':         'urea',
    'B-URETANE4':         'urea',
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
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _desalt_smiles(smiles):
    """
    Remove salt/counter-ion fragments and neutralise common charged species.

    1. Keep the largest fragment by heavy-atom count (removes Cl, Li+, Na+,
       etc.).
    2. Neutralise carboxylate [O-] → OH (from lithium/sodium salt forms)
       so that acid-based SMIRKS can match.

    Returns the raw fragment substring (not re-canonicalised) when no
    neutralisation is needed, matching ``extract-enamine-building-blocks``
    behaviour.
    """
    if not smiles or str(smiles).strip() in ('', 'nan'):
        return smiles
    smi = str(smiles).strip()
    frags = smi.split('.')
    if len(frags) == 1:
        best = smi
    else:
        best, best_size = frags[0], -1
        for frag in frags:
            mol = Chem.MolFromSmiles(frag)
            if mol is None:
                continue
            n = mol.GetNumHeavyAtoms()
            if n > best_size:
                best, best_size = frag, n

    # Neutralise carboxylate [O-] → OH
    mol = Chem.MolFromSmiles(best)
    if mol is None:
        return best
    changed = False
    for atom in mol.GetAtoms():
        # Carboxylate oxygen: [O-] bonded to C(=O)
        if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == -1:
            for nb in atom.GetNeighbors():
                if nb.GetSymbol() == 'C':
                    for nb2 in nb.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(nb.GetIdx(),
                                                       nb2.GetIdx())
                        if (nb2.GetIdx() != atom.GetIdx()
                                and nb2.GetSymbol() == 'O'
                                and bond is not None
                                and bond.GetBondTypeAsDouble() == 2.0):
                            atom.SetFormalCharge(0)
                            atom.SetNumExplicitHs(1)
                            changed = True
                            break
    if changed:
        try:
            Chem.SanitizeMol(mol)
            return Chem.MolToSmiles(mol)
        except Exception:
            return best
    return best


def _canonical_smiles(smiles):
    """Return canonical SMILES or None if unparseable."""
    if not smiles or str(smiles).strip() in ('', 'nan'):
        return None
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def _safe_str(val):
    """Return cleaned string or empty string for NaN / None."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ''
    return str(val).strip()


def get_smirks_for_reaction(reaction_name):
    """
    Look up the SMIRKS for an Enamine reaction name.
    Returns (reaction_class, smirks_string).
    """
    rxn_class = ENAMINE_REACTION_TO_CLASS.get(reaction_name, '')
    if not rxn_class:
        return ('', '')
    smirks = REACTION_CLASS_SMIRKS.get(rxn_class, '')
    return (rxn_class, smirks)


# ── Reaction class refinement based on actual reactants ──────────────────

REFINED_CLASS_DISPLAY = {
    'stille': 'STILLE',
    'urea_carbamoyl_chloride': 'UREA_CARBAMOYL_CHLORIDE',
}

_CARBAMOYL_CHLORIDE_PAT = Chem.MolFromSmarts('[#7]-C(=O)-[Cl]')


def _refine_reaction_class(rxn_class, r1_smiles, r2_smiles):
    """
    Detect reactant-specific sub-classes that share an Enamine reaction
    name but require different SMIRKS.

    Returns the refined reaction class (unchanged if no refinement needed).
    """
    if rxn_class == 'urea_isocyanate':
        for smi in (r1_smiles, r2_smiles):
            if not smi or smi == 'nan':
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol and mol.HasSubstructMatch(_CARBAMOYL_CHLORIDE_PAT):
                return 'urea_carbamoyl_chloride'

    if rxn_class == 'suzuki':
        for smi in (r1_smiles, r2_smiles):
            if not smi or smi == 'nan':
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol and any(a.GetAtomicNum() == 50 for a in mol.GetAtoms()):
                return 'stille'

    return rxn_class


# ── Deprotection transforms for SMIRKS validation ────────────────────────

_DEPROTECTION_RXNS = None


def _get_deprotection_rxns():
    """Lazily initialise deprotection reaction objects."""
    global _DEPROTECTION_RXNS
    if _DEPROTECTION_RXNS is None:
        _DEPROTECTION_RXNS = []
        # Order matters: Boc removal before tBu ester hydrolysis
        smirks_list = [
            # Boc removal from aliphatic nitrogen: N-C(=O)-O-C(C)(C)C → NH
            '[N:1]C(=O)OC(C)(C)C>>[N:1]',
            # Boc removal from aromatic nitrogen (indole, pyrazole, etc.)
            '[n:1]C(=O)OC(C)(C)C>>[nH:1]',
            # MOM (methoxymethyl) removal from aromatic nitrogen
            '[n:1]COC>>[nH:1]',
            # tert-Butyl ester hydrolysis → free acid
            '[C:1](=[O:2])OC(C)(C)C>>[C:1](=[O:2])O',
            # Methyl ester hydrolysis → free acid
            '[C:1](=[O:2])O[CH3]>>[C:1](=[O:2])O',
        ]
        for s in smirks_list:
            rxn = AllChem.ReactionFromSmarts(s)
            if rxn:
                _DEPROTECTION_RXNS.append(rxn)
    return _DEPROTECTION_RXNS


def _deprotect_mol(mol):
    """Apply common deprotection transforms exhaustively to a molecule."""
    if mol is None:
        return None
    for rxn in _get_deprotection_rxns():
        for _ in range(10):  # max iterations per transform
            try:
                products = rxn.RunReactants((mol,))
            except Exception:
                break
            if not products:
                break
            try:
                new_mol = products[0][0]
                Chem.SanitizeMol(new_mol)
                mol = new_mol
            except Exception:
                break
    return mol


def _validate_smirks(smirks, r1, r2, target_smi, r3=None, r4=None,
                     found_smi=None):
    """
    Validate that the SMIRKS produces the expected target from reactants.

    For 2-component reactions: checks product == target (exact or after
    common deprotection transforms for Boc, tBu ester, methyl ester).
    Also checks against ``found_smi`` (Enamine's canonical SMILES) which
    may differ from target_smi due to stereochemistry normalisation.
    For multi-component: checks that step-1 fires on at least one
    reactant pair.

    Returns (status, product_smiles) where status is one of:
      'pass'              - exact product match
      'pass_found'        - matched Enamine found_smiles (stereo variant)
      'pass_deprotected'  - match after deprotection
      'mismatch'          - products generated but none match target
      'no_products'       - SMIRKS did not fire on any reactant pair
      'no_smirks'         - no SMIRKS available for this reaction
      'multi_step1_pass'  - multi-component: step-1 fired
      'multi_step1_fail'  - multi-component: step-1 did not fire
      'error'             - bad input (unparseable SMILES, etc.)
    """
    if not smirks:
        return ('no_smirks', '')

    rxn = AllChem.ReactionFromSmarts(smirks)
    if rxn is None:
        return ('error', '')

    is_multi = bool(r3 or r4)

    # Collect parseable reactant molecules
    smi_list = [s for s in [r1, r2] if s and s not in ('', 'nan')]
    if is_multi:
        for s in [r3, r4]:
            if s and s not in ('', 'nan'):
                smi_list.append(s)

    mols = {}
    for smi in smi_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mols[smi] = mol

    if len(mols) < 2:
        return ('error', '')

    # Run reaction on all reactant pairs (both orders)
    all_products = set()
    smi_keys = list(mols.keys())
    for i, smi_a in enumerate(smi_keys):
        for j, smi_b in enumerate(smi_keys):
            if i == j:
                continue
            if rxn.GetNumReactantTemplates() != 2:
                continue
            try:
                outcomes = rxn.RunReactants((mols[smi_a], mols[smi_b]))
                for product_set in outcomes:
                    for prod in product_set:
                        try:
                            Chem.SanitizeMol(prod)
                            all_products.add(Chem.MolToSmiles(prod))
                        except Exception:
                            pass
            except Exception:
                pass

    if not all_products:
        if is_multi:
            return ('multi_step1_fail', '')
        return ('no_products', '')

    if is_multi:
        return ('multi_step1_pass', sorted(all_products)[0])

    # ── 2-component: compare product to target ────────────────────────
    target_canon = _canonical_smiles(target_smi)
    if target_canon is None:
        return ('error', '')

    # Exact match
    if target_canon in all_products:
        return ('pass', target_canon)

    # Match against Enamine's found_smiles (may differ in stereochemistry)
    if found_smi:
        found_canon = _canonical_smiles(found_smi)
        if found_canon and found_canon in all_products:
            return ('pass_found', found_canon)

    # Match after deprotection
    for prod_smi in all_products:
        prod_mol = Chem.MolFromSmiles(prod_smi)
        deprot = _deprotect_mol(prod_mol)
        if deprot is not None:
            deprot_smi = Chem.MolToSmiles(deprot)
            if deprot_smi == target_canon:
                return ('pass_deprotected', prod_smi)
            # Also try deprotected product vs found_smiles
            if found_smi:
                found_canon = _canonical_smiles(found_smi)
                if found_canon and deprot_smi == found_canon:
                    return ('pass_deprotected', prod_smi)

    # ── Fuzzy fallbacks: tautomer-canonical and stereo-insensitive ────
    target_mol = Chem.MolFromSmiles(target_canon)
    if target_mol is None:
        return ('mismatch', sorted(all_products)[0])

    # Tautomer-canonical comparison (resolves pyrazole [nH]n vs n[nH] etc.)
    try:
        target_taut = Chem.MolToSmiles(_TAUT_ENUMERATOR.Canonicalize(target_mol))
    except Exception:
        target_taut = None

    if target_taut:
        for prod_smi in all_products:
            try:
                pm = Chem.MolFromSmiles(prod_smi)
                if pm and Chem.MolToSmiles(_TAUT_ENUMERATOR.Canonicalize(pm)) == target_taut:
                    return ('pass_tautomer', prod_smi)
                # Also try after deprotection + tautomer
                dm = _deprotect_mol(pm)
                if dm and Chem.MolToSmiles(_TAUT_ENUMERATOR.Canonicalize(dm)) == target_taut:
                    return ('pass_deprotected', prod_smi)
            except Exception:
                continue

    # Stereo-insensitive comparison (resolves E/Z unspecified vs specified)
    try:
        t_nostereo = Chem.RWMol(target_mol)
        Chem.RemoveStereochemistry(t_nostereo)
        target_nostereo = Chem.MolToSmiles(t_nostereo)
    except Exception:
        target_nostereo = None

    if target_nostereo:
        for prod_smi in all_products:
            try:
                pm = Chem.MolFromSmiles(prod_smi)
                if pm:
                    p_ns = Chem.RWMol(pm)
                    Chem.RemoveStereochemistry(p_ns)
                    if Chem.MolToSmiles(p_ns) == target_nostereo:
                        return ('pass_stereo', prod_smi)
                    dm = _deprotect_mol(pm)
                    if dm:
                        d_ns = Chem.RWMol(dm)
                        Chem.RemoveStereochemistry(d_ns)
                        if Chem.MolToSmiles(d_ns) == target_nostereo:
                            return ('pass_deprotected', prod_smi)
            except Exception:
                continue

    return ('mismatch', sorted(all_products)[0])


# ═══════════════════════════════════════════════════════════════════════════
# Core: merge auto + enamine → manual CSV
# ═══════════════════════════════════════════════════════════════════════════

def merge_to_manual(auto_df, enamine_df):
    """
    Merge a Syndirella auto input DataFrame with an Enamine REAL Tools
    output DataFrame.  Returns a DataFrame in Syndirella manual CSV
    format, plus a set of (reaction_name, smirks) used.

    Matching is done by canonical SMILES.  For each target in the auto
    CSV, the first FOUND route from Enamine (method_no=1, or first
    available) is selected.
    """
    # Canonicalise SMILES for matching
    auto_df = auto_df.copy()
    auto_df['_canon'] = auto_df['smiles'].apply(_canonical_smiles)

    enamine_found = enamine_df[enamine_df['status'] == 'FOUND'].copy()
    enamine_found['_canon'] = enamine_found['target_smiles'].apply(
        _canonical_smiles
    )

    # Pick the first route per target (lowest method_no, then first row)
    if 'method_no' in enamine_found.columns:
        enamine_found = enamine_found.sort_values('method_no')
    enamine_first = enamine_found.drop_duplicates(subset='_canon', keep='first')

    # Build output rows
    manual_rows = []
    used_reactions = {}  # reaction_name → smirks
    unmatched = []
    reclassified = []    # (original_name, new_name, refined_class, smiles)

    for _, auto_row in auto_df.iterrows():
        canon = auto_row['_canon']
        if canon is None:
            unmatched.append(_safe_str(auto_row['smiles']))
            continue

        match = enamine_first[enamine_first['_canon'] == canon]
        if match.empty:
            unmatched.append(_safe_str(auto_row['smiles']))
            continue

        en = match.iloc[0]
        rxn_name = _safe_str(en.get('reaction_name', ''))
        r1 = _desalt_smiles(_safe_str(en.get('reactant_1_smiles', '')))
        r2 = _desalt_smiles(_safe_str(en.get('reactant_2_smiles', '')))
        r3 = _desalt_smiles(_safe_str(en.get('reactant_3_smiles', '')))
        r4 = _desalt_smiles(_safe_str(en.get('reactant_4_smiles', '')))

        # Look up SMIRKS; refine class based on actual reactants
        rxn_class, smirks = get_smirks_for_reaction(rxn_name)
        refined_class = _refine_reaction_class(rxn_class, r1, r2)
        if refined_class != rxn_class:
            smirks = REACTION_CLASS_SMIRKS.get(refined_class, smirks)
            effective_name = REFINED_CLASS_DISPLAY.get(
                refined_class, rxn_name)
            reclassified.append(
                (rxn_name, effective_name, refined_class,
                 _safe_str(auto_row['smiles'])))
            rxn_name = effective_name
        if rxn_name and smirks:
            used_reactions[rxn_name] = smirks

        # For single-step (2-component): populate step1 only.
        # For multi-component: step1 gets first 2 reactants, step2 third,
        # step3 fourth — but product columns are left empty because
        # intermediates are not in the Enamine output.
        row = {
            'smiles': _safe_str(auto_row['smiles']),
            'reaction_name_step1': rxn_name,
            'reactant_step1': r1,
            'reactant2_step1': r2,
            'product_step1': '',
            'reaction_name_step2': '',
            'reactant_step2': r3 if r3 else '',
            'product_step2': '',
            'reaction_name_step3': '',
            'reactant_step3': r4 if r4 else '',
            'hit1': _safe_str(auto_row.get('hit1', '')),
            'hit2': _safe_str(auto_row.get('hit2', '')),
            'hit3': _safe_str(auto_row.get('hit3', '')),
            'template': _safe_str(auto_row.get('template', '')),
            'compound_set': 'syndirella_manual',
        }

        # ── Validate SMIRKS against target ─────────────────────────────
        val_status, val_product = _validate_smirks(
            smirks, r1, r2,
            _safe_str(auto_row['smiles']),
            r3=r3 if r3 else None,
            r4=r4 if r4 else None,
            found_smi=_safe_str(en.get('found_smiles', '')),
        )
        row['smirks_validated'] = val_status
        row['smirks_product'] = val_product

        manual_rows.append(row)

    manual_df = pd.DataFrame(manual_rows, columns=[
        'smiles',
        'reaction_name_step1', 'reactant_step1', 'reactant2_step1',
        'product_step1',
        'reaction_name_step2', 'reactant_step2', 'product_step2',
        'reaction_name_step3', 'reactant_step3',
        'hit1', 'hit2', 'hit3', 'template', 'compound_set',
        'smirks_validated', 'smirks_product',
    ])

    return manual_df, used_reactions, unmatched, reclassified


# ═══════════════════════════════════════════════════════════════════════════
# SMIRKS JSON output
# ═══════════════════════════════════════════════════════════════════════════

def build_smirks_json(used_reactions, reclassified=None):
    """
    Build a dict suitable for merging into Syndirella's
    RXN_SMIRKS_CONSTANTS.json.  Keyed by Enamine reaction name.
    """
    # Build reverse map: refined_name → original Enamine name
    reclassed_origins = {}
    if reclassified:
        for orig, new, _cls, _smi in reclassified:
            reclassed_origins.setdefault(new, orig)

    output = {}
    for rxn_name, smirks in sorted(used_reactions.items()):
        entry = {
            'smirks': smirks,
            'source': 'enamine_real_tools',
            'type': 'parent',
        }
        if rxn_name in reclassed_origins:
            entry['original_enamine_name'] = reclassed_origins[rxn_name]
        output[rxn_name] = entry
    return output


# ═══════════════════════════════════════════════════════════════════════════
# Enamine REAL Tools API integration
# ═══════════════════════════════════════════════════════════════════════════

def call_enamine_api(smiles_list, reactions_file=None,
                     batch_size=1000, reaction_types=None,
                     fetch_synthons=True):
    """
    Call the Enamine REAL Tools API for a list of SMILES and return a
    DataFrame with the same column schema as ``enamine_real_tools.py``
    produces (target_smiles, status, reaction_name, reactant_*_smiles, …).

    Imports are intentionally lazy so the module can be used in offline
    mode without the ``ENAMINE_REAL_TOOLS_API_KEY`` env-var being set.
    """
    # Lazy imports — only needed when the API is actually called
    try:
        from enamine_real_tools import (
            process_targets_batch, collect_all_synthon_ids,
            fetch_synthon_smarts_batch, create_output_dataframe,
            load_reaction_mapping,
        )
    except ImportError:
        from .enamine_real_tools import (
            process_targets_batch, collect_all_synthon_ids,
            fetch_synthon_smarts_batch, create_output_dataframe,
            load_reaction_mapping,
        )

    # Default reactions file is in the same directory as this script
    if reactions_file is None:
        reactions_file = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'enamine-real-tools-rxns.csv',
        )

    reaction_df = load_reaction_mapping(reactions_file)

    # Process in batches
    all_results = []
    for i in range(0, len(smiles_list), batch_size):
        batch = smiles_list[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(smiles_list) + batch_size - 1) // batch_size
        print(f"  API batch {batch_num}/{total_batches} "
              f"({len(batch)} targets) …")
        batch_results = process_targets_batch(
            batch, reaction_df, reaction_types,
        )
        all_results.extend(batch_results)

    # Fetch synthon SMARTS (optional — only needed by extract-enamine-building-blocks)
    synthon_smarts_map = {}
    if fetch_synthons:
        all_synthon_ids = collect_all_synthon_ids(all_results)
        if all_synthon_ids:
            print(f"  Fetching SMARTS for {len(all_synthon_ids)} synthons …")
            synthon_smarts_map = fetch_synthon_smarts_batch(all_synthon_ids)
    else:
        print("  Skipping synthon SMARTS fetch (--no-synthons)")

    # Build output DataFrame (same schema as enamine_real_tools.py)
    output_df = create_output_dataframe(
        all_results, reaction_df, synthon_smarts_map,
    )
    return output_df


# ═══════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════

def print_summary(auto_df, manual_df, used_reactions, unmatched,
                  reclassified=None):
    n_auto = len(auto_df)
    n_manual = len(manual_df)
    n_unmatched = len(unmatched)

    print(f"\n{'='*70}")
    print(f"PREPARE SYNDIRELLA INPUT — SUMMARY")
    print(f"{'='*70}")
    print(f"Syndirella auto input  : {n_auto} compounds")
    print(f"Matched to Enamine     : {n_manual}")
    print(f"Unmatched (no route)   : {n_unmatched}")
    if n_auto > 0:
        print(f"Match rate             : {n_manual/n_auto*100:.1f}%")

    if n_manual > 0:
        # Reaction type breakdown
        rxn_counts = manual_df['reaction_name_step1'].value_counts()
        print(f"\n--- Enamine Reaction Types ---")
        for rxn, count in rxn_counts.items():
            if rxn:
                print(f"  {rxn:40s}: {count:4d}")

        # Multi-step rows
        has_step2 = (manual_df['reactant_step2'].astype(str).str.strip() != '')
        n_multi = has_step2.sum()
        if n_multi > 0:
            print(f"\n  Multi-component reactions (3+ reactants): {n_multi}")

    # Reclassified reactions
    if reclassified:
        from collections import Counter
        reclass_counts = Counter(
            (orig, new) for orig, new, _cls, _smi in reclassified)
        print(f"\n--- Reclassified Reactions ---")
        for (orig, new), count in reclass_counts.most_common():
            print(f"  {orig:30s} → {new:30s}: {count:4d}")

    # SMIRKS validation results
    if 'smirks_validated' in manual_df.columns and n_manual > 0:
        val_counts = manual_df['smirks_validated'].value_counts()
        print(f"\n--- SMIRKS Validation ---")
        for status, count in val_counts.items():
            print(f"  {status:25s}: {count:4d}")
        n_pass = (val_counts.get('pass', 0)
                  + val_counts.get('pass_found', 0)
                  + val_counts.get('pass_deprotected', 0)
                  + val_counts.get('pass_tautomer', 0)
                  + val_counts.get('pass_stereo', 0))
        n_testable = n_pass + val_counts.get('mismatch', 0) + val_counts.get('no_products', 0)
        if n_testable > 0:
            print(f"  2-comp pass rate       : {n_pass}/{n_testable} "
                  f"({n_pass/n_testable*100:.1f}%)")

    print(f"\n--- SMIRKS for Syndirella JSON ---")
    print(f"  Unique reaction names with SMIRKS: {len(used_reactions)}")

    if unmatched:
        print(f"\n--- Unmatched SMILES (first 10) ---")
        for smi in unmatched[:10]:
            print(f"  {smi[:80]}")

    print(f"{'='*70}")


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description='Prepare Syndirella manual input CSV from a Syndirella '
                    'auto input.  Calls the Enamine REAL Tools API for '
                    'retrosynthesis unless --enamine-csv is given.'
    )
    parser.add_argument(
        'auto_csv',
        help='Syndirella automatic input CSV '
             '(columns: smiles, hit1, hit2, hit3, template, compound_set)'
    )
    # ── Offline shortcut ─────────────────────────────────────────────────
    parser.add_argument(
        '--enamine-csv', default=None, metavar='FILE',
        help='Skip the API call and use this pre-existing Enamine REAL '
             'Tools output CSV instead.'
    )
    # ── Output paths ─────────────────────────────────────────────────────
    parser.add_argument(
        '-e', '--enamine-output', default=None, metavar='FILE',
        help='Path for the Enamine REAL Tools CSV output '
             '(default: <auto_basename>_enamine_real_tools.csv)'
    )
    parser.add_argument(
        '-o', '--output', default=None, metavar='FILE',
        help='Output Syndirella manual CSV path '
             '(default: <auto_basename>_syndirella_manual.csv)'
    )
    parser.add_argument(
        '-j', '--json', default=None, metavar='FILE',
        help='Output SMIRKS JSON path '
             '(default: <auto_basename>_smirks.json)'
    )
    # ── API options ──────────────────────────────────────────────────────
    parser.add_argument(
        '--batch-size', type=int, default=1000,
        help='Batch size for API requests (default: 1000)'
    )
    parser.add_argument(
        '--reactions-file', default=None, metavar='FILE',
        help='Reaction mapping CSV '
             '(default: enamine-real-tools-rxns.csv in script directory)'
    )
    parser.add_argument(
        '--reaction-types', nargs='+', type=int, default=None,
        help='Reaction type codes to search (default: 0 1050)'
    )
    parser.add_argument(
        '--no-synthons', action='store_true', default=False,
        help='Skip fetching synthon SMARTS from the API (faster; '
             'only needed for extract-enamine-building-blocks)'
    )
    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.auto_csv):
        print(f"ERROR: File not found: {args.auto_csv}")
        sys.exit(1)
    if args.enamine_csv and not os.path.exists(args.enamine_csv):
        print(f"ERROR: File not found: {args.enamine_csv}")
        sys.exit(1)

    # Default output paths
    base = os.path.splitext(args.auto_csv)[0]
    output_csv = args.output or f"{base}_syndirella_manual.csv"
    output_json = args.json or f"{base}_smirks.json"
    enamine_output = args.enamine_output or f"{base}_enamine_real_tools.csv"

    # ── Read auto input ──────────────────────────────────────────────────
    auto_df = pd.read_csv(args.auto_csv)
    print(f"Syndirella auto input : {os.path.abspath(args.auto_csv)}")
    print(f"  {len(auto_df)} compounds")

    # ── Obtain Enamine REAL Tools results ────────────────────────────────
    if args.enamine_csv:
        # Offline mode
        enamine_df = pd.read_csv(args.enamine_csv)
        print(f"Enamine CSV (offline) : {os.path.abspath(args.enamine_csv)}")
        print(f"  {len(enamine_df)} rows")
    else:
        # API mode
        if not os.environ.get('ENAMINE_REAL_TOOLS_API_KEY'):
            print("ERROR: ENAMINE_REAL_TOOLS_API_KEY env var not set.\n"
                  "Set it or supply --enamine-csv to use offline mode.")
            sys.exit(1)

        smiles_list = auto_df['smiles'].dropna().unique().tolist()
        print(f"\nCalling Enamine REAL Tools API for "
              f"{len(smiles_list)} unique SMILES …")
        enamine_df = call_enamine_api(
            smiles_list,
            reactions_file=args.reactions_file,
            batch_size=args.batch_size,
            reaction_types=args.reaction_types,
            fetch_synthons=not args.no_synthons,
        )
        # Save the REAL Tools output
        enamine_df.to_csv(enamine_output, index=False)
        print(f"\nEnamine CSV written   : {os.path.abspath(enamine_output)}")
        print(f"  {len(enamine_df)} rows")

    found_count = (
        (enamine_df['status'] == 'FOUND').sum()
        if 'status' in enamine_df.columns else 0
    )
    print(f"  FOUND rows          : {found_count}")

    # ── Merge and produce outputs ────────────────────────────────────────
    manual_df, used_reactions, unmatched, reclassified = merge_to_manual(
        auto_df, enamine_df)
    print_summary(auto_df, manual_df, used_reactions, unmatched, reclassified)

    # Write manual CSV
    manual_df.to_csv(output_csv, index=False)
    print(f"\nManual CSV written    : {os.path.abspath(output_csv)}")
    print(f"  {len(manual_df)} rows")

    # Write SMIRKS JSON
    smirks_dict = build_smirks_json(used_reactions, reclassified)
    with open(output_json, 'w') as f:
        json.dump(smirks_dict, f, indent=2)
    print(f"SMIRKS JSON written   : {os.path.abspath(output_json)}")
    print(f"  {len(smirks_dict)} reaction entries")


if __name__ == '__main__':
    main()
