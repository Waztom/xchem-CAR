# Encoded recipes

encoded_recipes = {
    "Amidation": [
        {
            "name": "add",
            "content": {
                "action_no": 1,
                "material": {
                    "SMARTS": ["[CX3](=O)[OX2H1]"],
                    "SMILES": None,
                    "quantity": {"value": 1.0, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 2,
                "material": {
                    "SMARTS": None,
                    "SMILES": "CCCP1(=O)OP(=O)(OP(=O)(O1)CCC)CCC",
                    "quantity": {"value": 1.2, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 3,
                "material": {
                    "SMARTS": None,
                    "SMILES": "CCN(C(C)C)C(C)C",
                    "quantity": {"value": 3.5, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 4,
                "material": {
                    "SMARTS": ["[NX3;H1,H2;!$(NC=O)]"],
                    "SMILES": None,
                    "quantity": {"value": 1.1, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "stir",
            "content": {
                "action_no": 5,
                "temperature": 25,  # degrees celcius
                "duration": {"value": 12, "unit": "hours"},  # in hours
            },
        },
    ],
    "Amide schotten - baumann": [
        {
            "name": "add",
            "content": {
                "action_no": 1,
                "material": {
                    "SMARTS": ["[NX3;H1,H2;!$(NC=O)]"],
                    "SMILES": None,
                    "quantity": {"value": 1.1, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 2,
                "material": {
                    "SMARTS": ["[CX3](=[OX1])[Cl]"],
                    "SMILES": None,
                    "quantity": {"value": 1.0, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 3,
                "material": {
                    "SMARTS": None,
                    "SMILES": "CCN(C(C)C)C(C)C",
                    "quantity": {"value": 3.5, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "stir",
            "content": {
                "action_no": 4,
                "temperature": 25,  # degrees celcius
                "duration": {"value": 12, "unit": "hours"},  # in hours
            },
        },
    ],
    "Reductive amination": [
        {
            "name": "add",
            "content": {
                "action_no": 1,
                "material": {
                    "SMARTS": ["[#6][CX3](=O)[#6]", "[CX3H1](=O)[#6]"],  # check if works
                    "SMILES": None,
                    "quantity": {"value": 1.0, "unit": "moleq"},
                    "solvent": "ACN",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 2,
                "material": {
                    "SMARTS": ["[NX3;H1,H2;!$(NC=O)]"],
                    "SMILES": None,
                    "quantity": {"value": 1.0, "unit": "moleq"},
                    "solvent": "ACN",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 3,
                "material": {
                    "SMARTS": None,
                    "SMILES": "[Na+].CC(=O)O[BH-](OC(C)=O)OC(C)=O",
                    "quantity": {"value": 1.4, "unit": "moleq"},
                    "solvent": "MeOH",
                    "concentration": 0.25,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 4,
                "material": {
                    "SMARTS": None,
                    "SMILES": "CC(=O)O",
                    "quantity": {"value": 1.0, "unit": "moleq"},
                    "solvent": "ACN",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "stir",
            "content": {
                "action_no": 5,
                "temperature": 25,  # degrees celcius
                "duration": {"value": 12, "unit": "hours"},  # in hours
            },
        },
    ],
    "N-nucelophilic aromatic substitution": [
        {
            "name": "add",
            "content": {
                "action_no": 1,
                "material": {
                    "SMARTS": ["[c][F,Cl,Br,I]", "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"],
                    # includes halides & NO2 not specifically attached to aromatic c as unlikely a Nu has NO2 attached
                    "SMILES": None,  # leaving group question eg OTs
                    "quantity": {"value": 1, "unit": "moleq"},
                    "solvent": "MeOH",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "add",
            "content": {
                "action_no": 2,
                "material": {
                    "SMARTS": [
                        "[NX3;H2,H1;!$(NC=O)]"
                    ],  # allowing both 1' & 2' amine to be Nu but not 3' due to sterics
                    "SMILES": None,
                    "quantity": {"value": 1.2, "unit": "moleq"},
                    "solvent": "MeOH",
                    "concentration": 0.5,
                },
            },
        },
        {
            "name": "stir",
            "content": {
                "action_no": 3,
                "temperature": 50,  # generic elevated temperature
                "duration": {"value": 12, "unit": "hours"},  # in hours
            },
        },
    ],
}
