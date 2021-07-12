encoded_recipes = {
    "Amidation": {
        "reactionSMARTS": ["[#6:1](=[#8:2])-[#8].[A#7:3]>>[#6:1](=[#8:2])-[A#7:3]"],
        "recipe": [
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
                        "solvent": None,
                        "density": 0.79,
                        "concentration": None,
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
    },
    "Amide schotten - baumann": {
        "reactionSMARTS": ["[#6:1](=[#8:2])-[#17].[A#7:3]>>[#6:1](=[#8:2])-[A#7:3]"],
        "recipe": [
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
                        "solvent": None,
                        "density": 0.79,
                        "concentration": None,
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
    },
    "Reductive amination": {
        "reactionSMARTS": [
            "[#7:1].[#6:2]=O>>[CH1:2]-[#7R0:1]",
            "[#7:1].[#6:2]=O>>[CH2R0:2]-[#7:1]",
            "[#7:1].[#6:2]=O>>[CH1R0:2]-[#7:1]",
            "[#7:1].[#6:2]=O>>[CH2:2]-[#7R0:1]",
            "[#7:1].[#6:2]=O>>[CH3:2]-[#7R0:1]",
        ],
        "recipe": [
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
    },
    "N-nucleophilic aromatic substitution": {
        "reactionSMARTS": ["[#6:3]-[#7:1].[c:2]Br>>[c:2][N:1][#6:3]"],
        "recipe": [
            {
                "name": "add",
                "content": {
                    "action_no": 1,
                    "material": {
                        # What about looking for other nucelophiles? Check with Harry -> SMARTS above also not very general -> only
                        # looking for amines
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
    },
}
