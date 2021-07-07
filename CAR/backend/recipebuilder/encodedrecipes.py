# Encoded recipes

encoded_recipes = {
    "Amidation": [
        {
            "name": "add",
            "content": {
                "action_no": 1,
                "material": {
                    "SMARTS": "[CX3](=O)[OX2H1]",
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
                    "SMARTS": "[NX3;H1,H2;!$(NC=O)]",
                    "SMILES": None,
                    "quantity": {"value": 1.1, "unit": "moleq"},
                    "solvent": "DMA",
                    "concentration": 0.5,
                },
            },
        },
    ]
}