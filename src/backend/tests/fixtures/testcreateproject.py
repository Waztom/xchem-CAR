import pandas as pd

uploaded_df = pd.DataFrame(
    columns=["targets", "amount-required-mg", "batch-tag"],
    data=[
        [
            "C[NH2+]CC(Nc1nccc(C[NH2+]CC[C@@H](C([O-])=O)NC(C2CCCC2)=O)c1)=O",
            "O=C(N1CCCCC1)c1ccnn1C",
        ],
        [10, 10],
        ["multi-step", "single-step"],
    ],
)
grouped_targets = uploaded_df.groupby("batch-tag")

routes = [
    {
        "molecules": [
            {
                "smiles": "CNCC(=O)Nc1cc(CNCCC(NC(=O)C2CCCC2)C(=O)O)ccn1",
                "isBuildingBlock": False,
                "catalogEntries": [],
            },
            {
                "smiles": "NCCC(NC(=O)C1CCCC1)C(=O)O",
                "isBuildingBlock": False,
                "catalogEntries": [],
            },
            {
                "smiles": "NCCC(N)C(=O)O",
                "isBuildingBlock": True,
                "catalogEntries": [
                    {
                        "catalogName": "emolecules",
                        "catalogId": "713319",
                        "smiles": "NCCC(N)C(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 4.0,
                            "bbPriceRange": "$100-500 / g",
                        },
                    },
                    {
                        "catalogName": "generic",
                        "catalogId": "",
                        "smiles": "NCCC(N)C(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "enamine_bb",
                        "catalogId": "EN300-296225",
                        "smiles": "NCCC(N)C(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 12.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "mcule",
                        "catalogId": "MCULE-1239764966",
                        "smiles": "NCCC(N)C(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "enamine_made",
                        "catalogId": "BBV-38371885",
                        "smiles": "NCCC(N)C(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": False,
                            "isScreening": True,
                            "scrLeadTimeWeeks": 6.0,
                            "scrPriceRange": "unknown",
                        },
                    },
                    {
                        "catalogName": "molport",
                        "catalogId": "MolPort-006-069-148",
                        "smiles": "NCCC(N)C(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                ],
            },
            {
                "smiles": "CC(C)(C)OC(=O)C1CCCC1",
                "isBuildingBlock": True,
                "catalogEntries": [
                    {
                        "catalogName": "emolecules",
                        "catalogId": "40102374",
                        "smiles": "CC(C)(C)OC(=O)C1CCCC1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.4,
                            "bbPriceRange": "$100-500 / g",
                        },
                    },
                    {
                        "catalogName": "mcule",
                        "catalogId": "MCULE-3526043917",
                        "smiles": "CC(C)(C)OC(=O)C1CCCC1",
                        "purchaseInfo": {
                            "isBuildingBlock": False,
                            "isScreening": True,
                            "scrLeadTimeWeeks": 4.0,
                            "scrPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "molport",
                        "catalogId": "MolPort-020-664-169",
                        "smiles": "CC(C)(C)OC(=O)C1CCCC1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.0,
                            "bbPriceRange": "$100-500 / g",
                        },
                    },
                ],
            },
            {
                "smiles": "CNCC(=O)Nc1cc(C=O)ccn1",
                "isBuildingBlock": False,
                "catalogEntries": [],
            },
            {
                "smiles": "Nc1cc(C=O)ccn1",
                "isBuildingBlock": True,
                "catalogEntries": [
                    {
                        "catalogName": "enamine_made",
                        "catalogId": "BBV-38325159",
                        "smiles": "Nc1cc(C=O)ccn1",
                        "purchaseInfo": {
                            "isBuildingBlock": False,
                            "isScreening": True,
                            "scrLeadTimeWeeks": 6.0,
                            "scrPriceRange": "unknown",
                        },
                    },
                    {
                        "catalogName": "emolecules",
                        "catalogId": "70144010",
                        "smiles": "Nc1cc(C=O)ccn1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 12.0,
                            "bbPriceRange": "$1k-10k / g",
                        },
                    },
                    {
                        "catalogName": "mcule",
                        "catalogId": "MCULE-8065423265",
                        "smiles": "Nc1cc(C=O)ccn1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 3.0,
                            "bbPriceRange": "$100-500 / g",
                        },
                    },
                    {
                        "catalogName": "enamine_bb",
                        "catalogId": "EN300-2987077",
                        "smiles": "Nc1cc(C=O)ccn1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 12.0,
                            "bbPriceRange": "$500-1k / g",
                        },
                    },
                ],
            },
            {
                "smiles": "CNCC(=O)O",
                "isBuildingBlock": True,
                "catalogEntries": [
                    {
                        "catalogName": "enamine_bb_EU-US",
                        "catalogId": "EN300-20732",
                        "smiles": "CNCC(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.4,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "emolecules",
                        "catalogId": "480611",
                        "smiles": "CNCC(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.4,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "mcule",
                        "catalogId": "MCULE-4484444573",
                        "smiles": "CNCC(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "generic",
                        "catalogId": "",
                        "smiles": "CNCC(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "molport",
                        "catalogId": "MolPort-001-785-681",
                        "smiles": "CNCC(=O)O",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                ],
            },
        ],
        "reactions": [
            {
                "name": "Reductive amination",
                "reactantSmiles": [
                    "NCCC(NC(=O)C1CCCC1)C(=O)O",
                    "CNCC(=O)Nc1cc(C=O)ccn1",
                ],
                "productSmiles": "CNCC(=O)Nc1cc(CNCCC(NC(=O)C2CCCC2)C(=O)O)ccn1",
            },
            {
                "name": "Ester amidation",
                "reactantSmiles": ["NCCC(N)C(=O)O", "CC(C)(C)OC(=O)C1CCCC1"],
                "productSmiles": "NCCC(NC(=O)C1CCCC1)C(=O)O",
            },
            {
                "name": "Amidation",
                "reactantSmiles": ["Nc1cc(C=O)ccn1", "CNCC(=O)O"],
                "productSmiles": "CNCC(=O)Nc1cc(C=O)ccn1",
            },
        ],
    },
    {
        "molecules": [
            {
                "smiles": "Cn1nccc1C(=O)N1CCCCC1",
                "isBuildingBlock": False,
                "catalogEntries": [],
            },
            {
                "smiles": "COC(=O)c1ccnn1C",
                "isBuildingBlock": True,
                "catalogEntries": [
                    {
                        "catalogName": "mcule",
                        "catalogId": "MCULE-1338013575",
                        "smiles": "COC(=O)c1ccnn1C",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "molport",
                        "catalogId": "MolPort-000-162-104",
                        "smiles": "COC(=O)c1ccnn1C",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "enamine_real",
                        "catalogId": "s_1458____2164982____1327596",
                        "smiles": "COC(=O)c1ccnn1C",
                        "purchaseInfo": {
                            "isBuildingBlock": False,
                            "isScreening": True,
                            "scrLeadTimeWeeks": 4.0,
                            "scrPriceRange": "unknown",
                        },
                    },
                    {
                        "catalogName": "enamine_made",
                        "catalogId": "BBV-36627156",
                        "smiles": "COC(=O)c1ccnn1C",
                        "purchaseInfo": {
                            "isBuildingBlock": False,
                            "isScreening": True,
                            "scrLeadTimeWeeks": 6.0,
                            "scrPriceRange": "unknown",
                        },
                    },
                    {
                        "catalogName": "emolecules",
                        "catalogId": "842571",
                        "smiles": "COC(=O)c1ccnn1C",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 4.0,
                            "bbPriceRange": "$500-1k / g",
                        },
                    },
                    {
                        "catalogName": "wuxi_bb_screening",
                        "catalogId": "LN00002790",
                        "smiles": "COC(=O)c1ccnn1C",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "generic",
                        "catalogId": "",
                        "smiles": "COC(=O)c1ccnn1C",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                ],
            },
            {
                "smiles": "C1CCNCC1",
                "isBuildingBlock": True,
                "catalogEntries": [
                    {
                        "catalogName": "generic",
                        "catalogId": "",
                        "smiles": "C1CCNCC1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "< $100 / g",
                        },
                    },
                    {
                        "catalogName": "emolecules",
                        "catalogId": "477538",
                        "smiles": "C1CCNCC1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 1.4,
                            "bbPriceRange": "$1k-10k / g",
                        },
                    },
                    {
                        "catalogName": "mcule",
                        "catalogId": "MCULE-4094049864",
                        "smiles": "C1CCNCC1",
                        "purchaseInfo": {
                            "isBuildingBlock": True,
                            "isScreening": False,
                            "bbLeadTimeWeeks": 2.0,
                            "bbPriceRange": "$100-500 / g",
                        },
                    },
                ],
            },
        ],
        "reactions": [
            {
                "name": "Ester amidation",
                "reactantSmiles": ["COC(=O)c1ccnn1C", "C1CCNCC1"],
                "productSmiles": "Cn1nccc1C(=O)N1CCCCC1",
            }
        ],
    },
]
