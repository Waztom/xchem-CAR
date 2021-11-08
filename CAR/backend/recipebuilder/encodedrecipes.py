encoded_recipes = {
    "Amidation": {
        "reactionSMARTS": "[#6:1](=[#8:2])-[#8].[#7;H3,H2,H1:3]>>[#6:1](=[#8:2])-[#7:3]",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#6](=[#8])-[#8]",
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
                                "density": 0.74,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": "[#7;H3,H2,H1]",
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
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
            "Intramolecular": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#6](=[#8])-[#8]",
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
                                "density": 0.74,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Amide schotten - baumann": {
        "reactionSMARTS": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7:3]",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#7;H2,H1:3]",
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
                                "SMARTS": "[#6:1](=[#8:2])-[#17]",
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
                                "density": 0.74,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 4,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
            "Intramolecular": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#7;H2,H1:3]",
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
                            "action_no": 3,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CCN(C(C)C)C(C)C",
                                "quantity": {"value": 3.5, "unit": "moleq"},
                                "solvent": None,
                                "density": 0.74,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 4,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Reductive amination": {
        "reactionSMARTS": "[#6:2](=[#8])(-[#6:1]).[#7;H3,H2,H1:3]>>[#6:2](-[#6:1])-[#7:3]",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#6:2](=[#8])(-[#6:1])",
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
                                "SMARTS": "[#7;H2,H1:3]",
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
                                "solvent": "DMA",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
            "Intramolecular": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#6:2](=[#8])(-[#6:1])",
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
                                "solvent": "DMA",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "N-nucleophilic aromatic substitution": {
        "reactionSMARTS": "[#6:3]-[#7;H3,H2,H1:2].[c:1]-[F,Cl,Br,I]>>[#6:3]-[#7:2]-[c:1]",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#6:3]-[#7;H2,H1:2]",
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": "[c:1]-[F,Cl,Br,I]",
                                "SMILES": None,
                                "quantity": {"value": 1.2, "unit": "moleq"},
                                "solvent": "NMP",
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
                                "quantity": {"value": 2.5, "unit": "moleq"},
                                "solvent": None,
                                "density": 0.74,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 4,
                            "temperature": {"value": 50, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Sp2-sp2 Suzuki coupling": {
        "reactionSMARTS": "[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[c:1]-[F,Cl,Br,I]",
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": "[#6:2]-[B]",
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "EtOH",
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
                                "SMILES": "CC(C)C1=CC(=C(C(=C1)C(C)C)C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4)C(C)C.CS(=O)(=O)[O-].C1=CC=C(C=C1)C2=CC=CC=C2N.[Pd]",
                                # Smiles for XPhosPdG3
                                "quantity": {"value": 0.1, "unit": "moleq"},  # 10mol% catalyst
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "C1CCN2CCCN=C2CC1",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 100, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    ##############################################################################################################
    ################ Reactions we will not use for, are not working on and still have to test ###########################################################
    "Sulfonamide schotten-baumann": {
        "reactionSMARTS": "[#16:5](=[#8])(=[#8:7])-[#17].[#6]-[#7;H2,H1:2]>>[#16:5](=[#8])(=[#8:7])-[#7:2]",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "MeOH",  # solvent for all can be IPA, EtOAc or MeOH will need to test
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": "[#6]-[#7;H2,H1:2]",
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
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
            "Intramolecular": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#16:5](=[#8])(=[#8:7])-[#17]",
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "MeOH",  # solvent for all can be IPA, EtOAc or MeOH will need to test
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 3,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Boc protection": {
        "reactionSMARTS": [
            "[#7:2].[#8:3]-[#6:1](=[#8:4])-[#8:5]>>[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])"
        ],
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": [
                                    "[#8:3]-[#6:1](=[#8:4])-[#8:5]-[#6](-[#6])(-[#6])(-[#6])"
                                ],  # BOC
                                "SMILES": None,
                                "quantity": {"value": 3, "unit": "moleq"},
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
                                "SMARTS": None,
                                "SMILES": "n1ccc(N(C)C)cc1",
                                "quantity": {"value": 1, "unit": "moleq"},  # DMAP
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
                                "SMARTS": ["[#7;H2:2]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},  # primary amine
                                "solvent": "ACN",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 4,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 0.5, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Boc deprotection": {  # this does not work yet, SMARTS needs rethinking
        "reactionSMARTS": [
            "[#6:9]-[#6:8]-[#7:2]-[#6](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6]).[#1]-[#17]>>[#6:9]-[#6:8]-[#7:2]"
        ],  # [#7:2].[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])]
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": [
                                    "[#6:9]-[#6:8]-[#7:2]-[#6](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])"
                                ],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "EtOAc",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[#1]-[#17]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "EtOAc",
                                "concentration": 6,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 3,
                            "temperature": {"value": 25, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Buchwald hartwig amination - NMP": {
        "reactionSMARTS": ["[c:2]-[F,Cl,Br,I].[#6:1]-[#7;H2,H1]>>[c:2]-[#7:3]-[#6:1]"],
        "recipes": {
            "Standard": {
                "reference": "https://pubs.acs.org/doi/pdf/10.1021/acscatal.9b00981",
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 3,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.1,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[#6:1]-[#7;H2,H1:3]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "C1CCN2CCCN=C2CC1",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 1.02,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Buchwald hartwig amination - EtOH": {
        "reactionSMARTS": ["[c:2]-[F,Cl,Br,I].[#6:1]-[#7;H2,H1]>>[c:2]-[#7:3]-[#6:1]"],
        "recipes": {
            "Standard": {
                "reference": "https://pubs.acs.org/doi/pdf/10.1021/acscatal.9b00981",
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 3,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.1,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[#6:1]-[#7;H2,H1:3]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "C1CCN2CCCN=C2CC1",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 1.02,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Buchwald hartwig thiolation": {
        "reactionSMARTS": ["[c:2]-[F,Cl,Br,I].[#6:1]-[#16;H1]>>[c:2]-[#16:3]-[#6:1]"],
        "recipes": {
            "Standard": {
                "reference": "https://doi.org/10.1002/ejoc.201001393",
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "NMP",  # try NMP but reference says to use 2-MeTHF
                                "concentration": 0.1,
                                # try same as Suzuki - similar mechanism
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 3,
                            "material": {
                                "SMARTS": ["[#6:1]-[#16;H1]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "C1CCN2CCCN=C2CC1",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 1.02,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Buchwald hartwig etherification - EtOH": {
        "reactionSMARTS": ["[c:2]-[F,Cl,Br,I].[#6:1]-[#8;H1:3]>>[c:2]-[#8:3]-[#6:1]"],
        "recipes": {
            "Standard": {
                "reference": "https://pubs.acs.org/doi/pdf/10.1021/ja984321a",
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 3,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "EtOH",  # try NMP but reference says to use 2-MeTHF
                                "concentration": 0.1,
                                # try same as Suzuki - similar mechanism
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[#6:1]-[#8;H1:3]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "C1CCN2CCCN=C2CC1",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 1.02,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Buchwald hartwig etherification - NMP": {
        "reactionSMARTS": ["[c:2]-[F,Cl,Br,I].[#6:1]-[#8;H1:3]>>[c:2]-[#8:3]-[#6:1]"],
        "recipes": {
            "Standard": {
                "reference": "https://pubs.acs.org/doi/pdf/10.1021/ja984321a",
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 3,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "NMP",  # try NMP but reference says to use 2-MeTHF
                                "concentration": 0.1,
                                # try same as Suzuki - similar mechanism
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[#6:1]-[#8;H1:3]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "C1CCN2CCCN=C2CC1",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 1.02,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Sonogashira coupling - EtOH": {
        "reactionSMARTS": ["[CH:1].[c:2]-[Cl,Br,I]>>[C:1]-[c:2]"],
        "recipes": {
            "Standard": {
                "reference": [
                    "https://pubs.rsc.org/en/content/articlelanding/2008/cc/b810928a#!",
                    "https://pubs.acs.org/doi/pdf/10.1021/ol035632f"
                    "https://www.sciencedirect.com/science/article/abs/pii/S1381116908003257",
                ],
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 3,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.1,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": ["[c:2]-[Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[CH:1]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "EtOH",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "[Cu]I",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "ACN",
                                "concentration": 0.1,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 5,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CCN(C(C)C)C(C)C",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 0.74,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Sonogashira coupling - NMP": {
        "reactionSMARTS": ["[CH:1].[c:2]-[Cl,Br,I]>>[C:1]-[c:2]"],
        "recipes": {
            "Standard": {
                "reference": [
                    "https://pubs.rsc.org/en/content/articlelanding/2008/cc/b810928a#!",
                    "https://pubs.acs.org/doi/pdf/10.1021/ol035632f"
                    "https://www.sciencedirect.com/science/article/abs/pii/S1381116908003257",
                ],
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 3,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.1,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": ["[c:2]-[Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[CH:1]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CCN(C(C)C)C(C)C",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 0.74,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Heck coupling": {
        "reactionSMARTS": [
            "[c:2]-[F,Cl,Br,I].[CX3;H2:1]>>[c:2]-[CX3;H1:1]"
        ],  # SMARTS for terminal alkenes only
        "recipes": {
            "Standard": {
                "reference": ["Platinum Metals Rev., 1999, 43, (4), 138"],
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": ["[c:2]-[F,Cl,Br,I]"],
                                "SMILES": None,
                                "quantity": {"value": 1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.5,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 2,
                            "material": {
                                "SMARTS": ["[CX3;H2:1]"],
                                "SMILES": None,
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": "NMP",
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
                                "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                                "quantity": {"value": 0.1, "unit": "moleq"},
                                "solvent": "NMP",
                                "concentration": 0.1,
                            },
                        },
                    },
                    {
                        "name": "add",
                        "content": {
                            "action_no": 4,
                            "material": {
                                "SMARTS": None,
                                "SMILES": "CCN(CC)CC",
                                "quantity": {"value": 2, "unit": "moleq"},
                                "solvent": None,
                                "density": 0.73,
                                "concentration": None,
                            },
                        },
                    },
                    {
                        "name": "stir",
                        "content": {
                            "action_no": 5,
                            "temperature": {"value": 80, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
}