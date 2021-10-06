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
        "reactionSMARTS": "[#7;H2,H1:3].[#6:1](=[#8:2])-[#17]>>[#6:1](=[#8:2])-[#7;H2,H1:3]",
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
    "Boc protection": {
        "reactionSMARTS": "[#8:3]-[#6:1](=[#8:4])-[#8:5]-[#6](-[#6])(-[#6])(-[#6]).[#7;H2:2]>>[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#8:3]-[#6:1](=[#8:4])-[#8:5]-[#6](-[#6])(-[#6])(-[#6])",  # BOC
                                "SMILES": None,
                                "quantity": {"value": 1.5, "unit": "moleq"},
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
                                "SMARTS": "[#7;H2]",
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
                            "temperature": {"value": 100, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
    "Boc deprotection": {  # this does not work yet, SMARTS needs rethinking
        "reactionSMARTS": "[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])]>>[#7:2].[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])]",
        "recipes": {
            "Standard": {
                "reference": None,
                "actions": [
                    {
                        "name": "add",
                        "content": {
                            "action_no": 1,
                            "material": {
                                "SMARTS": "[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])]",
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
                                "SMARTS": None,
                                "SMILES": "HCl",
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
                            "temperature": {"value": 100, "unit": "degC"},
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
                            "temperature": {"value": 100, "unit": "degC"},
                            "duration": {"value": 12, "unit": "hours"},
                        },
                    },
                ],
            },
        },
    },
}
