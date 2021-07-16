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
                    "temperature": {"value": 25, "unit": "degC"},  # degrees celcius
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
                        "density": 0.74,
                        "concentration": None,
                    },
                },
            },
            {
                "name": "stir",
                "content": {
                    "action_no": 4,
                    "temperature": {"value": 25, "unit": "degC"},  # degrees celcius
                    "duration": {"value": 12, "unit": "hours"},  # in hours
                },
            },
        ],
    },
    "Reductive amination": {
        "reactionSMARTS": ["[#7;H2:1].[#6:2](=[#8:3])>>[#7:1]-[#6:2]"],
        "recipe": [
            {
                "name": "add",
                "content": {
                    "action_no": 1,
                    "material": {
                        "SMARTS": ["[#6](=[#8])"], #This now searches for any carbonyl group (low specificty)
                         #[#6][CX3](=O)[#6]", "[CX3H1](=O)[#6] OLD code for reference
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
                        "SMARTS": ["[#6]-[NX3]"],
                        #[#6]-[#7;H2]
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
                    "temperature": {"value": 25, "unit": "degC"},  # degrees celcius
                    "duration": {"value": 12, "unit": "hours"},  # in hours
                },
            },
        ],
    },
    "N-nucleophilic aromatic substitution": {
        "reactionSMARTS": ["[c:1]-[F,Cl,Br,I].[c:3]-[N:2]>>[c:1]-[N:2]-[c:3]"], #HR 15/07/21: Changed Reaction SMARTS
        "recipe": [
            {
                "name": "add",
                "content": {
                    "action_no": 1,
                    "material": {
                        # What about looking for other nucelophiles? Check with Harry -> SMARTS above also not very general -> only
                        # looking for amines
                        "SMARTS": ["[c:1]-[F,Cl,Br,I]"], #HR 15/07/21 removing nitro group as identifier
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
                        "SMARTS": ["[c:3]-[N:2]"],  # allowing amines 
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
                    "temperature": {"value": 50, "unit": "degC"},  # generic elevated temperature
                    "duration": {"value": 12, "unit": "hours"},  # in hours
                },
            },
        ],
    },
    "sp2-sp2 suzuki coupling": {
        "reactionSMARTS": ["[c:1]-[F,Cl,Br,I].[#6:2]-[B]>>[c:1]-[#6:2]"],
        "recipe": [
            {
                "name": "add",
                "content":{
                    "action_no": 1,
                    "material":{
                        "SMARTS": ["[c:1]-[F,Cl,Br,I]"],
                        "SMILES": None,
                        "quantity": {"value": "Unsure yet", "unit": "moleq"},
                        "solvent": "EtOH",
                        "concentration": "Whatever this is",
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 2,
                    "material":{
                        "SMARTS": ["[#6:2]-[B]"],
                        "SMILES": None,
                        "quantity": {"value": "Unsure yet", "unit": "moleq"},
                        "solvent": "EtOH",
                        "concentration": "Whatever this is",
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 3,
                    "material":{
                        "SMARTS": None,
                        "SMILES": "[Fe].Cl[Pd]Cl.[CH]1[CH][CH][C]([CH]1)P(c2ccccc2)c3ccccc3.[CH]4[CH][CH][C]([CH]4)P(c5ccccc5)c6ccccc6",
                        #Smiles for the Pd-Ferrocene catalyst we have
                        "quantity": {"value": "Unsure yet", "unit": "moleq"},
                        "solvent": "EtOH",
                        "concentration": "Whatever this is",
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 4,
                    "material":{
                        "SMARTS": None,
                        "SMILES": "C1CCN2CCCN=C2CC1",#DBU used in Thompson paper
                        "quantity": {"value": "Unsure yet", "unit": "moleq"},
                        "solvent": "EtOH",
                        "concentration": "Whatever this is",
                    },
                },
            },
            {
                "name": "stir",
                "content":{
                    "action_no": 5,
                    "temperature": {"value": 100, "unit": "degC"},  # high temp required
                    "duration": {"value": 12, "unit": "hours"}, # in hours
                },
            },
        ],
    },
}



                
            
           
                   

                    
          
  

