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
                    "temperature": {"value": 25, "unit": "degC"},  
                    "duration": {"value": 12, "unit": "hours"},  
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
                    "temperature": {"value": 25, "unit": "degC"},  
                    "duration": {"value": 12, "unit": "hours"},  
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
                        "SMARTS": ["[#6](=[#8])"], 
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
                        "SMARTS": ["[#6]-[NX3]"],
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
    "N-nucleophilic aromatic substitution": {
        "reactionSMARTS": ["[c:1]-[F,Cl,Br,I].[c:3]-[N:2]>>[c:1]-[N:2]-[c:3]"], 
        "recipe": [
            {
                "name": "add",
                "content": {
                    "action_no": 1,
                    "material": {
                        # What about looking for other nucelophiles? Check with Harry -> SMARTS above also not very general -> only
                        # looking for amines HR: for the moment keeping with amines as nucleophile, alcohols can be used in a further
                        #iteration as base will probably be required
                        "SMARTS": ["[c:1]-[F,Cl,Br,I]"], 
                        "SMILES": None,  
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
                        "SMARTS": ["[c:3]-[N:2]"],  
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
                    "temperature": {"value": 50, "unit": "degC"},  
                    "duration": {"value": 12, "unit": "hours"}, 
                },
            },
        ],
    },
    "sulfonamide schotten-baumann": {
        "reactionSMARTS": ["[#6:1]-[#17:3].[#7:8]>>[#6:1]-[#7:8]"],
        "recipe": [
            {
                "name": "add",
                "content":{
                    "action_no": 1,
                    "material":{
                        "SMARTS":["[#6:1](=[#8])-[#17:3]"],
                        "SMILES": None,
                        "quantity": {"value": 1, "unit": "moleq"},
                        "solvent": "MeOH", #solvent for all can be IPA, EtOAc or MeOH will need to test
                        "concentration": 0.5,
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 2,
                    "material":{
                        "SMARTS":["[#7:8]-[#16:5](=[#8:6])(=[#8:7])-[#6:4]"],
                        "SMILES": None,
                        "quanity": {"value": 1.2, "unit": "moleq"},
                        "solvent": "MeOH",
                        "concentration": 0.5,
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 3, # using DBU as pKa of sulfonamide ~10, DBU pKa ~13
                    "material":{
                        "SMARTS": None,
                        "SMILES": "C1CCC2=NCCCN2CC1",
                        "quantity": {"value": 3.5, "unit": "moleq"},
                        "solvent": "MeOH",
                        "concentration": 0.5,
                    },
                },
            },
            {
                "name": "stir",
                "content":{
                    "action_no": 4,
                    "temperature": {"value": 100, "unit": "degC"},  
                    "duration": {"value": 12, "unit": "hours"}, 
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
                        "SMARTS": ["[c:1]-[F,Cl,Br,I].[c]-[!B]"], 
                        "SMILES": None,
                        "quantity": {"value": 1, "unit": "moleq"},
                        "solvent": "EtOH",
                        "concentration": 0.5,
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
                        "quantity": {"value": 2, "unit": "moleq"}, 
                        "solvent": "EtOH",
                        "concentration": 0.5,
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 3,
                    "material":{
                        "SMARTS": None,
                        "SMILES": "CC(C)C(C=C1C(C)C)=CC(C(C)C)=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4.NC5=C(C6=C([Pd]OS(C)(=O)=O)C=CC=C6)C=CC=C5",
                        #Smiles for XPhosPdG3
                        "quantity": {"value": 0.1, "unit": "moleq"}, #10mol% catalyst
                        "solvent": "EtOH",
                        "concentration": 0.5,
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 4,
                    "material":{
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
                "content":{
                    "action_no": 5,
                    "temperature": {"value": 100, "unit": "degC"},  
                    "duration": {"value": 12, "unit": "hours"}, 
                },
            },
        ],
    },
    "Boc protection": {
        "reactionSMARTS": ["[#7:2].[#8:3]-[#6:1](=[#8:4])-[#8:5]>>[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])"],
        "recipe": [
            {
                "name": "add",
                "content":{
                    "action_no": 1,
                    "material":{
                        "SMARTS": ["[#8:3]-[#6:1](=[#8:4])-[#8:5]-[#6](-[#6])(-[#6])(-[#6])"],#BOC
                        "SMILES": None,
                        "quantity": {"value": 1.5, "unit": "moleq"},
                        "solvent": "ACN",
                        "concentration": 0.5,
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 2,
                    "material":{
                        "SMARTS": None,
                        "SMILES": "n1ccc(N(C)C)cc1",
                        "quantity": {"value": 1, "unit": "moleq"},#DMAP
                        "solvent": "ACN",
                        "concentration": 0.5,
                    },
                },
            },
            {   
                "name": "add",
                "content":{
                    "action_no": 3,
                    "material":{
                        "SMARTS": ["[#7;H2]"],
                        "SMILES": None,
                        "quantity": {"value": 1, "unit": "moleq"},#primary amine
                        "solvent": "ACN",
                        "concentration": 0.5,
                    },
                },
            },
            {
                 "name": "stir",
                 "content":{
                    "action_no": 4,
                    "temperature": {"value": 100, "unit": "degC"},  
                    "duration": {"value": 12, "unit": "hours"},
                },
            },
        ],
    },
    "Boc deprotection": {#this does not work yet, SMARTS needs rethinking
        "reactionSMARTS":["[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])]>>[#7:2].[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])]"],
        "recipe":[
            {
                "name": "add",
                "content":{
                    "action_no": 1,
                    "material":{
                        "SMARTS": ["[#7:2]-[#6:1](=[#8])-[#8]-[#6](-[#6])(-[#6])(-[#6])]"],
                        "SMILES": None,
                        "quantity": {"value": 1, "unit": "moleq"},
                        "solvent": "EtOAc",
                        "concentration": 0.5,
                    },
                },
            },
            {
                "name": "add",
                "content":{
                    "action_no": 2,
                    "material":{
                        "SMARTS": None,
                        "SMILES": "HCl",
                        "quantity": {"value": 1, "unit": "moleq"},
                        "solvent": "EtOAc",
                        "concentration": 6
                    },
                },
            },
            {
                "name": "stir",
                 "content":{
                    "action_no": 3,
                    "temperature": {"value": 25, "unit": "degC"},  
                    "duration": {"value": 12, "unit": "hours"},
                },
            },
        ],
    },
}


                
            
           
                   

                    
          
  

