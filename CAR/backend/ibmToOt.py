from pathlib import Path
import pandas as pd
import math
import re
#import os,sys
#sys.path.insert(0,'..')
import ibmRead

import opentrons.otWrite as otWrite
import opentrons.otDeck as otDeck




class otSession():
    def __init__(self, name, actions, author=None, description=None, split=None):

        self.name = name
        self.namecheck()

        self.author=author
        self.description=description,
        

        self.actions = []

        self.deck = otDeck.Deck()
        
        self.outputpath = "None"
        self.pathnamecheck()
        
        # self.getactions(actions)
        self.actions = actions
        self.tipoptions = []
        self.tipsneeded = {}
        self.tipRackList = []

        self.pipettesneeded = []

        self.output = otWrite.otScript(filepath=self.outputpath, protocolName=self.name, author=self.author, description=self.description)
        self.output.setupScript()
        self.setupPlate()
        
        
        self.preselecttips()
        
        
        
        self.tipOutput()
        self.output.setupLabwear(self.deck.PlateList, self.tipRackList)
        self.setupPipettes()
        self.output.setupPipettes(self.deck.PipetteList)

        #if split == None:
        #    self.ittrActions()
        #else:
        #    self.splitittract(split)
        self.ittrActions()

    
    def namecheck(self):
        trialname = self.name
        if re.match("^[\w-]+$", trialname):
            pass
        else:
            trialname = re.sub(r"^[\w-]+$", '', trialname)
            if trialname == '':
                trialname = "unamaed"

        self.name = trialname
        return trialname


    def pathnamecheck(self):
        self.namecheck()
        trialfile = Path("../output/" + self.name + ".py")
        if trialfile.exists():
            suffixnumber = 1
            while trialfile.exists():
                trialfile = Path("../output/" + self.name + "(" + str(suffixnumber) + ")" + ".py")
                #print("Debug, trying: "+str(trialfile))
                suffixnumber += 1
            
        self.outputpath = trialfile
        return self.outputpath

    def getactions(self, reactionno, split=None):
        if type(reactionno) == 'list':
            if len(reactionno) == 0:
                pass
            elif len(reacitonno) == 1:
                pass
            else:
                pass

        allactions = ibmRead.getactions()
        reactionsactions = ibmRead.getReactionActions(allactions, reactionno)
        self.actions = reactionsactions
        return reactionsactions

    def actionfilter(self, split):
        pass
    def setupPlate(self, incAddMaterials=True, incWashMaterials=True,incExtract=False):
        #print(self.actions)
        allmaterials = pd.DataFrame(columns=['material', 'materialquantity'])

        if incAddMaterials == True:
            addingSteps = self.actions.loc[(self.actions['actiontype']) == "add"]
            materials = addingSteps[['materialsmiles', 'materialquantity']]
            materials = materials.rename(columns={'materialsmiles':'material'})
            allmaterials = pd.concat([allmaterials, materials], ignore_index=True)

        if incWashMaterials==True:
            solventSteps = self.actions.loc[(self.actions['actiontype']) == "wash"]
            washmaterials = solventSteps[['material', 'materialquantity']]
            allmaterials = pd.concat([allmaterials, washmaterials], ignore_index=True)

        if incExtract == True:
            extractSteps = self.actions.loc[(self.actions['actiontype']) == "extract"]
            extractsolvents = solventSteps[['solvent', 'solventquantity']]
            extractsolvents = extractsolvents.rename(columns={'solvent':'material'})
            extractsolvents = extractsolvents.rename(columns={'solventquantity':'materialquantity'})
            allmaterials = pd.concat([allmaterials, extractsolvents], ignore_index=True)

        #print(allmaterials)

        self.maxtransfer = allmaterials['materialquantity'].max()
        self.totalvolume = allmaterials['materialquantity'].sum()

        self.orderPlate = self.deck.add("Plate", 1,12,500, "OrderPlate")

        plate = ['','','','','','','','','',''] #note: need to fix so not fixed posibles
        for i in allmaterials.index.values:
            wellnumber = self.orderPlate.nextfreewell()
            outcome = self.orderPlate.WellList[wellnumber].add(allmaterials.loc[i,'materialquantity'], smiles = allmaterials.loc[i,'material'])
            if outcome != False:
                #add otWrite setup plate here
                plate[wellnumber] = [allmaterials.loc[i,'materialquantity'], allmaterials.loc[i,'material']] 

        self.reactionPlate = self.deck.add("Plate", 2, 12, 500, "ReactionPlate")
        # for well in range(len(self.orderPlate)):
        #     print(self.orderPlate[well].StartSmiles)
        # print(self.orderPlate.smilesearch('C(Cl)Cl'))
        # print(self.orderPlate.smilesearch('brine'))

    def preselecttips(self):
        for actionindex in range(len(self.actions)):
            currentactionmask = self.actions['actionno'] == actionindex+1
            currentaction = self.actions[currentactionmask]
            currentactiontype = (currentaction['actiontype'].to_string(index=False)).strip()
            if currentactiontype == "add":
                self.choosetip(currentaction['materialquantity'].values[0], tipoptions = [])

    def choosetip(self, volume, tipoptions= []):
        if tipoptions == []:
            if self.tipoptions == []:
                #self.tipoptions = [10, 20, 200, 300, 1000]
                self.tipoptions = [300] # debuging
        else:
            pass
            # review wether to use tip options or self.tipoptions for rest of function
        for tip in self.tipoptions:
            if tip >= volume:
                if tip in self.tipsneeded:
                    self.tipsneeded[tip]+=1
                else:
                    self.tipsneeded[tip] = 1
                    #self.tipsneeded[1].append(1)
                return tip 
        return False
    
    #def setuptips

    def tipOutput(self):
        self.tipRackList = []
        for tip in self.tipsneeded :

            numplates = math.ceil(self.tipsneeded[tip]/96)

            while numplates >= 1:

                location = len(self.deck.PlateList)

                if tip == 10:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=10, platename = "opentrons_96_tiprack_10ul")
                    self.tipRackList.append(["opentrons_96_tiprack_10ul", self.deck.nextfreeplate(), tip])
                elif tip == 20:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=20, platename = "opentrons_96_tiprack_20ul")
                    self.tipRackList.append(["opentrons_96_tiprack_20ul", self.deck.nextfreeplate(), tip])
                elif tip == 200:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=200, platename = "opentrons_96_rtiprack_200ul")
                    self.tipRackList.append(["opentrons_96_tiprack_200ul", self.deck.nextfreeplate(), tip])
                elif tip == 300:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=300, platename = "opentrons_96_tiprack_300ul")
                    self.tipRackList.append(["opentrons_96_tiprack_300ul", self.deck.nextfreeplate(), tip])
                elif tip == 1000:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=1000, platename = "opentrons_96_tiprack_1000ul")
                    self.tipRackList.append(["opentrons_96_tiprack_1000ul", self.deck.nextfreeplate(), tip])

                numplates -= 1
        return self.tipRackList

    def setupPipettes (self):
        if len(self.tipsneeded) <= 2:
            if len(self.tipsneeded) > 0 :
                for pipette in self.tipsneeded:
                    self.pipettesneeded.append(pipette)
        mountnumber = 0
        for pipette in self.pipettesneeded:
            if mountnumber == 0:
                mount = "left"
                mountnumber += 1
            elif mountnumber == 1:
                mount = "right"
            else:
                break
            self.deck.addPipette(str(f"{mount}_{pipette}_pipette"), str("p"+str(pipette)+"_single"), mount, pipette)
            

    def ittrActions(self):
        # note: add checks to prevent errors if df of diffrent format or blank is presented
        # actionsCoppy = self.actions
        # actionsCoppy = actionsCoppy.sort_values(['id', 'actionno'], ascending=(True, False))
        print(self.actions)
        print(type(self.actions))
        # for actionitterator in range(len(self.actions)):
        #     print(actionitterator)
        #     currentaction = self.actions.loc[actionitterator]
        for currentaction in self.actions.items():
            # currentactionmask = self.actions['actionno'] == actionindex+1
            # currentaction = self.actions[currentactionmask]
            print(">>>>>>>>>>Debug>>>>>>>>>>>>>")
            print(currentaction)
            self.processAction(currentaction)

    def splitittract(self, split):
        print(split)
        print(type(split))
        splittype = False
        if type(split) == 'int':
            subset = self.actions[split]
            for actionindex in range(len(subset)+1):
                currentactionmask = subset['actionno'] == actionindex+1
                currentaction = subset[currentactionmask]                    
                #print("current: "+str(currentaction))
                self.processAction(currentaction)
            
        elif type(split) == "<class 'list'>":
            print(type(split[0]))
            if type(split[0]) == 'int':
                subset = self.actions[split[0]:split[1]]
                print("intlistsplit")
                for actionindex in range(len(subset)+1):
                    currentactionmask = subset['actionno'] == actionindex+1
                    currentaction = subset[currentactionmask]
                    #print("current: "+str(currentaction))
                    self.processAction(currentaction)

            elif type(split[0])=='str':
                print("stirlistsplit")
                subset = self.actions['actiontype' in split]
                for actionindex in range(len(subset)+1):
                    currentactionmask = subset['actionno'] == actionindex+1
                    currentaction = subset[currentactionmask]
                    #print("current: "+str(currentaction))
                    self.processAction(currentaction)
            else:
                print(type(split[0]))

        elif type(split) == 'str':
                subset = self.actions['actiontype' == split]
                for actionindex in range(len(subset)+1):
                    currentactionmask = subset['actionno'] == actionindex+1
                    currentaction = subset[currentactionmask]
                    #print("current: "+str(currentaction))
                    self.processAction(currentaction)
        else:
            print("debug1")


    def processAction(self, currentaction):
        print(f"currentaction {currentaction}")
        currentactiontype = currentaction['actiontype'].to_string(index=False)
        currentactiontype = currentactiontype.strip()
        if currentactiontype != 'Series([], )':
            print(f"#### {currentactiontype}")
            #print("current type: "+str(currentactiontype))
            numreps = currentaction['numberofrepetitions'].values
            if str(numreps) == "[nan]":
                numreps = 0
            elif str(numreps) == "[None]":
                numreps = 0
            elif str(numreps) == "[]":
                numreps = 0

            print("numrepsa: "+str(numreps))
            try:
                numreps = int(numreps)
            except:
                numreps = 0

            #print(currentaction)
            #print("numreps: "+str(numreps))


            repetitions = 0
            while repetitions  <= numreps:
                #print(currentaction)
                if currentactiontype == "add":
                    self.actionAdd(currentaction)
                elif currentactiontype == "collect-layer":
                    self.actionCollectLayer(currentaction)
                elif currentactiontype == "wash":
                    self.actionWash(currentaction)
                elif currentactiontype == "stir":
                    self.output.unsuportedAction("stir at a "+str(currentaction['stirringspeed'].values[0])+" speed at "+str(currentaction['temperature'].values[0])+" Celsius for "+str(currentaction['duration'].values[0])+" "+str(currentaction['durationunit'].values[0]))
                elif currentactiontype == "set-temperature":
                    if currentaction['duration'].values[0] in ['nan','NaN',0,'0']:
                        self.output.unsuportedAction(f"set temperature to {currentaction['temperature'].values[0]} Celsius for {currentaction['duration'].values[0]} {currentaction['durationunit'].values[0]}")
                    else:
                        self.output.unsuportedAction(f"set temperature to {currentaction['temperature'].values[0]} Celsius")

                elif currentactiontype == "store":
                    self.output.unsuportedAction("store product ("+str(currentaction['material'].values[0])+")")
                elif currentactiontype == "extract":
                    self.actionExtract(currentaction)
                elif currentactiontype == "concentrate":
                    self.actionConcentrate(currentaction)
                else:
                    print(f"unsupported{currentaction['actiontype'].values}")
                    self.output.unsuportedAction(f"{currentaction['actiontype'].values[0]} is not currently supported ")
                repetitions +=1
        

    def actionAdd(self, currentaction):
        print("add")
        tipvolume = self.choosetip(currentaction['materialquantity'].values[0])
        #print(tipvolume)
        pipetteName = (self.deck.findPippets(tipvolume)).name

        self.output.transferfluids( pipetteName, 
            (f"OrderPlate.wells(){self.deck['OrderPlate'].smilesearch(currentaction['materialsmiles'].values[0], start_smiles = True)[0]}"),
            (f"ReactionPlate.wells()[{self.deck['ReactionPlate'].activeWell}]"),
            currentaction['materialquantity'].values[0])

    def actionWash(self, currentaction):
        print("wash")
        tipvolume = self.choosetip(currentaction['materialquantity'].values[0])
        #print(tipvolume)
        pipetteName = self.deck.findTipRacks(tipvolume)
        pipetteName = pipetteName[0]

        self.output.transferfluids(pipetteName,
            (f"OrderPlate.wells()[{self.deck['OrderPlate'].smilesearch(currentaction['material'].values[0], start_smiles = True)[0]}]"),
            (f"ReactionPlate.wells()[{self.deck['ReactionPlate'].activeWell}]"),
            currentaction['materialquantity'].values[0])

    def actionCollectLayer(self, currentaction):
        print("collect")
        self.output.unsuportedAction("Collect-Layer not yet supported ")


    def actionExtract(self, currentaction):
        print("Extract")
        self.output.unsuportedAction(f"Extract using {currentaction['solventquantity'].values[0]} {currentaction['solventquantityunit'].values[0]} of {currentaction['solvent'].values[0]}")
        # for col in currentaction:
        #     print(currentaction[col].values[0])

    def actionConcentrate(self, currentaction):
        print("Concentrate")
        self.output.unsuportedAction("Concentrate not yet supported ")
    
allactions = ibmRead.getactions()
print(allactions)
print(allactions.columns)
#reactionsactions = ibmRead.getReactionActions(allactions, reactionno)
print(allactions['actiontype'].unique())

def docheck(row):
    if row['actiontype'] in ['add', 'wash', "extract"]:
        return True
    else:
        return False

def actionfilter(allactions,  actions=None, reactionset=None):
    if reactionset != None:
        subSetReactAct = allactions.loc[(allactions['reaction_id_id']).isin(reactionset)]
    else:
        subSetReactAct = allactions
    print(subSetReactAct)
    
    if actions != None:
        actionsfiltered = subSetReactAct.loc[(subSetReactAct['actiontype']).isin(actions)]
    else:
        actionsfiltered = subSetReactAct

    actionsfiltered['doable'] = actionsfiltered.apply (lambda row: docheck(row), axis=1)  
    # print(actionsfiltered)
    return actionsfiltered

def blockdefine(actionsfiltered):
    actionswithblocks =pd.DataFrame(columns=['id', 'reaction_id_id', 'actiontype', 'actionno', 'material',
        'materialsmiles', 'materialquantity', 'materialquantityunit',
        'dropwise', 'atmosphere', 'molecularweight', 'materialimage', 'layer',
        'solvent', 'solventquantity', 'solventquantityunit',
        'numberofrepetitions', 'temperature', 'duration', 'durationunit',
        'stirringspeed', 'doable', 'blocknum', 'blockbool']) 

    for reaction in actionsfiltered['reaction_id_id'].unique():
        actions = actionsfiltered.loc[actionsfiltered['reaction_id_id'] == reaction]
        currentblocknum = 0
        currentblockbool = False
        blocklist = {}
        blocklist = [[],[]]

        for index, row in actions.iterrows():
            # print(f"action {row}")
            # print(f"index #{index}#")
            # print(f"action: {row}")
            if currentblockbool != row.loc['doable']:
                currentblocknum+=1
                currentblockbool= row.loc['doable']
                #print(f"block {currentblocknum}, {currentblockbool})")
            blocklist[0].append(currentblocknum)
            blocklist[1].append(currentblockbool)

            #print(f"\t{currentblocknum}")
            #append block column
            #append block bool column
            # print(action['actiontype'])
        actions["blocknum"] = blocklist[0]
        actions["blockbool"] = blocklist[1]
        #print(actions)
        actionswithblocks = actionswithblocks.append(actions, ignore_index=True)
        #print(blocklist)
    #print(actionswithblocks.sort_values(by=['blocknum', 'reaction_id_id', 'actionno']))
    return actionswithblocks

actionsfiltered = actionfilter(allactions, actions=None, reactionset=[1,2,6])
actionsfiltered = blockdefine(actionsfiltered)



for blocknum in actionsfiltered['blocknum'].unique():
    #print(f"block num \n\t{blocknum}")
    block = actionsfiltered[actionsfiltered['blocknum'] == blocknum]
    #print(block)
    if block['blockbool'].values[0] == True:
        print("activeblock")
        b = otSession(f"block_{blocknum}", block, "Example Author", "example description")
    else:
        print("inactive block")

#print("test")
#a = otSession("test", 1)
# print(allactions['reaction_id_id'].unique())
# for reactionnumber in allactions['reaction_id_id'].unique():
#     print("genrating: "+"ittraexample"+str(reactionnumber))
#     b = otSession("ittraexample"+str(reactionnumber), int(reactionnumber), "Example Author", "example description", [0,3])
#print(a.deck)
#print(a.actions)
#print(a.outputpath)
#print(list(a.actions.columns))
#print(a.actions['material'])
#a.ittrActions()
#a.setupPlate()

#print(a.actions.materialquantity)
