from pathlib import Path
import pandas as pd
import math
import re
#import os,sys
#sys.path.insert(0,'..')
import ibmRead

import opentrons.otWrite as otWrite
import opentrons.otDeck as otDeck
import Ordering.OutputPlateTxt as OutputPlateTxt
import HumanRead.HumanRead as HumanRead




class otSession(): # WTOSCR: otsession could be renamed to otsessionblock or similar
    """A Class used to cordinate the conversion from the *actions* passed into a script to to exicute on the Robot
    

    :param name: A name to be used for this run
    :type name: str
    :param actions: the set of all actions to be exicuted in this run
    :type actions: Pandas DataFrame
    :param author: the name of the individual generating the script, defults to "None"
    :type author: str, optional
    :param description: a description of the protocol being genrated, defults to "None"
    :type description: str, optional
    :param split: NOT CURRENTLY FUNCTIONAL - a definition of a subset of the actions passed to implement, defults to "None"
    :type split: , optional

    """

    def __init__(self, name, actions, author=None, description=None, split=None, currentpipettesetup = [["Left", 300, "multi", "p300_multi_gen2"],["Right", 300, "single", "p300_single"]]):

        # setup name
        self.name = name # WTOSCR: currently using block name, would make sense to also include refrence to project name and block number
        self.namecheck()

        # setup protocol metadata
        self.author=author
        self.description=description
        
        # create blank actions list
        self.actions = actions

        #decalre defined pipettes if they exist
        # WTOSCR: could be improved by having currentpipettesetup defult to none
        if self.currentpipettesetup != None:
            self.definedpipettes = currentpipettesetup
        else:
            self.definedpipettes = None

        # Instantiate Deck to model physical constraints and locations involved in protocol
        self.deck = otDeck.Deck()

        #setup name for robot py file        
        self.outputpath = "None" # WTOSCR: create and save script in DJango model
        self.pathnamecheck() # generate unique name

        # Create and setup main robotic .py script to write output to
        self.output = otWrite.otScript(filepath=self.outputpath, protocolName=self.name, author=self.author, description=self.description)
        self.output.setupScript()
        
        # set create blank variables to be used for hardware setup
        self.tipoptions = []
        self.tipsneeded = {}
        self.tipRackList = []
        self.pipettesneeded = []
        
        # create list of plates to place plates to use in reactions
        self.setupPlate()
        
        # Generate idea of what tips will be used in protocol
        self.preselecttips()
        
        #
        self.tipOutput()
        
        # write hardwear setup to robotics .py script
        self.output.setupLabwear(self.deck.PlateList, self.tipRackList)

        self.setupPipettes()
        self.output.setupPipettes(self.deck.PipetteList)

        # legacy code to run only subset of selected actions form actions
        #if split == None:
        #    self.ittrActions()
        #else:
        #    self.splitittract(split)

        # loop through all listeded actions 
        self.ittrActions()


        # WTOSCR: add qc.text file at end of session once output plates are filled

    
    def namecheck(self):
        """checks if self.name is file name safe 

        :return: trialname: a name to use for the protocol
        :rtype: trialname: str
        """
        # creates trialname to hold the current idea of what the name should be
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
        """a function to ensure the output file name is unique (acheived though itterating suffix number till name is unique)

        :return: a unique filepaht for the genrated script
        :rtype: str
        """
        self.namecheck()

        #create trialfile to hold current idea of what the name should be
        trialfile = Path("../output/Opentrons/" + self.name + ".py")

        if trialfile.exists():
            # if the filename exists loop thorugh the suffix number increasing by one each time till the name is not matched
            suffixnumber = 1
            while trialfile.exists(): # WTOSCR: use unique filename/filepath in django model
                trialfile = Path("../output/Opentrons/" + self.name + "(" + str(suffixnumber) + ")" + ".py")
                suffixnumber += 1

        # set the filepath to trialfile and return filepath
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

    def combinestrings(self, row):
        print(f"{row['materialsmiles']}\t{row['solvent']})")
        return str(str(row['materialsmiles'])+str(row['solvent']))
    def setupPlate(self, inputnumwells = 24, inputwellvolume = 2500, incAddMaterials=True, incWashMaterials=True, incExtract=True):
        # WTOSCR: currently relying on user input of number of wells and input well volume, a function should be created along the line of the pipette selection process to automaticaly choose which plate should be used
        # WTOSCR: do we need autogenrated plate size selection - proberbly more for output plate but the same funciton should work for both
        #print(self.actions)
        
        allmaterials = pd.DataFrame(columns=['material', 'materialsmiles' 'materialquantity', 'solvent', 'materialKey']) # WTOSCR: should proberbly check about a Django model

        # each block getting a list of materials associated with that action type and concatinating them to the pandas Dataframe, All materials
        if incAddMaterials == True:
            addingSteps = self.actions.loc[(self.actions['actiontype']) == "add"] 
            materials = addingSteps[['material','materialsmiles', 'materialquantity', 'solvent']]
            key = addingSteps[['materialsmiles']]
            key = key.rename(columns={'materialsmiles':'materialKey'}) # material key is smiles
            materials = pd.concat([materials, key], axis = 1)
            print(materials)
            allmaterials = pd.concat([allmaterials, materials], ignore_index=True)

        if incWashMaterials==True:  # WTOSCR: review if properly implemented and bring inline with add 
            solventSteps = self.actions.loc[(self.actions['actiontype']) == "wash"]
            washmaterials = solventSteps[['material', 'materialquantitcombinestringsquantity':'materialquantity'})
            allmaterials = pd.concat([allmaterials, extractsolvents], ignore_index=True)
            
        allmaterials['matAndSolv'] = allmaterials.apply(lambda row: self.combinestrings(row), axis=1)  
        allmaterials = allmaterials.groupby(['matAndSolv']).aggregate({"materialquantity":"sum", "material":'first', 'solvent':'first', 'materialsmiles':'first'})
        allmaterials = allmaterials.sort_values(['solvent','materialquantity'], ascending = False)
        #experement with adding a working materials+sovlent collum to group by then drop that col
        
        
        self.maxvolume = allmaterials['materialquantity'].max() # not currently used
        self.totalvolume = allmaterials['materialquantity'].sum() # not currently used

        self.orderPlate = self.deck.add("Plate", 1, inputnumwells, inputwellvolume, "OrderPlate") #WTOSCR: index is hard coded

        # plate = ['','','','','','','','','',''] #note: need to fix so not fixed posibles
        for i in allmaterials.index.values:
            if allmaterials.loc[i, 'material'] == "" or allmaterials.loc[i, 'material'] == None or allmaterials.loc[i, 'material'] == 'NaN':
                matName = allmaterials.loc[i, 'materialsmiles'] # implement chemical name check for all reactants and reagents
            else:
                matName = allmaterials.loc[i, 'material']
            wellnumber = self.orderPlate.nextfreewell()
            outcome = self.orderPlate.WellList[wellnumber].add(allmaterials.loc[i,'materialquantity'], smiles = allmaterials.loc[i,'materialsmiles'], solvent = allmaterials.loc[i, 'solvent'], MaterialName=matName)
        #print(self.orderPlate.printplate())
        currentblocknum = self.actions['blocknum'].values[0]
        OutputPlateTxt.PlateTxt(self.orderPlate, self.name, self.author, currentblocknum )
        # WTOSCR: add qc.text file at end of session once output plates are filled
        # if well overvlows (outcome == False) needs to move to next weel
            # if outcome != False:
            #     plate[wellnumber] = [allmaterials.loc[i,'materialquantity'], allmaterials.loc[i,'material']] 
        # WTOSCR: add handleing materials spread across many wells in fluid handeling steps

        self.reactionPlate = self.deck.add("Plate", 2, 96, 500, "ReactionPlate") #WTOSCR: index is hard coded
        # for well in range(len(self.orderPlate)):
        #     print(self.orderPlate[well].StartSmiles)
        # print(self.orderPlate.smilesearch('C(Cl)Cl'))
        # print(self.orderPlate.smilesearch('brine'))


    def preselecttips(self):
        for index, row in self.actions.iterrows():
            currentactiontype = row["actiontype"]
            # currentactionmask = self.actions['actionno'] == actionindex+1
            # # print(currentactionmask)
            # currentaction = self.actions[currentactionmask]
            # print(currentaction)
            # currentactiontype = (currentaction['actiontype'].to_string(index=False)).strip()
            # print(currentactiontype)
            if currentactiontype == "add": # or wash or extract
                self.choosetip(row['materialquantity'], tipoptions = [])

    def choosetip(self, volume, tipoptions= [], overideselection= True):
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
            else:
                if overideselection == True:
                    if tip in self.tipsneeded:
                        self.tipsneeded[tip]+=1
                    else:
                        self.tipsneeded[tip] = 1
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
        print("setitng up pipettes")
        print(self.tipsneeded)
        if self.definedpipettes == None:
            if len(self.tipsneeded) <= 2:
                if len(self.tipsneeded) > 0 :
                    for pipette in self.tipsneeded:
                        self.pipettesneeded.append(pipette)
            else:
                print("pipetteDebug")
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
        else:
            for defined in self.definedpipettes:
                if defined[2] != "multi":
                    self.deck.addPipette(str(f"{defined[0]}_{defined[1]}_{defined[2]}_pipette"), defined[3], defined[0], defined[1])
                # WTOSCR: need to add handleing for multichannle pipettes


    def ittrActions(self):
        #print(self.actions)
        reactionids = self.actions['reaction_id_id'].unique().tolist() # WTOSCR: Django implementaiton
        for Index in range(len(self.actions.index.values)):
            columns  = self.actions.columns.values
            currentaction = pd.DataFrame(index = [Index], columns = columns)
            for col in columns:
                currentaction[col] = self.actions.iloc[Index][col]
            currentaction["outputwell"] = reactionids.index(currentaction["reaction_id_id"].values)

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
        currentactiontype = currentaction['actiontype'].to_string(index=False)
        currentactiontype = currentactiontype.strip()
        if currentactiontype != 'Series([], )':
            # WTOSCR: what to do if it is equal to series([],)
            # WTOSCR: check django compatability

            numreps = currentaction['numberofrepetitions'].values
            if str(numreps) == "[nan]":
                numreps = 1
            elif str(numreps) == "[None]":
                numreps = 1
            elif str(numreps) == "[]":
                numreps = 1
            elif numreps == 0:
                numreps = 1

            try:
                numreps = int(numreps)
            except:
                numreps = 1


            repetitions = 0
            while repetitions  < numreps:
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
        #print("add")
        tipvolume = self.choosetip(currentaction['materialquantity'].values[0])
        print(f"tipvolume {tipvolume}")
        pipetteName = (self.deck.findPippets(tipvolume)).name

        self.output.transferfluids( pipetteName, 
            (f"OrderPlate.wells(){self.deck['OrderPlate'].smilesearch(currentaction['materialsmiles'].values[0], start_smiles = True, start_solvent=currentaction['solvent'].values[0])[0]}"),
            (f"ReactionPlate.wells()[{currentaction['outputwell'].values[0]}]"),
            currentaction['materialquantity'].values[0])

    def actionWash(self, currentaction):
        print("wash")
        tipvolume = self.choosetip(currentaction['materialquantity'].values[0])
        print(tipvolume)
        pipetteName = (self.deck.findPippets(tipvolume)).name
        print("wash debug")
        print(pipetteName)
        pipetteName = pipetteName[0]

        self.output.transferfluids(pipetteName,
            (f"OrderPlate.wells()[{self.deck['OrderPlate'].smilesearch(currentaction['material'].values[0], start_smiles = True, start_solvent=currentaction['solvent'].values[0])[0]}]"),
            (f"ReactionPlate.wells()[{currentaction['outputwell'].values[0]}]"),
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


# WTOSCR: Actual start of script run, is it worth making another class here?
    
#allactions = ibmRead.getactions() #activate to enable workng with front end
allactions = pd.read_csv("../../debuging/for-Olivia-actions-final-test-3.csv", index_col=0, sep = ';')


# print(allactions.columns)
#reactionsactions = ibmRead.getReactionActions(allactions, reactionno)
#print(allactions['actiontype'].unique())

def docheck(row):
    if row['actiontype'] in ['add', 'wash', "extract"]:
        return True
    else:
        return False

def actionfilter(allactions,  actions=None, reactionset=None): # WTOSCR: needs splittting into filter and checking for doability
    if reactionset != None:
        subSetReactAct = allactions.loc[(allactions['reaction_id_id']).isin(reactionset)]
    else:
        subSetReactAct = allactions
    
    if actions != None:
        actionsfiltered = subSetReactAct.loc[(subSetReactAct['actiontype']).isin(actions)]
    else:
        actionsfiltered = subSetReactAct

    actionsfiltered['doable'] = actionsfiltered.apply(lambda row: docheck(row), axis=1)
    return actionsfiltered

def blockdefine(actionsfiltered):
    # WTOSCR: 1) add doable to action models,  2) retreive data from database

    actionswithblocks =pd.DataFrame(columns=['id', 'reaction_id_id', 'actiontype', 'actionno', 'material',
        'materialsmiles', 'materialquantity', 'materialquantityunit',
        'dropwise', 'atmosphere', 'molecularweight', 'materialimage', 'layer',
        'solvent', 'solventquantity', 'solventquantityunit',
        'numberofrepetitions', 'temperature', 'duration', 'durationunit',
        'stirringspeed', 'doable', 'blocknum', 'blockbool']) 

    for reaction in actionsfiltered['reaction_id_id'].unique(): #WTOSCR: check if two undoables produce 1 or two blocks
        actions = actionsfiltered.loc[actionsfiltered['reaction_id_id'] == reaction]
        currentblocknum = 0
        currentblockbool = False
        blocklist = {}
        blocklist = [[],[]]

        for index, row in actions.iterrows():
            if currentblockbool != row.loc['doable']:
                currentblocknum+=1
                currentblockbool= row.loc['doable']
            blocklist[0].append(currentblocknum) #WTOSCR: should be dictionary [[int,bool],[int,bool],[int,bool]]
            blocklist[1].append(currentblockbool)

        actions["blocknum"] = blocklist[0] #WTOSCR: check if can be done with django modles
        actions["blockbool"] = blocklist[1]
        actionswithblocks = actionswithblocks.append(actions, ignore_index=True)
    return actionswithblocks

actionsfiltered = actionfilter(allactions, actions=None)
actionsfiltered = blockdefine(actionsfiltered)


protocolOut = HumanRead.HumanReadable("../output/protocols/example.md")
protocolOut.setupDoc()

# loop through each block in actions filtered
for blocknum in actionsfiltered['blocknum'].unique():  
    
    print(f"block num \t{blocknum}")
    block = actionsfiltered[actionsfiltered['blocknum'] == blocknum]
    blockbool = None # WTOSCR: blockbool seems to be unused, should it be removed?
    if block['blockbool'].values[0] == True:
        blockbool = True # WTOSCR: blockbool seems to be unused, should it be removed?
        print("activeblock")
        blockSession = otSession(f"block_{blocknum}", block, "Example Author", "example description")
    else:
        blockbool = False # WTOSCR: blockbool seems to be unused, should it be removed?
        print("inactive block")
        blockSession = None


    protocolOut.newBlock(blocknum, block['blockbool'].values[0], blockSession) # human read stuff

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
