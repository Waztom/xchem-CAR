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
    def __init__(self, name, reactionnumber=1, author=None, description=None,):

        self.name = name
        self.namecheck()

        self.author=author
        self.description=description,
        

        self.actions = []

        self.deck = otDeck.Deck()
        
        self.outputpath = "None"
        self.pathnamecheck()
        
        self.getactions(reactionnumber)
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

    def getactions(self, reactionno):
        allactions = ibmRead.getactions()
        reactionsactions = ibmRead.getReactionActions(allactions, reactionno)
        self.actions = reactionsactions
        return reactionsactions

    def setupPlate(self):
        startingSteps = self.actions.loc[(self.actions['actiontype']) == "add"]
        materials = startingSteps[['materialsmiles', 'materialquantity']]

        self.maxtransfer = materials['materialquantity'].max()
        self.totalvolume = materials['materialquantity'].sum()

        self.orderPlate = self.deck.add("Plate", 1,12,500, "OrderPlate")

        plate = ['','','','','','','','','',''] #note: need to fix so not fixed posibles
        for i in materials.index.values:
            wellnumber = self.orderPlate.nextfreewell()
            outcome = self.orderPlate.WellList[wellnumber].add(materials.loc[i,'materialquantity'], smiles = materials.loc[i,'materialsmiles'])
            if outcome != False:
                #add otWrite setup plate here
                plate[wellnumber] = [materials.loc[i,'materialquantity'], materials.loc[i,'materialsmiles']] 

        self.reactionPlate = self.deck.add("Plate", 2, 12, 500, "ReactionPlate")

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
                self.tipoptions = [10, 20, 200, 300, 1000]
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
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=10, platename = "opentrons_96_filtertiprack_10ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_10ul", self.deck.nextfreeplate(), tip])
                elif tip == 20:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=20, platename = "opentrons_96_filtertiprack_20ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_20ul", self.deck.nextfreeplate(), tip])
                elif tip == 200:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=200, platename = "opentrons_96_filtertiprack_200ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_200ul", self.deck.nextfreeplate(), tip])
                elif tip == 300:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=300, platename = "opentrons_96_filtertiprack_300ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_300ul", self.deck.nextfreeplate(), tip])
                elif tip == 1000:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=1000, platename = "opentrons_96_filtertiprack_1000ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_1000ul", self.deck.nextfreeplate(), tip])

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
            self.deck.addPipette(str(str(mount)+"_"+str(pipette)), str("p"+str(pipette)+"_single"), mount, pipette)
            

    def ittrActions(self):
        # note: add checks to prevent errors if df of diffrent format or blank is presented
        # actionsCoppy = self.actions
        # actionsCoppy = actionsCoppy.sort_values(['id', 'actionno'], ascending=(True, False))
       # print(self.actions)
        for actionindex in range(len(self.actions)):
            currentactionmask = self.actions['actionno'] == actionindex+1
            currentaction = self.actions[currentactionmask]
            #print("current: "+str(currentaction))
            self.processAction(currentaction)

    def processAction(self, currentaction):
        currentactiontype = (currentaction['actiontype'].to_string(index=False)).strip()
        #print("current type: "+str(currentactiontype))
        numreps = currentaction['numberofrepetitions'].values
        if str(numreps) == "[nan]":
            numreps = 0
        elif str(numreps) == "[None]":
            numreps = 0
        elif str(numreps) == "[]":
            numreps = 0

        #print("numrepsa: "+str(numreps))
        numreps = int(numreps)

        #print(currentaction)
        #print("numreps: "+str(numreps))


        repetitions = 0
        while repetitions  <= numreps:
            if currentactiontype == "add":
                self.actionAdd(currentaction)
            if currentactiontype == "collect-layer":
                self.actionCollectLayer(currentaction)
            if currentactiontype == "wash":
                self.actionWash(currentaction)
            elif currentactiontype == "stir":
                self.output.unsuportedAction("stir at a "+str(currentaction['stirringspeed'].values[0])+" speed at "+str(currentaction['temperature'].values[0])+" Celsius for "+str(currentaction['duration'].values[0])+" "+str(currentaction['durationunit'].values[0]))
            elif currentactiontype == "set-temprature":
                self.output.unsuportedAction("set temprature to "+str(currentaction['temperature'].values[0])+" Celsius for"+str(currentaction['duration'].values[0])+" "+str(currentaction['durationunit'].values[0]))
            elif currentactiontype == "store":
                self.output.unsuportedAction("store product ("+str(currentaction['material'].values[0])+")")
            elif currentactiontype == "extract":
                self.output.unsuportedAction("extract ("+str(currentaction['material'].values[0])+")")
            elif currentactiontype == "concentrate":
                self.output.unsuportedAction("concentrate ("+str(currentaction['material'].values[0])+")")
            else:
                self.output.unsuportedAction(str(currentaction['actiontype'].values)+" is not currently supported ")
            repetitions +=1
        

    def actionAdd(self, currentaction):
        print(currentaction)
        tipvolume = self.choosetip(currentaction['materialquantity'].values[0])
        print(tipvolume)
        pipetteName = self.deck.findTipRacks(tipvolume)
        print(pipetteName)
        pipetteName = pipetteName[0]

        self.output.movefluids( pipetteName, 
            str("OrderPlate"+str(self.deck['OrderPlate'].smilesearch(currentaction['materialsmiles'].values[0], start_smiles = True)[0])+""),
            str("ReactionPlate["+str(self.deck['ReactionPlate'].activeWell)+"]"),
            currentaction['materialquantity'].values[0],
            dispenseVolume=None,
            writetoscript=True)

    def actionWash(self, currentaction):
        self.choosetip(currentaction['solventquantity'].values[0])
        self.output.movefluids(1,
            str("OrderPlate"+str(self.deck['OrderPlate'].smilesearch(currentaction['solvent'].values[0], start_smiles = True)[0])+""),
            str("ReactionPlate["+str(self.deck['ReactionPlate'].activeWell)+"]"),
            currentaction['solventquantity'].values[0],
            dispenseVolume=None,
            writetoscript=True)

    def actionCollectLayer(self, currentaction):
        self.output.unsuportedAction("collect layer not supported ")
allactions = ibmRead.getactions()
# print(allactions)
# print(allactions.columns)
#reactionsactions = ibmRead.getReactionActions(allactions, reactionno)
print(allactions['actiontype'].unique())

#print("test")
#a = otSession("test", 1)
# print(allactions['reaction_id_id'].unique())
for reactionnumber in allactions['reaction_id_id'].unique():
    print("genrating: "+"ittraexample"+str(reactionnumber))
    b = otSession("ittraexample"+str(reactionnumber), int(reactionnumber), "Example Author", "example description")
#print(a.deck)
#print(a.actions)
#print(a.outputpath)
#print(list(a.actions.columns))
#print(a.actions['material'])
#a.ittrActions()
#a.setupPlate()

#print(a.actions.materialquantity)
