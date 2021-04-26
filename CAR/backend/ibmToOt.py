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
    def __init__(self, name, reactionnumber=1):

        self.name = name
        self.namecheck()
        

        self.actions = []

        self.deck = otDeck.Deck()
        
        self.outputpath = "None"
        self.pathnamecheck()
        
        self.getactions(reactionnumber)
        self.tipoptions = []
        self.tipsneeded = {}
        self.tipRackList = []

        self.pipettesneeded = {}

        self.output = otWrite.otScript(filepath=self.outputpath, protocolName=self.name)
        self.output.setupScript()
        self.setupPlate()
        self.choosetip(190)
        self.choosetip(200)
        self.choosetip(201)
        self.tipOutput()
        self.output.setupLabwear(self.deck.PlateList, self.tipRackList)
        self.pipettesChoose()
        self.output.setupPipettes(self.pipettesneeded)
    
    def namecheck(self):
        trialname = self.name
        if re.match("^[\w-]+$", trialname):
            print("\""+str(trialname)+"\" is valid")
        else:
            print("\""+ str(trialname)+"\" is NOT valid")
            trialname = re.sub(r"^[\w-]+$", '', trialname)
            if trialname == '':
                trialname = "unamaed"
            print("new name: "+ str(trialname))

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

    def ittrActions(self):
        # note: add checks to prevent errors if df of diffrent format or blank is presented
        actionsCoppy = self.actions
        actionsCoppy = actionsCoppy.sort_values(['id', 'actionno'], ascending=(True, False))
        for actionindex in range(len(self.actions)):
            print(actionindex)
            currentactionmask = self.actions['actionno'] == actionindex+1
            currentaction = self.actions[currentactionmask]
            self.processAction(currentaction)

    def setupPlate(self):
        startingSteps = self.actions.loc[(self.actions['actiontype']) == "add"]
        materials = startingSteps[['materialsmiles', 'materialquantity']]
        print(materials)

        self.maxtransfer = materials['materialquantity'].max()
        self.totalvolume = materials['materialquantity'].sum()

        self.orderPlate = self.deck.add("Plate", 1,12,500, "OrderPlate")

        plate = ['','','','','','','','','',''] #note: need to fix so not fixed posibles
        for i in range(len(materials)):
            wellnumber = self.orderPlate.nextfreewell()
            outcome = self.orderPlate.WellList[wellnumber].add(materials.loc[i,'materialquantity'],
                smiles = materials.loc[i,'materialsmiles'])
            if outcome != False:
                #add otWrite setup plate here
                plate[wellnumber] = [materials.loc[i,'materialquantity'], materials.loc[i,'materialsmiles']] 
        print(plate)
        
        self.reactionPlate = self.deck.add("Plate", 2, 12, 500, "ReactionPlate")



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
        print("###tipsneeded###\n"+str(self.tipsneeded)+"\n######")
        for tip in self.tipsneeded :
            print("#tip##\n"+str(tip)+"\n###")

            numplates = math.ceil(self.tipsneeded[tip]/96)

            while numplates >= 1:

                location = len(self.deck.PlateList)

                if tip == 10:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=10, platename = "opentrons_96_filtertiprack_10ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_10ul", self.deck.nextfreeplate()])
                elif tip == 20:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=20, platename = "opentrons_96_filtertiprack_20ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_20ul", self.deck.nextfreeplate()])
                elif tip == 200:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=200, platename = "opentrons_96_filtertiprack_200ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_200ul", self.deck.nextfreeplate()])
                elif tip == 300:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=300, platename = "opentrons_96_filtertiprack_300ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_300ul", self.deck.nextfreeplate()])
                elif tip == 1000:
                    self.deck.add("TipRack", location, numwells=96, platewellVolume=1000, platename = "opentrons_96_filtertiprack_1000ul")
                    self.tipRackList.append(["opentrons_96_filtertiprack_1000ul", self.deck.nextfreeplate()])

                numplates -= 1
        print("###tipracklist###\n"+str(self.tipRackList)+"\n######")
        return self.tipRackList

    def pipettesChoose (self):
        if len(self.tipsneeded) <= 2:
            if len(self.tipsneeded) > 0 :
                for pipette in self.tipsneeded:
                    self.pipettesneeded[str(pipette)] = self.deck.findPippets(pipette)
        return self.pipettesneeded



    def processAction(self, currentaction):
        currentactiontype = (currentaction['actiontype'].to_string(index=False)).strip()
        print(currentactiontype)
        if currentactiontype == " add":
            print("\tadd")
            self.actionAdd(currentaction)

    def actionAdd(self, currentaction):
        well = self.deck.PlateList[1].nextfreewell()
        print(well)
        print(currentaction)

        #elf.output.movefluids(1, fromAdress, toAdress, volume, dispenseVolume=None, writetoscript=True)

#print("test")
a = otSession("test")
#print(a.deck)
a.getactions(1)
#print(a.actions)
#print(a.outputpath)
#print(list(a.actions.columns))
#print(a.actions['material'])
#a.ittrActions()
a.setupPlate()

#print(a.actions.materialquantity)
