from pathlib import Path
import pandas as pd
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
        self.tipsneeded = [[],[]]
        self.tipRackList = []

        self.output = otWrite.otScript(filepath=self.outputpath, protocolName=self.name)
        self.output.setupScript()
        self.setupPlate()
        self.output.setupLabwear(self.deck.PlateList, self.tipRackList)
    
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

        self.orderPlate = self.deck.addplate(1,12,500)

        plate = ['','','','','','','','','',''] #note: need to fix so not fixed posibles
        for i in range(len(materials)):
            wellnumber = self.orderPlate.nextfreewell()
            outcome = self.orderPlate.WellList[wellnumber].add(materials.loc[i,'materialquantity'],
                smiles = materials.loc[i,'materialsmiles'])
            if outcome != False:
                #add otWrite setup plate here
                plate[wellnumber] = [materials.loc[i,'materialquantity'], materials.loc[i,'materialsmiles']] 
        print(plate)
        
        self.reactionPlate = self.deck.addplate(2, 12, 500)

        self.orderPlate.printplate()
        print()
        self.reactionPlate.printplate()


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
                    neededindex = self.tipsneeded[0].index(tip)
                    self.tipsneeded[1][neededindex] += self.tipsneeded[1][neededindex]
                else:
                    self.tipsneeded[0].append(tip)
                    self.tipsneeded[1].append(1)
                return tip 
        return False
    
    #def setuptips

    def tipOutput(self):
        self.tipRackList = []
        for tip in self.tipsneeded :
            if tip == 10:
                self.tipRackList.append("opentrons_96_filtertiprack_10ul")
            elif tip == 20:
                self.tipRackList.append("opentrons_96_filtertiprack_20ul")
            elif tip == 200:
                self.tipRackList.append("opentrons_96_filtertiprack_200ul")
            elif tip == 300:
                self.tipRackList.append("opentrons_96_tiprack_300ul")
            elif tip == 1000:
                self.tipRackList.append("opentrons_96_filtertiprack_1000ul")
        return self.tipRackList




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

print("test")
a = otSession("test")
#print(a.deck)
a.getactions(1)
#print(a.actions)
#print(a.outputpath)
print(list(a.actions.columns))
#print(a.actions['material'])
#a.ittrActions()
a.setupPlate()

#print(a.actions.materialquantity)
