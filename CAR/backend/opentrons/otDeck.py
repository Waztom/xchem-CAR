import math
import pandas as pd
import numpy as np
import re
import ibmRead

class Deck ():

    def __init__(self, index=None):
        self.deckindex = index
        self.PlateList = []
        self.inputlist = []
        self.outputList = []
        self.PipetteList = []
        


    def __repr__(self):
        return 'D{}'.format(self.deckindex)

    def __str__(self):
        return 'Deck: {}'.format(self.deckindex)

    def __len__(self):
        return len(self.PlateList)

    def __getitem__(self, position):
        if type(position) == int:
            return self.PlateList[position]
        elif type(position) == str:
            for i in range(len(self.PlateList)):
                if (self.PlateList[i].plateName == position) or (self.PlateList[i].plateTypeName == position):
                    return self.PlateList[i]
        print("'"+str(position)+"' not found")
    
    def nextfreeplate(self):
        # freeplates = []
        # for plate in self.PlateList:
        #     if plate.isempty():
        #         freeplates.append(plate.plateIndex)
        # if len(freeplates)>=1:
        #     nextfreeplate= int(freeplates[0])
        # else:
        #     nextfreeplate = "no free wells"
        nextfreeplateindex = len(self.PlateList)
        #nextfreeplateindex = 1
        return nextfreeplateindex

    def add(self, Type, location, numwells=None, platewellVolume=None, platename = "",  inputPlate=False, outputPlate=False):
        # print(inputPlate)
        # print(outputPlate)
        index = self.nextfreeplate()
        DeckObject(self, Type, index = index, numwells=numwells,platewellVolume=platewellVolume, platename=platename,  inputPlate=inputPlate, outputPlate=outputPlate)
        return self.PlateList[index]

    def addDefinedPlate(self, Type, location, plateframe, numwells=None, platewellVolume=None, platename = "definedplate",  inputPlate=False, outputPlate=False):
        pass
        # index = self.nextfreeplate()
        # DeckObject(self, Type, index=index, numwells=numwells,platewellVolume=platewellVolume, platename=platename,  inputPlate=inputPlate, outputPlate=outputPlate)
        # plate = self.PlateList[index]
        # for i in plateframe.index.values:
        #     plate.WellList[plateframe.loc[i, 'WellIndex']].add(Ammount=plateframe.loc[i, 'Volume'], smiles=plateframe.loc[i, 'Contents'], solvent =plateframe.loc[i, 'Solvent'], MaterialName=plateframe.loc[i, 'ReactionID'])
        #     plate.WellList[plateframe.loc[i, 'WellIndex']].changereactionid(plateframe.loc[i, 'ReactionID'])
        #     plate.WellList[plateframe.loc[i, 'WellIndex']].changename(ibmRead.idToProduct(plate.ReactionID))
        
        # return plate

    def addDefinedFromFile(self, path, inputPlate=True, outputPlate=False):
        print(path)
        plateframe = pd.read_csv(path, sep= ';')
        print(plateframe)

        splitpath = path.split('/')[-1].split('.')[0].split('-')
        properties = {}
        for plateproperty in splitpath:
            properties[plateproperty[0]] = plateproperty[1:]
        print(f"splitpath: {splitpath}")
        print(f"properties: {properties}")

        index = self.nextfreeplate()
        print(self.PlateList)
        print(index)
        DeckObject(self, "Plate", index=index, numwells=int(properties['W']), platewellVolume=float(properties['V']), platename=properties['N'],  inputPlate=inputPlate, outputPlate=outputPlate)
        print(self.PlateList)
        plate = self.PlateList[index]
        print(plateframe)
        for i in plateframe.index.values:
            plate.WellList[plateframe.loc[i, 'WellIndex']].add(Ammount=plateframe.loc[i, 'Volume'], smiles=plateframe.loc[i, 'Contents'], solvent =plateframe.loc[i, 'Solvent'], MaterialName=plateframe.loc[i, 'ReactionID'])
            plate.WellList[plateframe.loc[i, 'WellIndex']].changereactionid(plateframe.loc[i, 'ReactionID'])
            plate.WellList[plateframe.loc[i, 'WellIndex']].changename(ibmRead.idToProduct(plateframe.loc[i, 'ReactionID'], spoof=True))
        
        return plate


    def findPippets (self, volume):
        for pipette in self.PipetteList:
            if int(pipette.volume) == int(volume):
                return pipette
            else:
                return False

    def findTipRacks (self, volume):
        print(f"searching for tips: {volume} in {self.PipetteList}")
        releventracks = []
        for plate in self.PlateList:
            if plate.platetype == "TipRack":
                if int(plate.tipVolume) == int(volume):
                    releventracks.append(plate.plateName)
        return releventracks

    def addPipette (self, name, model, mount, volume):
        print(f"adding pipette: name:{name}, model:{model}, mount:{mount}, volume{volume}")
        #self.PipetteList.append(self, pipette(len(self.PipetteList), name, model, mount, volume))
        self.PipetteList.append(Pipette(self, len(self.PipetteList), name, model, mount, volume))

    def decksearch(self, smiles, solvent, inputsearch = False, outputsearch = False):
        results = []
        if inputsearch == True:
            for plate in self.inputlist:
                results.append([plate.plateName, plate.smilesearch(smiles, start_smiles=True, start_solvent=solvent)[0]])
        # print(f"results for {smiles} and {solvent}: {results}")
        return results
    def nextplatewithspace(self, searchtype = "allplates", minvol=None):
        if searchtype == "allplates":
            for plate in self.PlateList:
                if plate.platetype == "Plate":
                    result = plate.nextfreewell( minvol=minvol)
                    if result is not False:
                        return plate
        elif searchtype == "inputs":
            for inplate in self.inputList:
                result = inplate.nextfreewell (minvol=minvol)
                if result is not False:
                    return inputplate
        elif searchtype == "outputs":
            print("searching {self.outputList}")
            for outplate in self.outputList:
                print(outplate)
                result = outplate.nextfreewell( minvol=minvol)
                print(f"result {result}")
                if result is not False:
                    return outplate
        return False

                


class DeckObject ():
    def __init__(self, Deck, Type, index=None, numwells=None, platewellVolume=None, platename="", numTips=None, tipVolume=None, inputPlate=False, outputPlate=False):
        currentplates = []
        for plate in Deck.PlateList:
            currentplates.append(plate.plateName)

        trialname = platename
        if trialname in currentplates:
            suffixnumber = 1
            while trialname in currentplates:
                trialname = (f"{trialname}_{suffixnumber}")
                #print("Debug, trying: "+str(trialname))
                suffixnumber += 1

        platename = trialname

        if Type == "Plate":
            print(f"adding plate {index}")
            plate = Plate(Deck, index, numwells, platewellVolume, platename, inputPlate, outputPlate)
            Deck.PlateList.append(plate)
                

        elif Type == "TipRack":
            Deck.PlateList.append(TipRack(Deck, index, numTips=numwells, tipVolume=platewellVolume, platename=platename))

class Plate ():

    def __init__(self, Deck, index=None, numwells=None, platewellVolume=None, platename="",  inputPlate=False, outputPlate=False):
        self.platetype = "Plate"
        self.deck = Deck
        self.plateIndex = index
        self.numwells = numwells
        self.platewellVolume = platewellVolume
        self.WellList = None
        self.setupwells()
        if inputPlate is True:
             self.deck.inputlist.append(self)
        if outputPlate is True:
            self.deck.outputList.append(self)


        #print(f"{self.numwells},{self.platewellVolume}")
        if numwells == 24:
            self.plateTypeName= "24_reservoir_2500ul"
        if numwells == 96:
            if platewellVolume == 2500:
                self.plateTypeName= "plateone_96_wellplate_2500ul"
                self.platewellVolume = 2500
            elif platewellVolume == 500:
                self.plateTypeName= "plateone_96_wellplate_500ul"
                self.platewellVolume = 500
        # print(self.plateTypeName)
        self.plateName = platename
        self.activeWell = self.nextfreewell()

    def __repr__(self):
        return 'D{}P{}'.format(self.deck.deckindex, self.plateIndex)

    def __str__(self):
        return 'Deck: {}, Plate number: {}'.format(self.deck.deckindex, self.plateIndex)

    def __len__(self):
        return len(self.WellList)

    def __getitem__(self, position):
        return self.WellList[position]

    def printself(self):
        print(f"{self.deck} - {self.plateIndex} ({self.plateTypeName})")
        print(f"\t{self.platetype}")
        print(f"\t{self.numwells} wells of volume {self.platewellVolume}ul")
        # print(f"\twell list:    {self.WellList}")
        print(f"\tactive well: {self.activeWell}")
        print(f"\tPlate name:  {self.plateName}")


    def dimentionShift(self, index, RowLetters = True, numrows = 8):
        if self.numwells == 96:
            numrows = 8
        elif self.numwells == 24:
            numrows = 4



        row = ((index)%numrows)+1
        #col = math.ceil((index+2)*10/self.numwells) #BUG: not correctly calculating column number
        col = math.floor((1/numrows)*index)+1

        if RowLetters == True:
            characters = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','T','U','V','W','X','Y','Z']
            row = characters[(row-1)%numrows]
        # print(f"{row},{col}")
        return [row, col]
        
    def printplate(self):
        if self.numwells == 96:
            welnum = 1
            for well in self.WellList:
                if welnum % 8 == 0:
                    print(f"{well.StartSmiles} {well.GoalSmiles} {well.CurrentVolume}/{well.VolumeAvailable}", end="\n")
                else:
                    print(f"{well.StartSmiles} {well.GoalSmiles} {well.CurrentVolume}/{well.VolumeAvailable}", end="\t")
                welnum += 1
        else:
            for well in self.WellList:
                print(f"{well.StartSmiles} {well.GoalSmiles} {well.CurrentVolume}/{well.VolumeAvailable}\t")

    def setupwells(self):
        self.WellList = []
        for well in range(self.numwells):
            self.WellList.append(Well(self, wellindex=well, VolumeAvailable=self.platewellVolume, CurrentVolume=0, StartSmiles="", GoalSmiles="", MaterialName=""))
        return self.WellList

    def addplacehoder(self):
        pass
        #self.WelList(self.nextfreewell()) = ""
        


    def nextfreewell(self, minvol=None):
        freewells = None
        freewells = []
        for well in self.WellList:
            if well.isempty():
                if well.reactionid == None:
                    if minvol == None:
                        freewells.append(well.wellindex)
                    else:
                        if round(well.platewellVolume, 3) >= round(minvol, 3):
                            freewells.append(well.wellindex)
        if len(freewells)>=1:
            nextfreewell = int(freewells[0])
        else:
            nextfreewell = False
        return nextfreewell

    def isempty(self):
        isempty = True
        for well in self.WellList:
            if not well.isempty():
                isempty = False
            elif well.reactionid != None:
                isempty = False
        return isempty

    def smilesearch(self, smiles, start_smiles=None, start_solvent= None, goal_smiles=None):
        if (start_smiles == None) & (goal_smiles == None):
            start_smiles = True
            goal_smiles = True
        startinstances = []
        goalinstances = []
        for well in self.WellList:
            if start_smiles == True:
                if well.StartSmiles == smiles:
                    if well.StartSolvent == start_solvent or start_solvent == '*' or (well.StartSolvent == None and str(start_solvent) == "nan"):
                        startinstances.append({'index':well.wellindex,'volume':well.CurrentVolume})
            if goal_smiles == True:
                if well.StartSmiles == smiles:
                    startinstances.append(well.wellindex)
        return[startinstances, goalinstances]

    def saveplate(self, saving=True):
        plateFrame = pd.DataFrame(columns = ['WellIndex', 'ReactionID', 'Name', 'Contents', 'Volume', 'VolumemUnits', 'Solvent'])
        iterlist = self.WellList
        for well in iterlist:
            plateFrame = plateFrame.append({'WellIndex':well.wellindex, 'ReactionID':well.reactionid,'Name':well.MaterialName, 'Contents':well.StartSmiles, 'Volume':well.CurrentVolume, 'VolumemUnits':'ul', 'Solvent':well.StartSolvent}, ignore_index=True)
            if saving == True:
                path = f"temp/plates/P{self.plateIndex}-T{self.plateTypeName}-W{self.numwells}-V{self.platewellVolume}-N{self.plateName}.csv"
                plateFrame.to_csv(path, sep=';', index=False) 
        
        if saving==True:
            return path
        else:
            return plateFrame




    
class TipRack ():
    
    def __init__(self, Deck, index=None, numTips=None, tipVolume=None, platename=""):
        self.platetype = "TipRack"
        self.deck = Deck
        self.plateIndex = index
        self.numTips = numTips
        self.tipVolume = tipVolume
        self.TipList = []
        self.setupTips()
        self.plateTypeName= platename
        self.plateName = f"tips_{self.numTips}_{self.tipVolume}_{self.plateIndex}"

    def __repr__(self):
        return 'D{}T{}'.format(self.deck.deckindex, self.plateIndex)

    def __str__(self):
        return self.plateName

    def __len__(self):
        return len(self.TipList)

    def __getitem__(self, position):
        return self.TipList[position]
    
    def setupTips (self):
        for i in range(self.numTips):
            self.TipList.append(True)
        return self.TipList

    def nextFreeTip (self):
        for i in range(len(self.TipList)):
            if self.TipList[i] == True:
                return i
        return False

    def useTip (self, index):
        self.TipList[index] = False

class Well ():
    
    def __init__(self, Plate, wellindex=None, VolumeAvailable=None, CurrentVolume=0, StartSmiles=None, StartSolvent = None, GoalSmiles=None, MaterialName = None, reactionid = None):
        self.plate = Plate
        self.wellindex = wellindex
        self.VolumeAvailable = float(VolumeAvailable)
        self.CurrentVolume = float(CurrentVolume)
        self.StartSmiles = StartSmiles
        self.StartSolvent = StartSolvent
        self.GoalSmiles = GoalSmiles
        self.MaterialName = "debug"#MaterialName
        self.SafetyMargin = 5 # percentage of the well's volume to be keept empty to prevent overvlow
        self.reactionid = reactionid

        
    def __repr__(self):
        return 'D{}P{}W{}'.format(self.plate.deck.deckindex, self.plate.plateIndex, self.wellindex)

    def __str__(self):
        return 'Deck: {}, Plate: {}, Well: {}'.format(self.plate.deck.deckindex, self.plate.plateIndex, self.wellindex)
    
    def add(self, Ammount, smiles="", solvent = "", MaterialName=""):
        WorkingCurrentVolume = self.CurrentVolume+Ammount
        if WorkingCurrentVolume >= float(self.VolumeAvailable)*(1-(self.SafetyMargin/100)):
            raise NameError("the resultant value of adding "+str(Ammount)+"ul  to "+str(self.CurrentVolume)+"ul would be "+str(WorkingCurrentVolume)+"ul exceeding the well's volume of "+str(self.VolumeAvailable)+"ul (with a saftey margin of "+str(self.SafetyMargin)+"%)")
            self.CurrentVolume = False
        elif WorkingCurrentVolume < 0:
            raise NameError("the resultant value of adding "+str(Ammount)+"ul to "+str(self.CurrentVolume)+"ul would be "+str(WorkingCurrentVolume)+"ul")
            self.CurrentVolume = False
        else:
            self.CurrentVolume = WorkingCurrentVolume
            if smiles != "":
                self.changesmiles(start_smiles=smiles)
                self.changesolvent(start_solvent=solvent)
                self.changename(name=MaterialName)
        return self.CurrentVolume
    
    def isempty(self):
        if self.CurrentVolume <= 0:
            isempty = True
        else:
            isempty = False
        return isempty


    def changesmiles(self, start_smiles=None, goal_smiles=None):
        if type(start_smiles) != None:
            if type(start_smiles) == str:
                self.StartSmiles = start_smiles
        if type(goal_smiles) != None:
            if type(goal_smiles) == str:
                self.GoalSmiles = goal_smiles
        return [self.StartSmiles, self.GoalSmiles]

    def changesolvent(self, start_solvent):
        if type(start_solvent) != None:
            if type(start_solvent) == str:
                self.StartSolvent = start_solvent
    
    def changename(self, name):
        if type(name) != None:
            print(f"old name:{self.MaterialName}, new name: {name}")
            self.MaterialName = name
            print(f"new name:{self.MaterialName}, new name: {name}")
        
    def changereactionid(self, rId):
        if type(rId) != None:
            print(f"old name:{self.reactionid}, new name: {rId}")
            self.reactionid = rId
            print(f"new name:{self.reactionid}, new name: {rId}")
    
    def changevolume(self, volumeincrease):
        if type(volumeincrease) != None:
            if type(volumeincrease) == int or type(volumeincrease) == float or isinstance(volumeincrease, np.float64) is True:
                newvolume = self.CurrentVolume + volumeincrease
                if newvolume >= float(self.VolumeAvailable)*(1-(self.SafetyMargin/100)):
                    raise NameError("the resultant value of adding "+str(volumeincrease)+"ul  to "+str(self.CurrentVolume)+"ul would be "+str(newvolume)+"ul exceeding the well's volume of "+str(self.VolumeAvailable)+"ul (with a saftey margin of "+str(SafetyMargin)+"%)")
                    return True
                elif newvolume < 0:
                    if round(newvolume,2) < 0:
                        raise NameError("the resultant value of adding "+str(volumeincrease)+"ul to "+str(self.CurrentVolume)+"ul would be "+str(newvolume)+"ul")
                        return False
                    else:
                        self.CurrentVolume = 0
                        return True
                else:
                    if round(newvolume, 4) == 0:
                        newvolume = 0 
                    self.CurrentVolume = newvolume
                    return True
    


class Pipette():
    def __init__(self, deck, pipetteIndex, name, model, mount, volume):
        self.deck = deck
        self.pipetteIndex = pipetteIndex
        self.name = name
        self.model = model
        self.mount = mount
        self.volume = volume
        self.tipRacks = self.deck.findTipRacks(self.volume)
