class Deck ():

    def __init__(self, index=None):
        self.deckindex = index
        self.PlateList = []

    def __repr__(self):
        return 'D{}'.format(self.deckindex)

    def __str__(self):
        return 'Deck {}'.format(self.deckindex)

    def __len__(self):
        return len(self.PlateList)

    def __getitem__(self, position):
        return self.PlateList[position]
    
    def nextfreeplate(self):
        freeplates = []
        for plate in self.PlateList:
            if plate.isempty():
                freeplates.append(plate.plateIndex)
        if len(freeplates)>=1:
            nextfreeplate= int(freeplates[0])
        else:
            nextfreeplate = "no free wells"
        return nextfreeplate

    def addplate(self, location, numwells=None, platewellVolume=None):
        self.PlateList.append(Plate(self, index=len(self.PlateList)+1, numwells=numwells,platewellVolume=platewellVolume))
        return self.PlateList[-1]



class Plate ():

    def __init__(self, Deck, index=None, numwells=None, platewellVolume=None):
        self.deck = Deck
        self.plateIndex = index
        self.numwells = numwells
        self.platewellVolume = platewellVolume
        self.WellList = None
        self.setupwells()
        self.strName = "plate_"+str(numwells)

    def __repr__(self):
        return 'D{}P{}'.format(self.deck.deckindex, self.plateIndex)

    def __str__(self):
        return 'Deck {}, Plate number {}'.format(self.deck.deckindex, self.plateIndex)

    def __len__(self):
        return len(self.WellList)

    def __getitem__(self, position):
        return self.WellList[position]
    
    def printplate(self):
        for well in self.WellList:
            print(str(well.StartSmiles) + " " + str(well.GoalSmiles) + " " + str(well.VolumeUsed)+ "/" + str(well.Volume))

    def setupwells(self):
        self.WellList = []
        for well in range(self.numwells):
            self.WellList.append(Well(self, wellindex=well, Volume=self.platewellVolume, VolumeUsed=0, StartSmiles="", GoalSmiles=""))
        return self.WellList


    def nextfreewell(self):
        freewells = []
        for well in self.WellList:
            if well.isempty():
                freewells.append(well.wellindex)
        if len(freewells)>=1:
            nextfreewell = int(freewells[0])
        else:
            nextfreewell = "no free wells"
        return nextfreewell

    def isempty(self):
        isempty = True
        for well in self.WellList:
            if not well.isempty():
                isempty = False
        return isempty

    def smilesearch(self, smiles, start_smiles=None, goal_smiles=None):
        if start_smiles == None & goal_smiles == None:
            start_smiles = True
            goal_smiles = True
        startinstances = []
        goalinstances = []
        if start_smiles == True:
            pass
            #start_smiles.append([i,j])
        if goal_smiles == True:
            pass
            #goalinstances.append([i,j])
        return[startinstances, goalinstances]


class Well ():
    
    def __init__(self, Plate, wellindex=None, Volume=None, VolumeUsed=0, StartSmiles=None, GoalSmiles=None):
        self.plate = Plate
        self.wellindex = wellindex
        self.Volume = float(Volume)
        self.VolumeUsed = float(VolumeUsed)
        self.StartSmiles = StartSmiles
        self.GoalSmiles = GoalSmiles
        
    def __repr__(self):
        return 'D{}P{}W{}'.format(self.plate.deck.deckindex, self.plate.plateIndex, self.wellindex)

    def __str__(self):
        return 'Deck {}, Plate {}, Well {}'.format(self.plate.deck.deckindex, self.plate.plateIndex, self.wellindex)
    
    def add(self, Ammount, smiles=""):
        SafetyMargin = 5 # percentage of the well's volume to be keept empty to prevent overvlow
        WorkingVolumeused = self.VolumeUsed+Ammount
        if WorkingVolumeused >= float(self.Volume)*(1-(SafetyMargin/100)):
            raise NameError("the resultant value of adding "+str(Ammount)+"ul  to "+str(self.VolumeUsed)+"ul would be "+str(WorkingVolumeused)+"ul exceeding the well's volume of "+str(self.Volume)+"ul (with a saftey margin of "+str(SafetyMargin)+"%)")
            self.VolumeUsed = False
        elif WorkingVolumeused < 0:
            raise NameError("the resultant value of adding "+str(Ammount)+"ul to "+str(self.VolumeUsed)+"ul would be "+str(WorkingVolumeused)+"ul")
            self.VolumeUsed = False
        else:
            self.VolumeUsed = WorkingVolumeused
            if smiles != "":
                self.changesmiles(start_smiles=smiles)
        return self.VolumeUsed
    
    def isempty(self):
        if self.VolumeUsed <= 0:
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