class Deck ():

    def __init__(self, index=None):
        self.deckindex = index
        self.PlateList = []
    
    def nextfreeplate(self):
        pass

    def addplate(self, location, numwells=None, platewellVolume=None):
        self.PlateList.append(Plate(self, index=len(self.PlateList), numwells=numwells,platewellVolume=platewellVolume))



class Plate ():

    def __init__(self, Deck, index=None, numwells=None, platewellVolume=None):
        self.deck = Deck
        self.plateindex = index
        self.numwells = numwells
        self.platewellVolume = None
        #self.WellList = None
        self.setupwells()

    def __str__(self):
        return 'Deck {}, Plate number {}'.format(self.deck.deckindex, self.plateindex)

    def setupwells(self):
        self.setofwells = []
        for well in range(self.numwells):
            self.setofwells.append(Well(self, wellindex=well, Volume=self.platewellVolume, VolumeUsed=0, StartSmiles="", GoalSmiles=""))
        return self.setofwells


    def nextfreewell(self):
        for i in self.setofwells:
            for j in i:
                if self.setofwells[i,j].VolumeUsed == 0:
                    print(str(i)+","+str(j))
                    return [i,j]

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
        self.Volume = Volume
        self.VolumeUsed = VolumeUsed
        self.StartSmiles = StartSmiles
        self.GoalSmiles = GoalSmiles
        
    def __repr__(self):
        return 'D{}P{}W{}'.format(self.plate.deck.deckindex, self.plate.plateindex, self.wellindex)

    def __str__(self):
        return 'Deck {}, Plate {}, Well {}'.format(self.plate.deck.deckindex, self.plate.plateindex, self.wellindex)
    
    def add(self, Ammount):
        SafetyMargin = 10
        WorkingVolumeused = self.VolumeUsed+Ammount
        if WorkingVolumeused >= self.Volume*(1-(SafetyMargin/100)):
            raise NameError("the resultant value of adding "+str(Ammount)+"ul  to "+str(self.VolumeUsed)+"ul would be "+str(WorkingVolumeused)+"ul exceeding the well's volume of "+str(self.Volume)+"ul (with a saftey margin of "+str(SafetyMargin)+"%)")
        elif WorkingVolumeused < 0:
            raise NameError("the resultant value of adding "+str(Ammount)+"ul  to "+str(self.VolumeUsed)+"ul would be "+str(WorkingVolumeused)+"ul")
        else:
            self.VolumeUsed = WorkingVolumeused
        return self.VolumeUsed

    def changesmiles(self, start_smiles=None, goal_smiles=None):
        if type(start_smiles) != None:
            if type(start_smiles) == str:
                self.StartSmiles = start_smiles
        if type(goal_smiles) != None:
            if type(goal_smiles) == str:
                self.GoalSmiles = goal_smiles
        return [self.StartSmiles, self.GoalSmiles]
