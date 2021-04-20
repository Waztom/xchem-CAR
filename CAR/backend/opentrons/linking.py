import backend.simpleread
import backend.opentrons.otWrite
import backend.opentrons.otDeck



class session():
    def __init__(actionsDF): #outputpath
        #self.deck = otDeck.Deck()
        #self.outputScript = otWrite.otScript(outputpath)
        self.actionsDF = actionsDF
        
    def loopactions():
        for action in self.actionsDF[actiontype]:
            print(action)

allactions = simpleread.getactions()
reactionActions = simpleread.getReactionActions(a,1)
session(reactionActions)