#currently to run from CAR
from pathlib import Path
import re
#import os,sys
#sys.path.insert(0,'..')
import ibmRead

import opentrons.otWrite as otWrite
import opentrons.otDeck as otDeck




class otSession():
    def __init__(self, name):
        
        self.name = name
        self.namecheck()
        
        self.deck = otDeck.Deck()
        
        self.outputpath = "None"
        self.pathnamecheck()
        
        self.output = otWrite.otScript(filepath=self.outputpath, protocolName=self.name)
        self.output.setup()
    
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
                print("Debug, trying: "+str(trialfile))
                suffixnumber += 1
            
        self.outputpath = trialfile
        return self.outputpath




print("test")
a = otSession("test")
print(a.deck)
print(a.outputpath)