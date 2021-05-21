# this is vunrable to python injection by the lack of checking of metadata inputs
# this opens and closes files frequently, could be improved by creating string to hold the file data before writing to file once

import os

class PlateTxt():
    def __init__(self, Plate, protocolName=None, author=None, block=None):
        self.plate = Plate
        self.protocolName = protocolName    
        self.author = author
        self.block = block
        self.filepath = f"../output/OrderPlates/{self.plate.plateName}_{protocolName}_{block}.txt"
        self.setupScript()

    def dirsetup(self, path="../output/OrderPlates"):
        if not os.path.exists(path):
            os.makedirs(path)

    def setupScript(self, indextype = "1D"):
        self.dirsetup()
        """This is vunrable to injection atacks """
        if self.protocolName == None:
                self.protocolName = input("Please name the Protocol Name: \t")
        if self.author == None:
                self.author = input("Please name the author Name: \t")

        script = open(self.filepath, "w")
        script.write(f"# plate to order for {self.protocolName}, plate: {self.plate.plateName}\n")
        script.write("#")
        if self.block != None:
            script.write(f"block: {self.block}, ")
        script.write(f"deck position: {self.plate.plateIndex}")
        
        script.write("## "+str(self.protocolName)+" for \""+str(self.author)+str("\" produced by XChem Car (https://car.xchem.diamond.ac.uk)"))
        script.write(f"\n# for a {self.plate.numwells} well plate with wells of volume: {self.plate.platewellVolume} ({self.plate.plateTypeName})")
        
        welnum = 1

        script.write("\n\n####################################################\n")
        for well in self.plate.WellList:
            if well.VolumeUsed > 0:
                if indextype == "1D":
                    script.write(f"{well.wellindex}\t")
                script.write(f"{well.StartSmiles}\t")
                script.write(f"\t{well.VolumeUsed}ul of {well.Volume}ul\t")
            
                script.write("\n####################################################\n")
            welnum += 1


        script.close()
    
   