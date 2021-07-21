import os
import pandas as pd

# can be used in ibmtoot by importing 
# import ordering.OutputPlateCSV as OutputPlateCSV
# and calling:
# OutputPlateTxt.PlateCSV(self.orderplate, self.name, self.author, currentblocknum)


class PlateCSV():
    def __init__(self, Plate, protocolName=None, author=None, block=None, indextype = "2D", volumeMultiplyer = 1.1):

        """ setup for, then run, script to generate an ordering csv for a specified plate
        
        :param indextype: 2 character description of plate indexing ("2D" takes form [letter][number] as printed on plate, "1D" counts wells increasing by on takes form [number]), defults to "2D"
        :type indextype: str
        :param volumeMultiplyer: the volume required is multiplied by the volume mulitplyer to enure the full vloume is avalible (avoides some surface tension and evaporation issues, defults to 1.1
        :type volumeMultiplyer: float
        ...

        ...
        :return: Pandas dataframe wih columns: "mculeid", "platename", "well","concentration", "solvent", "amount-ul"
        :rtype: Pandas.DataFrame
        
        """
        # WTOSCR: index type and volume multiplyer to be moved to _init_ arguments
        # WTOSCR: move volume multiplyer to IBMTOT aggrigation of all materials
        self.plate = Plate
        self.protocolName = protocolName    
        self.author = author
        self.block = block
        self.indextype = indextype
        self.volumeMultiplyer = volumeMultiplyer
        self.filepath = f"../output/OrderPlates/{self.plate.plateName}_{protocolName}_{block}.CSV"
        output = self.setupScript()
        return output

    def dirsetup(self, path="../output/OrderPlates"):
        if not os.path.exists(path):
            os.makedirs(path)

    def setupScript(self):

        """ loop through wells in plate identifying the relevent properties and write to .CSV File
        ...

        ...
        :return: Pandas dataframe wih columns: "mculeid", "platename", "well","concentration", "solvent", "amount-ul"
        :rtype: Pandas.DataFrame
        
        """

        self.dirsetup()     

        outputdf = pd.dataframe(columns=["mculeid", "platename", "well","concentration", "solvent", "amount-ul"])
        
        welnum = 1
        for well in self.plate.WellList:
            if well.VolumeUsed > 0:
                wellproperties = {
                    "mculeid":      self.getmculeid(well),
                    "platename":    "placeholder plate name",
                    "well":         self.getindex(well, self.indextype),
                    "concentration":self.getconcentration(),
                    "solvent":      self.getsolvent(well),
                    "amount-ul":    self.getvolume(well, self.volumeMultiplyer)
                }
            outputdf.append(wellproperties, ignore_index=True)
            welnum += 1

        outputdf.to_csv(self.filepath)
        print(outputdf)


    def getmculeid(self, well):
        """ function to pass obtain the relevent mcule ID
        can pass material name or start smiles to mculeid converter
        
        :param well: a well from the ot deck model that is to be reperesented in the order CSV
        :type well: CAR.Backend.OTDeck.well
        ...
        ...
        :return mculeid: Mcule ID
        :rtype mculeid: str
        """

        if well.MaterialName != None and well.MaterialName != "":
            materialname = well.MaterialName
            materialsmiles = well.StartSmiles
            mculeid = self.mculelookup(materialname) # or:  = self.mculelookup(materialsmiles)
            return mculeid

    def mculelookup(self, materialname):
        """ PLACEHOLDER to convert from a material to it's mcule id
        
        :param materialname: the name or smiles of a material to be ordered
        :type materialname: str
        ...
        ...
        :return mculeid: Mcule ID
        :rtype mculeid: str
        """
        # this should instead be held in a database with the material name proberbly
        mculeid = materialname # >> placehodler << 
        return mculeid

    def getindex(self, well):
        """ obtain the well number (1D) or alphanumeric index (2D)

        :param well: a well from the ot deck model that is to be reperesented in the order CSV
        :type well: CAR.Backend.OTDeck.well
        ...
        ...
        :return wellindex: the well's unique identifyer on the plate
        :rtype wellindex: str or (int if 1D)
        """
        if self.indextype == "1D":
            wellindex = well.wellindex
        elif self.indextype == "2D":
            alphanumericindex = self.plate.dimentionShift(well.wellindex)
            wellindex = (f"{alphanumericindex[0]}{alphanumericindex[1]}")
        return wellindex

    
    def getconcentration(self):
        """ PLACEHOLDER to obtain concentration required

        ...
        ...
        :return concentration: the concentration required
        :rtype concentration: float
        """
        concentration = 1.0 # >> placeholder <<
        return concentration

    def getsolvent(self, well):
        """ function to pass obtain solvent required
        
        :param well: a well from the ot deck model that is to be reperesented in the order CSV
        :type well: CAR.Backend.OTDeck.well
        ...
        ...
        :return solvent: name of solvent required
        :rtype solvent: str
        """
        solvent = well.StartSolvent
        return solvent

    def getvolume(self, well):
        """ function to pass obtain the volume of soloution required
        
        :param well: a well from the ot deck model that is to be reperesented in the order CSV
        :type well: CAR.Backend.OTDeck.well
        ...
        ...
        :return volumeneeded: volume from deck model, multiplied by volume multiplyer 
        :rtype volumeneeded: float

        """
        volumeneeded = well.VolumeUsed*self.volumeMultiplyer
        return volumeneeded

    
   
