import os
import pandas as pd


class PlateCSV:
    def __init__(
        self,
        platelist,
        protocolname=None,
        author=None,
        block=None,
        indextype="2D",
        volumeMultiplyer=1.1,
    ):

        """setup for, then run, script to generate an ordering csv for a specified plate

        :param indextype: 2 character description of plate indexing ("2D" takes form [letter][number] as printed on plate, "1D" counts wells increasing by on takes form [number]), defults to "2D"
        :type indextype: str
        :param volumeMultiplyer: the volume required is multiplied by the volume mulitplyer to enure the full vloume is avalible (avoides some surface tension and evaporation issues, defults to 1.1
        :type volumeMultiplyer: float
        ...

        ...
        :return: Pandas dataframe wih columns: "mculeid", "platename", "well","concentration", "solvent", "amount-ul"
        :rtype: Pandas.DataFrame

        """
        self.platelist = platelist
        self.protocolName = protocolname
        self.author = author
        self.block = block
        self.indextype = indextype
        self.volumeMultiplyer = volumeMultiplyer
        self.filepath = f"../output/OrderPlates/{self.plate.plateName}_{protocolname}_{block}.CSV"
        self.setupScript()

    def dirsetup(self, path="../output/OrderPlates"):
        if not os.path.exists(path):
            os.makedirs(path)

    def setupScript(self):

        """loop through wells in plate identifying the relevent properties and write to .CSV File
        ...

        ...
        :return: Pandas dataframe wih columns: "mculeid", "platename", "well","concentration", "solvent", "amount-ul"
        :rtype: Pandas.DataFrame

        """

        self.dirsetup()

        outputdf = pd.DataFrame(
            columns=["mculeid", "platename", "well", "concentration", "solvent", "amount-ul"]
        )
        for plate in self.platelist:
            for well in self.plate.WellList:
                if well.VolumeUsed > 0:
                    wellproperties = {
                        "mculeid": well.mculeid,
                        "platename": plate.plateName,
                        "well": self.getindex(well),
                        "concentration": well.concentration,
                        "solvent": well.StartSolvent,
                        "amount-ul": well.VolumeUsed * self.volumeMultiplyer,
                    }
                    outputdf = outputdf.append(wellproperties, ignore_index=True)

        outputdf.to_csv(self.filepath)

    def getindex(self, well):
        """obtain the well number (1D) or alphanumeric index (2D)

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
            wellindex = f"{alphanumericindex[0]}{alphanumericindex[1]}"
        return wellindex
