"""Create OT session"""
from __future__ import annotations
from pathlib import Path
from django.db.models.expressions import Col
import pandas as pd
import math
import re

from otwrite import otScript
from otdeck import Deck
from platesavailable import labware_plates

from ..models import (
    IBMAddAction,
    Product,
    Project,
    Target,
    Method,
    Reaction,
    OTSession,
    Deck,
    Plate,
    Well,
)

import mcule.outputplatecsv as OutputPlateCSV

# Create session
# Create deck
# create plate for order materials
# fill wells in order plate
# AND search plates for exisiting products to use in next step
# create order csv


class SortSessions(object):
    def __init__(self, projectid: int):
        self.projectid = projectid
        self.materialstoorder = []

    def getMethods(self, targetid):
        methodqueryset = Method.objects.filter(target_id=targetid)
        return methodqueryset

    def getReactions(self, methodid):
        reactionqueryset = Reaction.objects.filter(method_id=methodid)
        return reactionqueryset

    def getAddActions(self, reactionid):
        addactionqueryset = IBMAddAction.objects.filter(reaction_id=reactionid)
        return addactionqueryset

    def checkSessions(self):
        targetqueryset = Target.objects.filter(project_id=self.projectid)

        for target in targetqueryset:
            methodqueryset = self.getMethods(targetid=target.id)
            for method in methodqueryset:
                reactionqueryset = self.getReactions(methodid=method.id)


class CreateOTSession(object):
    """
    Creates a CreateOTSession object for generating a protocol
    from actions
    """

    def __init__(self, projectid: int, reactionqueryset: DjangoModel, platequeryset: list = None):
        self.projectid = projectid
        self.reactionqueryset = reactionqueryset
        self.otsessionid = self.createOTSessionModel()
        self.deckobj = self.createDeckModel()
        self.platequeryset = platequeryset
        if not platequeryset:
            self.createOrderPlate()

    def createOTSessionModel(self):
        otsessionobj = OTSession()
        project_obj = Project.objects.get(id=self.projectid)
        otsessionobj.project_id = project_obj
        otsessionobj.save()
        return otsessionobj.id

    def createDeckModel(self):
        deckobj = Deck()
        otsession_obj = OTSession.objects.get(id=self.otsessionid)
        deckobj.otsession_id = otsession_obj
        deckobj.save()
        return deckobj

    def updateDeckSlotAvailable(self):
        testslotavailable = self.deckobj.indexslotavailable + 1
        if testslotavailable <= self.deckobj.numberslots:
            self.deckobj.indexslotavailable = testslotavailable
            self.deckobj.save()
            return self.deckobj.indexslotavailable
        else:
            self.deckobj.slotavailable = False
            self.deckobj.save()
            return False

    def createPlateModel(self, platename, platetype, numberwells):
        testindexslot = self.updateDeckSlotAvailable()
        if testindexslot:
            plateindex = testindexslot
            plateobj = Plate()
            plateobj.deck_id = self.deck
            plateobj.platename = "{}-{}".format(platename, testindexslot)
            plateobj.plateindex = plateindex
            plateobj.labware = platetype
            plateobj.numberwells = numberwells
            plateobj.save()
            return plateobj
        else:
            print("No more deck slots available")

    def updatePlateWellsAvailable(self, plateobj):
        testwellavailable = plateobj.indexwellavailable + 1
        if testwellavailable <= plateobj.numberwells:
            plateobj.indexwellavailable = testwellavailable
            plateobj.save()
            return plateobj.indexwellavailable
        else:
            plateobj.wellavailable = False
            plateobj.save()
            return False

    def createWellModel(self, plateobj, reactionobj, addactionobj, volume):
        testwellslot = self.updatePlateWellsAvailable(plateobj=plateobj)
        if testwellslot:
            wellobj = Well()
            wellobj.plate_id = plateobj
            wellobj.reaction_id = reactionobj
            wellobj.addaction_id = addactionobj
            wellobj.wellindex = testwellslot
            wellobj.volume = volume
            wellobj.save()
            return wellobj
        else:
            print("No more well slots available")
            return False

    def getAddActions(self, reactionobj):
        addactionqueryset = IBMAddAction.objects.filter(reaction_id=reactionobj.id)
        return addactionqueryset

    def createOrderPlate(self):
        allactionslistdf = []
        firstreactionobjs = [reaction[0] for reaction in self.reactionqueryset]
        alladdactionqueryset = [
            self.getAddActions(reactionobj) for reactionobj in firstreactionobjs
        ]
        orderdf = pd.DataFrame(
            columns=["mculeid", "platename", "well", "concentration", "solvent", "amount-ul"]
        )

        for addactionqueryset in alladdactionqueryset:
            actions_to_add_df = pd.DataFrame(list(addactionqueryset.values()))
            if not actions_to_add_df.empty:
                allactionslistdf.append(actions_to_add_df)
        allactionslistdf = pd.concat(allactionslistdf)

        startingmaterialsdf = allactionslistdf.groupby(["id"]).aggregate(
            {
                "materialquantity": "sum",
                "material": "first",
                "solvent": "first",
                "materialsmiles": "first",
                "mculeid": "first",
                "concentration": "first",
            }
        )
        startingmaterialsdf = startingmaterialsdf.sort_values(
            ["solvent", "materialquantity"], ascending=False
        )

        plateobj = self.createPlateModel(
            platename="Orderplate", labware="24_reservoir_2500ul", numberwells=24
        )

        maxwellvolume = float(plateobj.platetype.split("_")[-1].strip("ul"))

        for i in startingmaterialsdf.index.values:
            totalvolume = startingmaterialsdf.at[i, "materialquantity"]
            if totalvolume > maxwellvolume:
                nowellsneeded = -(-totalvolume // maxwellvolume)
                for index in range(nowellsneeded):
                    if not plateobj.wellavailable:
                        plateobj = self.createPlateModel(
                            platename="Orderplate", platetype="24_reservoir_2500ul", numberwells=24
                        )
                    self.createWellModel(
                        plateobj=plateobj, reactionobj=None, addactionobj=None, volume=maxwellvolume
                    )
                    orderdf.at[i, "mculeid"] = (startingmaterialsdf.at[i, "mculeid"],)
                    orderdf.at[i, "platename"] = plateobj.platename
                    orderdf.at[i, "well"] = wellobj.wellindex
                    orderdf.at[i, "concentration"] = (startingmaterialsdf.at[i, "concentration"],)
                    orderdf.at[i, "solvent"] = (startingmaterialsdf.at[i, "solvent"],)
                    orderdf.at[i, "amount-ul"] = (startingmaterialsdf.at[i, "materialquantity"],)

            else:
                if not plateobj.wellavailable:
                    plateobj = self.createPlateModel(
                        platename="Orderplate", platetype="24_reservoir_2500ul", numberwells=24
                    )
                wellobj = self.createWellModel(
                    plateobj=plateobj, reactionobj=None, addactionobj=None, volume=totalvolume
                )
                orderdf.at[i, "mculeid"] = (startingmaterialsdf.at[i, "mculeid"],)
                orderdf.at[i, "platename"] = plateobj.platename
                orderdf.at[i, "well"] = wellobj.wellindex
                orderdf.at[i, "concentration"] = (startingmaterialsdf.at[i, "concentration"],)
                orderdf.at[i, "solvent"] = (startingmaterialsdf.at[i, "solvent"],)
                orderdf.at[i, "amount-ul"] = (startingmaterialsdf.at[i, "materialquantity"],)


class CollectActions(object):
    """
    Creates a CollectActions object for extracting action steps from the
    DB
    """

    def __init__(self, projectid: int):
        """
        ValidateFile constructor
        Args:
            project (int): Project id for genrating automated protocols
        """
        self.projectid = projectid
        self.actionmodels = [
            backend.models.IBMAddAction,
            backend.models.IBMCollectLayerAction,
            backend.models.IBMConcentrateAction,
            backend.models.IBMDegasAction,
            backend.models.IBMDrySolidAction,
            backend.models.IBMDrySolutionAction,
            backend.models.IBMExtractAction,
            backend.models.IBMFilterAction,
            backend.models.IBMMakeSolutionAction,
            backend.models.IBMPartitionAction,
            backend.models.IBMpHAction,
            backend.models.IBMPhaseSeparationAction,
            backend.models.IBMQuenchAction,
            backend.models.IBMRefluxAction,
            backend.models.IBMSetTemperatureAction,
            backend.models.IBMStirAction,
            backend.models.IBMStoreAction,
            backend.models.IBMWaitAction,
            backend.models.IBMWashAction,
        ]

        self.allactions_df = self.getActions()
        self.actionsfiltered = self.actionfilter()
        self.blockdefine()
        self.startingmaterials = self.getStartingMaterials()

    def getActions(self):
        allactions_list_df = []
        targets = backend.models.Target.objects.filter(project_id=self.projectid)
        methods = [backend.models.Method.objects.filter(target_id=target.id) for target in targets]
        reactions = [
            backend.models.Reaction.objects.filter(method_id=method[0].id) for method in methods
        ]

        for actionmodel in self.actionmodels:
            for reaction in reactions:
                actions_to_add_df = pd.DataFrame(
                    list(actionmodel.objects.filter(reaction_id=reaction[0].id).values())
                )
                if not actions_to_add_df.empty:
                    allactions_list_df.append(actions_to_add_df)
        allactions_list_df = pd.concat(allactions_list_df)

        return allactions_list_df.sort_values(["reaction_id_id", "actionno"])

    def docheck(self, row):
        if row["actiontype"] in ["add", "wash", "extract"]:
            return True
        else:
            return False

    def actionfilter(
        self, actions=None, reactionset=None
    ):  # WTOSCR: needs splittting into filter and checking for doability

        if reactionset != None:
            subSetReactAct = self.allactions_df.loc[
                (self.allactions_df["reaction_id_id"]).isin(reactionset)
            ]
        else:
            subSetReactAct = self.allactions_df

        if actions != None:
            actionsfiltered = subSetReactAct.loc[(subSetReactAct["actiontype"]).isin(actions)]
        else:
            actionsfiltered = subSetReactAct

        actionsfiltered["doable"] = actionsfiltered.apply(lambda row: self.docheck(row), axis=1)

        return actionsfiltered

    def blockdefine(self):
        # WTOSCR: 1) add doable to action models,  2) retreive data from database

        actionswithblocks = pd.DataFrame(
            columns=[
                "id",
                "reaction_id_id",
                "actiontype",
                "actionno",
                "material",
                "materialsmiles",
                "materialquantity",
                "materialquantityunit",
                "dropwise",
                "atmosphere",
                "molecularweight",
                "materialimage",
                "layer",
                "solvent",
                "solventquantity",
                "solventquantityunit",
                "numberofrepetitions",
                "temperature",
                "duration",
                "durationunit",
                "stirringspeed",
                "doable",
                "blocknum",
                "blockbool",
            ]
        )

        for reaction in self.actionsfiltered[
            "reaction_id_id"
        ].unique():  # WTOSCR: check if two undoables produce 1 or two blocks
            actions = self.actionsfiltered.loc[self.actionsfiltered["reaction_id_id"] == reaction]
            currentblocknum = 0
            currentblockbool = False
            blocklist = {}
            blocklist = [[], []]

            for index, row in actions.iterrows():
                if currentblockbool != row.loc["doable"]:
                    currentblocknum += 1
                    currentblockbool = row.loc["doable"]
                blocklist[0].append(
                    currentblocknum
                )  # WTOSCR: should be dictionary [[int,bool],[int,bool],[int,bool]]
                blocklist[1].append(currentblockbool)

            actions["blocknum"] = blocklist[0]
            actions["blockbool"] = blocklist[1]
            actionswithblocks = actionswithblocks.append(actions, ignore_index=True)
        self.actionsfiltered = actionswithblocks

    def getStartingMaterials(
        self,
        incAddMaterials=True,
    ):
        """ "
        Groups starting materials for easier sorting into orderplates
        """
        startingmaterials = pd.DataFrame(
            columns=[
                "mculeid",
                "concentration",
                "material",
                "materialsmiles",
                "materialquantity",
                "solvent",
                "materialKey",
            ]
        )

        if incAddMaterials == True:
            addingSteps = self.actionsfiltered.loc[(self.actionsfiltered["actiontype"]) == "add"]
            materials = addingSteps[
                [
                    "mculeid",
                    "concentration",
                    "material",
                    "materialsmiles",
                    "materialquantity",
                    "solvent",
                ]
            ]
            key = addingSteps[["materialsmiles"]]
            key = key.rename(columns={"materialsmiles": "materialKey"})  # material key is smiles
            materials = pd.concat([materials, key], axis=1)
            startingmaterials = pd.concat([startingmaterials, materials], ignore_index=True)

        startingmaterials["matAndSolv"] = startingmaterials.apply(
            lambda row: self.combinestrings(row), axis=1
        )
        startingmaterials = startingmaterials.groupby(["matAndSolv"]).aggregate(
            {
                "materialquantity": "sum",
                "material": "first",
                "solvent": "first",
                "materialsmiles": "first",
                "mculeid": "first",
                "concentration": "first",
            }
        )
        startingmaterials = startingmaterials.sort_values(
            ["solvent", "materialquantity"], ascending=False
        )

        return startingmaterials


class otSession:  # WTOSCR: otsession could be renamed to otsessionblock or similar
    """A Class used to cordinate the conversion from the *actions* passed into a script to to exicute on the Robot

    :param projectid: Project id from database
    :type name: int
    :param name: A name to be used for this run
    :type name: str
    :param actions: the set of all actions to be exicuted in this run
    :type actions: Pandas DataFrame
    :param author: the name of the individual generating the script, defults to "None"
    :type author: str, optional
    :param description: a description of the protocol being genrated, defults to "None"
    :type description: str, optional
    :param split: NOT CURRENTLY FUNCTIONAL - a definition of a subset of the actions passed to implement, defults to "None"
    :type split: , optional

    """

    def __init__(
        self,
        name,
        actions,
        startingmaterials,
        author=None,
        description=None,
        currentpipettesetup=[
            ["Left", 300, "multi", "p300_multi_gen2"],
            ["Right", 300, "single", "p300_single"],
        ],
    ):

        # setup name
        self.name = name  # WTOSCR: currently using block name, would make sense to also include refrence to project name and block number
        self.namecheck()

        # setup protocol metadata
        self.author = author
        self.description = description

        # create blank actions list
        self.actions = actions
        self.startingmaterials = startingmaterials
        self.currentpipettesetup = currentpipettesetup

        # decalre defined pipettes if they exist
        # WTOSCR: could be improved by having currentpipettesetup defult to none
        if self.currentpipettesetup != None:
            self.definedpipettes = currentpipettesetup
        else:
            self.definedpipettes = None

        # Instantiate Deck to model physical constraints and locations involved in protocol
        self.deck = Deck()

        # setup name for robot py file
        self.outputpath = "None"  # WTOSCR: create and save script in DJango model
        self.pathnamecheck()  # generate unique name

        # Create and setup main robotic .py script to write output to
        self.output = otScript(
            filepath=self.outputpath,
            protocolName=self.name,
            author=self.author,
            description=self.description,
        )
        self.output.setupScript()

        # set create blank variables to be used for hardware setup
        self.tipoptions = []
        self.tipsneeded = {}
        self.tipRackList = []
        self.pipettesneeded = []

        # create list of plates to place plates to use in reactions
        self.setupOrderPlate()

        # Generate idea of what tips will be used in protocol
        self.preselecttips()

        #
        self.tipOutput()

        # write hardwear setup to robotics .py script
        self.output.setupLabware(self.deck.PlateList, self.tipRackList)

        self.setupPipettes()
        self.output.setupPipettes(self.deck.PipetteList)
        self.ittrActions()

        # WTOSCR: add qc.text file at end of session once output plates are filled

    def startProtocol(self):
        for blocknum in self.actions["blocknum"].unique():
            actionsblock = self.actions[self.actions["blocknum"] == blocknum]
            if actionsblock["blockbool"].values[0] == True:
                otSession(
                    name=f"block_{blocknum}",
                    actions=actionsblock,
                    startingmaterials=self.startingmaterials,
                    author="Example Author",
                    description="example description",
                )

    def namecheck(self):
        """checks if self.name is file name safe

        :return: trialname: a name to use for the protocol
        :rtype: trialname: str
        """
        # creates trialname to hold the current idea of what the name should be
        trialname = self.name

        if re.match("^[\w-]+$", trialname):
            pass
        else:
            trialname = re.sub(r"^[\w-]+$", "", trialname)
            if trialname == "":
                trialname = "unamaed"

        self.name = trialname
        return trialname

    def pathnamecheck(self):
        """a function to ensure the output file name is unique (acheived though itterating suffix number till name is unique)

        :return: a unique filepaht for the genrated script
        :rtype: str
        """
        self.namecheck()

        # create trialfile to hold current idea of what the name should be
        trialfile = Path("../output/Opentrons/" + self.name + ".py")

        if trialfile.exists():
            # if the filename exists loop thorugh the suffix number increasing by one each time till the name is not matched
            suffixnumber = 1
            while trialfile.exists():  # WTOSCR: use unique filename/filepath in django model
                trialfile = Path(
                    "../output/Opentrons/" + self.name + "(" + str(suffixnumber) + ")" + ".py"
                )
                suffixnumber += 1

        # set the filepath to trialfile and return filepath
        self.outputpath = trialfile
        return self.outputpath

    def combinestrings(self, row):
        print(f"{row['materialsmiles']}\t{row['solvent']})")
        return str(str(row["materialsmiles"]) + str(row["solvent"]))

    def addPlate(self, platetype: str):
        if platetype == "Order":
            name = "OrderPlate"
            numwells = 24
            platewellvolume = 2500
        if platetype == "Reaction":
            name = "ReactionPlate"
            numwells = 96
            platewellvolume = 2500

        next_free_plate_index = self.deck.nextfreeplate()
        plate_name = "{}-{}".format(name, next_free_plate_index)

        self.deck.add(
            Type="Plate",
            numwells=numwells,
            platewellVolume=platewellvolume,
            platename=plate_name,
        )

        plate = [plate for plate in self.deck.PlateList if plate.plateName == plate_name][0]
        return plate

    def setupOrderPlate(
        self,
    ):
        orderplate = self.addPlate(platetype="Order")

        for i in self.startingmaterials.index.values:
            if (
                self.startingmaterials.at[i, "material"] == ""
                or self.startingmaterials.at[i, "material"] == None
                or self.startingmaterials.at[i, "material"] == "NaN"
            ):
                matName = self.startingmaterials.at[
                    i, "materialsmiles"
                ]  # implement chemical name check for all reactants and reagents
            else:
                matName = self.startingmaterials.at[i, "material"]

            orderplate.nextfreewell()

            if orderplate.isplatefull:
                orderplate = self.addPlate(platetype="Order")

            if not self.orderplate.isplatefull:
                wellnumber = orderplate.nextfreewellindex

                self.orderplate.WellList[wellnumber].add(
                    amount=self.startingmaterials.at[i, "materialquantity"],
                    smiles=self.startingmaterials.at[i, "materialsmiles"],
                    solvent=self.startingmaterials.at[i, "solvent"],
                    materialname=matName,
                    mculeid=self.startingmaterials.at[i, "mculeid"],
                    concentration=self.startingmaterials.at[i, "concentration"],
                )

        currentblocknum = self.actions["blocknum"].values[0]

        OutputPlateCSV.PlateCSV(
            platelist=self.deck.PlateList,
            protocolname=self.name,
            author=self.author,
            block=currentblocknum,
        )

    def preselecttips(self):
        for index, row in self.actions.iterrows():
            currentactiontype = row["actiontype"]
            if currentactiontype == "add":  # or wash or extract
                self.choosetip(row["materialquantity"], tipoptions=[])

    def choosetip(self, volume, tipoptions=[], overideselection=True):
        if tipoptions == []:
            if self.tipoptions == []:
                # self.tipoptions = [10, 20, 200, 300, 1000]
                self.tipoptions = [300]  # debuging
        else:
            pass
            # review wether to use tip options or self.tipoptions for rest of function
        for tip in self.tipoptions:
            if tip >= volume:
                if tip in self.tipsneeded:
                    self.tipsneeded[tip] += 1
                else:
                    self.tipsneeded[tip] = 1
                return tip
            else:
                if overideselection == True:
                    if tip in self.tipsneeded:
                        self.tipsneeded[tip] += 1
                    else:
                        self.tipsneeded[tip] = 1
                return tip
        return False

    def tipOutput(self):
        self.tipRackList = []
        for tip in self.tipsneeded:

            numplates = math.ceil(self.tipsneeded[tip] / 96)

            while numplates >= 1:

                if tip == 10:
                    self.deck.add(
                        "TipRack",
                        numwells=96,
                        platewellVolume=10,
                        platename="opentrons_96_tiprack_10ul",
                    )
                    self.tipRackList.append(
                        ["opentrons_96_tiprack_10ul", self.deck.nextfreeplate(), tip]
                    )
                elif tip == 20:
                    self.deck.add(
                        "TipRack",
                        numwells=96,
                        platewellVolume=20,
                        platename="opentrons_96_tiprack_20ul",
                    )
                    self.tipRackList.append(
                        ["opentrons_96_tiprack_20ul", self.deck.nextfreeplate(), tip]
                    )
                elif tip == 200:
                    self.deck.add(
                        "TipRack",
                        numwells=96,
                        platewellVolume=200,
                        platename="opentrons_96_rtiprack_200ul",
                    )
                    self.tipRackList.append(
                        ["opentrons_96_tiprack_200ul", self.deck.nextfreeplate(), tip]
                    )
                elif tip == 300:
                    self.deck.add(
                        Type="TipRack",
                        numwells=96,
                        platewellVolume=300,
                        platename="opentrons_96_tiprack_300ul",
                    )
                    self.tipRackList.append(
                        ["opentrons_96_tiprack_300ul", self.deck.nextfreeplate(), tip]
                    )
                elif tip == 1000:
                    self.deck.add(
                        "TipRack",
                        numwells=96,
                        platewellVolume=1000,
                        platename="opentrons_96_tiprack_1000ul",
                    )
                    self.tipRackList.append(
                        ["opentrons_96_tiprack_1000ul", self.deck.nextfreeplate(), tip]
                    )

                numplates -= 1
        return self.tipRackList

    def setupPipettes(self):
        print("setitng up pipettes")
        print(self.tipsneeded)
        if self.definedpipettes == None:
            if len(self.tipsneeded) <= 2:
                if len(self.tipsneeded) > 0:
                    for pipette in self.tipsneeded:
                        self.pipettesneeded.append(pipette)
            else:
                print("pipetteDebug")
            mountnumber = 0
            for pipette in self.pipettesneeded:
                if mountnumber == 0:
                    mount = "left"
                    mountnumber += 1
                elif mountnumber == 1:
                    mount = "right"
                else:
                    break
                self.deck.addPipette(
                    str(f"{mount}_{pipette}_pipette"),
                    str("p" + str(pipette) + "_single"),
                    mount,
                    pipette,
                )
        else:
            for defined in self.definedpipettes:
                if defined[2] != "multi":
                    self.deck.addPipette(
                        str(f"{defined[0]}_{defined[1]}_{defined[2]}_pipette"),
                        defined[3],
                        defined[0],
                        defined[1],
                    )
                # WTOSCR: need to add handleing for multichannle pipettes

    def ittrActions(self):
        reactionids = self.actions["reaction_id_id"].unique().tolist()
        reactionplate = self.addPlate(platetype="Reaction")

        for Index in range(len(self.actions.index.values)):
            columns = self.actions.columns.values
            currentaction = pd.DataFrame(index=[Index], columns=columns)
            for col in columns:
                currentaction[col] = self.actions.iloc[Index][col]
            currentaction["outputwell"] = reactionids.index(currentaction["reaction_id_id"].values)
            print()
            if currentaction["outputwell"].values[0] > len(reactionids):
                reactionplate = self.addPlate(platetype="Reaction")
            self.processAction(currentaction)

    def processAction(self, currentaction):
        currentactiontype = currentaction["actiontype"].to_string(index=False)
        currentactiontype = currentactiontype.strip()
        if currentactiontype != "Series([], )":
            # WTOSCR: what to do if it is equal to series([],)
            # WTOSCR: check django compatability

            numreps = currentaction["numberofrepetitions"].values
            if str(numreps) == "[nan]":
                numreps = 1
            elif str(numreps) == "[None]":
                numreps = 1
            elif str(numreps) == "[]":
                numreps = 1
            elif numreps == 0:
                numreps = 1

            try:
                numreps = int(numreps)
            except:
                numreps = 1

            repetitions = 0
            while repetitions < numreps:
                if currentactiontype == "add":
                    self.actionAdd(currentaction)
                elif currentactiontype == "collect-layer":
                    self.actionCollectLayer(currentaction)
                elif currentactiontype == "wash":
                    self.actionWash(currentaction)
                elif currentactiontype == "stir":
                    self.output.unsuportedAction(
                        "stir at a "
                        + str(currentaction["stirringspeed"].values[0])
                        + " speed at "
                        + str(currentaction["temperature"].values[0])
                        + " Celsius for "
                        + str(currentaction["duration"].values[0])
                        + " "
                        + str(currentaction["durationunit"].values[0])
                    )
                elif currentactiontype == "set-temperature":
                    if currentaction["duration"].values[0] in ["nan", "NaN", 0, "0"]:
                        self.output.unsuportedAction(
                            f"set temperature to {currentaction['temperature'].values[0]} Celsius for {currentaction['duration'].values[0]} {currentaction['durationunit'].values[0]}"
                        )
                    else:
                        self.output.unsuportedAction(
                            f"set temperature to {currentaction['temperature'].values[0]} Celsius"
                        )

                elif currentactiontype == "store":
                    self.output.unsuportedAction(
                        "store product (" + str(currentaction["material"].values[0]) + ")"
                    )
                elif currentactiontype == "extract":
                    self.actionExtract(currentaction)
                elif currentactiontype == "concentrate":
                    self.actionConcentrate(currentaction)
                else:
                    print(f"unsupported{currentaction['actiontype'].values}")
                    self.output.unsuportedAction(
                        f"{currentaction['actiontype'].values[0]} is not currently supported "
                    )
                repetitions += 1

    def actionAdd(self, currentaction):
        # handle multiple input plates!!!!! How?????
        tipvolume = self.choosetip(currentaction["materialquantity"].values[0])
        print(f"tipvolume {tipvolume}")
        pipetteName = (self.deck.findPippets(tipvolume)).name

        self.output.transferfluids(
            pipetteName,
            (
                f"OrderPlate.wells(){self.deck['OrderPlate'].smilesearch(currentaction['materialsmiles'].values[0], start_smiles = True, start_solvent=currentaction['solvent'].values[0])[0]}"
            ),
            (f"ReactionPlate.wells()[{currentaction['outputwell'].values[0]}]"),
            currentaction["materialquantity"].values[0],
        )

    def actionWash(self, currentaction):
        print("wash")
        tipvolume = self.choosetip(currentaction["materialquantity"].values[0])
        print(tipvolume)
        pipetteName = (self.deck.findPippets(tipvolume)).name
        print("wash debug")
        print(pipetteName)
        pipetteName = pipetteName[0]

        self.output.transferfluids(
            pipetteName,
            (
                f"OrderPlate.wells()[{self.deck['OrderPlate'].smilesearch(currentaction['material'].values[0], start_smiles = True, start_solvent=currentaction['solvent'].values[0])[0]}]"
            ),
            (f"ReactionPlate.wells()[{currentaction['outputwell'].values[0]}]"),
            currentaction["materialquantity"].values[0],
        )

    def actionCollectLayer(self, currentaction):
        print("collect")
        self.output.unsuportedAction("Collect-Layer not yet supported ")

    def actionExtract(self, currentaction):
        print("Extract")
        self.output.unsuportedAction(
            f"Extract using {currentaction['solventquantity'].values[0]} {currentaction['solventquantityunit'].values[0]} of {currentaction['solvent'].values[0]}"
        )

    def actionConcentrate(self, currentaction):
        print("Concentrate")
        self.output.unsuportedAction("Concentrate not yet supported ")


collected_actions = CollectActions(projectid=252)