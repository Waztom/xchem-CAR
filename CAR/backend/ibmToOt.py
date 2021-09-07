from pathlib import Path
from django.db.models.expressions import Col
import pandas as pd
import math
import re

import opentrons.otWrite as otWrite
import opentrons.otDeck as otDeck
import mcule.outputplatetxt as OutputPlateTxt

import backend.models

# can be used in ibmtoot by importing
import mcule.outputplatecsv as OutputPlateCSV

# and calling:
# OutputPlateTxt.PlateCSV(self.orderplate, self.name, self.author, currentblocknum)


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

    def getActions(self):
        # THIS NEEDS TO BE FIXED!!!!!!
        allactions_list_df = []
        targets = backend.models.Target.objects.filter(project_id=self.projectid)
        target_ids = targets.values_list("pk", flat=True)
        methods = [
            backend.models.Method.objects.filter(target_id=target_id) for target_id in target_ids
        ]
        method_ids_queryset = [method.values_list("pk", flat=True) for method in methods]
        method_ids = []
        for method_id_queryset in method_ids_queryset:
            ids_to_add = [query for query in method_id_queryset]
            method_ids.append(ids_to_add)
        method_ids = [item for sublist in method_ids for item in sublist]
        reactions = [
            backend.models.Reaction.objects.filter(method_id=method_id) for method_id in method_ids
        ]

        for actionmodel in self.actionmodels:
            for reaction in reactions:
                actions_to_add_df = pd.DataFrame(
                    list(actionmodel.objects.filter(reaction_id=reaction[0].id).values())
                )
                if not actions_to_add_df.empty:
                    allactions_list_df.append(actions_to_add_df)
        allactions_list_df = pd.concat(allactions_list_df)

        self.allactions_df = allactions_list_df.sort_values(["reaction_id_id", "actionno"])

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
            self.actionsfiltered = subSetReactAct.loc[(subSetReactAct["actiontype"]).isin(actions)]
        else:
            self.actionsfiltered = subSetReactAct

        self.actionsfiltered["doable"] = self.actionsfiltered.apply(
            lambda row: self.docheck(row), axis=1
        )

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

    def startProtocol(self):
        for blocknum in self.actionsfiltered["blocknum"].unique():
            actionsblock = self.actionsfiltered[self.actionsfiltered["blocknum"] == blocknum]
            if actionsblock["blockbool"].values[0] == True:
                otSession(
                    name=f"block_{blocknum}",
                    actions=actionsblock,
                    author="Example Author",
                    description="example description",
                )


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
        self.currentpipettesetup = currentpipettesetup

        # decalre defined pipettes if they exist
        # WTOSCR: could be improved by having currentpipettesetup defult to none
        if self.currentpipettesetup != None:
            self.definedpipettes = currentpipettesetup
        else:
            self.definedpipettes = None

        # Instantiate Deck to model physical constraints and locations involved in protocol
        self.deck = otDeck.Deck()

        # setup name for robot py file
        self.outputpath = "None"  # WTOSCR: create and save script in DJango model
        self.pathnamecheck()  # generate unique name

        # Create and setup main robotic .py script to write output to
        self.output = otWrite.otScript(
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
        self.setupPlate()

        # Generate idea of what tips will be used in protocol
        self.preselecttips()

        #
        self.tipOutput()

        # write hardwear setup to robotics .py script
        self.output.setupLabwear(self.deck.PlateList, self.tipRackList)

        self.setupPipettes()
        self.output.setupPipettes(self.deck.PipetteList)
        self.ittrActions()

        # WTOSCR: add qc.text file at end of session once output plates are filled

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

    def setupPlate(
        self,
        inputnumwells=24,
        inputwellvolume=2500,
        incAddMaterials=True,
        incWashMaterials=False,
    ):
        # WTOSCR: currently relying on user input of number of wells and input well volume, a function should be created along the line of the pipette selection process to automaticaly choose which plate should be used
        # WTOSCR: do we need autogenrated plate size selection - proberbly more for output plate but the same funciton should work for both
        # print(self.actions)

        allmaterials = pd.DataFrame(
            columns=[
                "mculeid",
                "concentration",
                "material",
                "materialsmiles",
                "materialquantity",
                "solvent",
                "materialKey",
            ]
        )  # WTOSCR: should proberbly check about a Django model

        # each block getting a list of materials associated with that action type and concatinating them to the pandas Dataframe, All materials
        if incAddMaterials == True:
            addingSteps = self.actions.loc[(self.actions["actiontype"]) == "add"]
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
            allmaterials = pd.concat([allmaterials, materials], ignore_index=True)

        if (
            incWashMaterials == True
        ):  # WTOSCR: review if properly implemented and bring inline with add
            solventSteps = self.actions.loc[(self.actions["actiontype"]) == "wash"]
            washmaterials = solventSteps[
                "material", "materialquantitcombinestringsquantity", "materialquantity"
            ]
            allmaterials = pd.concat([allmaterials, extractsolvents], ignore_index=True)

        allmaterials["matAndSolv"] = allmaterials.apply(
            lambda row: self.combinestrings(row), axis=1
        )
        allmaterials = allmaterials.groupby(["matAndSolv"]).aggregate(
            {
                "materialquantity": "sum",
                "material": "first",
                "solvent": "first",
                "materialsmiles": "first",
                "mculeid": "first",
                "concentration": "first",
            }
        )
        allmaterials = allmaterials.sort_values(["solvent", "materialquantity"], ascending=False)
        # experement with adding a working materials+sovlent collum to group by then drop that col
        self.maxvolume = allmaterials["materialquantity"].max()  # not currently used
        self.totalvolume = allmaterials["materialquantity"].sum()  # not currently used

        # Plate objects are stored in Deck object as a list in PlateList
        # Create and add order plate to deck
        self.deck.add(
            Type="Plate",
            numwells=inputnumwells,
            platewellVolume=inputwellvolume,
            platename="OrderPlate",
        )  # WTOSCR: index is hard coded

        # Get the order plate form the PlateList
        self.orderplate = [
            plate for plate in self.deck.PlateList if plate.plateName == "OrderPlate"
        ][0]

        for i in allmaterials.index.values:
            if (
                allmaterials.at[i, "material"] == ""
                or allmaterials.at[i, "material"] == None
                or allmaterials.at[i, "material"] == "NaN"
            ):
                matName = allmaterials.at[
                    i, "materialsmiles"
                ]  # implement chemical name check for all reactants and reagents
            else:
                matName = allmaterials.at[i, "material"]

            # Check if free well
            self.orderplate.nextfreewell()

            if not self.orderplate.isplatefull:
                wellnumber = self.orderplate.nextfreewellindex

                self.orderplate.WellList[wellnumber].add(
                    amount=allmaterials.at[i, "materialquantity"],
                    smiles=allmaterials.at[i, "materialsmiles"],
                    solvent=allmaterials.at[i, "solvent"],
                    materialname=matName,
                    mculeid=allmaterials.at[i, "mculeid"],
                    concentration=allmaterials.at[i, "concentration"],
                )

        currentblocknum = self.actions["blocknum"].values[0]
        OutputPlateCSV.PlateCSV(self.orderplate, self.name, self.author, currentblocknum)
        # OutputPlateTxt.PlateTxt(self.orderplate, self.name, self.author, currentblocknum)
        # WTOSCR: add qc.text file at end of session once output plates are filled
        # if well overvlows (outcome == False) needs to move to next weel
        # if outcome != False:
        #     plate[wellnumber] = [allmaterials.loc[i,'materialquantity'], allmaterials.loc[i,'material']]
        # WTOSCR: add handleing materials spread across many wells in fluid handeling steps

        # Add reaction plate
        # Harcoded plate values!!!! Must fix
        #######################################
        # Change reaction plate type
        #######################################
        self.deck.add(
            Type="Plate",
            numwells=96,
            platewellVolume=2500,
            platename="ReactionPlate",
        )

        self.reactionPlate = [
            plate for plate in self.deck.PlateList if plate.plateName == "ReactionPlate"
        ][0]
        # WTOSCR: index is hard coded

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
        reactionids = (
            self.actions["reaction_id_id"].unique().tolist()
        )  # WTOSCR: Django implementaiton
        for Index in range(len(self.actions.index.values)):
            columns = self.actions.columns.values
            currentaction = pd.DataFrame(index=[Index], columns=columns)
            for col in columns:
                currentaction[col] = self.actions.iloc[Index][col]
            currentaction["outputwell"] = reactionids.index(currentaction["reaction_id_id"].values)

            self.processAction(currentaction)

    def splitittract(self, split):
        if type(split) == "int":
            subset = self.actions[split]
            for actionindex in range(len(subset) + 1):
                currentactionmask = subset["actionno"] == actionindex + 1
                currentaction = subset[currentactionmask]
                self.processAction(currentaction)

        elif type(split) == "<class 'list'>":
            print(type(split[0]))
            if type(split[0]) == "int":
                subset = self.actions[split[0] : split[1]]
                print("intlistsplit")
                for actionindex in range(len(subset) + 1):
                    currentactionmask = subset["actionno"] == actionindex + 1
                    currentaction = subset[currentactionmask]
                    self.processAction(currentaction)

            elif type(split[0]) == "str":
                print("stirlistsplit")
                subset = self.actions["actiontype" in split]
                for actionindex in range(len(subset) + 1):
                    currentactionmask = subset["actionno"] == actionindex + 1
                    currentaction = subset[currentactionmask]
                    self.processAction(currentaction)
            else:
                print(type(split[0]))

        elif type(split) == "str":
            subset = self.actions["actiontype" == split]
            for actionindex in range(len(subset) + 1):
                currentactionmask = subset["actionno"] == actionindex + 1
                currentaction = subset[currentactionmask]
                self.processAction(currentaction)
        else:
            print("debug1")

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


collected_actions = CollectActions(projectid=138)
collected_actions.getActions()
collected_actions.actionfilter()
collected_actions.blockdefine()
collected_actions.startProtocol()
