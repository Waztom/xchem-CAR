"""Create OT session"""
from __future__ import annotations
from os import name
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from itertools import groupby
from statistics import mode

import pandas as pd
from otwrite import otWrite

from backend.models import (
    CompoundOrder,
    IBMAddAction,
    Pipette,
    TipRack,
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


class CreateOTSession(object):
    """
    Creates a StartOTSession object for generating a protocol
    from actions for reactionqueryset
    """

    def __init__(
        self,
        starterplatetype: str,
        reactionplatetype: str,
        projectid: int,
        reactiongroupqueryset: list,
        inputplatequeryset: list = None,
    ):
        self.starterplatetype = starterplatetype
        self.reactionplatetype = reactionplatetype
        self.projectid = projectid
        self.reactiongroupqueryset = reactiongroupqueryset
        self.alladdactionqueryset = [
            self.getAddActions(reactionobj) for reactionobj in self.reactiongroupqueryset
        ]
        self.alladdactionquerysetflat = [
            item for sublist in self.alladdactionqueryset for item in sublist
        ]

        self.modevolume = self.getModeVolumeAdded()
        self.numbertips = self.getNumberTips()
        self.tipracktype = self.getTipRackType()
        self.pipettetype = self.getPipetteType()

        self.productqueryset = [
            self.getProduct(reactionobj) for reactionobj in self.reactiongroupqueryset
        ]
        self.otsessionobj = self.createOTSessionModel()
        self.deckobj = self.createDeckModel()
        self.inputplatequeryset = inputplatequeryset
        if self.inputplatequeryset:
            self.cloneInputPlate()
        self.createPipetteModel()
        self.createStartingPlate()
        self.createReactionPlate()
        self.createTipRacks()
        self.startingreactionplatequeryset = self.getStartingReactionPlateQuerySet()

    def getReaction(self, reactionid):
        reactionobj = Reaction.objects.filter(id=reactionid)[0]
        return reactionobj

    def getProduct(self, reactionobj):
        productqueryset = Product.objects.filter(reaction_id=reactionobj.id).order_by("id")
        return productqueryset

    def getAddActions(self, reactionobj):
        addactionqueryset = IBMAddAction.objects.filter(reaction_id=reactionobj.id).order_by("id")
        return addactionqueryset

    def getStartingReactionPlateQuerySet(self):
        reactionplatequeryset = Plate.objects.filter(platename="Reactionplate").order_by("id")
        return reactionplatequeryset

    def getModeVolumeAdded(self):
        roundedvolumes = [
            round(addactionobj.materialquantity) for addactionobj in self.alladdactionquerysetflat
        ]
        modevolume = mode(roundedvolumes)
        return modevolume

    def getTipRackType(self):
        tipsavailable = {
            300: "opentrons_96_tiprack_300ul",
            10: "opentrons_96_tiprack_20ul",
        }
        tipkey = min(tipsavailable, key=lambda x: abs(x - self.modevolume))
        tipracktype = tipsavailable[tipkey]
        return tipracktype

    def createOTSessionModel(self):
        otsessionobj = OTSession()
        project_obj = Project.objects.get(id=self.projectid)
        otsessionobj.project_id = project_obj
        otsessionobj.save()
        return otsessionobj

    def createDeckModel(self):
        deckobj = Deck()
        deckobj.otsession_id = self.otsessionobj
        deckobj.save()
        self.deckobj = deckobj
        return deckobj

    def createPipetteModel(self):
        pipetteobj = Pipette()
        pipetteobj.otsession_id = self.otsessionobj
        pipetteobj.position = self.pipettetype["position"]
        pipetteobj.maxvolume = self.pipettetype["maxvolume"]
        pipetteobj.type = self.pipettetype["type"]
        pipetteobj.pipettename = "{}_{}".format(
            self.pipettetype["position"], self.pipettetype["labware"]
        )
        pipetteobj.labware = self.pipettetype["labware"]
        pipetteobj.save()

    def createTiprackModel(self, name):
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            tiprackindex = indexslot
            tiprackobj = TipRack()
            tiprackobj.otsession_id = self.otsessionobj
            tiprackobj.deck_id = self.deckobj
            tiprackobj.tiprackname = "{}_{}".format(name, indexslot)
            tiprackobj.tiprackindex = tiprackindex
            tiprackobj.labware = name
            tiprackobj.save()
        else:
            print("No more deck slots available")

    def createPlateModel(self, platename, labware, numberwells):
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            plateindex = indexslot
            plateobj = Plate()
            plateobj.otsession_id = self.otsessionobj
            plateobj.deck_id = self.deckobj
            plateobj.platename = "{}_{}".format(platename, indexslot)
            plateobj.plateindex = plateindex
            plateobj.labware = labware
            plateobj.numberwells = numberwells
            plateobj.save()
            return plateobj
        else:
            print("No more deck slots available")

    def createWellModel(
        self,
        plateobj,
        reactionobj,
        wellindex,
        volume,
        smiles,
        concentration,
        solvent=None,
        mculeid=None,
    ):
        wellobj = Well()
        wellobj.otsession_id = self.otsessionobj
        wellobj.plate_id = plateobj
        wellobj.reaction_id = reactionobj
        wellobj.wellindex = wellindex
        wellobj.volume = volume
        wellobj.smiles = smiles
        wellobj.concentration = concentration
        if solvent:
            wellobj.solvent = solvent
        if mculeid:
            wellobj.mculeid = mculeid
        wellobj.save()
        return wellobj

    def createCompoundOrderModel(self, orderdf):
        compoundorderobj = CompoundOrder()
        project_obj = Project.objects.get(id=self.projectid)
        compoundorderobj.project_id = project_obj
        csvdata = orderdf.to_csv(encoding="utf-8", index=False)
        ordercsv = default_storage.save(
            "mculeorders/"
            + "starterplate-for-project-{}-sessionid-{}".format(
                project_obj.name, str(self.otsessionobj.id)
            )
            + ".csv",
            ContentFile(csvdata),
        )
        compoundorderobj.ordercsv = ordercsv
        compoundorderobj.save()

    def getNumberTips(self):
        numbertips = len(self.alladdactionquerysetflat)
        return numbertips

    def createTipRacks(self):
        numbertips = self.getNumberTips()
        numberacks = int(-(-numbertips // 96))
        self.tipracktype = self.getTipRackType()
        for rack in range(numberacks):
            self.createTiprackModel(name=self.tipracktype)

    def getPipetteType(self):
        pipettesavailable = {
            10: {"labware": "p10_single", "position": "right", "type": "single", "maxvolume": 10},
            300: {
                "labware": "p300_single",
                "position": "right",
                "type": "single",
                "maxvolume": 300,
            },
        }
        pipettekey = min(pipettesavailable, key=lambda x: abs(x - self.modevolume))
        pipettetype = pipettesavailable[pipettekey]
        return pipettetype

    def checkDeckSlotAvailable(self):
        testslotavailable = self.deckobj.indexslotavailable + 1
        if testslotavailable <= self.deckobj.numberslots:
            self.deckobj.indexslotavailable = testslotavailable
            self.deckobj.save()
            return self.deckobj.indexslotavailable
        else:
            self.deckobj.slotavailable = False
            self.deckobj.save()
            return False

    def checkPlateWellsAvailable(self, plateobj):
        wellavailable = plateobj.indexswellavailable + 1
        numberwells = plateobj.numberwells
        if wellavailable <= numberwells:
            plateobj.indexswellavailable = wellavailable
            plateobj.save()
            return plateobj.indexswellavailable
        else:
            plateobj.wellavailable = False
            plateobj.save()
            return False

    def checkProductExists(self, smiles):
        testproduct = Product.objects.filter(smiles=smiles)
        if testproduct:
            return True
        else:
            return False

    def combinestrings(self, row):
        return str(str(row["materialsmiles"]) + str(row["solvent"]) + str(row["concentration"]))

    def getOrderAddActions(self):
        addactionslistdf = []

        for addactionqueryset in self.alladdactionqueryset:
            addactionsdf = pd.DataFrame(list(addactionqueryset.values()))
            if not addactionsdf.empty:
                addactionslistdf.append(addactionsdf)
        addactionsdf = pd.concat(addactionslistdf)

        addactionsdf["uniquesolution"] = addactionsdf.apply(
            lambda row: self.combinestrings(row), axis=1
        )

        return addactionsdf

    def createStartingPlate(self):
        addactionsdf = self.getOrderAddActions()

        startingmaterialsdf = addactionsdf.groupby(["uniquesolution"]).agg(
            {
                "materialquantity": "sum",
                "reaction_id_id": "first",
                "material": "first",
                "solvent": "first",
                "materialsmiles": "first",
                "mculeid": "first",
                "concentration": "first",
            }
        )

        startingmaterialsdf["productexists"] = startingmaterialsdf.apply(
            lambda row: self.checkProductExists(row["materialsmiles"]), axis=1
        )

        startingmaterialsdf = startingmaterialsdf[~startingmaterialsdf["productexists"]]

        startingmaterialsdf = startingmaterialsdf.sort_values(
            ["solvent", "materialquantity"], ascending=False
        )

        plateobj = self.createPlateModel(
            platename="Startingplate", labware=self.starterplatetype, numberwells=24
        )

        maxwellvolume = float(plateobj.labware.split("_")[-1].strip("ul"))

        orderdictslist = []

        for i in startingmaterialsdf.index.values:
            totalvolume = startingmaterialsdf.at[i, "materialquantity"]
            if totalvolume > maxwellvolume:
                nowellsneeded = int(-(-totalvolume // maxwellvolume))
                volumetoadd = totalvolume / nowellsneeded
                for index in range(nowellsneeded):
                    indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
                    if not indexwellavailable:
                        plateobj = self.createPlateModel(
                            platename="Startingplate", labware=self.starterplatetype, numberwells=24
                        )
                        indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
                    wellobj = self.createWellModel(
                        plateobj=plateobj,
                        reactionobj=self.getReaction(
                            reactionid=startingmaterialsdf.at[i, "reaction_id_id"]
                        ),
                        wellindex=indexwellavailable,
                        volume=volumetoadd,
                        smiles=startingmaterialsdf.at[i, "materialsmiles"],
                        concentration=startingmaterialsdf.at[i, "concentration"],
                        solvent=startingmaterialsdf.at[i, "solvent"],
                        mculeid=startingmaterialsdf.at[i, "mculeid"],
                    )

                    orderdictslist.append(
                        {
                            "material": startingmaterialsdf.at[i, "material"],
                            "mculeid": startingmaterialsdf.at[i, "mculeid"],
                            "platename": plateobj.platename,
                            "well": wellobj.wellindex,
                            "concentration": startingmaterialsdf.at[i, "concentration"],
                            "solvent": startingmaterialsdf.at[i, "solvent"],
                            "amount-ul": volumetoadd,
                        }
                    )

            else:
                indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
                if not indexwellavailable:
                    plateobj = self.createPlateModel(
                        platename="Startingplate", labware=self.starterplatetype, numberwells=24
                    )
                    indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)

                wellobj = self.createWellModel(
                    plateobj=plateobj,
                    reactionobj=self.getReaction(
                        reactionid=startingmaterialsdf.at[i, "reaction_id_id"]
                    ),
                    wellindex=indexwellavailable,
                    volume=totalvolume,
                    smiles=startingmaterialsdf.at[i, "materialsmiles"],
                    concentration=startingmaterialsdf.at[i, "concentration"],
                    solvent=startingmaterialsdf.at[i, "solvent"],
                    mculeid=startingmaterialsdf.at[i, "mculeid"],
                )

                orderdictslist.append(
                    {
                        "material": startingmaterialsdf.at[i, "material"],
                        "mculeid": startingmaterialsdf.at[i, "mculeid"],
                        "platename": plateobj.platename,
                        "well": wellobj.wellindex,
                        "concentration": startingmaterialsdf.at[i, "concentration"],
                        "solvent": startingmaterialsdf.at[i, "solvent"],
                        "amount-ul": startingmaterialsdf.at[i, "materialquantity"],
                    }
                )

        orderdf = pd.DataFrame(orderdictslist)

        self.createCompoundOrderModel(orderdf=orderdf)

    def createReactionPlate(self):
        plateobj = self.createPlateModel(
            platename="Reactionplate", labware=self.reactionplatetype, numberwells=96
        )

        for reactionobj in self.reactiongroupqueryset:
            productqueryset = self.getProduct(reactionobj)
            productobj = productqueryset[0]
            indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
            if not indexwellavailable:
                plateobj = self.createPlateModel(
                    platename="Reactionplate", labware=self.reactionplatetype, numberwells=24
                )

                indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)

            self.createWellModel(
                plateobj=plateobj,
                reactionobj=reactionobj,
                wellindex=indexwellavailable,
                volume=None,
                smiles=productobj.smiles,
                concentration=None,
                solvent=None,
                mculeid=None,
            )

    def cloneInputPlate(self):
        for plateobj in self.inputplatequeryset:
            plateobj.deck_id = self.deckobj
            plateobj.save()
            self.cloneInputWells(plateobj)

    def getCloneWells(self, plateobj):
        clonewellqueryset = Well.objects.filter(plate_id=plateobj.id)
        return clonewellqueryset

    def cloneInputWells(self, plateobj):
        clonewellqueryset = self.getCloneWells(plateobj=plateobj)
        for clonewellobj in clonewellqueryset:
            clonewellobj.plate_id = plateobj
            clonewellobj.save()


def getMethods(targetid):
    methodqueryset = Method.objects.filter(target_id=targetid).order_by("id")
    return methodqueryset


def getReactions(methodid):
    reactionqueryset = Reaction.objects.filter(method_id=methodid).order_by("id")
    return reactionqueryset


def getProjectReactions(projectid):
    targetqueryset = Target.objects.filter(project_id=projectid).order_by("id")
    allreactionquerysets = []
    for target in targetqueryset:
        methodqueryset = getMethods(targetid=target.id)
        for method in methodqueryset:
            reactionqueryset = getReactions(methodid=method.id)
            allreactionquerysets.append(reactionqueryset)
    return allreactionquerysets


def findnoallreactionsteps(allreactionquerysets: list):
    """ "
    Function to get all possible number of reactions steps for
    multiple methods
    """
    allnumberofsteps = set([len(x) for x in allreactionquerysets])
    return allnumberofsteps


def findmaxlist(allreactionquerysets: list):
    maxlength = max([len(i) for i in allreactionquerysets])
    return maxlength


def groupReactions(allreactionquerysets: list, maxsteps: int):
    """
    Groups reactionqueries into first reactions, second reactions and so on
    """
    groupedreactionquerysets = []
    for i in range(maxsteps):
        reactiongroup = [
            reactionqueryset[i]
            for reactionqueryset in allreactionquerysets
            if i <= len(reactionqueryset) - 1
        ]
        groupedreactionquerysets.append(reactiongroup)
    return groupedreactionquerysets


projectid = 3

allreactionquerysets = getProjectReactions(projectid=projectid)
maxsteps = findmaxlist(allreactionquerysets=allreactionquerysets)
groupedreactionquerysets = groupReactions(
    allreactionquerysets=allreactionquerysets, maxsteps=maxsteps
)

platequeryset = []

for index, reactiongroup in enumerate(groupedreactionquerysets):
    if index == 0:
        otsession = CreateOTSession(
            starterplatetype="24_resovoir_2500ul",
            reactionplatetype="96_labcyte_2500ul",
            projectid=projectid,
            reactiongroupqueryset=reactiongroup,
        )

        otsessionobj = otsession.otsessionobj
        alladdactionsquerysetflat = otsession.alladdactionquerysetflat
        reactionplatequeryset = otsession.startingreactionplatequeryset

        otWrite(
            otsessionobj=otsessionobj,
            alladdactionsquerysetflat=alladdactionsquerysetflat,
        )
    if index > 0:
        otsession = CreateOTSession(
            starterplatetype="24_resovoir_2500ul",
            reactionplatetype="96_labcyte_2500ul",
            projectid=projectid,
            reactiongroupqueryset=reactiongroup,
            inputplatequeryset=reactionplatequeryset,
        )

        otsessionobj = otsession.otsessionobj
        alladdactionsquerysetflat = otsession.alladdactionquerysetflat
        reactionplatequeryset = otsession.startingreactionplatequeryset

        otWrite(
            otsessionobj=otsessionobj,
            alladdactionsquerysetflat=alladdactionsquerysetflat,
        )
