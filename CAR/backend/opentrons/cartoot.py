"""Create OT session"""
from __future__ import annotations
from os import name
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile

from django.forms.models import model_to_dict

from statistics import mean, median, mode

import pandas as pd
from pandas.core.frame import DataFrame
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

import math
from itertools import groupby
from labwareavailable import labware_plates


class CreateOTSession(object):
    """
    Creates a StartOTSession object for generating a protocol
    from actions for reactionqueryset
    """

    def __init__(
        self,
        projectid: int,
        reactiongroupqueryset: list,
        inputplatequeryset: list = None,
    ):
        self.projectid = projectid
        self.reactiongroupqueryset = reactiongroupqueryset
        self.groupedtemperaturereactionobjs = self.getGroupedTemperatureReactions()
        self.inputplatequeryset = inputplatequeryset
        self.alladdactionqueryset = [
            self.getAddActions(reactionobj) for reactionobj in self.reactiongroupqueryset
        ]
        self.alladdactionquerysetflat = [
            item for sublist in self.alladdactionqueryset for item in sublist
        ]
        self.roundedvolumes = self.getRoundedVolumes(
            addactionqueryset=self.alladdactionquerysetflat
        )
        self.numbertips = self.getNumberTips()
        self.tipracktype = self.getTipRackType()
        self.pipettetype = self.getPipetteType()

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
        productobj = Product.objects.filter(reaction_id=reactionobj.id).order_by("id")[0]
        return productobj

    def getAddActions(self, reactionobj):
        addactionqueryset = IBMAddAction.objects.filter(reaction_id=reactionobj.id).order_by("id")
        return addactionqueryset

    def getStartingReactionPlateQuerySet(self):
        reactionplatequeryset = Plate.objects.filter(
            otsession_id=self.otsessionobj.id, platename__contains="Reactionplate"
        ).order_by("id")
        return reactionplatequeryset

    def getRoundedVolumes(self, addactionqueryset):
        roundedvolumes = [
            round(addactionobj.materialquantity) for addactionobj in addactionqueryset
        ]
        return roundedvolumes

    def getTipRackType(self):
        tipsavailable = {
            300: "opentrons_96_tiprack_300ul",
            10: "opentrons_96_tiprack_20ul",
        }
        tipkey = min(tipsavailable, key=lambda x: self.getNumberTransfers(pipettevolume=x))
        tipracktype = tipsavailable[tipkey]
        return tipracktype

    def getNumberTips(self):
        numbertips = len(self.alladdactionquerysetflat)
        return numbertips

    def getNumberTransfers(self, pipettevolume):
        numbertransfers = sum(
            [
                round(volume / pipettevolume) if pipettevolume < volume else 1
                for volume in self.roundedvolumes
            ]
        )
        return numbertransfers

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
        pipettekey = min(pipettesavailable, key=lambda x: self.getNumberTransfers(pipettevolume=x))
        pipettetype = pipettesavailable[pipettekey]
        return pipettetype

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

    def getMaxWellVolume(self, plateobj):
        maxwellvolume = plateobj.maxwellvolume
        return maxwellvolume

    def getDeadVolume(self, maxwellvolume):
        deadvolume = maxwellvolume * 0.05
        return deadvolume

    def getCloneWells(self, plateobj):
        clonewellqueryset = Well.objects.filter(plate_id=plateobj.id)
        return clonewellqueryset

    def getUniqueTemperatures(self):
        temperatures = sorted(
            set([reactionobj.reactiontemperature for reactionobj in self.reactiongroupqueryset])
        )
        return temperatures

    def getGroupedTemperatureReactions(self):
        temperatures = self.getUniqueTemperatures()
        groupedtemperaturereactionobjs = []

        for temperature in temperatures:
            temperaturereactiongroup = [
                reactionobj
                for reactionobj in self.reactiongroupqueryset
                if reactionobj.reactiontemperature == temperature
            ]
            groupedtemperaturereactionobjs.append(temperaturereactiongroup)

        return groupedtemperaturereactionobjs

    def getMedianValue(self, values):
        medianvalue = median(values)
        return medianvalue

    def getMaxValue(self, values):
        maxvalue = max(values)
        return maxvalue

    def getSumValue(self, values):
        sumvalue = sum(values)
        return sumvalue

    def getNumberObjs(self, queryset: list):
        numberobjs = len(queryset)
        return numberobjs

    def getReactionLabwarePlateType(self, grouptemperaturereactionobjs):
        numberreactions = self.getNumberObjs(queryset=grouptemperaturereactionobjs)
        reactionvolumes = []

        for reactionobj in grouptemperaturereactionobjs:
            addactionqueryset = self.getAddActions(reactionobj=reactionobj)
            roundedvolumes = self.getRoundedVolumes(addactionqueryset=addactionqueryset)
            sumvolume = self.getSumValue(values=roundedvolumes)
            reactionvolumes.append(sumvolume)
        maxvolume = self.getMaxValue(values=reactionvolumes)
        medianvolume = self.getMedianValue(values=reactionvolumes)
        headspacevolume = maxvolume + (maxvolume * 0.2)
        reactiontemperature = grouptemperaturereactionobjs[0].reactiontemperature

        if reactiontemperature > 25:
            labwareplatetypes = [
                labware_plate
                for labware_plate in labware_plates
                if labware_plates[labware_plate]["reflux"] == True
                and labware_plates[labware_plate]["volume_well"] > headspacevolume
                and labware_plates[labware_plate]["no_wells"] >= numberreactions
            ]
        else:
            labwareplatetypes = [
                labware_plate
                for labware_plate in labware_plates
                if not labware_plates[labware_plate]["reflux"]
                and labware_plates[labware_plate]["volume_well"] > headspacevolume
                and labware_plates[labware_plate]["no_wells"] >= numberreactions
            ]

        if len(labwareplatetypes) > 1:
            volumewells = [
                labware_plates[labware_plate]["volume_well"] for labware_plate in labwareplatetypes
            ]
            indexclosestvalue = min(range(len(volumewells)), key=lambda x: abs(x - medianvolume))
            labwareplatetype = labwareplatetypes[indexclosestvalue]
        else:
            labwareplatetype = labwareplatetypes[0]

        return labwareplatetype

    def getNumberVials(self, maxvolumevial, volumematerial):
        if maxvolumevial > volumematerial:
            novialsneeded = 1
        else:
            volumestoadd = []
            deadvolume = self.getDeadVolume(maxwellvolume=maxvolumevial)
            novialsneededratio = volumematerial / (maxvolumevial - deadvolume)
            frac, whole = math.modf(novialsneededratio)
            volumestoadd = [maxvolumevial for i in range(int(whole))]
            volumestoadd.append(frac * maxvolumevial + deadvolume)
            novialsneeded = sum(volumestoadd)
        return novialsneeded

    def getStarterPlateType(self, startingmaterialsdf: DataFrame):
        labwareplatetypes = [
            labware_plate
            for labware_plate in labware_plates
            if labware_plates[labware_plate]["starting_plate"] == True
        ]

        vialcomparedict = {}

        for labwareplate in labwareplatetypes:
            maxvolumevial = labware_plates[labwareplate]["volume_well"]
            noplatevials = labware_plates[labwareplate]["no_wells"]

            vialsneeded = startingmaterialsdf.apply(
                lambda row: self.getNumberVials(
                    maxvolumevial=maxvolumevial, volumematerial=row["materialquantity"]
                ),
                axis=1,
            )

            totalvialsneeded = sum(vialsneeded)
            platesneeded = int(math.ceil(totalvialsneeded / noplatevials))

            vialcomparedict[maxvolumevial] = platesneeded

        minimumnovialsvolume = min(vialcomparedict, key=vialcomparedict.get)

        labwareplatetype = [
            labware_plate
            for labware_plate in labwareplatetypes
            if labware_plates[labware_plate]["volume_well"] == minimumnovialsvolume
        ][0]

        return labwareplatetype

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

    def createPlateModel(self, platename, labwaretype):
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            plateindex = indexslot
            maxwellvolume = labware_plates[labwaretype]["volume_well"]
            numberwells = labware_plates[labwaretype]["no_wells"]
            plateobj = Plate()
            plateobj.otsession_id = self.otsessionobj
            plateobj.deck_id = self.deckobj
            plateobj.platename = "{}_{}".format(platename, indexslot)
            plateobj.plateindex = plateindex
            plateobj.labware = labwaretype
            plateobj.maxwellvolume = maxwellvolume
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

    def createTipRacks(self):
        numberacks = int(-(-self.numbertips // 96))
        self.tipracktype = self.getTipRackType()
        for rack in range(numberacks):
            self.createTiprackModel(name=self.tipracktype)

    def checkDeckSlotAvailable(self):
        testslotavailable = self.deckobj.indexslotavailable
        if testslotavailable <= self.deckobj.numberslots:
            self.deckobj.indexslotavailable = testslotavailable + 1
            self.deckobj.save()
            return testslotavailable
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

    def createStartingPlate(self):
        # Can we do this in Django queries and aggregation instead?
        addactionsdf = self.getOrderAddActions()

        startingmaterialsdf = addactionsdf.groupby(["uniquesolution"]).agg(
            {
                "reaction_id_id": "first",
                "material": "first",
                "materialsmiles": "first",
                "materialquantity": "sum",
                "solvent": "first",
                "concentration": "first",
                "mculeid": "first",
            }
        )

        startingmaterialsdf["productexists"] = startingmaterialsdf.apply(
            lambda row: self.checkProductExists(
                reaction_id=row["reaction_id_id"], smiles=row["materialsmiles"]
            ),
            axis=1,
        )

        startingmaterialsdf = startingmaterialsdf[~startingmaterialsdf["productexists"]]

        startingmaterialsdf = startingmaterialsdf.sort_values(
            ["solvent", "materialquantity"], ascending=False
        )

        startinglabwareplatetype = self.getStarterPlateType(startingmaterialsdf=startingmaterialsdf)

        plateobj = self.createPlateModel(
            platename="Startingplate", labwaretype=startinglabwareplatetype
        )

        maxwellvolume = self.getMaxWellVolume(plateobj=plateobj)
        deadvolume = self.getDeadVolume(maxwellvolume=maxwellvolume)

        orderdictslist = []

        for i in startingmaterialsdf.index.values:
            totalvolume = startingmaterialsdf.at[i, "materialquantity"]
            if totalvolume > maxwellvolume:
                nowellsneededratio = totalvolume / (maxwellvolume - deadvolume)

                frac, whole = math.modf(nowellsneededratio)
                volumestoadd = [maxwellvolume for i in range(int(whole))]
                volumestoadd.append(frac * maxwellvolume + deadvolume)

                for volumetoadd in volumestoadd:
                    indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
                    if not indexwellavailable:

                        plateobj = self.createPlateModel(
                            platename="Startingplate", labwaretype=startinglabwareplatetype
                        )

                        indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)

                    wellobj = self.createWellModel(
                        plateobj=plateobj,
                        reactionobj=self.getReaction(
                            reactionid=startingmaterialsdf.at[i, "reaction_id_id"]
                        ),
                        wellindex=indexwellavailable - 1,
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
                volumetoadd = totalvolume + deadvolume

                if not indexwellavailable:
                    plateobj = self.createPlateModel(
                        platename="Startingplate", labwaretype=startinglabwareplatetype
                    )
                    indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)

                wellobj = self.createWellModel(
                    plateobj=plateobj,
                    reactionobj=self.getReaction(
                        reactionid=startingmaterialsdf.at[i, "reaction_id_id"]
                    ),
                    wellindex=indexwellavailable - 1,
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

        orderdf = pd.DataFrame(orderdictslist)

        self.createCompoundOrderModel(orderdf=orderdf)

    def createReactionPlate(self):
        for grouptemperaturereactionobjs in self.groupedtemperaturereactionobjs:
            labwareplatetype = self.getReactionLabwarePlateType(
                grouptemperaturereactionobjs=grouptemperaturereactionobjs
            )

            plateobj = self.createPlateModel(
                platename="Reactionplate", labwaretype=labwareplatetype
            )

            for reactionobj in grouptemperaturereactionobjs:
                productobj = self.getProduct(reactionobj)
                indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)
                if not indexwellavailable:
                    plateobj = self.createPlateModel(
                        platename="Reactionplate", labwaretype=labwareplatetype
                    )

                    indexwellavailable = self.checkPlateWellsAvailable(plateobj=plateobj)

                self.createWellModel(
                    plateobj=plateobj,
                    reactionobj=reactionobj,
                    wellindex=indexwellavailable - 1,
                    volume=None,
                    smiles=productobj.smiles,
                    concentration=None,
                    solvent=None,
                    mculeid=None,
                )

    def checkProductExists(self, reaction_id, smiles):
        previousreactionid = reaction_id - 1
        testproduct = Product.objects.filter(reaction_id=previousreactionid, smiles=smiles)
        if testproduct:
            return True
        else:
            return False

    def combinestrings(self, row):
        return (
            str(row["materialsmiles"]) + "-" + str(row["solvent"]) + "-" + str(row["concentration"])
        )

    def cloneInputPlate(self):
        for plateobj in self.inputplatequeryset:
            indexslot = self.checkDeckSlotAvailable()
            if indexslot:
                clonewellqueryset = self.getCloneWells(plateobj=plateobj)
                plateindex = indexslot
                platename = "Startingplate"
                plateobj.pk = None
                plateobj.deck_id = self.deckobj
                plateobj.otsession_id = self.otsessionobj
                plateobj.plateindex = plateindex
                plateobj.platename = "{}_{}".format(platename, indexslot)
                plateobj.save()
                self.cloneInputWells(clonewellqueryset, plateobj)
            else:
                print("No more deck slots available")

    def cloneInputWells(self, clonewellqueryset, plateobj):
        for clonewellobj in clonewellqueryset:
            clonewellobj.pk = None
            clonewellobj.plate_id = plateobj
            clonewellobj.otsession_id = self.otsessionobj
            clonewellobj.save()


def getTargets(projectid):
    targetqueryset = Target.objects.filter(project_id=projectid).order_by("id")
    return targetqueryset


def getMethods(targetid):
    methodqueryset = Method.objects.filter(target_id=targetid).order_by("id")
    return methodqueryset


def getReactions(methodid):
    reactionqueryset = Reaction.objects.filter(method_id=methodid).order_by("id")
    return reactionqueryset


def getProjectReactions(projectid):
    targetqueryset = getTargets(projectid=projectid)
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


projectid = 20

allreactionquerysets = getProjectReactions(projectid=projectid)
maxsteps = findmaxlist(allreactionquerysets=allreactionquerysets)
groupedreactionquerysets = groupReactions(
    allreactionquerysets=allreactionquerysets, maxsteps=maxsteps
)

platequeryset = []

for index, reactiongroup in enumerate(groupedreactionquerysets):
    if index == 0:
        otsession = CreateOTSession(
            projectid=projectid,
            reactiongroupqueryset=reactiongroup,
        )

        otsessionobj = otsession.otsessionobj
        alladdactionsquerysetflat = otsession.alladdactionquerysetflat
        startingreactionplatequeryset = otsession.startingreactionplatequeryset

        otWrite(
            otsessionobj=otsessionobj,
            alladdactionsquerysetflat=alladdactionsquerysetflat,
        )
    if index > 0:
        otsession = CreateOTSession(
            projectid=projectid,
            reactiongroupqueryset=reactiongroup,
            inputplatequeryset=startingreactionplatequeryset,
        )

        otsessionobj = otsession.otsessionobj
        alladdactionsquerysetflat = otsession.alladdactionquerysetflat
        reactionplatequeryset = otsession.startingreactionplatequeryset

        otWrite(
            otsessionobj=otsessionobj,
            alladdactionsquerysetflat=alladdactionsquerysetflat,
        )
