"""Create OT session"""

from __future__ import annotations
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.db.models import QuerySet, Q
import logging

logger = logging.getLogger(__name__)

from statistics import median

from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
from pandas.core.frame import DataFrame

from ..utils import (
    canonSmiles,
    getProductSmiles,
    checkPreviousReactionProducts,
    getReactionQuerySet,
    getProduct,
    getReaction,
    getChemicalName,
    getInchiKey,
    wellIndexToWellName,
    stripSalts,
)

from ..models import (
    ActionSession,
    Batch,
    Column,
    ExtractAction,
    OTBatchProtocol,
    Reaction,
    Product,
    AddAction,
    OTSession,
    Deck,
    Plate,
    SolventPrep,
    Well,
    Pipette,
    TipRack,
    CompoundOrder,
)

import math
from .labwareavailable import labware_plates


class CreateOTSession(object):
    """
    Creates a StartOTSession object for generating a protocol
    from actions
    """

    def __init__(
        self,
        reactionstep: int,
        otbatchprotocolobj: OTBatchProtocol,
        actionsessionqueryset: QuerySet[ActionSession],
        customSMcsvpath: str = None,
    ):
        """Initiates a CreateOTSession

        Parameters
        ----------
        reactionstep: int
            The reaction step that the protocol is being created for
        otbatchprotocolobj: OTBatchProtocol
            The Django OTBactchProtocol model object, collects protocols
            for a group, that all OT sessions are linked to
        actionsession: QuerySet[ActionSession]
            The action sessions being executed on the OT for a type eg. reaction,
            stir, analayse
        groupreactionqueryset: QuerySet[Reaction]
            The reactions the session needs to execute the actions for
        customSMcsvpath: str (Optional)
            The path to the custom starting material plate csv file
        """
        self.reactionstep = reactionstep
        self.otbatchprotocolobj = otbatchprotocolobj
        self.actionsessionqueryset = actionsessionqueryset
        self.customSMcsvpath = customSMcsvpath
        self.actionsessionnumber = actionsessionqueryset.values_list(
            "sessionnumber", flat=True
        )[0]
        self.actionsession_ids = actionsessionqueryset.values_list("id", flat=True)
        self.reaction_ids = [
            actionsession_obj.reaction_id.id
            for actionsession_obj in self.actionsessionqueryset
        ]
        self.groupreactionqueryset = getReactionQuerySet(reaction_ids=self.reaction_ids)
        self.otsessionqueryset = self.otbatchprotocolobj.otsessions.all()
        self.batchobj = Batch.objects.get(id=otbatchprotocolobj.batch_id_id)
        self.actionsessiontype = self.getActionSessionType()
        self.otsessionobj = self.createOTSessionModel()
        self.createActionSession()

    def createActionSession(self):
        """Calls the functions to create the appropriate action session"""
        actionSessionTypes = {
            "reaction": self.createReactionSession,
            "workup": self.createWorkUpSession,
            "analyse": self.createAnalyseSession,
        }

        if self.actionsessiontype in actionSessionTypes:
            actionSessionTypes[self.actionsessiontype]()

    def createReactionSession(self):
        """Creates a reaction OT session"""
        self.groupedreactiontemperaturequerysets = self.getGroupedTemperatureReactions(
            reactionqueryset=self.groupreactionqueryset
        )
        self.addactionqueryset = self.getAddActionQuerySet(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids,
        )
        self.addactionsdf = self.getAddActionsDataFrame(
            addactionqueryset=self.addactionqueryset
        )
        self.extractactionqueryset = self.getExtractActionQuerySet(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids,
        )
        self.roundedaddvolumes = self.getRoundedAddActionVolumes(
            addactionqueryset=self.addactionqueryset
        )
        self.roundedextractvolumes = self.getRoundedExtractActionVolumes(
            extractactionqueryset=self.extractactionqueryset
        )
        self.roundedvolumes = self.roundedaddvolumes + self.roundedextractvolumes
        self.deckobj = self.createDeckModel()
        self.tipracktype = self.getTipRackType(roundedvolumes=self.roundedvolumes)
        self.createTipRacks(tipracktype=self.tipracktype)
        self.pipettetype = self.getPipetteType(roundedvolumes=self.roundedvolumes)
        continuationactionsessions = self.actionsessionqueryset.filter(
            continuation=True
        )
        noncontinuationactionsessions = self.actionsessionqueryset.filter(
            continuation=False
        )
        if continuationactionsessions.exists():
            searchsmiles = getProductSmiles(reaction_ids=self.reaction_ids)
            searchsmiles += list(
                self.addactionqueryset.values_list("smiles", flat=True)
            )
        if noncontinuationactionsessions.exists():
            searchsmiles = self.addactionqueryset.values_list("smiles", flat=True)
            self.createReactionPlate(platetype="reaction")
        inputplatequeryset = self.getInputPlatesNeeded(searchsmiles=searchsmiles)
        self.updatePlateDeckOTSessionIDs(platequeryset=inputplatequeryset)
        self.createPipetteModel()
        if self.customSMcsvpath:
            self.createStartingMaterialPlatesFromCSV(csv_path=self.customSMcsvpath)
        self.createReactionStartingPlate()
        if self.reactionstep > 1:
            self.solventmaterialsdf = self.getAddActionsMaterialDataFrame(
                productexists=True
            )
            self.createSolventPlate(materialsdf=self.solventmaterialsdf)

    def createWorkUpSession(self):
        """Creates a workup OT session"""
        self.roundedvolumes = []
        self.addactionqueryset = self.getAddActionQuerySet(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids,
        )
        if self.addactionqueryset:
            self.roundedaddvolumes = self.getRoundedAddActionVolumes(
                addactionqueryset=self.addactionqueryset
            )
            self.roundedvolumes = self.roundedvolumes + self.roundedaddvolumes
        self.extractactionqueryset = self.getExtractActionQuerySet(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids,
        )
        if self.extractactionqueryset:
            self.roundedextractvolumes = self.getRoundedExtractActionVolumes(
                extractactionqueryset=self.extractactionqueryset
            )
            self.roundedvolumes = self.roundedvolumes + self.roundedextractvolumes
        self.deckobj = self.createDeckModel()
        self.tipracktype = self.getTipRackType(roundedvolumes=self.roundedvolumes)
        self.createTipRacks(tipracktype=self.tipracktype)
        searchsmiles = self.getProductSmiles(reaction_ids=self.reaction_ids)
        inputplatequeryset = self.getInputPlatesNeeded(searchsmiles=searchsmiles)
        self.updatePlateDeckOTSessionIDs(platequeryset=inputplatequeryset)
        self.pipettetype = self.getPipetteType(
            roundedvolumes=self.roundedvolumes, channeltype="multi"
        )
        self.addactionsdf = self.getAddActionsDataFrame(
            addactionqueryset=self.addactionqueryset
        )

        self.createPipetteModel()
        self.solventmaterialsdf = self.getAddActionsMaterialDataFrame(
            productexists=False
        )
        print("Creating solvent plates for workup")
        self.createSolventPlate(materialsdf=self.solventmaterialsdf)
        self.workupplatesneeded = self.getUniqueToPlates(
            actionsessionqueryset=self.actionsessionqueryset,
            platetypes=["workup1", "workup2", "workup3", "spefilter"],
        )
        for workuplateneeded in self.workupplatesneeded:
            self.createWorkUpPlate(platetype=workuplateneeded)

    def createAnalyseSession(self):
        """Creates an analyse OT session"""
        self.roundedvolumes = []
        self.addactionqueryset = self.getAddActionQuerySet(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids,
        )
        if self.addactionqueryset:
            self.roundedaddvolumes = self.getRoundedAddActionVolumes(
                addactionqueryset=self.addactionqueryset
            )
            self.roundedvolumes = self.roundedvolumes + self.roundedaddvolumes
        self.extractactionqueryset = self.getExtractActionQuerySet(
            reaction_ids=self.reaction_ids,
            actionsession_ids=self.actionsession_ids,
        )
        if self.extractactionqueryset:
            self.roundedextractvolumes = self.getRoundedExtractActionVolumes(
                extractactionqueryset=self.extractactionqueryset
            )
            self.roundedvolumes = self.roundedvolumes + self.roundedextractvolumes

        self.deckobj = self.createDeckModel()
        self.tipracktype = self.getTipRackType(roundedvolumes=self.roundedvolumes)
        self.createTipRacks(tipracktype=self.tipracktype)
        searchsmiles = self.getProductSmiles(reaction_ids=self.reaction_ids)
        inputplatequeryset = self.getInputPlatesNeeded(searchsmiles=searchsmiles)
        self.updatePlateDeckOTSessionIDs(platequeryset=inputplatequeryset)
        self.pipettetype = self.getPipetteType(
            roundedvolumes=self.roundedvolumes, channeltype="multi"
        )
        self.addactionsdf = self.getAddActionsDataFrame(
            addactionqueryset=self.addactionqueryset
        )
        self.createPipetteModel()
        self.solventmaterialsdf = self.getAddActionsMaterialDataFrame(
            productexists=False
        )
        self.createSolventPlate(materialsdf=self.solventmaterialsdf)
        self.analyseplatesneeded = self.getUniqueToPlates(
            actionsessionqueryset=self.actionsessionqueryset,
            platetypes=["lcms", "xchem"],
        )
        for analyseplateneeded in self.analyseplatesneeded:
            self.createAnalysePlate(platetype=analyseplateneeded)

    def getInputPlatesNeeded(
        self, searchsmiles: list, reaction_ids: list = None
    ) -> list[Plate]:
        """Gets plates, created in previous reaction and workup
        sessions with reaction products that are required as
        reactants in current reaction session

        Parameters
        ----------
        searchsmiles: list
            The list of SMILES that are required from previous
            reaction plate wells
        reaction_ids: list
            The optional reaction ids to match wells and plates with.

        Returns
        -------
        inputplatesneeded: list
            The list of previous OT session reaction plates in
            an OT batch protocol that have products needed as
            reactants for current reaction OT session
        """
        inputplatesneeded = []
        otbatchprotocolplatequeryset = self.getAllOTBatchProtocolPlates(
            otbatchprotocol_id=self.otbatchprotocolobj
        )
        if not reaction_ids:
            methodids = [
                reactionobj.method_id for reactionobj in self.groupreactionqueryset
            ]
            criterion1 = Q(method_id__in=methodids)
        if reaction_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)

        criterion2 = Q(reactantfornextstep=True)
        criterion3 = Q(smiles__in=searchsmiles)
        criterion4 = Q(
            type__in=["reaction", "workup1", "workup2", "workup3", "spefilter"]
        )
        if otbatchprotocolplatequeryset:
            for plateobj in otbatchprotocolplatequeryset:
                wellmatchqueryset = plateobj.well_set.all().filter(
                    criterion1 & criterion2 & criterion3 & criterion4
                )
                if wellmatchqueryset:
                    inputplatesneeded.append(plateobj)
        return inputplatesneeded

    def getAllOTBatchProtocolPlates(
        self, otbatchprotocol_id: OTBatchProtocol
    ) -> QuerySet[Plate]:
        """Get all input reaction plates used for an OT batch protocol

        Parameters
        ----------
        otbatchprotocol_id: OTBatchProtocol
            All OT batch protocol to find all matching plates for
        Returns
        -------
        otbatchprotocolplatequeryset: QuerySet[Plate]
            The plates used for all previous reaction and workup
            sessions
        status: False
            The status if no plates were found
        """
        criterion1 = Q(otbatchprotocol_id=otbatchprotocol_id)
        criterion2 = Q(
            type__in=["reaction", "workup1", "workup2", "workup3", "spefilter"]
        )

        otbatchprotocolplatequeryset = Plate.objects.filter(criterion1 & criterion2)
        return otbatchprotocolplatequeryset

    def updatePlateDeckOTSessionIDs(self, platequeryset: QuerySet[Plate]):
        """Updates the plates to link to the current Deck and
        OT session

        Parameters
        ----------
        platequeryset: QuerySet[Plate]
            The plates to update
        """
        for plateobj in platequeryset:
            indexslot = self.checkDeckSlotAvailable()
            if indexslot:
                wellqueryset = self.getPlateWells(plateobj=plateobj)
                columnqueryset = self.getPlateColumns(plateobj=plateobj)
                previoustype = plateobj.type
                plateobj.deck_id = self.deckobj
                plateobj.otsession_id = self.otsessionobj
                if previoustype == "spefilter":
                    plateobj.labware = "plateone_96_wellplate_2500ul"
                plateobj.index = indexslot
                plateobj.save()
                self.updateColumnOTSessionIDs(
                    columnqueryset=columnqueryset, plateobj=plateobj
                )
                self.updateWellOTSessionIDs(
                    wellqueryset=wellqueryset, plateobj=plateobj
                )
            else:
                print("cloneInputPlate")
                print("No more deck slots available")
        for plateobj in platequeryset:
            plateobj.deck_id = self.deckobj
            plateobj.otsession_id = self.otsessionobj
            plateobj.save()

    def getActionSessionByPlateType(self, platetype: str) -> QuerySet[ActionSession]:
        criterion1 = Q(id__in=self.actionsessionqueryset)
        criterion2 = Q(addaction__toplatetype=platetype)
        criterion3 = Q(extractaction__toplatetype=platetype)
        actionsessionqueryset = ActionSession.objects.filter(
            criterion1 & (criterion2 | criterion3)
        ).distinct()
        return actionsessionqueryset

    def getUniqueToPlates(
        self, actionsessionqueryset: QuerySet[ActionSession], platetypes: list
    ) -> list:
        """Gets the distinct to plate types to plates an action session
        queryset

        Parameters
        ----------
        actionsessionqueryset: QuerySet[ActionSession]
            The action queryset to get the to plates for
        platetypes: list
            The plate types to try and find in an action session
            queryset eg. ["reaction", "workup1"]

        Returns
        -------
        toplatetypes: list
            The to plates needed for an action session
        """
        criterion1 = Q(actionsession_id__in=actionsessionqueryset)
        criterion2 = Q(toplatetype__in=platetypes)
        toaddplates = (
            AddAction.objects.filter(criterion1 & criterion2)
            .values_list("toplatetype", flat=True)
            .distinct()
        )
        toextractplates = (
            ExtractAction.objects.filter(criterion1 & criterion2)
            .values_list("toplatetype", flat=True)
            .distinct()
        )
        toplatetypes = set(list(toaddplates) + list(toextractplates))
        return toplatetypes

    def getProductSmiles(self, reaction_ids: list) -> list:
        """Get product smiles of reactions

        Parameters
        ----------
        reaction_ids: list
            The reactions to get product smiles for

        Returns
        -------
        productsmiles: list
            The list of product smiles
        """

        productsmiles = Product.objects.filter(
            reaction_id__in=reaction_ids
        ).values_list("smiles", flat=True)
        if productsmiles:
            return list(productsmiles)
        else:
            return None

    def getAddActionQuerySet(
        self,
        reaction_ids: list,
        actionsession_ids: list = None,
        actionsessiontype: str = None,
    ) -> QuerySet[AddAction]:
        """Get add actions queryset for reaction_id

        Parameters
        ----------
        reaction_ids: list
            The reactions to search for related add actions
        actionsession_ids: list
            Optional action session ids to match add actions with
        actionsessiontype: str
            The optional action session type to look for add actions for

        Returns
        -------
        addactionqueryset: QuerySet[AddAction]
            The add actions related to the reaction
        """

        if actionsession_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__in=actionsession_ids)
            addactionqueryset = AddAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return addactionqueryset
        if actionsessiontype:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            addactionqueryset = AddAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return addactionqueryset

    def getExtractActionQuerySet(
        self,
        reaction_ids: list,
        actionsession_ids: list = None,
        actionsessiontype: str = None,
    ) -> QuerySet[ExtractAction]:
        """Get extract actions queryset for reaction_id

        Parameters
        ----------
        reaction_ids: list
            The reactions to search for related add actions
        actionsession_ids: list
            Optional action session ids to match add actions with
        actionsessiontype: str
            The optional action session type to look for add actions for

        Returns
        -------
        extractactionqueryset: QuerySet[ExtractAction]
            The extract actions related to the reaction
        """

        if actionsession_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__in=actionsession_ids)
            extractactionqueryset = ExtractAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")
            return extractactionqueryset
        if actionsessiontype:
            criterion1 = Q(reaction_id__in=reaction_ids)
            extractactionqueryset = ExtractAction.objects.filter(
                criterion1, actionsession_id__type=actionsessiontype
            ).order_by("id")
            return extractactionqueryset

    def getRoundedAddActionVolumes(
        self, addactionqueryset: QuerySet[AddAction]
    ) -> list:
        """Gets the total rounded volume (ul) for all the add actions

        Parameters
        ----------
        addactionqueryset: QuerySet[AddAction]
            The add actions to calculate the rounded volumes (ul)

        Returns
        -------
        roundedvolumes: list
            The list of rounded volumes for the add actions
        """
        roundedvolumes = [
            round(addactionobj.volume) for addactionobj in addactionqueryset
        ]
        return roundedvolumes

    def getRoundedExtractActionVolumes(
        self, extractactionqueryset: QuerySet[ExtractAction]
    ) -> list:
        """Gets the total rounded volume (ul) for all the extract actions

        Parameters
        ----------
        extractactionqueryset: QuerySet[ExtractAction]
            The extract actions to calculate the rounded volumes (ul)

        Returns
        -------
        roundedvolumes: list
            The list of rounded volumes for the extract actions
        """
        roundedvolumes = [
            round(extractactionobj.volume) for extractactionobj in extractactionqueryset
        ]
        return roundedvolumes

    def getRoundedReactionVolumes(
        self, groupedreactiontemperaturequeryset: QuerySet[Reaction]
    ) -> list:
        """Gets the total rounded volume (ul) for all the reactions done at the
           same temperature

        Parameters
        ----------
        groupedreactiontemperaturequerysets: QuerySet[Reaction]
            The reactions executed at the same temperature to
            calculate the rounded volumes (ul)

        Returns
        -------
        roundedvolumes: list
            The list of rounded volumes (ul) for the add actions
        """
        reactionvolumes = []
        for reactionobj in groupedreactiontemperaturequeryset:
            addactionqueryset = self.getAddActionQuerySet(
                reaction_ids=[reactionobj.id], actionsession_ids=self.actionsession_ids
            )
            roundedvolumes = self.getRoundedAddActionVolumes(
                addactionqueryset=addactionqueryset
            )
            sumvolume = self.getSumValue(values=roundedvolumes)
            reactionvolumes.append(sumvolume)
        return reactionvolumes

    def getTipRackType(self, roundedvolumes: list) -> str:
        """Gets OT tiprack best suited for transferring volumes (ul)
           that minimises the number of transfers required

        Parameters
        ----------
        roundedvolumes: list
            The list of rounded volumes needed for a set of actions

        Returns
        -------
        tipracktype: str
            The most suitable tiprack type
        """
        tipsavailable = {
            300: "opentrons_96_tiprack_300ul",
            10: "opentrons_96_tiprack_20ul",
        }
        tipkey = min(
            tipsavailable,
            key=lambda x: self.getNumberTransfers(
                pipettevolume=x, roundedvolumes=roundedvolumes
            ),
        )
        tipracktype = tipsavailable[tipkey]
        return tipracktype

    def getPipetteType(self, roundedvolumes: list, channeltype: str = "single") -> str:
        """Gets the type of pippete that minmises the number of transfers
           needed for transferring volumes (ul).

        Parameters
        ----------
        roundedvolumes: list
            The list of rounded volumes that need to be transferred
        channeltype: str
            The type of channel eg. single or multi. Default set to single

        Returns
        -------
        pipettetype: str
            The pipette type needed
        """
        pipettesavailable = [
            {
                "labware": "p10_single",
                "position": "right",
                "type": "single",
                "maxvolume": 10,
            },
            {
                "labware": "p300_single",
                "position": "right",
                "type": "single",
                "maxvolume": 300,
            },
            {
                "labware": "p300_multi_gen2",
                "position": "left",
                "type": "multi",
                "maxvolume": 300,
            },
            {
                "labware": "p10_multi",
                "position": "left",
                "type": "multi",
                "maxvolume": 10,
            },
        ]

        mediantransfervolume = self.getMedianValue(values=roundedvolumes)
        pipettecomparedict = {}

        for pipette in pipettesavailable:
            pipettevolume = pipette["maxvolume"]
            pipettetype = pipette["type"]
            pipettelabware = pipette["labware"]
            volumedifference = pipettevolume - mediantransfervolume
            if volumedifference < 0 or pipettetype != channeltype:
                continue
            numbertransfers = self.getNumberTransfers(
                pipettevolume=pipettevolume, roundedvolumes=roundedvolumes
            )
            pipettecomparedict[pipettelabware] = {
                "notransfers": numbertransfers,
                "volumedifference": volumedifference,
            }

        minimumnotransfers = min(
            [pipette["notransfers"] for pipette in pipettecomparedict.values()]
        )
        pipettetypes = [
            pipettelabware
            for pipettelabware in pipettecomparedict
            if pipettecomparedict[pipettelabware]["notransfers"] == minimumnotransfers
        ]
        if len(pipettetypes) > 1:
            minimumvolumedifference = min(
                [
                    pipettecomparedict[pipettelabware]["volumedifference"]
                    for pipettelabware in pipettetypes
                ]
            )
            pipettetypes = [
                pipettelabware
                for pipettelabware in pipettecomparedict
                if pipettecomparedict[pipettelabware]["volumedifference"]
                == minimumvolumedifference
            ]
        pipettetype = next(
            (
                pipette
                for pipette in pipettesavailable
                if pipette["labware"] == pipettetypes[0]
            ),
            None,
        )

        return pipettetype

    def getNumberTransfers(self, pipettevolume: int, roundedvolumes: list) -> int:
        """Gets the number of transfers required for transferring
           a list of rounded volumes

        Parameters
        ----------
        pipettevolume: int
            The pippette's maximum transfer volume (ul)
        roundedvolumes: list
            The list fo rounded volumes that need to be transferred

        Returns
        -------
        numbertransfers: int
            The number of tranfers required for the pipette type used
        """
        numbertransfers = sum(
            [
                round(volume / pipettevolume) if pipettevolume < volume else 1
                for volume in roundedvolumes
            ]
        )
        return numbertransfers

    def getAddActionsDataFrame(
        self, addactionqueryset: QuerySet[AddAction]
    ) -> DataFrame:
        """Creates a Pandas dataframe from the add actions.
           Current version of code uses a dataframe to create plates and wells.

        Parameters
        ----------
        addactionqueryset: QuerySet[AddAction]
            The add action queryset top convert to a dataframe

        Returns
        -------
        addactionsdf: DataFrame
            The dataframe of the add actions
        """
        # Optimise -> https://stackoverflow.com/questions/11697887/converting-django-queryset-to-pandas-dataframe
        addactionsdf = pd.DataFrame(list(addactionqueryset.values()))
        addactionsdf["uniquesolution"] = addactionsdf.apply(
            lambda row: self.combinestrings(row), axis=1
        )

        return addactionsdf

    def getMaxWellVolume(self, plateobj: Plate) -> float:
        """Get max well volume of a well plate

        Parameters
        ----------
        plateobj: Plate
            The plate to get the max well volume of

        Returns
        -------
        maxwellvolume: float
            The maximum well volume of a well plate
        """
        maxwellvolume = plateobj.maxwellvolume
        return maxwellvolume

    def getDeadVolume(self, maxwellvolume: float) -> float:
        """Calculates the dead volume (5%) of a well

        Parameters
        ----------
        maxwellvolume: float
            The well's maximum volume

        Returns
        -------
        deadvolume: float
            The dead volume of the well
        """
        deadvolume = maxwellvolume * 0.05
        return deadvolume

    def getPlateWells(self, plateobj: Plate) -> QuerySet[Well]:
        """Retrieves the wells for a plate

        Parameters
        ----------
        plateobj: Plate
            The plate to get all the related wells to

        Returns
        -------
        wellqueryset: QuerySet[Well]
            The plates wells
        """
        wellqueryset = Well.objects.filter(plate_id=plateobj.id)
        return wellqueryset

    def getPlateColumns(self, plateobj: Plate) -> QuerySet[Column]:
        """Retrieves the columns for a plate

        Parameters
        ----------
        plateobj: Plate
            The plate to get all the related wells to

        Returns
        -------
        columnqueryset: QuerySet[Column]
            The plates columns
        """
        columnqueryset = Column.objects.filter(plate_id=plateobj.id)
        return columnqueryset

    def getActionSessionType(self) -> str:
        """Get the action session type

        Returns
        ------
        actionsessiontype: list
            The set of action session types
        """
        actionsessiontype = self.actionsessionqueryset.values_list(
            "type", flat=True
        ).distinct()[0]
        return actionsessiontype

    def getUniqueTemperatures(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """Set of reaction temperatures

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to get unique list of reaction classes for

        Returns
        ------
        temperatures: list
            The set of reaction temperatures
        """
        temperatures = (
            self.groupreactionqueryset.values_list("temperature", flat=True)
            .order_by("temperature")
            .distinct()
        )
        return temperatures

    def getUniqueReactionClasses(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """Set of unique reaction classes

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to get unique list of reaction classes for

        Returns
        ------
        reactionclasses: list
            The set of reaction classes
        """
        reactionclasses = (
            reactionqueryset.values_list("reactionclass", flat=True)
            .order_by("reactionclass")
            .distinct()
        )
        return reactionclasses

    def getUniqueReactionRecipes(
        self, reactionclass: str, reactionqueryset: QuerySet[Reaction]
    ) -> list:
        """Set of unique reaction recipes

        Parameters
        ----------
        reactionclass: str
            The reaction type to find linked recipes for
        reactionqueryset: QuerySet[Reaction]
            The reactions to get unique list of reaction recipes for

        Returns
        ------
        reactionrecipes: list
            The set of reaction receipes
        """
        reactionrecipes = (
            reactionqueryset.filter(reactionclass=reactionclass)
            .values_list("recipe", flat=True)
            .order_by("recipe")
            .distinct()
        )
        return reactionrecipes

    def getGroupedTemperatureReactions(
        self, reactionqueryset: QuerySet[Reaction]
    ) -> list:
        """Group reactions done at the same temperature

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to to group by reaction class

        Returns
        -------
        groupedreactiontemperaturequerysets: list
            A list of sublists of reaction querysets grouped by reaction
            temperature
        """
        temperatures = self.getUniqueTemperatures(reactionqueryset=reactionqueryset)
        groupedreactiontemperaturequerysets = []
        for temperature in temperatures:
            reactiontemperaturequeryset = (
                reactionqueryset.filter(temperature=temperature)
                .distinct()
                .order_by("id")
            )
            if reactiontemperaturequeryset:
                groupedreactiontemperaturequerysets.append(reactiontemperaturequeryset)
        return groupedreactiontemperaturequerysets

    def getGroupedReactionByClassRecipe(
        self, reactionqueryset: QuerySet[Reaction]
    ) -> list:
        """Group reactions by reaction class and recipe type

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to group by reaction class

        Returns
        -------
        groupedreactionbyclassquerysets: list
            The list of sublists of reaction querysets grouped by reaction class
        """
        reactionclasses = self.getUniqueReactionClasses(
            reactionqueryset=reactionqueryset
        )
        groupedreactionquerysets = []

        for reactionclass in reactionclasses:
            reactionclassqueryset = reactionqueryset.filter(reactionclass=reactionclass)
            if reactionclassqueryset:
                reactionrecipes = (
                    reactionclassqueryset.values_list("recipe", flat=True)
                    .distinct()
                    .order_by("recipe")
                )
                if len(reactionrecipes) == 1:
                    groupedreactionquerysets.append(reactionclassqueryset)
                if len(reactionrecipes) > 1:
                    notgroupbycolumnreactionqueryset = reactionclassqueryset.filter(
                        groupbycolumn=False
                    )
                    if notgroupbycolumnreactionqueryset:
                        groupedreactionquerysets.append(
                            notgroupbycolumnreactionqueryset
                        )

                    groupbycolumnreactionqueryset = reactionclassqueryset.filter(
                        groupbycolumn=True
                    )
                    if groupbycolumnreactionqueryset:
                        for reactionrecipe in reactionrecipes:
                            reactionbyrecipequeryset = (
                                groupbycolumnreactionqueryset.filter(
                                    recipe=reactionrecipe,
                                    groupbycolumn=True,
                                )
                            )
                            if reactionbyrecipequeryset:
                                groupedreactionquerysets.append(
                                    reactionbyrecipequeryset
                                )

        return groupedreactionquerysets

    def getMedianValue(self, values: list) -> float:
        medianvalue = median(values)
        return medianvalue

    def getSumValue(self, values: list) -> float:
        sumvalue = sum(values)
        return sumvalue

    def getNumberVials(self, maxvolumevial: float, volumematerial: float) -> int:
        """Gets the total number of vials needed to prepare a starter plate

        Parameters
        ----------
        maxvolumevial: float
            The maximum volume of a vial
        volumematerial: float
            The volume of the material that needs to be stored in a vial

        Returns
        -------
        novialneeded: int
            The number of vials needed to store the material
        """
        if maxvolumevial > volumematerial:
            novialsneeded = 1
        else:
            volumestoadd = []
            deadvolume = self.getDeadVolume(maxwellvolume=maxvolumevial)
            novialsneededratio = volumematerial / (maxvolumevial - deadvolume)
            frac, whole = math.modf(novialsneededratio)
            volumestoadd = [maxvolumevial for maxvolumevial in range(int(whole))]
            volumestoadd.append(frac * maxvolumevial + deadvolume)
            novialsneeded = sum(volumestoadd)
        return novialsneeded

    def getPlateType(
        self,
        platetype: str,
        volumes: list,
        temperature: int = 25,
        wellsneeded: int = None,
    ):
        """Gets best suited plate
        Optimises for (in order of decreasing preference) minimum:
        number of plates, volume difference and mumber of vials

        Parameters
        ---------
        platetype: str
            The plateype eg. reaction, starting plate
        volumes: list
            The volumes that the plate and wells need to accomodate
        temperature: int = 25
            The temperature (degC) the plate will be used at
        wellsneeded: int = None
            Optional to specify if a specific number of wells are needed
        """
        platetype = platetype.lower()
        medianvolume = self.getMedianValue(values=volumes)

        possiblelabwareplatetypes = [
            labware_plate
            for labware_plate in labware_plates
            if platetype in labware_plates[labware_plate]["type"]
            and labware_plates[labware_plate]["max_temp"] >= temperature
        ]

        vialcomparedict = {}

        for labwareplate in possiblelabwareplatetypes:
            maxtemp = labware_plates[labwareplate]["max_temp"]
            maxvolumevial = labware_plates[labwareplate]["volume_well"]
            noplatevials = labware_plates[labwareplate]["no_wells"]
            if not wellsneeded:
                wellsneeded = sum(
                    [
                        self.getNumberVials(
                            maxvolumevial=maxvolumevial, volumematerial=volume
                        )
                        for volume in volumes
                    ]
                )
            noplatesneeded = int(math.ceil(wellsneeded / noplatevials))
            volumedifference = maxvolumevial - medianvolume
            tempdifference = maxtemp - temperature
            if platetype == "reaction":  # Reaction plates can only fill one well
                maxvolumeexceededtest = all(
                    [False if maxvolumevial - vol <= 0 else True for vol in volumes]
                )
                if (
                    volumedifference < 0
                    or tempdifference < 0
                    or not maxvolumeexceededtest
                ):
                    continue
            if platetype == "starting":  # Starting plates can fill multiple wells
                if volumedifference < 0 or tempdifference < 0:
                    continue
            vialcomparedict[labwareplate] = {
                "noplatesneeded": noplatesneeded,
                "volumedifference": volumedifference,
                "novialsneeded": wellsneeded,
                "tempdifference": tempdifference,
            }
        mininumnoplatesneeeded = min(
            (d["noplatesneeded"] for d in vialcomparedict.values())
        )
        mininumtempdifference = min(
            (d["tempdifference"] for d in vialcomparedict.values())
        )

        labwareplatetypes = [
            labwareplate
            for labwareplate in vialcomparedict
            if vialcomparedict[labwareplate]["noplatesneeded"] == mininumnoplatesneeeded
            and vialcomparedict[labwareplate]["tempdifference"] == mininumtempdifference
        ]
        if len(labwareplatetypes) > 1:
            mininumvolumedifference = min(
                (d["volumedifference"] for d in vialcomparedict.values())
            )
            labwareplatetypes = [
                labwareplate
                for labwareplate in vialcomparedict
                if vialcomparedict[labwareplate]["volumedifference"]
                == mininumvolumedifference
            ]
            if len(labwareplatetypes) > 1:
                mininumnovialsneeded = min(
                    (d["novialsneeded"] for d in vialcomparedict.values())
                )
                labwareplatetypes = [
                    labwareplate
                    for labwareplate in vialcomparedict
                    if vialcomparedict[labwareplate]["novialsneeded"]
                    == mininumnovialsneeded
                ]

        return labwareplatetypes[0]

    def getAddActionsMaterialDataFrame(self, productexists: bool) -> DataFrame:
        """Aggregates all add actions materials and sums up volume requires using solvent type and
        concentration. Checks existing materials and only requests additional volume if needed.

        Parameters
        ----------
        productexists: bool
            Set to true to check if the add action material needed is a product from
            a previous reaction

        Returns
        -------
        materialsdf: DataFrame
            The add action material as dataframe grouping materials by SMILES, concentration
            and solvent, with volumes adjusted based on existing materials
        """
        try:
            materialsdf = self.addactionsdf.groupby(["uniquesolution"]).agg(
                {
                    "reaction_id_id": "first",
                    "smiles": "first",
                    "volume": "sum",
                    "solvent": "first",
                    "concentration": "first",
                    "molecularweight": "first",
                }
            )

            # Apply canonicalization to SMILES for consistent comparison
            materialsdf["smiles"] = materialsdf["smiles"].apply(canonSmiles)

            # Check if materials are products from previous reactions
            materialsdf["productexists"] = materialsdf.apply(
                lambda row: checkPreviousReactionProducts(
                    reaction_id=row["reaction_id_id"], smiles=row["smiles"]
                ),
                axis=1,
            )

            # Check existing materials and get remaining volume needed
            material_status = []
            for idx, row in materialsdf.iterrows():
                (
                    exists,
                    matching_wells,
                    plate,
                    remaining_volume,
                ) = self.checkStartingMaterialExists(
                    smiles=row["smiles"],
                    volume=row["volume"],
                    concentration=row["concentration"],
                    solvent=row["solvent"],
                )

                material_status.append(
                    {
                        "exists": exists,  # True if enough volume exists
                        "remaining_volume": remaining_volume,  # 0 if exists=True, otherwise the additional volume needed
                    }
                )

            # Add material status information to materialsdf
            status_df = pd.DataFrame(material_status, index=materialsdf.index)
            materialsdf = pd.concat([materialsdf, status_df], axis=1)

            # Store original volume for reference
            materialsdf["original_volume"] = materialsdf["volume"]

            # Update volume to only request what's additionally needed
            materialsdf["volume"] = materialsdf["remaining_volume"]

            # Mark if material exists with enough volume
            materialsdf["materialexists"] = materialsdf["exists"]

            # Filter out materials based on existence and product status
            if productexists:
                # Keep only products from previous reactions
                materialsdf = materialsdf[materialsdf["productexists"]]
            else:
                # Keep only non-products
                materialsdf = materialsdf[~materialsdf["productexists"]]

            print("Materials df after filtering: {}".format(materialsdf))

            # Filter out materials that have enough volume (remaining_volume = 0)
            materialsdf = materialsdf[materialsdf["remaining_volume"] > 0]

            print("The materials df is: {}".format(materialsdf))

            # Log materials that need additional volume
            for idx, row in materialsdf.iterrows():
                if row["original_volume"] > row["volume"]:
                    logger.info(
                        f"Partial material available for {row['smiles']}. "
                        f"Total needed: {row['original_volume']}µL, "
                        f"Adding {row['volume']}µL more"
                    )

            materialsdf = materialsdf.sort_values(
                ["solvent", "volume"], ascending=False
            )
            return materialsdf

        except Exception as e:
            logger.error(f"Error in getAddActionsMaterialDataFrame: {str(e)}")
            return None

    def createOTSessionModel(self):
        """Create an OT Session object"""
        otsessionobj = OTSession()
        otsessionobj.otbatchprotocol_id = self.otbatchprotocolobj
        otsessionobj.reactionstep = self.reactionstep
        otsessionobj.sessiontype = self.actionsessiontype
        otsessionobj.save()
        return otsessionobj

    def createDeckModel(self):
        """Create a deck object"""
        deckobj = Deck()
        deckobj.otsession_id = self.otsessionobj
        deckobj.numberslots = 11
        deckobj.save()
        return deckobj

    def createPipetteModel(self):
        """Create a pipette object"""
        pipetteobj = Pipette()
        pipetteobj.otsession_id = self.otsessionobj
        pipetteobj.position = self.pipettetype["position"]
        pipetteobj.maxvolume = self.pipettetype["maxvolume"]
        pipetteobj.type = self.pipettetype["type"]
        pipetteobj.name = "{}_{}".format(
            self.pipettetype["position"], self.pipettetype["labware"]
        )
        pipetteobj.labware = self.pipettetype["labware"]
        pipetteobj.save()

    def createTiprackModel(self, name: str):
        """Creates TipRack object"""
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            index = indexslot
            tiprackobj = TipRack()
            tiprackobj.otsession_id = self.otsessionobj
            tiprackobj.deck_id = self.deckobj
            tiprackobj.name = "{}_{}".format(name, indexslot)
            tiprackobj.index = index
            tiprackobj.labware = name
            tiprackobj.save()
        else:
            print("createTiprackModel")
            print("No more deck slots available")

    def createPlateModel(self, platetype: str, platename: str, labwaretype: str):
        """Creates Plate model if Deck index is available"""
        indexslot = self.checkDeckSlotAvailable()
        if indexslot:
            plateindex = indexslot
            numberwellsincolumn = labware_plates[labwaretype]["no_wells_in_column"]
            maxwellvolume = labware_plates[labwaretype]["volume_well"]
            numberwells = labware_plates[labwaretype]["no_wells"]
            numbercolumns = labware_plates[labwaretype]["no_columns"]
            plateobj = Plate()
            plateobj.otbatchprotocol_id = self.otbatchprotocolobj
            plateobj.otsession_id = self.otsessionobj
            plateobj.deck_id = self.deckobj
            plateobj.labware = labwaretype
            plateobj.index = plateindex
            plateobj.type = platetype
            plateobj.maxwellvolume = maxwellvolume
            plateobj.numberwells = numberwells
            plateobj.numberwellsincolumn = numberwellsincolumn
            plateobj.numbercolumns = numbercolumns
            plateobj.save()
            plateobj.name = "Reaction_step_{}_{}_{}".format(
                self.reactionstep, platename, plateobj.id
            )
            plateobj.save()
            return plateobj
        else:
            print("CreatePlateModel")
            print("No more deck slots available")

    def createColumnModel(
        self,
        plateobj: Plate,
        columnindex: int,
        columntype: str,
        reactionclass: str,
    ) -> Column:
        """Creates a column object

        Parameters
        ----------
        plateobj: Plate
            The plate that the column is linked to
        columnindex: int
            The index of the well in the plate
        columntype: str
            The type of plate the column is used on for eg.
            reaction, workup1, workup2, lcms
        reactionclass: str
            The reaction class occupying the column eg. amidation.
            Only one type of reaction class can occupy a column

        Returns
        -------
        columnobj: Column
            The column created
        """
        columnobj = Column()
        columnobj.otsession_id = self.otsessionobj
        columnobj.plate_id = plateobj
        columnobj.index = columnindex
        columnobj.type = columntype
        columnobj.reactionclass = reactionclass
        columnobj.save()
        return columnobj

    def createWellModel(
        self,
        plateobj: Plate,
        welltype: str,
        wellindex: int,
        volume: float = None,
        reactionobj: Reaction = None,
        columnobj: Column = None,
        smiles: str = None,
        concentration: float = None,
        solvent: str = None,
        reactantfornextstep: bool = False,
    ) -> Well:
        """Creates a well object

        Parameters
        ----------
        plateobj: Plate
            The plate that the well is linked to
        welltype: str
            The well type eg. reaction, analyse
        wellindex: int
            The index of the well in the plate eg. 0, 1, 2, 3 etc
        wellname: str
            The name of the well eg. A1, B1, C1 etc
        volume: float = None
            The optional volume of the well contents
        reactionobj: Reaction = None
            The optional reaction the well is linked to
        columnobj: Column = None
            The optional column object the well is linked to
        smiles: str = None
            The optional contents of the well
        concentration: float = None
            The optional cocentration of the well contents
        solvent: str = None
            The optional solvent used to prepare the content of the well
        reactantfornextstep: bool = False
            The optional setting if the contents of the well are
            used in any proceeding reactions

        Returns
        -------
        wellobj: Well
            The well created
        """
        wellobj = Well()
        wellobj.otsession_id = self.otsessionobj
        wellobj.plate_id = plateobj
        if reactionobj:
            wellobj.reaction_id = reactionobj
            wellobj.method_id = reactionobj.method_id
        if columnobj:
            wellobj.column_id = columnobj
        wellobj.type = welltype
        wellobj.index = wellindex
        wellobj.name = wellIndexToWellName(
            wellindex=wellindex, platesize=plateobj.numberwells
        )
        wellobj.volume = volume
        wellobj.smiles = smiles
        wellobj.concentration = concentration
        wellobj.solvent = solvent
        wellobj.reactantfornextstep = reactantfornextstep
        wellobj.save()
        return wellobj

    def createCompoundOrderModel(
        self, orderdf: DataFrame, is_custom_starter_plate: bool = False
    ):
        """Creates a compound order object

        Parameters
        ----------
        orderdf: DataFrame
            DataFrame containing the compound order data
        is_custom_starter: bool, optional
            Flag to indicate this is a custom starter plate (default: False)
        """
        compoundorderobj = CompoundOrder()
        compoundorderobj.otsession_id = self.otsessionobj

        # Create a prefix based on whether this is a custom starter plate
        prefix = (
            f"custom-SMs-{self.actionsessiontype}-session-starterplate"
            if is_custom_starter_plate
            else f"{self.actionsessiontype}-session-starterplate"
        )

        csvdata = orderdf.to_csv(encoding="utf-8", index=False)
        ordercsv = default_storage.save(
            "compoundorders/"
            + f"{prefix}-for-batch-{self.batchobj.batchtag}-reactionstep-{self.reactionstep}-sessionid-{str(self.otsessionobj.id)}"
            + ".csv",
            ContentFile(csvdata),
        )

        compoundorderobj.ordercsv = ordercsv
        compoundorderobj.iscustomSMplate = is_custom_starter_plate  # This requires adding this field to the CompoundOrder model
        compoundorderobj.save()

        if is_custom_starter_plate:
            logger.info(
                f"Created compound order model for custom starting material plate"
            )
        else:
            logger.info(f"Created standard compound order model")

    def createSolventPrepModel(self, solventdf: DataFrame):
        """Creates a Django solvent prep object - a solvent prep file

        Parameters
        ----------
        solventdf: DataFrame
            The solvent dataframe grouped by type of solvent
            and contains the: platenames, well index, solvent
            and volume required
        """
        solventprepobj = SolventPrep()
        solventprepobj.otsession_id = self.otsessionobj
        csvdata = solventdf.to_csv(encoding="utf-8", index=False)
        ordercsv = default_storage.save(
            "solventprep/"
            + "{}-session-solventplate-for-batch-{}-reactionstep-{}-sessionid-{}".format(
                self.actionsessiontype,
                self.batchobj.batchtag,
                self.reactionstep,
                str(self.otsessionobj.id),
            )
            + ".csv",
            ContentFile(csvdata),
        )
        solventprepobj.solventprepcsv = ordercsv
        solventprepobj.save()

    def createTipRacks(self, tipracktype: str, numbertipracks: int = 3):
        """Creates three tipracks

        Parameters
        ----------
        tipracktype: str
            The type of tiprack needed

        numbertipracks: int
            The number of tip racks to create. Default is three.
        """
        for _ in range(numbertipracks):
            self.createTiprackModel(name=tipracktype)

    def calcMass(self, row) -> float:
        """Calculates the mass of material (mg) from the
           concentration (mol/L) and volume (ul) needed

        Parameters
        ----------
        row: DataFrame row
            The row from the dataframe containing the
            concentration and volume information

        Retruns
        -------
        massmg: float
            The mass of the material needed
        """
        mols = row["concentration"] * row["amount-uL"] * 1e-6
        smiles = row["SMILES"]
        mw = Descriptors.MolWt(Chem.MolFromSmiles(smiles))
        massmg = mols * mw * 1e3
        return round(massmg, 2)

    def checkDeckSlotAvailable(self) -> int:
        """Check if a deck slot is available

        Returns
        -------
        testslotavailable: int
            The index of the deck slot available
        status: False
            Returns false if no deck slot/index is available
        """
        testslotavailable = self.deckobj.indexslotavailable
        if testslotavailable <= self.deckobj.numberslots:
            self.deckobj.indexslotavailable = testslotavailable + 1
            self.deckobj.save()
            return testslotavailable
        else:
            self.deckobj.slotavailable = False
            self.deckobj.save()
            return False

    def checkStartingMaterialExists(
        self, smiles: str, volume: float, concentration: float, solvent: str
    ) -> tuple:
        """Checks if starting material exists with enough total volume across wells in current OT batch protocol.

        Parameters
        ----------
        smiles: str
            The SMILES string of the starting material to check
        volume: float
            The total volume needed in microliters
        concentration: float
            The concentration in moles per liter
        solvent: str
            The solvent used

        Returns
        -------
        tuple
            (exists: bool, matching_wells: list[Well], plate: Plate, remaining_volume_needed: float)
            Returns if material exists with enough volume, list of matching wells, plate containing wells,
            and remaining volume needed (0 if all volume is available)
        """
        try:
            # Canonicalize the SMILES for consistent comparison
            canonical_smiles = canonSmiles(smiles)

            # Get all plates in current OT batch protocol
            plates = Plate.objects.filter(
                otbatchprotocol_id=self.otbatchprotocolobj, type="startingmaterial"
            )

            # Track total volume and matching wells across all plates
            total_volume = 0
            all_matching_wells = []
            containing_plate = None
            remaining_volume_needed = volume  # Initialize with full volume needed

            for plate in plates:
                # Find wells with matching material properties
                matching_wells = Well.objects.filter(
                    plate_id=plate,
                    smiles=canonical_smiles,
                    concentration=concentration,
                    solvent=solvent,
                    type="startingmaterial",
                ).order_by(
                    "volume"
                )  # Order by volume for optimal usage

                if matching_wells.exists():
                    containing_plate = plate

                    # Sum up volumes until we have enough
                    for well in matching_wells:
                        all_matching_wells.append(well)
                        total_volume += well.volume
                        remaining_volume_needed = max(0, volume - total_volume)

                        if total_volume >= volume:
                            logger.info(
                                f"Found enough material: {canonical_smiles} across {len(all_matching_wells)} wells "
                                f"in plate {plate.name}. Required: {volume}µL, Available: {total_volume}µL"
                            )
                            return (
                                True,
                                all_matching_wells,
                                containing_plate,
                                0,
                            )  # No additional volume needed

            if all_matching_wells:
                logger.warning(
                    f"Found material {canonical_smiles} but insufficient volume. "
                    f"Required: {volume}µL, Available: {total_volume}µL, "
                    f"Still need: {remaining_volume_needed}µL"
                )
            else:
                logger.info(f"No existing material found for: {canonical_smiles}")

            return False, all_matching_wells, containing_plate, remaining_volume_needed

        except Exception as e:
            logger.error(f"Error checking starting material existence: {str(e)}")
            return False, None, None, volume  # Return full volume needed on error

    def getPlateCurrentColumnIndex(self, plateobj: Plate) -> int:
        """Check if any columns available on a plate

        Parameters
        ----------
        plateobj: Plate
            The plate to search for a well available

        Returns
        -------
        plateobj.indexcolumnavailable: int
            The index of the column available on a plate
        status: False
            Returns false if no column is available
        """
        indexcolumnavailable = plateobj.indexcolumnavailable
        numbercolumns = plateobj.numbercolumns
        if indexcolumnavailable + 1 <= numbercolumns:
            return indexcolumnavailable
        else:
            return False

    def getPlateWellIndexAvailable(self, plateobj: Plate) -> int:
        """Check if any wells available on a plate
        Parameters
        ----------
        plateobj: Plate
            The plate to search for a well available
        Returns
        -------
        plateobj.indexswellavailable: int
            The index of the well available on a plate
        status: False
            Returns false if no well is available
        """
        indexwellavailable = plateobj.indexswellavailable
        numberwells = plateobj.numberwells
        if indexwellavailable + 1 <= numberwells:
            return indexwellavailable
        else:
            return False

    def updatePlateWellIndex(self, plateobj: Plate, wellindexupdate: int):
        """Updates the plates well index used

        Parameters
        ----------
        plateobj: Plate
            The plate to update the next available well index
        wellindexupdate: int
            The well index to update on the plate
        """
        plateobj.indexswellavailable = wellindexupdate
        plateobj.save()

    def updatePlateColumnIndexAvailable(self, plateobj: Plate, columnindexupdate: int):
        """Updates the column index used

        Parameters
        ----------
        plateobj: Plate
            The plate to update the next available column index
        columnindexupdate: int
            The column index to update on the plate
        """
        plateobj.indexcolumnavailable = columnindexupdate
        plateobj.save()

    def combinestrings(self, row):
        return (
            str(row["smiles"])
            + "-"
            + str(row["solvent"])
            + "-"
            + str(row["concentration"])
        )

    def updateWellOTSessionIDs(self, wellqueryset: QuerySet[Well], plateobj: Plate):
        """Updates well to link to current OT session

        Parameters
        ----------
        wellqueryset: QuerySet[Well]
            The wells to be updated
        plateobj: Plate
            The plate object related to the cloned well
        """
        for wellobj in wellqueryset:
            wellobj.plate_id = plateobj
            wellobj.otsession_id = self.otsessionobj
            wellobj.save()

    def updateColumnOTSessionIDs(self, columnqueryset: QuerySet[Well], plateobj: Plate):
        """Updates column to link to current OT session

        Parameters
        ----------
        columnqueryset: QuerySet[Column]
            The columns to be updated
        plateobj: Plate
            The plate object related to the updated column
        """
        for columnobj in columnqueryset:
            columnobj.plate_id = plateobj
            columnobj.otsession_id = self.otsessionobj
            columnobj.save()

    def createStartingMaterialPlatesFromCSV(self, csv_path: str) -> dict:
        """Creates starting material plates from a CSV file, supporting multiple plates identified by plate-ID."""
        try:
            df = pd.read_csv(csv_path)
            required_cols = ["plate-ID", "labware-type", "well-index", "SMILES", "amount-uL"]
            if not all(col in df.columns for col in required_cols):
                missing_cols = [col for col in required_cols if col not in df.columns]
                raise ValueError(f"CSV must contain columns: {required_cols}. Missing: {missing_cols}")
            
            # Keep the original SMILES for use in the compound records
            df["original_smiles"] = df["SMILES"].copy()
            
            # Canonicalize both original and desalted versions for consistent comparison
            df["canonical_smiles"] = df["SMILES"].apply(canonSmiles)
            df["desalted_smiles"] = df["canonical_smiles"].apply(stripSalts)
            
            # Log cases where salts were removed
            for i, (orig, desalted) in enumerate(zip(df["canonical_smiles"], df["desalted_smiles"])):
                if orig != desalted:
                    logger.info(f"Removed salts from CSV entry {i}: {orig} -> {desalted}")

            # Get required SMILES from the current reaction session's add actions
            if not hasattr(self, "addactionsdf") or self.addactionsdf.empty:
                logger.warning("No add actions found in current session, skipping custom starting materials")
                return {}

            # Create dictionary to store add action SMILES in both original and desalted forms
            required_smiles_map = {}
            
            for _, row in self.addactionsdf.iterrows():
                original = canonSmiles(row["smiles"])
                desalted = stripSalts(original)
                
                # Store both forms with original as the value - we'll use the original in the well
                required_smiles_map[original] = original
                required_smiles_map[desalted] = original
            
            logger.info(f"Required unique SMILES for reaction session: {len(set(required_smiles_map.values()))}")
            
            # Create a map to track which SM SMILES corresponds to which add action SMILES
            matching_smiles_map = {}
            matched_rows = []
            
            # Check all combinations for matches
            for i, row in df.iterrows():
                csv_canonical = row["canonical_smiles"]
                csv_desalted = row["desalted_smiles"]
                
                # Check if any form matches any required SMILES
                matched = False
                matched_add_action_smiles = None
                
                # Check original CSV SMILES against all required SMILES
                if csv_canonical in required_smiles_map:
                    matched = True
                    matched_add_action_smiles = required_smiles_map[csv_canonical]
                    logger.info(f"CSV SMILES {csv_canonical} matched add action SMILES directly")
                
                # Check desalted CSV SMILES against all required SMILES
                elif csv_desalted in required_smiles_map:
                    matched = True
                    matched_add_action_smiles = required_smiles_map[csv_desalted]
                    logger.info(f"Desalted CSV SMILES {csv_desalted} matched add action SMILES")
                
                if matched:
                    matching_smiles_map[i] = matched_add_action_smiles
                    matched_rows.append(i)
                    
            # Create filtered dataframe with only matching rows
            if not matched_rows:
                logger.warning("No matching SMILES found in custom starting materials CSV")
                return {}
                
            filtered_df = df.loc[matched_rows].copy()
            logger.info(f"Found {len(filtered_df)} matching entries in custom starting materials CSV")
            
            # Add the matched add action SMILES to use for wells
            filtered_df["add_action_smiles"] = filtered_df.index.map(matching_smiles_map)

            # Dictionary to store created plates
            created_plates = {}

            # List to collect data for compound order
            custom_plate_entries = []

            # First group by plate-ID, then by labware_type
            for plate_id, plate_df in filtered_df.groupby("plate-ID"):
                logger.info(f"Processing plate ID: {plate_id} with {len(plate_df)} entries")
                
                for labware_type, group_df in plate_df.groupby("labware-type"):
                    if labware_type not in labware_plates:
                        raise ValueError(f"Invalid labware type: {labware_type}")
                    
                    plate_key = f"{plate_id}_{labware_type}"
                    
                    plateobj = self.createPlateModel(
                        platetype="startingmaterial",
                        platename=f"Startingplate_{plate_id}",  # Include plate ID in name
                        labwaretype=labware_type,
                    )

                    created_plates[plate_key] = plateobj
                    logger.info(f"Created plate {plateobj.name} for plate ID: {plate_id}, labware type: {labware_type}")

                    for _, row in group_df.iterrows():
                        indexwellavailable = self.getPlateWellIndexAvailable(plateobj=plateobj)
                        if type(indexwellavailable) == bool:
                            plateobj = self.createPlateModel(
                                platetype="startingmaterial",
                                platename=f"Startingplate_{plate_id}",
                                labwaretype=labware_type,
                            )
                            indexwellavailable = self.getPlateWellIndexAvailable(plateobj=plateobj)

                        # Create the well using the ORIGINAL ADD ACTION SMILES for proper matching
                        wellobj = self.createWellModel(
                            plateobj=plateobj,
                            welltype="startingmaterial",
                            wellindex=int(row["well-index"]),
                            volume=float(row["amount-uL"]),
                            smiles=str(row["add_action_smiles"]),  # Use add action SMILES for well
                            concentration=float(row["concentration"]) if "concentration" in row else None,
                            solvent=str(row["solvent"]) if "solvent" in row else None,
                        )

                        # Update plate well index
                        self.updatePlateWellIndex(
                            plateobj=plateobj, wellindexupdate=indexwellavailable + 1
                        )

                        # Add entry for compound order using the ORIGINAL SMILES (with salts)
                        custom_plate_entries.append(
                            {
                                "SMILES": row["original_smiles"],  # Use original CSV SMILES with salts
                                "plate-ID": row["plate-ID"],  # Include plate ID in compound order
                                "plate-name": plateobj.name,
                                "labware": labware_type,
                                "well-index": row["well-index"],
                                "well-name": wellIndexToWellName(
                                    wellindex=row["well-index"],
                                    platesize=labware_plates[labware_type]["no_wells"],
                                ),
                                "concentration": row["concentration"] if "concentration" in row else None,
                                "solvent": row["solvent"] if "solvent" in row else None,
                                "amount-uL": float(row["amount-uL"]),
                                "molecularweight": None,
                            }
                        )

                    logger.info(f"Added {len(group_df)} wells to plate {plateobj.name}")

            # Create DataFrame for compound order
            if custom_plate_entries:
                customplatedf = pd.DataFrame(custom_plate_entries)
                
                # Calculate molecular weight and add to DataFrame
                customplatedf["molecularweight"] = customplatedf["SMILES"].apply(
                    lambda smiles: Descriptors.MolWt(Chem.MolFromSmiles(smiles))
                    if Chem.MolFromSmiles(smiles)
                    else None
                )

                # Calculate mass
                customplatedf["mass-mg"] = customplatedf.apply(
                    lambda row: self.calcMass(row), axis=1
                )

                # Add inchikey
                customplatedf["inchikey"] = customplatedf.apply(
                    lambda row: getInchiKey(row["SMILES"]), axis=1
                )

                # Add compound name
                compoundnames = customplatedf.apply(
                    lambda row: getChemicalName(row["inchikey"]), axis=1
                )
                customplatedf.insert(1, column="compound-name", value=compoundnames)

                # Create compound order
                self.createCompoundOrderModel(
                    orderdf=customplatedf, is_custom_starter_plate=True
                )
                
                logger.info(f"Created compound order for {len(customplatedf)} custom starting materials across {len(customplatedf['plate-ID'].unique())} plates")
                
            return created_plates

        except Exception as e:
            logger.error(f"Error creating starting material plates from CSV: {str(e)}")
            raise

    def createReactionStartingPlate(self):
        """Creates the starting material plate/s for executing a reaction's add actions"""
        startingmaterialsdf = self.getAddActionsMaterialDataFrame(productexists=False)
        if startingmaterialsdf is not None and not startingmaterialsdf.empty:

            logger.info("Creating starting material plates for new materials.")

            startinglabwareplatetype = self.getPlateType(
                platetype="startingmaterial", volumes=startingmaterialsdf["volume"]
            )

            plateobj = self.createPlateModel(
                platetype="startingmaterial",
                platename="Startingplate",
                labwaretype=startinglabwareplatetype,
            )
            maxwellvolume = self.getMaxWellVolume(plateobj=plateobj)
            deadvolume = self.getDeadVolume(maxwellvolume=maxwellvolume)
            orderdictslist = []
            for i in startingmaterialsdf.index.values:
                extraerrorvolume = startingmaterialsdf.at[i, "volume"] * 0.05
                totalvolume = startingmaterialsdf.at[i, "volume"] + extraerrorvolume
                if totalvolume > maxwellvolume:
                    nowellsneededratio = totalvolume / (maxwellvolume - deadvolume)

                    frac, whole = math.modf(nowellsneededratio)
                    volumestoadd = [maxwellvolume for i in range(int(whole))]
                    volumestoadd.append(frac * maxwellvolume + deadvolume)

                    for volumetoadd in volumestoadd:
                        indexwellavailable = self.getPlateWellIndexAvailable(
                            plateobj=plateobj
                        )
                        if type(indexwellavailable) == bool:
                            plateobj = self.createPlateModel(
                                platetype="startingmaterial",
                                platename="Startingplate",
                                labwaretype=startinglabwareplatetype,
                            )
                            indexwellavailable = self.getPlateWellIndexAvailable(
                                plateobj=plateobj
                            )

                        wellobj = self.createWellModel(
                            plateobj=plateobj,
                            welltype="startingmaterial",
                            wellindex=indexwellavailable,
                            volume=volumetoadd,
                            smiles=startingmaterialsdf.at[i, "smiles"],
                            concentration=startingmaterialsdf.at[i, "concentration"],
                            solvent=startingmaterialsdf.at[i, "solvent"],
                        )
                        self.updatePlateWellIndex(
                            plateobj=plateobj, wellindexupdate=indexwellavailable + 1
                        )

                        orderdictslist.append(
                            {
                                "SMILES": startingmaterialsdf.at[i, "smiles"],
                                "plate-name": plateobj.name,
                                "labware": plateobj.labware,
                                "well-index": wellobj.index,
                                "well-name": wellobj.name,
                                "concentration": startingmaterialsdf.at[
                                    i, "concentration"
                                ],
                                "solvent": startingmaterialsdf.at[i, "solvent"],
                                "molecularweight": startingmaterialsdf.at[
                                    i, "molecularweight"
                                ],
                                "amount-uL": round(volumetoadd, 2),
                            }
                        )

                else:
                    indexwellavailable = self.getPlateWellIndexAvailable(
                        plateobj=plateobj
                    )
                    volumetoadd = totalvolume + deadvolume
                    if type(indexwellavailable) == bool:
                        plateobj = self.createPlateModel(
                            platetype="startingmaterial",
                            platename="Startingplate",
                            labwaretype=startinglabwareplatetype,
                        )
                        indexwellavailable = self.getPlateWellIndexAvailable(
                            plateobj=plateobj
                        )

                    wellobj = self.createWellModel(
                        plateobj=plateobj,
                        reactionobj=getReaction(
                            reaction_id=startingmaterialsdf.at[i, "reaction_id_id"]
                        ),
                        welltype="startingmaterial",
                        wellindex=indexwellavailable,
                        volume=volumetoadd,
                        smiles=startingmaterialsdf.at[i, "smiles"],
                        concentration=startingmaterialsdf.at[i, "concentration"],
                        solvent=startingmaterialsdf.at[i, "solvent"],
                    )
                    self.updatePlateWellIndex(
                        plateobj=plateobj, wellindexupdate=indexwellavailable + 1
                    )
                    orderdictslist.append(
                        {
                            "SMILES": startingmaterialsdf.at[i, "smiles"],
                            "plate-name": plateobj.name,
                            "labware": plateobj.labware,
                            "well-index": wellobj.index,
                            "well-name": wellobj.name,
                            "concentration": startingmaterialsdf.at[i, "concentration"],
                            "solvent": startingmaterialsdf.at[i, "solvent"],
                            "molecularweight": startingmaterialsdf.at[
                                i, "molecularweight"
                            ],
                            "amount-uL": round(volumetoadd, 2),
                        }
                    )

            orderdf = pd.DataFrame(orderdictslist)
            orderdf["mass-mg"] = orderdf.apply(lambda row: self.calcMass(row), axis=1)
            orderdf["inchikey"] = orderdf.apply(
                lambda row: getInchiKey(row["SMILES"]), axis=1
            )
            compoundnames = orderdf.apply(
                lambda row: getChemicalName(row["inchikey"]), axis=1
            )
            orderdf.insert(1, column="compound-name", value=compoundnames)
            self.createCompoundOrderModel(orderdf=orderdf)

    def getNewColumnAndWellIndexAvailable(self, plateobj: Plate) -> tuple:
        """Checks if a new column is available and updates the plates
        column and well index
           eg. Current plate column index = 1 and well index = 5,
               then set column index to 2 and well index to the new column's
               starting well index
               New well index = 1 * number wells in column (8) = 8
               New column index = 2

        Parameters
        ----------
        plateobj: Plate
            The plate to search for a column available

        Returns
        -------
        newcolumnindex: int
            The index of the new column
        newcolumnwellindex : int
            The well index that matches a new column
        """
        wellindexcorrection = plateobj.numberwellsincolumn
        indexcolumnavailable = plateobj.indexcolumnavailable
        if indexcolumnavailable + 1 <= plateobj.numbercolumns:
            newwellindex = indexcolumnavailable * wellindexcorrection
            return (indexcolumnavailable, newwellindex)
        else:
            return False

    def checkIndexWellIsNewColumn(self, plateobj: Plate):
        """Checks if current available well index on plate is the beginning
           of a new plate column

        Parameters
        ----------
        plateobj: Plate
            The plate containing the well to check

        Returns
        -------
        isnewcolumn: bool
            True id the well is at the start of a new plate column
        """
        indexwellavailable = plateobj.indexswellavailable
        numberwellsincolumn = plateobj.numberwellsincolumn
        if indexwellavailable == 0:
            return True
        if indexwellavailable != 0:
            if (indexwellavailable % numberwellsincolumn) == 0:
                return True
            if (indexwellavailable % numberwellsincolumn) != 0:
                return False

    def createPlateByReactionClassRecipe(
        self,
        reactionqueryset: QuerySet[Reaction],
        labwareplatetype: str,
        platetype: str,
    ):
        """Creates plates based on reaction class where columns in a plate
           are occupied by the same reaction

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to group by reaction class
        labwareplatetype: str
            The type of labware to be created eg. "plateone_96_wellplate_2500ul"
        platetype: str
            The type of plate being created eg. reaction, workup1, workup2, workup3, analyse
        """
        platename = "{}_plate".format(platetype.capitalize())
        groupedreactionclassrecipequerysets = self.getGroupedReactionByClassRecipe(
            reactionqueryset=reactionqueryset
        )
        plateobj = self.createPlateModel(
            platetype=platetype,
            platename=platename,
            labwaretype=labwareplatetype,
        )
        for groupreactionclassrecipequeryset in groupedreactionclassrecipequerysets:
            reactionclass = groupreactionclassrecipequeryset.values_list(
                "reactionclass", flat=True
            ).distinct()[0]
            indexwellavailable = self.getPlateWellIndexAvailable(plateobj=plateobj)
            if type(indexwellavailable) == bool:
                plateobj = self.createPlateModel(
                    platetype=platetype,
                    platename=platename,
                    labwaretype=labwareplatetype,
                )
                indexcurrentcolumn = self.getPlateCurrentColumnIndex(plateobj=plateobj)
                columnobj = self.createColumnModel(
                    plateobj=plateobj,
                    columnindex=indexcurrentcolumn,
                    columntype=platetype,
                    reactionclass=reactionclass,
                )
                self.updatePlateColumnIndexAvailable(
                    plateobj=plateobj, columnindexupdate=indexcurrentcolumn + 1
                )
            if type(indexwellavailable) == int:
                wellindexisnewcolumn = self.checkIndexWellIsNewColumn(plateobj=plateobj)
                if wellindexisnewcolumn:
                    indexcurrentcolumn = self.getPlateCurrentColumnIndex(
                        plateobj=plateobj
                    )
                    columnobj = self.createColumnModel(
                        plateobj=plateobj,
                        columnindex=indexcurrentcolumn,
                        columntype=platetype,
                        reactionclass=reactionclass,
                    )
                    self.updatePlateColumnIndexAvailable(
                        plateobj=plateobj, columnindexupdate=indexcurrentcolumn + 1
                    )
                if not wellindexisnewcolumn:
                    indexnewcolumnwellavailable = (
                        self.getNewColumnAndWellIndexAvailable(plateobj=plateobj)
                    )
                    if type(indexnewcolumnwellavailable) == bool:
                        plateobj = self.createPlateModel(
                            platetype=platetype,
                            platename=platename,
                            labwaretype=labwareplatetype,
                        )
                        indexcurrentcolumn = self.getPlateCurrentColumnIndex(
                            plateobj=plateobj
                        )
                        columnobj = self.createColumnModel(
                            plateobj=plateobj,
                            columnindex=indexcurrentcolumn,
                            columntype=platetype,
                            reactionclass=reactionclass,
                        )
                        self.updatePlateColumnIndexAvailable(
                            plateobj=plateobj, columnindexupdate=indexcurrentcolumn + 1
                        )
                    if indexnewcolumnwellavailable:
                        (
                            indexcurrentcolumn,
                            indexwellavailable,
                        ) = indexnewcolumnwellavailable
                        columnobj = self.createColumnModel(
                            plateobj=plateobj,
                            columnindex=indexcurrentcolumn,
                            columntype=platetype,
                            reactionclass=reactionclass,
                        )
                        self.updatePlateColumnIndexAvailable(
                            plateobj=plateobj, columnindexupdate=indexcurrentcolumn + 1
                        )
                        self.updatePlateWellIndex(
                            plateobj=plateobj, wellindexupdate=indexwellavailable
                        )

            for index, reactionobj in enumerate(groupreactionclassrecipequeryset):
                productobj = getProduct(reaction_id=reactionobj.id)
                indexwellavailable = self.getPlateWellIndexAvailable(plateobj=plateobj)
                if type(indexwellavailable) == bool:
                    plateobj = self.createPlateModel(
                        platetype=platetype,
                        platename=platename,
                        labwaretype=labwareplatetype,
                    )
                    indexcurrentcolumn = self.getPlateCurrentColumnIndex(
                        plateobj=plateobj
                    )
                    columnobj = self.createColumnModel(
                        plateobj=plateobj,
                        columnindex=indexcurrentcolumn,
                        columntype=platetype,
                        reactionclass=reactionclass,
                    )
                    self.updatePlateColumnIndexAvailable(
                        plateobj=plateobj, columnindexupdate=indexcurrentcolumn + 1
                    )
                if index != 0 and type(indexwellavailable) == int:
                    wellindexisnewcolumn = self.checkIndexWellIsNewColumn(
                        plateobj=plateobj
                    )
                    if wellindexisnewcolumn:
                        indexcurrentcolumn = self.getPlateCurrentColumnIndex(
                            plateobj=plateobj
                        )
                        columnobj = self.createColumnModel(
                            plateobj=plateobj,
                            columnindex=indexcurrentcolumn,
                            columntype=platetype,
                            reactionclass=reactionclass,
                        )
                        self.updatePlateColumnIndexAvailable(
                            plateobj=plateobj, columnindexupdate=indexcurrentcolumn + 1
                        )
                indexwellavailable = self.getPlateWellIndexAvailable(plateobj=plateobj)
                if platetype == "reaction":
                    reactantfornextstep = True
                else:
                    reactantfornextstep = False
                self.createWellModel(
                    plateobj=plateobj,
                    reactionobj=reactionobj,
                    columnobj=columnobj,
                    welltype=platetype,
                    wellindex=indexwellavailable,
                    smiles=productobj.smiles,
                    reactantfornextstep=reactantfornextstep,
                )
                self.updatePlateWellIndex(
                    plateobj=plateobj, wellindexupdate=indexwellavailable + 1
                )

    def createReactionPlate(self, platetype: str):
        """Creates reaction plate/s for executing reaction's add actions"""
        for (
            groupreactiontemperaturequeryset
        ) in self.groupedreactiontemperaturequerysets:
            reactiontemperature = groupreactiontemperaturequeryset[0].temperature
            volumes = self.getRoundedReactionVolumes(
                groupedreactiontemperaturequeryset=groupreactiontemperaturequeryset
            )
            labwareplatetype = self.getPlateType(
                platetype=platetype, temperature=reactiontemperature, volumes=volumes
            )
            self.createPlateByReactionClassRecipe(
                reactionqueryset=groupreactiontemperaturequeryset,
                labwareplatetype=labwareplatetype,
                platetype=platetype,
            )

    def createWorkUpPlate(self, platetype: str):
        """Creates workup plate/s for executing work up actions"""

        actionsessionqueryset = self.getActionSessionByPlateType(platetype=platetype)
        wellsneeded = len(actionsessionqueryset)
        plateroundedvolumes = []
        addactionqueryset = AddAction.objects.filter(
            actionsession_id__in=actionsessionqueryset, toplatetype=platetype
        )
        if addactionqueryset:
            plateroundedvolumes += self.getRoundedAddActionVolumes(
                addactionqueryset=addactionqueryset
            )
        extractactionqueryset = ExtractAction.objects.filter(
            actionsession_id__in=actionsessionqueryset, toplatetype=platetype
        )
        if extractactionqueryset:
            plateroundedvolumes += self.getRoundedExtractActionVolumes(
                extractactionqueryset=extractactionqueryset
            )
        labwareplatetype = self.getPlateType(
            platetype=platetype,
            volumes=plateroundedvolumes,
            wellsneeded=wellsneeded,
        )
        reaction_ids = actionsessionqueryset.values_list(
            "reaction_id", flat=True
        ).order_by("reaction_id")
        reactionqueryset = getReactionQuerySet(reaction_ids=reaction_ids)
        self.createPlateByReactionClassRecipe(
            reactionqueryset=reactionqueryset,
            labwareplatetype=labwareplatetype,
            platetype=platetype,
        )

    def createAnalysePlate(self, platetype: str):
        """Creates analyse plate/s for executing analyse actions"""
        actionsessionqueryset = self.getActionSessionByPlateType(platetype=platetype)
        wellsneeded = len(actionsessionqueryset)
        plateroundedvolumes = []
        addactionqueryset = AddAction.objects.filter(
            actionsession_id__in=actionsessionqueryset, toplatetype=platetype
        )
        if addactionqueryset:
            plateroundedvolumes += self.getRoundedAddActionVolumes(
                addactionqueryset=addactionqueryset
            )
        extractactionqueryset = ExtractAction.objects.filter(
            actionsession_id__in=actionsessionqueryset, toplatetype=platetype
        )
        if extractactionqueryset:
            plateroundedvolumes += self.getRoundedExtractActionVolumes(
                extractactionqueryset=extractactionqueryset
            )
        labwareplatetype = self.getPlateType(
            platetype=platetype,
            volumes=plateroundedvolumes,
            wellsneeded=wellsneeded,
        )
        reaction_ids = actionsessionqueryset.values_list(
            "reaction_id", flat=True
        ).order_by("reaction_id")
        reactionqueryset = getReactionQuerySet(reaction_ids=reaction_ids)
        self.createPlateByReactionClassRecipe(
            reactionqueryset=reactionqueryset,
            labwareplatetype=labwareplatetype,
            platetype=platetype,
        )

    def createSolventPlate(self, materialsdf: DataFrame):
        """Creates solvent plate/s for diluting reactants for reactions or analysis."""
        if not materialsdf.empty:
            solventdictslist = []
            materialsdf = materialsdf.groupby(["solvent"])["volume"].sum().to_frame()
            startinglabwareplatetype = self.getPlateType(
                platetype="solvent", volumes=materialsdf["volume"]
            )

            plateobj = self.createPlateModel(
                platetype="solvent",
                platename="Solventplate",
                labwaretype=startinglabwareplatetype,
            )

            maxwellvolume = self.getMaxWellVolume(plateobj=plateobj)
            deadvolume = self.getDeadVolume(maxwellvolume=maxwellvolume)
            for solventgroup in materialsdf.index.values:
                totalvolume = materialsdf.at[solventgroup, "volume"]
                if totalvolume > maxwellvolume:
                    nowellsneededratio = totalvolume / (maxwellvolume - deadvolume)
                    frac, whole = math.modf(nowellsneededratio)
                    volumestoadd = [maxwellvolume for i in range(int(whole))]
                    volumestoadd.append(frac * maxwellvolume + deadvolume)
                    for volumetoadd in volumestoadd:
                        indexwellavailable = self.getPlateWellIndexAvailable(
                            plateobj=plateobj
                        )
                        if type(indexwellavailable) == bool:
                            plateobj = self.createPlateModel(
                                platetype="solvent",
                                platename="Solventplate",
                                labwaretype=startinglabwareplatetype,
                            )
                            indexwellavailable = self.getPlateWellIndexAvailable(
                                plateobj=plateobj
                            )

                        wellobj = self.createWellModel(
                            plateobj=plateobj,
                            welltype="solvent",
                            wellindex=indexwellavailable,
                            volume=volumetoadd,
                            solvent=solventgroup,
                        )
                        self.updatePlateWellIndex(
                            plateobj=plateobj, wellindexupdate=indexwellavailable + 1
                        )
                        solventdictslist.append(
                            {
                                "name": plateobj.name,
                                "labware": plateobj.labware,
                                "well-index": wellobj.index,
                                "well-name": wellobj.name,
                                "solvent": solventgroup,
                                "amount-uL": volumetoadd,
                            }
                        )
                else:
                    indexwellavailable = self.getPlateWellIndexAvailable(
                        plateobj=plateobj
                    )
                    volumetoadd = totalvolume + deadvolume
                    if type(indexwellavailable) == bool:
                        plateobj = self.createPlateModel(
                            platetype="solvent",
                            platename="Solventplate",
                            labwaretype=startinglabwareplatetype,
                        )
                        indexwellavailable = self.getPlateWellIndexAvailable(
                            plateobj=plateobj
                        )
                    wellobj = self.createWellModel(
                        plateobj=plateobj,
                        welltype="solvent",
                        wellindex=indexwellavailable,
                        volume=volumetoadd,
                        solvent=solventgroup,
                    )
                    self.updatePlateWellIndex(
                        plateobj=plateobj, wellindexupdate=indexwellavailable + 1
                    )
                    solventdictslist.append(
                        {
                            "name": plateobj.name,
                            "labware": plateobj.labware,
                            "well-index": wellobj.index,
                            "well-name": wellobj.name,
                            "solvent": solventgroup,
                            "amount-uL": volumetoadd,
                        }
                    )
            solventdf = pd.DataFrame(solventdictslist)
            self.createSolventPrepModel(solventdf=solventdf)
