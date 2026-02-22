"""
Query Service for OpenTrons protocols.

This module provides methods for querying database objects needed
for protocol generation.
"""

import logging
from typing import List, Any
from django.db.models import QuerySet, Q

from backend.models import (
    ActionSession,
    AddAction,
    Column,
    ExtractAction,
    MixAction,
    Reaction,
    Pipette,
    TipRack,
    Plate,
    Well,
)
from backend.db_utils import (
    get_product_smiles,
    get_previous_reaction_query_sets,
    get_reaction_query_set,
)

logger = logging.getLogger(__name__)


class QueryService:
    """
    Provides database query functionality for the script generator.

    This class centralizes all database access needed for protocol generation.
    """

    def __init__(self, script_generator):
        """
        Initialize the query service.

        Parameters
        ----------
        script_generator : ScriptGenerator
            Parent script generator
        """
        self.script_generator = script_generator
        self.otsession_id = script_generator.otsession_id
        self.actionsession_ids = script_generator.actionsession_ids
        logger.info(f"QueryService initialized for OT session ID {self.otsession_id}")

    def get_action_session_query_set(self) -> QuerySet[ActionSession]:
        """Get action session queryset for actionsession_ids

        Returns
        -------
        QuerySet[ActionSession]
            The action sessions related to the provided IDs
        """
        logger.info(f"Querying action sessions for IDs: {self.actionsession_ids}")
        criterion = Q(id__in=self.actionsession_ids)

        actionsessionqueryset = ActionSession.objects.filter(criterion).order_by("id")
        count = actionsessionqueryset.count()
        logger.info(f"Found {count} action session(s)")
        return actionsessionqueryset

    def get_add_action_query_set(
        self,
        reaction_ids: List[int],
        actionsessiontype: str = None,
        actionnumber: int = None,
    ) -> QuerySet[AddAction]:
        """Get add actions queryset for reaction_ids

        Parameters
        ----------
        reaction_ids: List[int]
            The reactions to search for related add actions
        actionsessiontype: str
            The optional action session type to look for add actions for
        actionnumber: int
            The optional number the action is executed in a process

        Returns
        -------
        QuerySet[AddAction]
            The add actions related to the reaction
        """
        if isinstance(reaction_ids, QuerySet):
            reaction_id_list = list(reaction_ids.values_list("id", flat=True))
        else:
            reaction_id_list = reaction_ids

        logger.info(f"Querying add actions for {len(reaction_id_list)} reaction(s)")
        if actionsessiontype:
            logger.info(f"Filtering for action session type: {actionsessiontype}")
        if actionnumber:
            logger.info(f"Filtering for action number: {actionnumber}")

        if actionsessiontype and not actionnumber:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            addactionqueryset = AddAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")

            count = addactionqueryset.count()
            logger.info(
                f"Found {count} add action(s) for session type {actionsessiontype}"
            )
            return addactionqueryset

        if actionnumber and actionsessiontype:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            criterion3 = Q(number=actionnumber)
            addactionqueryset = AddAction.objects.filter(
                criterion1 & criterion2 & criterion3
            ).order_by("id")

            count = addactionqueryset.count()
            logger.info(
                f"Found {count} add action(s) for session type {actionsessiontype}, number {actionnumber}"
            )
            return addactionqueryset

        # Default case if no conditions are met
        criterion = Q(reaction_id__in=reaction_ids)
        addactionqueryset = AddAction.objects.filter(criterion).order_by("id")

        count = addactionqueryset.count()
        logger.info(f"Found {count} add action(s) with no additional filters")
        return addactionqueryset

    def get_extract_action_query_set(
        self,
        reaction_ids: List[int],
        actionsessiontype: str = None,
        actionnumber: int = None,
    ) -> QuerySet[ExtractAction]:
        """Get extract actions queryset for reaction_ids

        Parameters
        ----------
        reaction_ids: List[int]
            The reactions to search for related extract actions
        actionsessiontype: str
            The optional action session type to look for extract actions for
        actionnumber: int
            The optional number the action is executed in a process

        Returns
        -------
        QuerySet[ExtractAction]
            The extract actions related to the reaction
        """
        if isinstance(reaction_ids, QuerySet):
            reaction_id_list = list(reaction_ids.values_list("id", flat=True))
        else:
            reaction_id_list = reaction_ids

        logger.info(f"Querying extract actions for {len(reaction_id_list)} reaction(s)")
        if actionsessiontype:
            logger.info(f"Filtering for action session type: {actionsessiontype}")
        if actionnumber:
            logger.info(f"Filtering for action number: {actionnumber}")

        if actionsessiontype and not actionnumber:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            extractactionqueryset = ExtractAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")

            count = extractactionqueryset.count()
            logger.info(
                f"Found {count} extract action(s) for session type {actionsessiontype}"
            )
            return extractactionqueryset

        if actionnumber and actionsessiontype:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            criterion3 = Q(number=actionnumber)
            extractactionqueryset = ExtractAction.objects.filter(
                criterion1 & criterion2 & criterion3
            ).order_by("id")

            count = extractactionqueryset.count()
            logger.info(
                f"Found {count} extract action(s) for session type {actionsessiontype}, number {actionnumber}"
            )
            return extractactionqueryset

        # Default case if no conditions are met
        criterion = Q(reaction_id__in=reaction_ids)
        extractactionqueryset = ExtractAction.objects.filter(criterion).order_by("id")

        count = extractactionqueryset.count()
        logger.info(f"Found {count} extract action(s) with no additional filters")
        return extractactionqueryset

    def get_mix_action_query_set(
        self,
        reaction_ids: List[int],
        actionsession_ids: List[int] = None,
        actionsessiontype: str = None,
        actionnumber: int = None,
    ) -> QuerySet[MixAction]:
        """Get mix actions queryset for reaction_ids

        Parameters
        ----------
        reaction_ids: List[int]
            The reactions to search for related mix actions
        actionsession_ids: List[int]
            Optional action session ids to match mix actions with
        actionsessiontype: str
            The optional action session type to look for extract actions for
        actionnumber: int
            The optional number the action is executed in a process

        Returns
        -------
        QuerySet[MixAction]
            The mix actions related to the reaction
        """
        if isinstance(reaction_ids, QuerySet):
            reaction_id_list = list(reaction_ids.values_list("id", flat=True))
        else:
            reaction_id_list = reaction_ids

        logger.info(f"Querying mix actions for {len(reaction_id_list)} reaction(s)")
        if actionsession_ids:
            logger.info(f"Filtering for {len(actionsession_ids)} action session ID(s)")
        if actionsessiontype:
            logger.info(f"Filtering for action session type: {actionsessiontype}")
        if actionnumber:
            logger.info(f"Filtering for action number: {actionnumber}")

        if actionsession_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__in=actionsession_ids)
            mixactionqueryset = MixAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")

            count = mixactionqueryset.count()
            logger.info(f"Found {count} mix action(s) for specified action sessions")
            return mixactionqueryset

        if actionsessiontype and not actionnumber:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            mixactionqueryset = MixAction.objects.filter(
                criterion1 & criterion2
            ).order_by("id")

            count = mixactionqueryset.count()
            logger.info(
                f"Found {count} mix action(s) for session type {actionsessiontype}"
            )
            return mixactionqueryset

        if actionsessiontype and actionnumber:
            criterion1 = Q(reaction_id__in=reaction_ids)
            criterion2 = Q(actionsession_id__type=actionsessiontype)
            criterion3 = Q(number=actionnumber)
            mixactionqueryset = MixAction.objects.filter(
                criterion1 & criterion2 & criterion3
            ).order_by("id")

            count = mixactionqueryset.count()
            logger.info(
                f"Found {count} mix action(s) for session type {actionsessiontype}, number {actionnumber}"
            )
            return mixactionqueryset

        # Default case if no conditions are met
        criterion = Q(reaction_id__in=reaction_ids)
        mixactionqueryset = MixAction.objects.filter(criterion).order_by("id")

        count = mixactionqueryset.count()
        logger.info(f"Found {count} mix action(s) with no additional filters")
        return mixactionqueryset

    def get_plates(self) -> QuerySet[Plate]:
        """Gets plates for an OT session

        Returns
        -------
        QuerySet[Plate]
            The plates linked to the OT session
        """
        logger.info(f"Querying plates for OT session {self.otsession_id}")
        platequeryset = Plate.objects.filter(otsession_id=self.otsession_id).order_by(
            "id"
        )

        count = platequeryset.count()
        if count == 0:
            logger.warning(f"No plates found for OT session {self.otsession_id}")
        else:
            plate_names = ", ".join([plate.name for plate in platequeryset])
            logger.info(f"Found {count} plate(s): {plate_names}")

        return platequeryset

    def get_plate_by_id(self, plateid: int) -> Plate:
        """Gets the plate object

        Parameters
        ----------
        plateid: int
            The id of the plate to search for

        Returns
        -------
        Plate
            The plate object
        """
        logger.info(f"Retrieving plate with ID {plateid}")
        plate = Plate.objects.filter(id=plateid).first()

        if plate:
            logger.info(f"Found plate: {plate.name} (role: {plate.role}, role_index: {plate.role_index})")
        else:
            logger.warning(f"No plate found with ID {plateid}")

        return plate

    def get_tip_racks(self) -> QuerySet[TipRack]:
        """Get the tip racks for an OT Session

        Returns
        -------
        QuerySet[TipRack]
            The tipracks linked to an OT Session
        """
        logger.info(f"Querying tip racks for OT session {self.otsession_id}")
        tipracksqueryset = TipRack.objects.filter(
            otsession_id=self.otsession_id
        ).order_by("id")

        count = tipracksqueryset.count()
        if count == 0:
            logger.warning(f"No tip racks found for OT session {self.otsession_id}")
        else:
            tiprack_info = [f"{tr.name} (slot {tr.index})" for tr in tipracksqueryset]
            logger.info(f"Found {count} tip rack(s): {', '.join(tiprack_info)}")

        return tipracksqueryset

    def get_pipette(self) -> Pipette:
        """Get the pipette for an OT session

        Returns
        -------
        Pipette
            The pipette used for an OT session
        """
        logger.info(f"Retrieving pipette for OT session {self.otsession_id}")
        try:
            pipetteobj = Pipette.objects.get(otsession_id=self.otsession_id)
            logger.info(
                f"Found pipette: {pipetteobj.name} ({pipetteobj.type}) on mount {pipetteobj.position}"
            )
            return pipetteobj
        except Pipette.DoesNotExist:
            logger.error(f"No pipette found for OT session {self.otsession_id}")
            raise
        except Pipette.MultipleObjectsReturned:
            logger.warning(
                f"Multiple pipettes found for OT session {self.otsession_id}, using first (single) one"
            )
            return Pipette.objects.filter(
                otsession_id=self.otsession_id, type="single"
            ).first() or Pipette.objects.filter(otsession_id=self.otsession_id).first()

    def get_multichannel_pipette(self):
        """Get the multi-channel pipette for an OT session, if one exists.

        Returns
        -------
        Pipette or None
            The multi-channel pipette, or None if only single-channel is set up.
        """
        try:
            pipetteobj = Pipette.objects.get(
                otsession_id=self.otsession_id, type="multi"
            )
            logger.info(
                f"Found multi-channel pipette: {pipetteobj.name} on mount {pipetteobj.position}"
            )
            return pipetteobj
        except Pipette.DoesNotExist:
            logger.info("No multi-channel pipette found for this session")
            return None
        except Pipette.MultipleObjectsReturned:
            return Pipette.objects.filter(
                otsession_id=self.otsession_id, type="multi"
            ).first()

    def get_multichannel_source_wells(self):
        """Get all multichannel-tagged starting material wells for this session.

        Returns
        -------
        QuerySet[Well]
            Wells with ``transfer_type='multichannel'`` ordered by index.
        """
        qs = Well.objects.filter(
            otsession_id=self.otsession_id,
            role="startingmaterial",
            transfer_type="multichannel",
            available=True,
        ).order_by("index")
        logger.info(f"Found {qs.count()} multichannel source well(s)")
        return qs

    def get_multichannel_add_actions(
        self,
        actionsession_queryset,
        smiles: str,
        solvent: str = None,
        concentration: float = None,
    ):
        """Get AddActions matching a multichannel source material.

        Parameters
        ----------
        actionsession_queryset : QuerySet
            The action sessions to search within.
        smiles : str
            SMILES of the material.
        solvent : str, optional
            Solvent filter.
        concentration : float, optional
            Concentration filter.

        Returns
        -------
        QuerySet[AddAction]
            Matching add actions, ordered by id.
        """
        qs = AddAction.objects.filter(
            actionsession_id__in=actionsession_queryset,
            smiles=smiles,
        )
        if solvent:
            qs = qs.filter(solvent=solvent)
        if concentration is not None:
            qs = qs.filter(concentration=concentration)
        qs = qs.order_by("id")
        logger.info(
            f"Found {qs.count()} AddAction(s) for MC material "
            f"{smiles[:30]}… (solvent={solvent}, conc={concentration})"
        )
        return qs

    def get_well_by_plate_and_index(self, plate_id, well_index: int):
        """Look up a single well by plate and positional index.

        Parameters
        ----------
        plate_id : int or Plate
            The plate (or its id) containing the well.
        well_index : int
            The well index on the plate.

        Returns
        -------
        Well or None
        """
        return Well.objects.filter(
            plate_id=plate_id, index=well_index,
        ).first()

    def get_column_query_set(
        self, role: str, role_index: int, reactionclass: str
    ) -> QuerySet[Column]:
        """Get column queryset for role/index and reactionclass

        Parameters
        ----------
        role: str
            The plate role (e.g., "reaction", "workup").
        role_index: int
            The plate role index.
        reactionclass: str
            The reaction class that occupy the columns

        Returns
        -------
        QuerySet[Column]
            The columns related to the role and reaction class
        """
        logger.info(
            f"Querying columns with role='{role}', index={role_index} for reaction class '{reactionclass}'"
        )
        criterion1 = Q(otsession_id=self.otsession_id)
        criterion2 = Q(role=role, role_index=role_index)
        criterion3 = Q(reactionclass=reactionclass)

        columnqueryset = Column.objects.filter(
            criterion1 & criterion2 & criterion3
        ).order_by("id")

        count = columnqueryset.count()
        if count == 0:
            logger.warning(
                f"No columns found for role='{role}', index={role_index} and reaction class '{reactionclass}'"
            )
        else:
            logger.info(
                f"Found {count} column(s) for role='{role}', index={role_index} and reaction class '{reactionclass}'"
            )

        return columnqueryset

    def get_well_by_reaction_id(
        self, reaction_id: int, role: str, role_index: int
    ) -> Well:
        """Find the reaction plate well

        Parameters
        ----------
        reaction_id: int
            The reaction's id linked to the well on the reaction plate
        role: str
            The plate role (e.g., "reaction", "workup").
        role_index: int
            The plate role index.

        Returns
        -------
        Well
            The well used in the reaction
        """
        logger.info(f"Finding well for reaction ID {reaction_id}, role='{role}', index={role_index}")
        try:
            productsmiles = get_product_smiles(reaction_ids=[reaction_id])[0]
            logger.info(
                f"Found product SMILES for reaction {reaction_id}: {productsmiles[:20]}..."
            )

            wellobj = Well.objects.get(
                otsession_id=self.otsession_id,
                reaction_id=reaction_id,
                role=role,
                role_index=role_index,
                smiles=productsmiles,
            )
            logger.info(
                f"Found well: plate {wellobj.plate_id.name}, well {wellobj.index}"
            )
            return wellobj

        except IndexError:
            logger.error(f"No product SMILES found for reaction {reaction_id}")
            raise
        except Well.DoesNotExist:
            logger.error(f"No well found for reaction {reaction_id}, role='{role}'")
            raise
        except Well.MultipleObjectsReturned:
            logger.warning(
                f"Multiple wells found for reaction {reaction_id}, role='{role}', using first one"
            )
            wells = Well.objects.filter(
                otsession_id=self.otsession_id,
                reaction_id=reaction_id,
                role=role,
                role_index=role_index,
                smiles=productsmiles,
            )
            return wells.first()

    def find_solvent_plate_well_obj(
        self, solvent: str, transfervolume: float
    ) -> List[List]:
        """Finds solvent well for diluting a previous reaction steps product. If
        solvent available in well, updates well volume by removing transfer volume from
        available volume

        Parameters
        ----------
        solvent: str
            The solvent needed for diluting the contents of a well
        transfervolume: float
            The volume of solvent required for dilution

        Returns
        -------
        List[List]
            List of lists containing well objects and transfer volumes
        """
        logger.info(
            f"Finding solvent wells for {solvent}, volume needed: {transfervolume} µL"
        )
        wellinfo = []

        try:
            solventplatequeryset = Plate.objects.filter(
                otsession_id=self.otsession_id, role="solvent"
            )

            plate_count = solventplatequeryset.count()
            if plate_count == 0:
                logger.warning(
                    f"No solvent plates found for OT session {self.otsession_id}"
                )
                return wellinfo

            logger.info(f"Found {plate_count} solvent plate(s)")

            if solventplatequeryset:
                wellobjs = []

                for solventplate in solventplatequeryset:
                    logger.info(f"Searching plate {solventplate.name} for {solvent}")
                    wellqueryset = solventplate.well_set.all().filter(
                        solvent=solvent,
                        available=True,
                        role="solvent",
                    )

                    well_count = wellqueryset.count()
                    if well_count == 0:
                        logger.info(
                            f"No wells with {solvent} found in plate {solventplate.name}"
                        )
                    else:
                        logger.info(
                            f"Found {well_count} wells with {solvent} in plate {solventplate.name}"
                        )
                        wellobjs.extend(list(wellqueryset))

                original_transfer_volume = transfervolume
                for wellobj in wellobjs:
                    volume_manager = self.script_generator.volume_manager

                    areclose = volume_manager.check_volumes_close(
                        volume1=transfervolume, volume2=0.00
                    )

                    if areclose:
                        logger.info("Required volume has been fulfilled")
                        break

                    wellvolumeavailable = volume_manager.get_well_volume_available(
                        wellobj=wellobj
                    )

                    if wellvolumeavailable > 0:
                        if wellvolumeavailable >= transfervolume:
                            logger.info(
                                f"Taking full volume {transfervolume} µL from well {wellobj.index}"
                            )
                            volume_manager.update_well_volume(
                                wellobj=wellobj, transfervolume=transfervolume
                            )
                            wellinfo.append([wellobj, transfervolume])
                            transfervolume = 0.00

                        if wellvolumeavailable < transfervolume:
                            logger.info(
                                f"Taking partial volume {wellvolumeavailable} µL from well {wellobj.index}"
                            )
                            volume_manager.update_well_volume(
                                wellobj=wellobj, transfervolume=wellvolumeavailable
                            )
                            wellinfo.append([wellobj, wellvolumeavailable])
                            transfervolume = transfervolume - wellvolumeavailable
                            logger.info(f"Still need {transfervolume} µL")
                    else:
                        logger.info(f"Well {wellobj.index} has insufficient volume")

                if transfervolume > 0:
                    logger.warning(
                        f"Could not find enough {solvent}: needed {original_transfer_volume} µL, found {original_transfer_volume - transfervolume} µL"
                    )

        except Exception as e:
            logger.error(f"Error finding solvent plate well: {e}")

        if not wellinfo:
            logger.warning(f"No solvent well info found for {solvent}!")
        else:
            logger.info(f"Found {len(wellinfo)} solvent well(s) for transfer")

        return wellinfo

    def find_starting_plate_well_obj(
        self,
        reaction_step_no: int,
        reaction_id: int,
        smiles: str,
        solvent: str,
        concentration: float,
        transfervolume: float,
    ) -> List[List]:
        """Finds starting plate well for executing an add action

        Parameters
        ----------
        reaction_step_no: int
            The reaction step number eg. Step 2
        reaction_id: int
            The reaction's id used in the transfer
        smiles: str
            The SMILES of the starting material needed to be transferred
        solvent: str
            The solvent used to prepare the starting material
        concentration: float
            The concentration of the starting material
        transfervolume: float
            The volume of starting material needed for the transfer

        Returns
        -------
        List[List]
            List of lists containing previous reaction querysets, well objects and transfer volumes
        """
        logger.info(
            f"Finding starting wells for reaction {reaction_id}, step {reaction_step_no}"
        )
        logger.info(
            f"Material: SMILES={smiles[:20]}..., transfer volume={transfervolume} µL"
        )
        if solvent:
            logger.info(f"Solvent: {solvent}, concentration: {concentration}")

        previousreactionqueryset = get_previous_reaction_query_sets(
            reaction_id=reaction_id, smiles=smiles
        )

        if previousreactionqueryset:
            prev_count = len(previousreactionqueryset)
            logger.info(f"Found {prev_count} previous reaction(s) for this material")

        wellinfo = []
        volume_manager = self.script_generator.volume_manager

        try:
            if reaction_step_no == 1:
                # Find wells for first reaction step
                logger.info("Processing first reaction step")
                if not concentration and not solvent:
                    logger.info(
                        "Searching for starting material wells without solvent/concentration constraints"
                    )
                    wellobjects = Well.objects.filter(
                        otsession_id=self.otsession_id,
                        smiles=smiles,
                        available=True,
                        role="startingmaterial",
                    ).order_by("id")
                else:
                    logger.info(
                        f"Searching for starting material wells with solvent={solvent}, concentration={concentration}"
                    )
                    wellobjects = Well.objects.filter(
                        otsession_id=self.otsession_id,
                        smiles=smiles,
                        solvent=solvent,
                        concentration=concentration,
                        available=True,
                        role="startingmaterial",
                    ).order_by("id")

                well_count = wellobjects.count()
                if well_count == 0:
                    logger.warning(
                        "No starting material wells found with specified criteria"
                    )
                else:
                    logger.info(f"Found {well_count} starting material well(s)")

                self._process_well_objects(
                    wellobjects,
                    previousreactionqueryset,
                    transfervolume,
                    volume_manager,
                    wellinfo,
                )

            elif reaction_step_no > 1:
                # Find wells for subsequent reaction steps
                logger.info(f"Processing reaction step {reaction_step_no}")

                if not concentration and not solvent:
                    logger.info(
                        "Searching for starting material wells without solvent/concentration constraints"
                    )
                    wellobjects = Well.objects.filter(
                        otsession_id=self.otsession_id,
                        smiles=smiles,
                        available=True,
                        role="startingmaterial",
                    ).order_by("id")
                else:
                    logger.info(
                        f"Searching for starting material wells with solvent={solvent}, concentration={concentration}"
                    )
                    wellobjects = Well.objects.filter(
                        otsession_id=self.otsession_id,
                        smiles=smiles,
                        solvent=solvent,
                        concentration=concentration,
                        available=True,
                        role="startingmaterial",
                    ).order_by("id")

                well_count = wellobjects.count()
                if well_count == 0:
                    logger.warning(
                        "No starting material wells found with specified criteria"
                    )
                else:
                    logger.info(f"Found {well_count} starting material well(s)")

                self._process_well_objects(
                    wellobjects,
                    previousreactionqueryset,
                    transfervolume,
                    volume_manager,
                    wellinfo,
                )

                # If no wells found, try finding from reaction and workup plates
                if not wellinfo:
                    logger.info(
                        "No starting material wells found, checking reaction and workup plates"
                    )
                    try:
                        roles_to_search = ["reaction", "workup", "spefilter"]
                        logger.info(f"Searching in plate roles: {', '.join(roles_to_search)}")

                        criterion1 = Q(otsession_id=self.otsession_id)
                        criterion2 = Q(reaction_id=reaction_id)
                        criterion3 = Q(reaction_id__in=previousreactionqueryset)
                        criterion4 = Q(smiles=smiles)
                        criterion5 = Q(role__in=["reaction", "workup", "spefilter"])
                        criterion6 = Q(reactantfornextstep=True)

                        wellqueryset = Well.objects.filter(
                            criterion1
                            & (criterion2 | criterion3)
                            & criterion3
                            & criterion4
                            & criterion5
                            & criterion6
                        )

                        well_count = wellqueryset.count()
                        if well_count == 0:
                            logger.warning("No wells found in reaction/workup plates")
                        else:
                            logger.info(
                                f"Found {well_count} well(s) in reaction/workup plates"
                            )
                            wellobj = wellqueryset[0]
                            logger.info(
                                f"Using well {wellobj.index} in plate {wellobj.plate_id.name}"
                            )
                            wellinfo.append(
                                [previousreactionqueryset, wellobj, transfervolume]
                            )
                    except Exception as e:
                        logger.error(
                            f"Error finding well in reaction/workup plates: {e}"
                        )

        except Exception as e:
            logger.error(f"Error finding starting plate well: {e}")

        if not wellinfo:
            logger.warning(
                f"No starting plate well info found for material with SMILES: {smiles[:20]}..."
            )

            starterwellsmilesavailable = (
                Well.objects.values_list("smiles", flat=True)
                .filter(
                    otsession_id=self.otsession_id,
                    role="startingmaterial",
                    available=True,
                )
                .distinct()
            )

            available_count = len(starterwellsmilesavailable)
            logger.warning(
                f"Available starting materials: {available_count} different compounds"
            )

            return None
        else:
            logger.info(f"Found {len(wellinfo)} well(s) for material transfer")

        return wellinfo

    def _process_well_objects(
        self,
        wellobjects,
        previousreactionqueryset,
        transfervolume,
        volume_manager,
        wellinfo,
    ):
        """Helper method to process well objects for the find_starting_plate_well_obj method"""
        logger.info(
            f"Processing {wellobjects.count()} well objects for transfer of {transfervolume} µL"
        )
        remaining_volume = transfervolume
        wells_processed = 0

        for wellobj in wellobjects:
            areclose = volume_manager.check_volumes_close(
                volume1=remaining_volume, volume2=0.00
            )

            if areclose:
                logger.info("Required volume has been fulfilled")
                break

            wells_processed += 1
            wellvolumeavailable = volume_manager.get_well_volume_available(
                wellobj=wellobj
            )

            if wellvolumeavailable > 0:
                plate_name = wellobj.plate_id.name if wellobj.plate_id else "unknown"

                if wellvolumeavailable >= remaining_volume:
                    logger.info(
                        f"Taking {remaining_volume} µL from well {wellobj.index} in plate {plate_name}"
                    )
                    volume_manager.update_well_volume(
                        wellobj=wellobj, transfervolume=remaining_volume
                    )
                    wellinfo.append(
                        [previousreactionqueryset, wellobj, remaining_volume]
                    )
                    remaining_volume = 0.00

                elif wellvolumeavailable < remaining_volume:
                    logger.info(
                        f"Taking partial volume {wellvolumeavailable} µL from well {wellobj.index} in plate {plate_name}"
                    )
                    volume_manager.update_well_volume(
                        wellobj=wellobj, transfervolume=wellvolumeavailable
                    )
                    wellinfo.append(
                        [previousreactionqueryset, wellobj, wellvolumeavailable]
                    )
                    remaining_volume = remaining_volume - wellvolumeavailable
                    logger.info(f"Still need {remaining_volume} µL")
            else:
                plate_name = wellobj.plate_id.name if wellobj.plate_id else "unknown"
                logger.info(
                    f"Well {wellobj.index} in plate {plate_name} has insufficient volume"
                )

        if wells_processed > 0:
            if remaining_volume > 0:
                logger.warning(
                    f"Only found {transfervolume - remaining_volume} µL out of {transfervolume} µL needed"
                )
            logger.info(f"Processed {wells_processed} well(s) for transfer")
        else:
            logger.warning("No wells processed for transfer")

    def get_unique_reaction_classes(
        self, reactionqueryset: QuerySet[Reaction]
    ) -> List[str]:
        """Set of unique reaction classes

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to get unique list of reaction classes for

        Returns
        ------
        List[str]
            The set of reactionclasses
        """
        logger.info("Getting unique reaction classes")
        reactionclasses = (
            reactionqueryset.values_list("reactionclass", flat=True)
            .order_by("reactionclass")
            .distinct()
        )

        class_count = len(reactionclasses)
        if class_count == 0:
            logger.warning("No reaction classes found")
        else:
            logger.info(
                f"Found {class_count} unique reaction class(es): {', '.join(reactionclasses)}"
            )

        return reactionclasses

    def get_unique_reaction_recipes(
        self, reactionclass: str, reactionqueryset: QuerySet[Reaction]
    ) -> List[str]:
        """Set of unique reaction recipes

        Parameters
        ----------
        reactionclass: str
            The reaction type to find linked recipes for
        reactionqueryset: QuerySet[Reaction]
            The reactions to get unique list of reaction recipes for

        Returns
        ------
        List[str]
            The set of reaction recipes
        """
        logger.info(f"Getting unique recipes for reaction class '{reactionclass}'")
        reactionrecipes = (
            reactionqueryset.filter(reactionclass=reactionclass)
            .values_list("recipe", flat=True)
            .order_by("recipe")
            .distinct()
        )

        recipe_count = len(reactionrecipes)
        if recipe_count == 0:
            logger.warning(f"No recipes found for reaction class '{reactionclass}'")
        else:
            logger.info(
                f"Found {recipe_count} unique recipe(s) for reaction class '{reactionclass}': {', '.join(reactionrecipes)}"
            )

        return reactionrecipes

    def get_grouped_reaction_by_class_recipe(
        self, reactionqueryset: QuerySet[Reaction]
    ) -> List[QuerySet[Reaction]]:
        """Group reactions by reaction class and recipe type

        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reactions to group by reaction class

        Returns
        -------
        List[QuerySet[Reaction]]
            The list of sublists of reaction querysets grouped by reaction class
        """
        logger.info("Grouping reactions by reaction class and recipe")
        reaction_count = reactionqueryset.count()
        logger.info(f"Processing {reaction_count} reactions")

        reactionclasses = self.get_unique_reaction_classes(
            reactionqueryset=reactionqueryset
        )
        groupedreactionquerysets = []

        for reactionclass in reactionclasses:
            logger.info(f"Processing reaction class: {reactionclass}")
            reactionclassqueryset = reactionqueryset.filter(reactionclass=reactionclass)
            class_count = reactionclassqueryset.count()

            if reactionclassqueryset:
                reactionrecipes = (
                    reactionclassqueryset.values_list("recipe", flat=True)
                    .distinct()
                    .order_by("recipe")
                )
                recipe_count = len(reactionrecipes)
                logger.info(
                    f"Found {recipe_count} unique recipe(s) for class '{reactionclass}'"
                )

                if len(reactionrecipes) == 1:
                    logger.info(
                        f"Adding {class_count} reactions with recipe '{reactionrecipes[0]}' as one group"
                    )
                    groupedreactionquerysets.append(reactionclassqueryset)

                if len(reactionrecipes) > 1:
                    notgroupbycolumnreactionqueryset = reactionclassqueryset.filter(
                        groupbycolumn=False
                    )
                    not_grouped_count = notgroupbycolumnreactionqueryset.count()

                    if notgroupbycolumnreactionqueryset:
                        logger.info(
                            f"Adding {not_grouped_count} non-column-grouped reactions as one group"
                        )
                        groupedreactionquerysets.append(
                            notgroupbycolumnreactionqueryset
                        )

                    groupbycolumnreactionqueryset = reactionclassqueryset.filter(
                        groupbycolumn=True
                    )
                    grouped_count = groupbycolumnreactionqueryset.count()

                    if groupbycolumnreactionqueryset:
                        logger.info(
                            f"Processing {grouped_count} column-grouped reactions by recipe"
                        )
                        for reactionrecipe in reactionrecipes:
                            reactionbyrecipequeryset = (
                                groupbycolumnreactionqueryset.filter(
                                    recipe=reactionrecipe,
                                    groupbycolumn=True,
                                )
                            )
                            recipe_group_count = reactionbyrecipequeryset.count()

                            if reactionbyrecipequeryset:
                                logger.info(
                                    f"Adding {recipe_group_count} reactions with recipe '{reactionrecipe}' as one group"
                                )
                                groupedreactionquerysets.append(
                                    reactionbyrecipequeryset
                                )

        group_count = len(groupedreactionquerysets)
        logger.info(f"Created {group_count} groups of reactions")
        return groupedreactionquerysets

    def get_next_obj_entries(self, queryset: QuerySet, obj: Any) -> QuerySet:
        """Finds all proceeding Django model object relative to the Django model
           object in a queryset

        Parameters
        ----------
        queryset: QuerySet
            The queryset to search for proceeding entries
        obj: Any
            The object that you want to find all proceeding object entries relative to

        Returns
        -------
        QuerySet
            The proceeding Django model objects as a queryset
        """
        logger.info(f"Finding proceeding entries after object with ID {obj.pk}")
        nextqueryset = queryset.filter(pk__gt=obj.pk).order_by("pk")

        count = nextqueryset.count()
        logger.info(f"Found {count} proceeding entries")
        return nextqueryset

    def check_next_reactions_add_actions(
        self, reactionobj: Reaction, productsmiles: str
    ) -> List[AddAction]:
        """Checks if there are any reaction objects following the reaction in a method.
           If there is, checks if any of the proceeding reaction add actions match
           the reaction product's SMILES

        Parameters
        ----------
        reactionobj: Reaction
            The Django reaction model object to search for it's product SMILES
            matching any add actions needing the product as a reactant in the
            proceeding reactions
        productsmiles: str
            The SMILES of the reaction's product

        Returns
        -------
        List[AddAction]
            The Django add action model objects that require the reaction product
            as an input reactant
        """
        logger.info(
            f"Checking if reaction {reactionobj.id} product is used in subsequent reactions"
        )
        logger.info(f"Product SMILES: {productsmiles[:20]}...")

        reactionqueryset = get_reaction_query_set(method_id=reactionobj.method_id.id)
        nextreactionqueryset = self.get_next_obj_entries(
            queryset=reactionqueryset, obj=reactionobj
        )

        addactionsmatches = []
        for next_reaction_obj in nextreactionqueryset:
            logger.info(
                f"Checking if reaction {next_reaction_obj.id} uses the product as input"
            )
            addactionmatch = self.get_add_action_query_set(
                reaction_ids=[next_reaction_obj.id],
                actionsessiontype="reaction",
            ).filter(smiles=productsmiles)

            if addactionmatch:
                logger.info(
                    f"Found add action match in reaction {next_reaction_obj.id}"
                )
                addactionsmatches.append(addactionmatch[0])

        match_count = len(addactionsmatches)
        if match_count > 0:
            logger.info(f"Product is used in {match_count} subsequent reaction(s)")
        else:
            logger.info("Product is not used in any subsequent reactions")

        return addactionsmatches
