"""
Well Finder for OpenTrons protocols.

This module provides methods for finding wells for various operations.
"""

import logging
from typing import List
from django.db.models import Q

from backend.models import Well, Plate
from backend.db_utils import get_previous_reaction_query_sets, get_product_smiles

logger = logging.getLogger(__name__)


class WellFinder:
    """
    Finds wells for various protocol operations.

    This class encapsulates the well-finding logic from the original OTWrite class
    in a more modular form.
    """

    def __init__(self, script_generator):
        """
        Initialize the well finder.

        Parameters
        ----------
        script_generator : ScriptGenerator
            Parent script generator
        """
        self.script_generator = script_generator
        logger.info("WellFinder initialized")

    def find_reaction_well(self, reaction_id: int, role: str, role_index: int = 1) -> Well:
        """
        Find the reaction plate well for a given reaction.

        Parameters
        ----------
        reaction_id : int
            The reaction's ID
        role : str
            The plate role (e.g., "reaction", "workup", "lcms")
        role_index : int
            The role index (default 1)

        Returns
        -------
        Well
            The well used in the reaction
        """
        logger.info(
            f"Finding reaction well for reaction ID {reaction_id}, role={role}, index={role_index}"
        )
        product_smiles = get_product_smiles(reaction_ids=[reaction_id])[0]
        logger.info(f"Looking for well with product SMILES: {product_smiles[:20]}...")

        try:
            well_obj = Well.objects.get(
                otsession_id=self.script_generator.otsession_id,
                reaction_id=reaction_id,
                role=role,
                role_index=role_index,
                smiles=product_smiles,
            )
            logger.info(
                f"Found well (ID: {well_obj.id}) for reaction {reaction_id} in plate {well_obj.plate_id.name}, well index {well_obj.index}"
            )
            return well_obj
        except Well.DoesNotExist:
            logger.error(
                f"Could not find well for reaction {reaction_id}, role={role}, index={role_index}"
            )
            logger.error(
                f"Searched with params: otsession_id={self.script_generator.otsession_id}, reaction_id={reaction_id}, role={role}, index={role_index}"
            )
            raise

    def find_solvent_wells(self, solvent: str, transfer_volume: float) -> List[List]:
        """
        Find solvent wells for diluting a previous reaction step's product.

        If solvent is available in wells, updates well volumes by removing
        transfer volume from available volume.

        Parameters
        ----------
        solvent : str
            The solvent needed for diluting the contents of a well
        transfer_volume : float
            The volume of solvent required for dilution

        Returns
        -------
        List[List]
            List of lists containing well objects and transfer volumes
        """
        logger.info(
            f"Finding solvent wells for {solvent}, volume needed: {transfer_volume} µL"
        )
        well_info = []
        volume_manager = self.script_generator.volume_manager

        try:
            solvent_plate_queryset = Plate.objects.filter(
                otsession_id=self.script_generator.otsession_id, role="solvent"
            )

            if not solvent_plate_queryset.exists():
                logger.warning(
                    f"No solvent plates found for session ID {self.script_generator.otsession_id}"
                )
                return well_info

            logger.info(f"Found {solvent_plate_queryset.count()} solvent plate(s)")
            well_objs = []

            for solvent_plate in solvent_plate_queryset:
                logger.info(f"Searching solvent plate: {solvent_plate.name}")
                well_queryset = solvent_plate.well_set.all().filter(
                    solvent=solvent,
                    available=True,
                    role="solvent",
                )

                if not well_queryset.exists():
                    logger.warning(
                        f"No wells with {solvent} found in plate {solvent_plate.name}"
                    )
                    continue

                logger.info(
                    f"Found {well_queryset.count()} wells with {solvent} in plate {solvent_plate.name}"
                )
                well_objs.extend(list(well_queryset))

            remaining_volume = transfer_volume
            for well_obj in well_objs:
                volumes_are_close = volume_manager.check_volumes_close(
                    volume1=remaining_volume, volume2=0.00
                )

                if volumes_are_close:
                    logger.info("Required volume has been fulfilled")
                    break

                well_volume_available = volume_manager.get_well_volume_available(
                    wellobj=well_obj
                )

                if well_volume_available > 0:
                    if well_volume_available >= remaining_volume:
                        logger.info(
                            f"Taking {remaining_volume} µL from well {well_obj.index} in plate {well_obj.plate_id.name}"
                        )
                        volume_manager.update_well_volume(
                            wellobj=well_obj, transfervolume=remaining_volume
                        )
                        well_info.append([well_obj, remaining_volume])
                        remaining_volume = 0.00

                    elif well_volume_available < remaining_volume:
                        logger.info(
                            f"Taking partial volume {well_volume_available} µL from well {well_obj.index} in plate {well_obj.plate_id.name}"
                        )
                        volume_manager.update_well_volume(
                            wellobj=well_obj, transfervolume=well_volume_available
                        )
                        well_info.append([well_obj, well_volume_available])
                        remaining_volume = remaining_volume - well_volume_available
                        logger.info(f"Still need {remaining_volume} µL")
                else:
                    logger.info(
                        f"Well {well_obj.index} has insufficient volume available"
                    )

        except Exception as e:
            logger.error(f"Error finding solvent plate well: {e}")

        if not well_info:
            logger.warning(f"No solvent well info found for solvent: {solvent}")
        else:
            logger.info(f"Found {len(well_info)} solvent well(s) for transfer")

        return well_info

    def find_starting_material_wells(
        self,
        reaction_step_no: int,
        reaction_id: int,
        smiles: str,
        solvent: str,
        concentration: float,
        transfer_volume: float,
        transfer_type: str = None,
    ) -> List[List]:
        """
        Find starting plate wells for executing an add action.

        Parameters
        ----------
        reaction_step_no : int
            The reaction step number
        reaction_id : int
            The reaction's ID used in the transfer
        smiles : str
            The SMILES of the starting material needed to be transferred
        solvent : str
            The solvent used to prepare the starting material
        concentration : float
            The concentration of the starting material
        transfer_volume : float
            The volume of starting material needed for the transfer
        transfer_type : str, optional
            When provided, restrict the search to wells with this
            ``transfer_type`` (e.g. ``'single'`` or ``'multichannel'``).
            ``None`` means no filtering (default).

        Returns
        -------
        List[List]
            List of lists containing previous reaction querysets, well objects and transfer volumes
        """
        logger.info(
            f"Finding starting material wells for reaction {reaction_id}, step {reaction_step_no}"
        )
        logger.info(f"Material: SMILES={smiles[:20]}..., volume={transfer_volume} µL")
        if solvent:
            logger.info(f"Solvent: {solvent}, concentration: {concentration}")

        previous_reaction_queryset = get_previous_reaction_query_sets(
            reaction_id=reaction_id, smiles=smiles
        )

        if previous_reaction_queryset:
            logger.info(
                f"Found {len(previous_reaction_queryset)} previous reaction(s) for this material"
            )

        well_info = []
        volume_manager = self.script_generator.volume_manager

        try:
            # First try to find wells from starting materials
            if not concentration and not solvent:
                logger.info(
                    "Searching for starting material wells without solvent/concentration constraints"
                )
                filter_kwargs = dict(
                    otsession_id=self.script_generator.otsession_id,
                    smiles=smiles,
                    available=True,
                    role="startingmaterial",
                )
                if transfer_type:
                    filter_kwargs["transfer_type"] = transfer_type
                well_objects = Well.objects.filter(
                    **filter_kwargs
                ).order_by("id")
            else:
                logger.info(
                    f"Searching for starting material wells with solvent={solvent}, concentration={concentration}"
                )
                filter_kwargs = dict(
                    otsession_id=self.script_generator.otsession_id,
                    smiles=smiles,
                    solvent=solvent,
                    concentration=concentration,
                    available=True,
                    role="startingmaterial",
                )
                if transfer_type:
                    filter_kwargs["transfer_type"] = transfer_type
                well_objects = Well.objects.filter(
                    **filter_kwargs
                ).order_by("id")

            if not well_objects.exists():
                logger.warning(
                    "No starting material wells found with initial search criteria"
                )
            else:
                logger.info(
                    f"Found {well_objects.count()} starting material well candidates"
                )

            self._process_well_objects(
                well_objects=well_objects,
                previous_reaction_queryset=previous_reaction_queryset,
                transfer_volume=transfer_volume,
                volume_manager=volume_manager,
                well_info=well_info,
            )

            # If this is a later reaction step and no wells were found, look in reaction/workup plates
            if reaction_step_no > 1 and not well_info:
                logger.info(
                    "This is a later reaction step - searching for wells in reaction/workup plates"
                )
                try:
                    criterion1 = Q(otsession_id=self.script_generator.otsession_id)
                    criterion2 = Q(reaction_id=reaction_id)
                    criterion3 = Q(reaction_id__in=previous_reaction_queryset)
                    criterion4 = Q(smiles=smiles)
                    criterion5 = Q(role__in=["reaction", "workup", "spefilter"])
                    criterion6 = Q(reactantfornextstep=True)

                    well_queryset = Well.objects.filter(
                        criterion1
                        & (criterion2 | criterion3)
                        & criterion3
                        & criterion4
                        & criterion5
                        & criterion6
                    )

                    if well_queryset.exists():
                        well_obj = well_queryset[0]
                        logger.info(
                            f"Found material in {well_obj.role} plate, well {well_obj.index}"
                        )
                        well_info.append(
                            [previous_reaction_queryset, well_obj, transfer_volume]
                        )
                    else:
                        logger.warning("No wells found in reaction/workup plates")
                except Exception as e:
                    logger.error(f"Error finding well in reaction/workup plates: {e}")

        except Exception as e:
            logger.error(f"Error finding starting material wells: {e}")

        if not well_info:
            logger.warning(f"No starting plate well info found for SMILES: {smiles}")

            starter_well_smiles_available = (
                Well.objects.values_list("smiles", flat=True)
                .filter(
                    otsession_id=self.script_generator.otsession_id,
                    role="startingmaterial",
                    available=True,
                )
                .distinct()
            )

            logger.warning(
                f"Available starting materials: {len(starter_well_smiles_available)} different compounds"
            )
        else:
            logger.info(
                f"Found {len(well_info)} well(s) for starting material transfer"
            )

        return well_info

    def _process_well_objects(
        self,
        well_objects,
        previous_reaction_queryset,
        transfer_volume,
        volume_manager,
        well_info,
    ):
        """Helper method to process well objects and update volume information."""
        remaining_volume = transfer_volume
        wells_processed = 0

        for well_obj in well_objects:
            volumes_are_close = volume_manager.check_volumes_close(
                volume1=remaining_volume, volume2=0.00
            )

            if volumes_are_close:
                logger.info("Required volume has been fulfilled")
                break

            wells_processed += 1
            well_volume_available = volume_manager.get_well_volume_available(
                wellobj=well_obj
            )

            if well_volume_available > 0:
                if well_volume_available >= remaining_volume:
                    logger.info(
                        f"Taking {remaining_volume} µL from well {well_obj.index} in plate {well_obj.plate_id.name}"
                    )
                    volume_manager.update_well_volume(
                        wellobj=well_obj, transfervolume=remaining_volume
                    )
                    well_info.append(
                        [previous_reaction_queryset, well_obj, remaining_volume]
                    )
                    remaining_volume = 0.00

                elif well_volume_available < remaining_volume:
                    logger.info(
                        f"Taking partial volume {well_volume_available} µL from well {well_obj.index} in plate {well_obj.plate_id.name}"
                    )
                    volume_manager.update_well_volume(
                        wellobj=well_obj, transfervolume=well_volume_available
                    )
                    well_info.append(
                        [previous_reaction_queryset, well_obj, well_volume_available]
                    )
                    remaining_volume = remaining_volume - well_volume_available
                    logger.info(f"Still need {remaining_volume} µL")
            else:
                logger.info(f"Well {well_obj.index} has insufficient volume available")

        if wells_processed > 0:
            logger.info(f"Processed {wells_processed} well(s) for transfer")
