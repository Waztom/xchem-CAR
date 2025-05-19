"""
Volume Manager for OpenTrons protocols.

This module provides methods for managing well volumes, checking availability,
calculating aspirate heights, and handling dead volumes.
"""

import math
import logging

from backend.models import Well
from backend.opentrons.labwareavailable import labware_plates

logger = logging.getLogger(__name__)


class VolumeManager:
    """
    Manages volume calculations and well volume tracking.

    This class handles all volume-related functionality including checking
    available volumes, updating wells after transfers, and calculating
    aspirate heights.
    """

    def __init__(self, script_generator):
        """
        Initialize the volume manager.

        Parameters
        ----------
        script_generator : ScriptGenerator
            Parent script generator
        """
        self.script_generator = script_generator
        logger.info("VolumeManager initialized")

    def check_volumes_close(self, volume1: float, volume2: float) -> bool:
        """
        Check if two volumes are almost the same value.

        Parameters
        ----------
        volume1 : float
            First volume to compare
        volume2 : float
            Second volume to compare

        Returns
        -------
        bool
            True if volumes are close, False otherwise
        """
        check_close = math.isclose(volume1, volume2, rel_tol=0.001)
        return check_close

    def get_well_volume_available(self, wellobj: Well) -> float:
        """
        Get volume of well available for transfer.

        Parameters
        ----------
        wellobj : Well
            The well to check for available volume

        Returns
        -------
        float
            Available volume in the well (accounting for dead volume)
        """
        plate_id = wellobj.plate_id.id
        plate_name = wellobj.plate_id.name
        logger.info(
            f"Checking available volume in well {wellobj.index} of plate {plate_name} (ID: {plate_id})"
        )

        max_well_volume = self.get_max_well_volume(plateid=plate_id)
        dead_volume = self.get_dead_volume(maxwellvolume=max_well_volume)
        well_volume = wellobj.volume
        well_volume_available = well_volume - dead_volume

        logger.info(f"Well contains {well_volume} µL, dead volume is {dead_volume} µL")
        logger.info(f"Available volume for transfer: {well_volume_available} µL")

        self.update_well_available(
            wellvolumeavailable=well_volume_available, wellobj=wellobj
        )

        return well_volume_available

    def update_well_available(self, wellvolumeavailable: float, wellobj: Well):
        """
        Update well availability status based on available volume.

        Parameters
        ----------
        wellvolumeavailable : float
            Available volume in the well
        wellobj : Well
            Well object to update
        """
        if wellvolumeavailable < 0:
            logger.warning(
                f"Well {wellobj.index} in plate {wellobj.plate_id.name} has insufficient volume"
            )
            logger.info(f"Marking well {wellobj.index} as unavailable")
            wellobj.available = False
            wellobj.save()

    def update_well_volume(self, wellobj: Well, transfervolume: float):
        """
        Update the volume available in a well after a transfer.

        Parameters
        ---------
        wellobj : Well
            Well to update volume for
        transfervolume : float
            The volume transferred from the well
        """
        original_volume = wellobj.volume
        new_volume = original_volume - transfervolume

        logger.info(
            f"Updating well {wellobj.index} volume: {original_volume} µL → {new_volume} µL (transfer: {transfervolume} µL)"
        )

        wellobj.volume = new_volume
        wellobj.save()

        if new_volume <= 0:
            logger.warning(f"Well {wellobj.index} is now empty after transfer")

    def get_max_well_volume(self, plateid: int) -> float:
        """
        Get the maximum volume of a plate's wells.

        Parameters
        ----------
        plateid : int
            The plate ID to get the maximum well volume of

        Returns
        -------
        float
            The plate's maximum well volume
        """
        query_service = self.script_generator.query_service
        plateobj = query_service.get_plate_by_id(plateid)
        max_well_volume = plateobj.maxwellvolume
        logger.info(
            f"Maximum well volume for plate {plateobj.name}: {max_well_volume} µL"
        )
        return max_well_volume

    def get_dead_volume(self, maxwellvolume: float) -> float:
        """
        Calculate the dead volume (minimum volume that should remain in well).

        Parameters
        ----------
        maxwellvolume : float
            Maximum well volume

        Returns
        -------
        float
            Dead volume (typically 5% of maximum well volume)
        """
        dead_volume = maxwellvolume * 0.05
        logger.info(
            f"Dead volume calculated as {dead_volume} µL (5% of {maxwellvolume} µL)"
        )
        return dead_volume

    def calculate_aspirate_height(
        self, bottomlayervolume: float, labware: str
    ) -> float:
        """
        Calculate the height at which to extract the top layer based on
        volume occupying the bottom layer.

        Parameters
        ----------
        bottomlayervolume : float
            The volume occupying the bottom layer
        labware : str
            The type of labware being used for the extraction

        Returns
        -------
        float
            The height (mm) from the bottom of the plate,
            the pipette tip should aspirate from
        """
        logger.info(
            f"Calculating aspirate height for {labware} with bottom layer volume {bottomlayervolume} µL"
        )

        try:
            aspirateheightconvesion_m = labware_plates[labware][
                "aspirateheightconversion-m"
            ]
            aspirateheightconvesion_c = labware_plates[labware][
                "aspirateheightconversion-c"
            ]

            bottomlayerheight = (
                aspirateheightconvesion_m * bottomlayervolume
            ) + aspirateheightconvesion_c

            aspirateheight = bottomlayerheight * 1.15

            logger.info(f"Calculated aspirate height: {aspirateheight} mm")
            return aspirateheight

        except KeyError:
            logger.error(
                f"Labware type '{labware}' not found in labware_plates dictionary"
            )
            logger.warning(f"Using default aspirate height calculation")
            # Fallback calculation if labware not found
            aspirateheight = (bottomlayervolume / 200) + 0.5
            return aspirateheight

    def get_optimal_transfer_volume(
        self, available_volume: float, requested_volume: float
    ) -> float:
        """
        Calculate optimal transfer volume based on available and requested volumes.

        Parameters
        ----------
        available_volume : float
            Volume available in the source well
        requested_volume : float
            Volume requested for transfer

        Returns
        -------
        float
            Volume that can actually be transferred
        """
        logger.info(
            f"Determining optimal transfer volume: requested={requested_volume} µL, available={available_volume} µL"
        )

        if available_volume >= requested_volume:
            logger.info(
                f"Sufficient volume available, transferring requested {requested_volume} µL"
            )
            return requested_volume

        # If not enough volume, return the available amount
        optimal_volume = max(0, available_volume)
        logger.warning(
            f"Insufficient volume! Requested: {requested_volume} µL, available: {available_volume} µL"
        )
        logger.info(
            f"Will transfer {optimal_volume} µL instead of requested {requested_volume} µL"
        )
        return optimal_volume

    def update_well_reactant_status(self, wellobj: Well, is_for_next_step: bool):
        """
        Update well reactant status for next step.

        Parameters
        ----------
        wellobj : Well
            The well to update
        is_for_next_step : bool
            Whether the well contains a reactant for the next step
        """
        status_text = "available" if is_for_next_step else "not available"
        logger.info(
            f"Updating well {wellobj.index} in plate {wellobj.plate_id.name}: marking as {status_text} for next step"
        )

        wellobj.reactantfornextstep = is_for_next_step
        wellobj.save()

    def update_column_reactant_status(self, columnobj, is_for_next_step: bool):
        """
        Update all wells in a column for reactant status.

        Parameters
        ----------
        columnobj : Column
            The column containing wells to update
        is_for_next_step : bool
            Whether the wells contain reactants for the next step
        """
        from backend.models import Well

        status_text = "available" if is_for_next_step else "not available"
        logger.info(
            f"Updating all wells in column {columnobj.index}: marking as {status_text} for next step"
        )

        wellqueryset = Well.objects.filter(column_id=columnobj)
        well_count = wellqueryset.count()

        if well_count:
            logger.info(
                f"Found {well_count} wells to update in column {columnobj.index}"
            )
            for well in wellqueryset:
                well.reactantfornextstep = is_for_next_step
                well.save()
        else:
            logger.warning(f"No wells found in column {columnobj.index}")
