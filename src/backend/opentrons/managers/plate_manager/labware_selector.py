import logging
import math

from ....models import Plate
from ...labwareavailable import labware_plates

logger = logging.getLogger(__name__)

class LabwareSelector :
    """Class to select labware based on the given criteria."""

    def __init__(self, session):
        """
        Initialize with a reference to the parent session.

        Parameters
        ----------
        session: BaseSession
            The parent session object
        """
        self.session = session


    def get_plate_type(
        self,
        platetype: str,
        volumes: list,
        temperature: int = 25,
        wellsneeded: int = None,
    ):
        """
        Gets best suited plate based on several criteria including temperature,
        volume, and number of wells needed.

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

        Returns
        -------
        labware_plate_type: str
            The best suited labware type for the requirements
        """
        platetype = platetype.lower()
        median_volume = self.session.data_manager.get_median_value(values=volumes)

        possible_labware_plate_types = [
            labware_plate
            for labware_plate in labware_plates
            if platetype in labware_plates[labware_plate]["type"]
            and labware_plates[labware_plate]["max_temp"] >= temperature
        ]

        if not possible_labware_plate_types:
            logger.warning(
                f"No labware types found for {platetype} at {temperature}°C for type {platetype} and volumes {platetype}"
            )
            return None 
        
        vial_compare_dict = {}

        for labware_plate in possible_labware_plate_types:
            max_temp = labware_plates[labware_plate]["max_temp"]
            max_volume_vial = labware_plates[labware_plate]["volume_well"]
            no_plate_vials = labware_plates[labware_plate]["no_wells"]

            if not wellsneeded:
                wellsneeded = sum(
                    [
                        self.session.material_manager.get_number_vials(
                            max_volume_vial=max_volume_vial, volume_material=volume
                        )
                        for volume in volumes
                    ]
                )

            no_plates_needed = int(math.ceil(wellsneeded / no_plate_vials))
            volume_difference = max_volume_vial - median_volume
            temp_difference = max_temp - temperature

            logger.debug(
                f"Labware: {labware_plate}, Max Temp: {max_temp}, "
                f"Max Volume: {max_volume_vial}, "
                f"Volume Difference: {volume_difference}, "
                f"Temperature Difference: {temp_difference}"
            )
            if platetype == "reaction":  # Reaction plates can only fill one well
                max_volume_exceeded_test = all(
                    [False if max_volume_vial - vol <= 0 else True for vol in volumes]
                )
                if (
                    volume_difference < 0
                    or temp_difference < 0
                    or not max_volume_exceeded_test
                ):
                    continue

            if platetype == "starting":  # Starting plates can fill multiple wells
                if volume_difference < 0 or temp_difference < 0:
                    continue

            vial_compare_dict[labware_plate] = {
                "noplatesneeded": no_plates_needed,
                "volumedifference": volume_difference,
                "novialsneeded": wellsneeded,
                "tempdifference": temp_difference,
            }
    
        logger.debug(f"Vial compare dict: {vial_compare_dict}")
        
        # Add check for empty dictionary
        if not vial_compare_dict:
            logger.warning(f"No suitable labware found for {platetype}-plate with volumes {volumes} at {temperature}°C")
            if platetype == "starting":
                # For starting materials, use the largest available plate that can handle the temperature
                fallback_plates = [
                    p for p in labware_plates 
                    if "starting" in labware_plates[p]["type"]
                    and labware_plates[p]["max_temp"] >= temperature
                ]
                if fallback_plates:
                    # Sort by volume capacity (descending)
                    fallback_plates.sort(
                        key=lambda x: labware_plates[x]["volume_well"], 
                        reverse=True
                    )
                    logger.info(f"Falling back to largest compatible plate: {fallback_plates[0]}")
                    return fallback_plates[0]
                else:
                    # Last resort fallback
                    logger.warning("No compatible labware found, defaulting to standard 96-well plate")
                    return "plateone_96_wellplate_2500ul"
            else:
                # For other plate types
                return "plateone_96_wellplate_2500ul"

        minimum_no_plates_needed = min(
            (d["noplatesneeded"] for d in vial_compare_dict.values())
        )
        minimum_temp_difference = min(
            (d["tempdifference"] for d in vial_compare_dict.values())
        )

        labware_plate_types = [
            labware_plate
            for labware_plate in vial_compare_dict
            if vial_compare_dict[labware_plate]["noplatesneeded"]
            == minimum_no_plates_needed
            and vial_compare_dict[labware_plate]["tempdifference"]
            == minimum_temp_difference
        ]

        if len(labware_plate_types) > 1:
            minimum_volume_difference = min(
                (vial_compare_dict[p]["volumedifference"] for p in labware_plate_types)
            )
            labware_plate_types = [
                labware_plate
                for labware_plate in labware_plate_types
                if vial_compare_dict[labware_plate]["volumedifference"]
                == minimum_volume_difference
            ]

            if len(labware_plate_types) > 1:
                minimum_no_vials_needed = min(
                    (vial_compare_dict[p]["novialsneeded"] for p in labware_plate_types)
                )
                labware_plate_types = [
                    labware_plate
                    for labware_plate in labware_plate_types
                    if vial_compare_dict[labware_plate]["novialsneeded"]
                    == minimum_no_vials_needed
                ]

        return labware_plate_types[0]

    def get_max_well_volume(self, plate_obj: Plate) -> float:
        """
        Get max well volume of a well plate.

        Parameters
        ----------
        plate_obj: Plate
            The plate to get the max well volume of

        Returns
        -------
        max_well_volume: float
            The maximum well volume of a well plate
        """
        max_well_volume = plate_obj.maxwellvolume
        return max_well_volume
    
    def get_dead_volume(self, max_well_volume: float) -> float:
        """
        Calculates the dead volume (5%) of a well.

        Parameters
        ----------
        max_well_volume: float
            The well's maximum volume

        Returns
        -------
        dead_volume: float
            The dead volume of the well
        """
        dead_volume = max_well_volume * 0.05
        return dead_volume

