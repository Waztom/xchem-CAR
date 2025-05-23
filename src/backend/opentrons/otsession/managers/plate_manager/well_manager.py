import logging
from django.db.models import QuerySet

from .....models import Plate, Well, Reaction, Column
from .....utils import wellIndexToWellName

logger = logging.getLogger(__name__)


class WellManager:
    """Manages wells in a plate."""

    def __init__(self, sesssion):
        """
        Initialize with a reference to the parent session.

        Parameters
        ----------
        session: BaseSession
            The parent session object
        """
        self.session = sesssion

    def get_plate_wells(self, plate_obj: Plate) -> QuerySet[Well]:
        """
        Retrieves the wells for a plate.

        Parameters
        ----------
        plate_obj: Plate
            The plate to get all the related wells for

        Returns
        -------
        well_queryset: QuerySet[Well]
            The plate's wells
        """
        well_queryset = Well.objects.filter(plate_id=plate_obj.id)
        return well_queryset

    def get_plate_well_index_available(self, plate_obj: Plate) -> int:
        """
        Check if any wells are available on a plate.

        Parameters
        ----------
        plate_obj: Plate
            The plate to search for an available well

        Returns
        -------
        index_well_available: int or False
            The index of the well available on a plate, or False if none available
        """
        index_well_available = plate_obj.indexswellavailable
        number_wells = plate_obj.numberwells

        if index_well_available + 1 <= number_wells:
            return index_well_available
        else:
            return False

    def update_plate_well_index(self, plate_obj: Plate, wellindexupdate: int):
        """
        Updates the plate's well index used.

        Parameters
        ----------
        plate_obj: Plate
            The plate to update the next available well index
        wellindexupdate: int
            The well index to update on the plate
        """
        plate_obj.indexswellavailable = wellindexupdate
        plate_obj.save()

    def update_well_ot_session_ids(
        self, well_queryset: QuerySet[Well], plate_obj: Plate
    ):
        """
        Updates wells to link to current OT session.

        Parameters
        ----------
        well_queryset: QuerySet[Well]
            The wells to be updated
        plate_obj: Plate
            The plate object related to the updated wells
        """
        for well_obj in well_queryset:
            well_obj.plate_id = plate_obj
            well_obj.otsession_id = self.session.otsessionobj
            well_obj.save()

    def create_well_model(
        self,
        plate_obj: Plate,
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
        """
        Creates a well object.

        Parameters
        ----------
        plate_obj: Plate
            The plate that the well is linked to
        welltype: str
            The well type eg. reaction, analyse
        wellindex: int
            The index of the well in the plate eg. 0, 1, 2, 3 etc
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
        well_obj: Well
            The created well object
        """
        well_obj = Well()
        well_obj.otsession_id = self.session.otsessionobj
        well_obj.plate_id = plate_obj

        if reactionobj:
            well_obj.reaction_id = reactionobj
            well_obj.method_id = reactionobj.method_id

        if columnobj:
            well_obj.column_id = columnobj

        well_obj.type = welltype
        well_obj.index = wellindex
        well_obj.name = wellIndexToWellName(
            wellindex=wellindex, platesize=plate_obj.numberwells
        )
        well_obj.volume = volume
        well_obj.smiles = smiles
        well_obj.concentration = concentration
        well_obj.solvent = solvent
        well_obj.reactantfornextstep = reactantfornextstep
        well_obj.save()

        return well_obj

    def check_index_well_is_new_column(self, plate_obj: Plate) -> bool:
        """
        Checks if current available well index on plate is the beginning
        of a new plate column.

        Parameters
        ----------
        plate_obj: Plate
            The plate containing the well to check

        Returns
        -------
        is_new_column: bool
            True if the well is at the start of a new plate column
        """
        index_well_available = plate_obj.indexswellavailable
        number_wells_in_column = plate_obj.numberwellsincolumn

        if index_well_available == 0:
            return True

        if index_well_available != 0:
            if (index_well_available % number_wells_in_column) == 0:
                return True
            if (index_well_available % number_wells_in_column) != 0:
                return False

    def get_new_column_and_well_index_available(self, plate_obj: Plate) -> tuple:
        """
        Checks if a new column is available and calculates the
        well index for the start of that column.

        Parameters
        ----------
        plate_obj: Plate
            The plate to search for a column available

        Returns
        -------
        tuple: (int, int) or False
            The index of the new column and the corresponding well index,
            or False if no new column is available
        """
        well_index_correction = plate_obj.numberwellsincolumn
        index_column_available = plate_obj.indexcolumnavailable

        if index_column_available + 1 <= plate_obj.numbercolumns:
            new_well_index = index_column_available * well_index_correction
            return (index_column_available, new_well_index)
        else:
            return False

    def get_max_well_volume(self, plate_obj):
        """
        Get the maximum well volume for a plate, adjusted for necessary headspace.

        Parameters
        ----------
        plate_obj : Plate
            The plate object to get the maximum well volume for

        Returns
        -------
        float
            The effective maximum well volume with headspace consideration
        """
        # Get the raw maximum volume
        max_volume = plate_obj.maxwellvolume

        # Get labware type
        labware_type = plate_obj.labware

        # Get max fill percentage from labware definition (default to 80% if not specified)
        from backend.opentrons.labwareavailable import labware_plates

        max_fill_percentage = labware_plates.get(labware_type, {}).get(
            "max_fill_percentage", 80
        )

        # Apply fill percentage limit
        effective_max_volume = max_volume * (max_fill_percentage / 100.0)

        logger.info(
            f"Using effective well volume of {effective_max_volume}µL ({max_fill_percentage}% of {max_volume}µL)"
        )
        return effective_max_volume
