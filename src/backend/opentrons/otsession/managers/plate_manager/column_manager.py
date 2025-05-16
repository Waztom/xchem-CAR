import logging
from django.db.models import QuerySet

from .....models import Plate, Column

logger = logging.getLogger(__name__)

class ColumnManager:
    """
    Class to manage columns in a plate.
    """

    def __init__(self, session):
        """
        Initialize with a reference to the parent session.

        Parameters
        ----------
        session: BaseSession
            The parent session object
        """
        self.session = session

    def get_plate_columns(self, plate_obj: Plate) -> QuerySet[Column]:
        """
        Retrieves the columns for a plate.

        Parameters
        ----------
        plate_obj: Plate
            The plate to get all the related columns for

        Returns
        -------
        column_queryset: QuerySet[Column]
            The plate's columns
        """
        column_queryset = Column.objects.filter(plate_id=plate_obj.id)
        return column_queryset
    
    def get_plate_current_column_index(self, plate_obj: Plate) -> int:
        """
        Check if any columns available on a plate.

        Parameters
        ----------
        plate_obj: Plate
            The plate to search for a column available

        Returns
        -------
        index_column_available: int or False
            The index of the column available on a plate, or False if none available
        """
        index_column_available = plate_obj.indexcolumnavailable
        number_columns = plate_obj.numbercolumns

        if index_column_available + 1 <= number_columns:
            return index_column_available
        else:
            return False
        
    def update_column_ot_session_ids(
        self, column_queryset: QuerySet[Column], plate_obj: Plate
    ):
        """
        Updates columns to link to current OT session.

        Parameters
        ----------
        column_queryset: QuerySet[Column]
            The columns to be updated
        plate_obj: Plate
            The plate object related to the updated columns
        """
        for column_obj in column_queryset:
            column_obj.plate_id = plate_obj
            column_obj.otsession_id = self.session.otsessionobj
            column_obj.save()

    def update_plate_column_index_available(
        self, plate_obj: Plate, columnindexupdate: int
    ):
        """
        Updates the column index used.

        Parameters
        ----------
        plate_obj: Plate
            The plate to update the next available column index
        columnindexupdate: int
            The column index to update on the plate
        """
        plate_obj.indexcolumnavailable = columnindexupdate
        plate_obj.save()
    
    def create_column_model(
        self, plate_obj: Plate, columnindex: int, columntype: str, reactionclass: str
    ) -> Column:
        """
        Creates a column object.

        Parameters
        ----------
        plate_obj: Plate
            The plate that the column is linked to
        columnindex: int
            The index of the column in the plate
        columntype: str
            The type of plate the column is used on
        reactionclass: str
            The reaction class occupying the column

        Returns
        -------
        column_obj: Column
            The created column object
        """
        column_obj = Column()
        column_obj.otsession_id = self.session.otsessionobj
        column_obj.plate_id = plate_obj
        column_obj.index = columnindex
        column_obj.type = columntype
        column_obj.reactionclass = reactionclass
        column_obj.save()

        return column_obj
