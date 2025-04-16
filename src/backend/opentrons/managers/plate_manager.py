"""Manages plate and well creation and organization for OpenTrons sessions."""

import logging
import math
import pandas as pd
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.db.models import QuerySet, Q
from rdkit import Chem
from rdkit.Chem import Descriptors

from ...models import (
    Plate, Column, Well, Reaction, 
    ActionSession, AddAction, ExtractAction, CompoundOrder
)
from ...utils import (
    wellIndexToWellName, getInchiKey, getChemicalName, 
    canonSmiles, stripSalts, getReaction, getProduct
)
from ..labwareavailable import labware_plates

logger = logging.getLogger(__name__)

class PlateManager:
    """
    Manages plate and well creation, organization, and related operations.
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
    
    def create_plate_model(self, platetype: str, platename: str, labwaretype: str):
        """
        Creates Plate model if Deck index is available.
        
        Parameters
        ----------
        platetype: str
            The type of plate (e.g., reaction, startingmaterial)
        platename: str
            The name of the plate
        labwaretype: str
            The labware type (e.g., plateone_96_wellplate_2500ul)
            
        Returns
        -------
        plate_obj: Plate or None
            The created plate object, or None if no deck slots available
        """
        index_slot = self.session.deck_manager.check_deck_slot_available()
        if index_slot:
            plate_index = index_slot
            number_wells_in_column = labware_plates[labwaretype]["no_wells_in_column"]
            max_well_volume = labware_plates[labwaretype]["volume_well"]
            number_wells = labware_plates[labwaretype]["no_wells"]
            number_columns = labware_plates[labwaretype]["no_columns"]
            
            plate_obj = Plate()
            plate_obj.otbatchprotocol_id = self.session.otbatchprotocolobj
            plate_obj.otsession_id = self.session.otsessionobj
            plate_obj.deck_id = self.session.deckobj
            plate_obj.labware = labwaretype
            plate_obj.index = plate_index
            plate_obj.type = platetype
            plate_obj.maxwellvolume = max_well_volume
            plate_obj.numberwells = number_wells
            plate_obj.numberwellsincolumn = number_wells_in_column
            plate_obj.numbercolumns = number_columns
            plate_obj.save()
            
            plate_obj.name = f"Reaction_step_{self.session.reactionstep}_{platename}_{plate_obj.id}"
            plate_obj.save()
            
            return plate_obj
        else:
            logger.warning("CreatePlateModel - No more deck slots available")
            return None
    
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
    
    def update_plate_column_index_available(self, plate_obj: Plate, columnindexupdate: int):
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
    
    def update_well_ot_session_ids(self, well_queryset: QuerySet[Well], plate_obj: Plate):
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
    
    def update_column_ot_session_ids(self, column_queryset: QuerySet[Column], plate_obj: Plate):
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
    
    def update_plate_deck_ot_session_ids(self, plate_queryset):
        """
        Updates plates to link to the current Deck and OT session.
        
        Parameters
        ----------
        plate_queryset: QuerySet[Plate]
            The plates to update
        """
        for plate_obj in plate_queryset:
            index_slot = self.session.deck_manager.check_deck_slot_available()
            if index_slot:
                well_queryset = self.get_plate_wells(plate_obj=plate_obj)
                column_queryset = self.get_plate_columns(plate_obj=plate_obj)
                previous_type = plate_obj.type
                plate_obj.deck_id = self.session.deckobj
                plate_obj.otsession_id = self.session.otsessionobj
                if previous_type == "spefilter":
                    plate_obj.labware = "plateone_96_wellplate_2500ul"
                plate_obj.index = index_slot
                plate_obj.save()
                self.update_column_ot_session_ids(
                    column_queryset=column_queryset, plate_obj=plate_obj
                )
                self.update_well_ot_session_ids(
                    well_queryset=well_queryset, plate_obj=plate_obj
                )
            else:
                logger.warning("cloneInputPlate - No more deck slots available")
                
        # Make sure all plates are linked to the current session and deck
        for plate_obj in plate_queryset:
            plate_obj.deck_id = self.session.deckobj
            plate_obj.otsession_id = self.session.otsessionobj
            plate_obj.save()
            
        return plate_queryset
    
    def get_custom_starter_plates(self, custom_compound_orders):
        """
        Get plates associated with custom starter materials.
        
        Parameters
        ----------
        custom_compound_orders: QuerySet[CompoundOrder]
            The compound orders marked as custom starter plates
            
        Returns
        -------
        plates: QuerySet[Plate]
            The plates used for custom starting materials
        """
        # Get session IDs from the compound orders
        custom_session_ids = custom_compound_orders.values_list(
            "otsession_id", flat=True
        )
        
        # Get only custom starting material plates
        plates = Plate.objects.filter(
            otbatchprotocol_id=self.session.otbatchprotocolobj,
            type="startingmaterial",
            otsession_id__in=custom_session_ids,
        )
        
        return plates
    
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
    
    def create_column_model(self, plate_obj: Plate, columnindex: int, columntype: str, reactionclass: str) -> Column:
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
    
    def get_action_session_by_plate_type(self, platetype: str) -> QuerySet[ActionSession]:
        """
        Get action sessions that use a specific plate type.
        
        Parameters
        ----------
        platetype: str
            The plate type to look for in action sessions
            
        Returns
        -------
        action_session_queryset: QuerySet[ActionSession]
            Action sessions using the specified plate type
        """
        criterion1 = Q(id__in=self.session.actionsessionqueryset)
        criterion2 = Q(addaction__toplatetype=platetype)
        criterion3 = Q(extractaction__toplatetype=platetype)
        
        action_session_queryset = ActionSession.objects.filter(
            criterion1 & (criterion2 | criterion3)
        ).distinct()
        
        return action_session_queryset
    
    def get_unique_to_plates(self, action_session_queryset: QuerySet[ActionSession], platetypes: list) -> list:
        """
        Gets the distinct plate types needed for an action session queryset.
        
        Parameters
        ----------
        action_session_queryset: QuerySet[ActionSession]
            The action queryset to get the to plates for
        platetypes: list
            The plate types to try and find in an action session
            queryset eg. ["reaction", "workup1"]
            
        Returns
        -------
        to_plate_types: list
            The plate types needed for an action session
        """
        criterion1 = Q(actionsession_id__in=action_session_queryset)
        criterion2 = Q(toplatetype__in=platetypes)
        
        to_add_plates = (
            AddAction.objects.filter(criterion1 & criterion2)
            .values_list("toplatetype", flat=True)
            .distinct()
        )
        
        to_extract_plates = (
            ExtractAction.objects.filter(criterion1 & criterion2)
            .values_list("toplatetype", flat=True)
            .distinct()
        )
        
        to_plate_types = set(list(to_add_plates) + list(to_extract_plates))
        return list(to_plate_types)
    
    def get_plate_type(self, platetype: str, volumes: list, temperature: int = 25, wellsneeded: int = None):
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
            
        minimum_no_plates_needed = min(
            (d["noplatesneeded"] for d in vial_compare_dict.values())
        )
        minimum_temp_difference = min(
            (d["tempdifference"] for d in vial_compare_dict.values())
        )

        labware_plate_types = [
            labware_plate
            for labware_plate in vial_compare_dict
            if vial_compare_dict[labware_plate]["noplatesneeded"] == minimum_no_plates_needed
            and vial_compare_dict[labware_plate]["tempdifference"] == minimum_temp_difference
        ]
        
        if len(labware_plate_types) > 1:
            minimum_volume_difference = min(
                (vial_compare_dict[p]["volumedifference"] for p in labware_plate_types)
            )
            labware_plate_types = [
                labware_plate
                for labware_plate in labware_plate_types
                if vial_compare_dict[labware_plate]["volumedifference"] == minimum_volume_difference
            ]
            
            if len(labware_plate_types) > 1:
                minimum_no_vials_needed = min(
                    (vial_compare_dict[p]["novialsneeded"] for p in labware_plate_types)
                )
                labware_plate_types = [
                    labware_plate
                    for labware_plate in labware_plate_types
                    if vial_compare_dict[labware_plate]["novialsneeded"] == minimum_no_vials_needed
                ]

        return labware_plate_types[0]
    
    def get_all_ot_batch_protocol_plates(self, otbatchprotocol_id):
        """
        Get all plates used for an OT batch protocol.
        
        Parameters
        ----------
        otbatchprotocol_id: OTBatchProtocol
            The OT batch protocol to find all matching plates for
            
        Returns
        -------
        otbatchprotocol_plate_queryset: QuerySet[Plate]
            The plates used for all previous reaction and workup sessions
        """
        criterion1 = Q(otbatchprotocol_id=otbatchprotocol_id)
        criterion2 = Q(
            type__in=["reaction", "workup1", "workup2", "workup3", "spefilter"]
        )

        otbatchprotocol_plate_queryset = Plate.objects.filter(criterion1 & criterion2)
        return otbatchprotocol_plate_queryset
    
    def create_plate_by_reaction_class_recipe(self, reaction_queryset, platetype: str):
        """
        Creates plates by reaction class and recipe.
        
        Parameters
        ----------
        reaction_queryset: QuerySet[Reaction]
            The reactions to create plates for
        platetype: str
            The type of plate to create
            
        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Get grouped reactions by class and recipe
        grouped_reaction_querysets = self.session.data_manager.get_grouped_reaction_by_class_recipe(
            reactionqueryset=reaction_queryset
        )
        
        # Store created plates
        plate_objs = []
        
        # Create a plate for each group
        for group_key, group_reactions in grouped_reaction_querysets.items():
            # Extract class and recipe from group key
            reaction_class, recipe = group_key.split("-", 1)
            
            # Generate plate name
            plate_name = f"{platetype}-{reaction_class}-{recipe}"
            
            # Get volumes for the reaction group
            reaction_volumes = self.session.data_manager.get_rounded_reaction_volumes(
                reactionqueryset=group_reactions
            )
            
            # Determine best labware type
            labware_type = self.get_plate_type(
                platetype="reaction", 
                volumes=reaction_volumes, 
                temperature=25
            )
            
            # Create the plate
            plate_obj = self.create_plate_model(
                platetype=platetype,
                platename=plate_name,
                labwaretype=labware_type
            )
            
            if plate_obj:
                plate_objs.append(plate_obj)
                
                # Create a column for the reaction class-recipe
                column_index = self.get_plate_current_column_index(plate_obj=plate_obj)
                if column_index is not False:
                    column_obj = self.create_column_model(
                        plate_obj=plate_obj,
                        columnindex=column_index,
                        columntype=platetype,
                        reactionclass=reaction_class
                    )
                    
                    # Update column index
                    self.update_plate_column_index_available(
                        plate_obj=plate_obj, columnindexupdate=column_index + 1
                    )
                    
                    # Add wells for each reaction in this group
                    for reaction in group_reactions:
                        # Check if well index is available
                        well_index = self.get_plate_well_index_available(plate_obj=plate_obj)
                        if well_index is not False:
                            # Create well
                            well_obj = self.create_well_model(
                                plate_obj=plate_obj,
                                welltype=platetype,
                                wellindex=well_index,
                                reactionobj=reaction,
                                columnobj=column_obj
                            )
                            
                            # Update well index
                            self.update_plate_well_index(
                                plate_obj=plate_obj, wellindexupdate=well_index + 1
                            )
                        else:
                            logger.warning(f"No wells available on plate {plate_obj.name}")
                            break
            else:
                logger.warning(f"Could not create plate for {group_key}")
        
        return plate_objs

    def create_reaction_plate(self, platetype: str = "reaction"):
        """
        Creates a reaction plate.
        
        Parameters
        ----------
        platetype: str
            The type of plate to create, default is "reaction"
            
        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Get grouped reactions by temperature
        grouped_reaction_querysets = self.session.data_manager.get_grouped_temperature_reactions(
            reactionqueryset=self.session.groupreactionqueryset
        )
        
        # Store created plates
        plate_objs = []
        
        # Create reaction plates based on temperature groups
        for temp, reaction_queryset in grouped_reaction_querysets.items():
            # Create plates by reaction class and recipe
            temp_plates = self.create_plate_by_reaction_class_recipe(
                reaction_queryset=reaction_queryset,
                platetype=platetype
            )
            
            plate_objs.extend(temp_plates)
        
        return plate_objs

    def create_workup_plate(self, platetype: str):
        """
        Creates a workup plate.
        
        Parameters
        ----------
        platetype: str
            The type of workup plate (e.g., "workup1", "workup2")
            
        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Check if we need this plate type for the current actions
        action_session_queryset = self.get_action_session_by_plate_type(platetype=platetype)
        
        if not action_session_queryset.exists():
            logger.info(f"No action sessions require {platetype} plate")
            return []
        
        # Get volumes for determining plate type
        add_action_queryset = self.session.data_manager.get_add_action_query_set(
            reaction_ids=self.session.reaction_ids,
            actionsession_ids=self.session.actionsession_ids
        )
        
        rounded_volumes = self.session.data_manager.get_rounded_add_action_volumes(
            addactionqueryset=add_action_queryset
        )
        
        # Determine best labware type
        labware_type = self.get_plate_type(
            platetype="reaction",  # Use reaction type for labware compatibility
            volumes=rounded_volumes,
            temperature=25
        )
        
        # Create the plate
        plate_obj = self.create_plate_model(
            platetype=platetype,
            platename=f"{platetype}-plate",
            labwaretype=labware_type
        )
        
        if plate_obj:
            return [plate_obj]
        else:
            logger.warning(f"Could not create {platetype} plate")
            return []

    def create_analyse_plate(self):
        """
        Creates an analysis plate.
        
        Returns
        -------
        plate_objs: list
            List of created plate objects
        """
        # Get volumes for determining plate type
        extract_action_queryset = self.session.data_manager.get_extract_action_query_set(
            reaction_ids=self.session.reaction_ids,
            actionsession_ids=self.session.actionsession_ids
        )
        
        rounded_volumes = self.session.data_manager.get_rounded_extract_action_volumes(
            extractactionqueryset=extract_action_queryset
        )
        
        # Determine best labware type - analysis plates often use smaller volumes
        labware_type = self.get_plate_type(
            platetype="analyse",
            volumes=rounded_volumes,
            temperature=25
        )
        
        # Create the plate
        plate_obj = self.create_plate_model(
            platetype="analyse",
            platename="analyse-plate",
            labwaretype=labware_type
        )
        
        if plate_obj:
            return [plate_obj]
        else:
            logger.warning("Could not create analyse plate")
            return []

    def create_solvent_plate(self, materials_df: pd.DataFrame):
        """
        Creates a solvent plate. Used to create solvent plates
        for workup/cleanup solvents.
        
        Parameters
        ----------
        materials_df: DataFrame
            Dataframe containing material information
            
        Returns
        -------
        plate_obj: Plate
            The created solvent plate
        """
        if materials_df.empty:
            logger.info("No materials provided for solvent plate")
            return None
        
        # Group materials by solvent
        solvent_groups = materials_df.groupby("solvent")
        
        solvent_dicts_list = []
        
        for solvent_group, group_df in solvent_groups:
            # Calculate total volume needed for this solvent
            total_volume = group_df["volume"].sum() * 1.1  # Add 10% safety margin
            
            # Create or use existing plate based on volume
            labware_type = self.get_plate_type(
                platetype="starting",
                volumes=[total_volume],
                temperature=25
            )
            
            # Create a solvent plate (or reuse existing one)
            plate_obj = None
            existing_solvent_plates = Plate.objects.filter(
                otsession_id=self.session.otsessionobj,
                type="solvent"
            )
            
            if existing_solvent_plates.exists():
                for p in existing_solvent_plates:
                    if p.labware == labware_type:
                        plate_obj = p
                        break
            
            if not plate_obj:
                plate_obj = self.create_plate_model(
                    platetype="solvent",
                    platename=f"solvent-{solvent_group}",
                    labwaretype=labware_type
                )
            
            if not plate_obj:
                logger.warning(f"Could not create solvent plate for {solvent_group}")
                continue
            
            # Check for available well
            index_well_available = self.get_plate_well_index_available(plate_obj=plate_obj)
            if index_well_available is False:
                logger.warning(f"No wells available on solvent plate {plate_obj.name}")
                continue
            
            # Create well for this solvent
            well_obj = self.create_well_model(
                plate_obj=plate_obj,
                welltype="solvent",
                wellindex=index_well_available,
                volume=total_volume,
                solvent=solvent_group
            )
            
            # Update well index
            self.update_plate_well_index(
                plate_obj=plate_obj, wellindexupdate=index_well_available + 1
            )
            
            # Add to solvent prep list
            solvent_dicts_list.append({
                "name": plate_obj.name,
                "labware": plate_obj.labware,
                "well-index": well_obj.index,
                "well-name": well_obj.name,
                "solvent": solvent_group,
                "amount-uL": total_volume,
            })
        
        # Create solvent prep model if we have solvents
        if solvent_dicts_list:
            solvent_df = pd.DataFrame(solvent_dicts_list)
            self.session.material_manager.create_solvent_prep_model(solvent_df=solvent_df)
        
        return plate_obj

    def create_starting_material_plates_from_csv(self, csv_path: str):
        """
        Creates starting material plates from a CSV file.
        
        Parameters
        ----------
        csv_path: str
            Path to the CSV file containing starting material information
            
        Returns
        -------
        plate_obj: Plate
            The created starting material plate
        """
        try:
            # Read the CSV file
            materials_df = pd.read_csv(csv_path)
            
            if materials_df.empty:
                logger.warning("CSV file is empty")
                return None
            
            # Make sure required columns exist
            required_columns = ["SMILES", "concentration", "solvent", "amount-uL"]
            missing_columns = [col for col in required_columns if col not in materials_df.columns]
            
            if missing_columns:
                logger.error(f"CSV is missing required columns: {', '.join(missing_columns)}")
                return None
            
            # Get volumes for plate sizing
            volumes = materials_df["amount-uL"].tolist()
            
            # Determine best labware type
            labware_type = self.get_plate_type(
                platetype="starting",
                volumes=volumes,
                temperature=25
            )
            
            # Create the plate
            plate_obj = self.create_plate_model(
                platetype="startingmaterial",
                platename="custom-starting-materials",
                labwaretype=labware_type
            )
            
            if not plate_obj:
                logger.warning("Could not create custom starting material plate")
                return None
            
            # Add each material to the plate
            for _, row in materials_df.iterrows():
                # Check for available well
                index_well_available = self.get_plate_well_index_available(plate_obj=plate_obj)
                if index_well_available is False:
                    logger.warning(f"No wells available on custom starting plate {plate_obj.name}")
                    break
                
                # Create well for this material
                well_obj = self.create_well_model(
                    plate_obj=plate_obj,
                    welltype="startingmaterial",
                    wellindex=index_well_available,
                    volume=row["amount-uL"],
                    smiles=row["SMILES"],
                    concentration=row["concentration"],
                    solvent=row["solvent"]
                )
                
                # Update well index
                self.update_plate_well_index(
                    plate_obj=plate_obj, wellindexupdate=index_well_available + 1
                )
            
            # Create compound order record
            self.session.data_manager.create_compound_order_model(
                order_df=materials_df, is_custom_starter_plate=True
            )
            
            return plate_obj
            
        except Exception as e:
            logger.error(f"Error creating starting material plates from CSV: {str(e)}")
            return None

    def create_reaction_starting_plate(self):
        """
        Creates a starting material plate for reaction materials.
        
        Returns
        -------
        plate_obj: Plate
            The created starting material plate
        """
        # Get materials that need to be prepared
        materials_df = self.session.material_manager.get_add_actions_material_dataframe(
            product_exists=True
        )
        
        if materials_df.empty:
            logger.info("No materials needed for starting plate")
            return None
        
        # Get total volumes needed for plate sizing
        volumes = materials_df["volume"].tolist()
        
        # Determine best labware type
        labware_type = self.get_plate_type(
            platetype="starting",
            volumes=volumes,
            temperature=25
        )
        
        # Create the plate
        plate_obj = self.create_plate_model(
            platetype="startingmaterial",
            platename="reaction-starting-materials",
            labwaretype=labware_type
        )
        
        if not plate_obj:
            logger.warning("Could not create starting material plate")
            return None
        
        # Add each material to the plate
        for _, row in materials_df.iterrows():
            # Check for available well
            index_well_available = self.get_plate_well_index_available(plate_obj=plate_obj)
            if index_well_available is False:
                logger.warning(f"No wells available on starting plate {plate_obj.name}")
                break
            
            # Create well for this material
            well_obj = self.create_well_model(
                plate_obj=plate_obj,
                welltype="startingmaterial",
                wellindex=index_well_available,
                volume=row["volume"],
                smiles=row["smiles"],
                concentration=row["concentration"],
                solvent=row["solvent"]
            )
            
            # Update well index
            self.update_plate_well_index(
                plate_obj=plate_obj, wellindexupdate=index_well_available + 1
            )
        
        return plate_obj

    def get_input_plates_needed(self, searchsmiles: list, otbatchprotocolobj, reaction_ids: list = None, groupreactionqueryset=None) -> list:
        """
        Gets plates created in previous reaction and workup sessions with reaction products 
        that are required as reactants in current reaction session.
        
        Parameters
        ----------
        searchsmiles: list
            The list of SMILES that are required from previous reaction plate wells
        otbatchprotocolobj: OTBatchProtocol
            The batch protocol object
        reaction_ids: list, optional
            The optional reaction ids to match wells and plates with
        groupreactionqueryset: QuerySet, optional
            The reaction queryset to extract method IDs from if reaction_ids not provided
            
        Returns
        -------
        inputplatesneeded: list
            List of previous OT session plates that have products needed for this reaction
        """
        input_plates_needed = []
        
        # Get all plates from the current batch protocol
        otbatchprotocol_plates = self.get_all_ot_batch_protocol_plates(
            otbatchprotocol_id=otbatchprotocolobj
        )
        
        if not otbatchprotocol_plates:
            return []
            
        # Build filter criteria
        if not reaction_ids and groupreactionqueryset:
            method_ids = [reaction.method_id for reaction in groupreactionqueryset]
            criterion1 = Q(method_id__in=method_ids)
        elif reaction_ids:
            criterion1 = Q(reaction_id__in=reaction_ids)
        else:
            logger.warning("Neither reaction_ids nor groupreactionqueryset provided")
            return []
            
        criterion2 = Q(reactantfornextstep=True)
        criterion3 = Q(smiles__in=searchsmiles)
        criterion4 = Q(type__in=["reaction", "workup1", "workup2", "workup3", "spefilter"])
        
        # Find plates with wells matching our criteria
        for plate in otbatchprotocol_plates:
            well_matches = plate.well_set.filter(
                criterion1 & criterion2 & criterion3 & criterion4
            )
            if well_matches.exists():
                input_plates_needed.append(plate)
                
        logger.info(f"Found {len(input_plates_needed)} input plates with required materials")
        return input_plates_needed