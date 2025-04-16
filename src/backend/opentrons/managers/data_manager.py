"""
Manages data query, processing, and analysis for OpenTrons sessions.
"""

import logging
import pandas as pd
import numpy as np
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.db.models import QuerySet, Q

from ...models import (
    AddAction, ExtractAction, Reaction, CompoundOrder
)
from ...utils import (
    getReactionTemperature, getReactionClass, getReactionRecipe,
    getReactionQuerySet
)

logger = logging.getLogger(__name__)

class DataManager:
    """
    Manages data retrieval, transformation, and analysis operations.
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
    
    def get_median_value(self, values: list) -> float:
        """
        Gets median value from a list of values.
        
        Parameters
        ----------
        values: list
            The list of values to get the median of
            
        Returns
        -------
        median: float
            The median value
        """
        if not values:
            return 0
            
        n_values = len(values)
        sorted_values = sorted(values)
        
        if n_values % 2 == 0:
            median = (sorted_values[n_values//2 - 1] + sorted_values[n_values//2]) / 2
        else:
            median = sorted_values[n_values//2]
            
        return median
    
    def get_sum_value(self, values: list) -> float:
        """
        Gets sum of a list of values.
        
        Parameters
        ----------
        values: list
            The list of values to get the sum of
            
        Returns
        -------
        sum: float
            The sum of values
        """
        return sum(values)
    
    def get_add_action_query_set(
        self,
        reaction_ids: list,
        actionsession_ids: list = None,
        actionsessiontype: str = None,
    ) -> QuerySet[AddAction]:
        """
        Get the add action queryset from a list of reactions
        and action session ids.
        
        Parameters
        ----------
        reaction_ids: list
            The list of reaction IDs to get the add actions for
        actionsession_ids: list = None
            Optional list of actionsession IDs to filter by
        actionsessiontype: str = None
            Optional action session type to filter by
            
        Returns
        -------
        add_action_queryset: QuerySet[AddAction]
            The queryset containing all add actions
        """
        if actionsession_ids is None:
            actionsession_ids = []
            
        criterion1 = Q(reaction_id__in=reaction_ids)
        
        if actionsession_ids:
            criterion2 = Q(actionsession_id__in=actionsession_ids)
        else:
            criterion2 = Q()
            
        if actionsessiontype:
            criterion3 = Q(actionsession_id__type=actionsessiontype)
        else:
            criterion3 = Q()
            
        add_action_queryset = AddAction.objects.filter(criterion1 & criterion2 & criterion3)
        return add_action_queryset
    
    def get_extract_action_query_set(
        self,
        reaction_ids: list,
        actionsession_ids: list = None,
        actionsessiontype: str = None,
    ) -> QuerySet[ExtractAction]:
        """
        Get the extract action queryset from a list of reactions
        and action session ids.
        
        Parameters
        ----------
        reaction_ids: list
            The list of reaction IDs to get the extract actions for
        actionsession_ids: list = None
            Optional list of actionsession IDs to filter by
        actionsessiontype: str = None
            Optional action session type to filter by
            
        Returns
        -------
        extract_action_queryset: QuerySet[ExtractAction]
            The queryset containing all extract actions
        """
        if actionsession_ids is None:
            actionsession_ids = []
            
        criterion1 = Q(reaction_id__in=reaction_ids)
        
        if actionsession_ids:
            criterion2 = Q(actionsession_id__in=actionsession_ids)
        else:
            criterion2 = Q()
            
        if actionsessiontype:
            criterion3 = Q(actionsession_id__type=actionsessiontype)
        else:
            criterion3 = Q()
            
        extract_action_queryset = ExtractAction.objects.filter(criterion1 & criterion2 & criterion3)
        return extract_action_queryset
    
    def get_rounded_add_action_volumes(self, addactionqueryset: QuerySet[AddAction]) -> list:
        """
        Get the rounded volumes from the add action queryset.
        
        Parameters
        ----------
        addactionqueryset: QuerySet[AddAction]
            The queryset containing add actions to get volumes from
            
        Returns
        -------
        rounded_volumes: list
            The list of rounded volumes
        """
        volumes = addactionqueryset.values_list("volume", flat=True)
        
        # Round volumes to nearest 5 μL for easier pipetting
        rounded_volumes = [5 * round(vol / 5) for vol in volumes]
        
        # Make sure no volumes are 0
        rounded_volumes = [max(5, vol) for vol in rounded_volumes]
        
        return rounded_volumes
    
    def get_rounded_extract_action_volumes(self, extractactionqueryset: QuerySet[ExtractAction]) -> list:
        """
        Get the rounded volumes from the extract action queryset.
        
        Parameters
        ----------
        extractactionqueryset: QuerySet[ExtractAction]
            The queryset containing extract actions to get volumes from
            
        Returns
        -------
        rounded_volumes: list
            The list of rounded volumes
        """
        volumes = extractactionqueryset.values_list("volume", flat=True)
        
        # Round volumes to nearest 5 μL for easier pipetting
        rounded_volumes = [5 * round(vol / 5) for vol in volumes]
        
        # Make sure no volumes are 0
        rounded_volumes = [max(5, vol) for vol in rounded_volumes]
        
        return rounded_volumes
    
    def get_rounded_reaction_volumes(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """
        Get rounded reaction volumes from reaction queryset.
        
        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reaction queryset to get volumes from
            
        Returns
        -------
        rounded_volumes: list
            The list of rounded volumes
        """
        volumes = reactionqueryset.values_list("totalvolume", flat=True)
        
        # Round volumes to nearest 5 μL for easier pipetting
        rounded_volumes = [5 * round(vol / 5) for vol in volumes]
        
        # Make sure no volumes are 0
        rounded_volumes = [max(5, vol) for vol in rounded_volumes]
        
        return rounded_volumes
    
    def get_add_actions_dataframe(self, addactionqueryset: QuerySet[AddAction]) -> pd.DataFrame:
        """
        Converts add action queryset to a dataframe.
        
        Parameters
        ----------
        addactionqueryset: QuerySet[AddAction]
            The queryset containing add actions
            
        Returns
        -------
        add_actions_df: DataFrame
            The dataframe of add actions
        """
        try:
            if not addactionqueryset.exists():
                return pd.DataFrame()
                
            # Create a list to store all the action data
            add_actions_list = []
            
            for addaction in addactionqueryset:
                action_data = {
                    "id": addaction.id,
                    "reaction_id_id": addaction.reaction_id.id,
                    "actionsession_id_id": addaction.actionsession_id.id,
                    "smiles": addaction.smiles,
                    "name": addaction.name,
                    "volume": addaction.volume,
                    "concentration": addaction.concentration,
                    "solvent": addaction.solvent,
                    "molecularweight": addaction.molecularweight,
                    "fromplatetype": addaction.fromplatetype,
                    "toplatetype": addaction.toplatetype,
                    "uniquesolution": f"{addaction.smiles}-{addaction.concentration}-{addaction.solvent}"
                }
                add_actions_list.append(action_data)
            
            # Convert to dataframe
            add_actions_df = pd.DataFrame(add_actions_list)
            return add_actions_df
            
        except Exception as e:
            logger.error(f"Error creating add actions dataframe: {str(e)}")
            return pd.DataFrame()
    
    def get_grouped_temperature_reactions(self, reactionqueryset: QuerySet[Reaction]) -> dict:
        """
        Groups reactions by temperature.
        
        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reaction queryset to group
            
        Returns
        -------
        grouped_reactions: dict
            Dictionary of reaction querysets grouped by temperature
        """
        unique_temps = self.get_unique_temperatures(reactionqueryset=reactionqueryset)
        
        grouped_reactions = {}
        for temp in unique_temps:
            reactions_at_temp = []
            for reaction in reactionqueryset:
                if getReactionTemperature(reaction_id=reaction.id) == temp:
                    reactions_at_temp.append(reaction.id)
            
            if reactions_at_temp:
                grouped_reactions[temp] = getReactionQuerySet(reaction_ids=reactions_at_temp)
        
        return grouped_reactions
    
    def get_grouped_reaction_by_class_recipe(self, reactionqueryset: QuerySet[Reaction]) -> dict:
        """
        Groups reactions by class and recipe.
        
        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reaction queryset to group
            
        Returns
        -------
        grouped_reactions: dict
            Dictionary of reaction querysets grouped by class-recipe combination
        """
        unique_classes = self.get_unique_reaction_classes(reactionqueryset=reactionqueryset)
        unique_recipes = self.get_unique_reaction_recipes(reactionqueryset=reactionqueryset)
        
        grouped_reactions = {}
        
        for reaction_class in unique_classes:
            for recipe in unique_recipes:
                reactions_in_group = []
                
                for reaction in reactionqueryset:
                    if (getReactionClass(reaction_id=reaction.id) == reaction_class and
                        getReactionRecipe(reaction_id=reaction.id) == recipe):
                        reactions_in_group.append(reaction.id)
                
                if reactions_in_group:
                    group_key = f"{reaction_class}-{recipe}"
                    grouped_reactions[group_key] = getReactionQuerySet(reaction_ids=reactions_in_group)
        
        return grouped_reactions
    
    def get_unique_temperatures(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """
        Gets the unique temperatures from reaction queryset.
        
        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reaction queryset to get unique temperatures from
            
        Returns
        -------
        unique_temperatures: list
            List of unique temperatures
        """
        unique_temperatures = []
        
        for reaction in reactionqueryset:
            temp = getReactionTemperature(reaction_id=reaction.id)
            if temp not in unique_temperatures:
                unique_temperatures.append(temp)
        
        return unique_temperatures
    
    def get_unique_reaction_classes(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """
        Gets the unique reaction classes from reaction queryset.
        
        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reaction queryset to get unique classes from
            
        Returns
        -------
        unique_classes: list
            List of unique reaction classes
        """
        unique_classes = []
        
        for reaction in reactionqueryset:
            reaction_class = getReactionClass(reaction_id=reaction.id)
            if reaction_class not in unique_classes:
                unique_classes.append(reaction_class)
        
        return unique_classes
    
    def get_unique_reaction_recipes(self, reactionqueryset: QuerySet[Reaction]) -> list:
        """
        Gets the unique reaction recipes from reaction queryset.
        
        Parameters
        ----------
        reactionqueryset: QuerySet[Reaction]
            The reaction queryset to get unique recipes from
            
        Returns
        -------
        unique_recipes: list
            List of unique reaction recipes
        """
        unique_recipes = []
        
        for reaction in reactionqueryset:
            recipe = getReactionRecipe(reaction_id=reaction.id)
            if recipe not in unique_recipes:
                unique_recipes.append(recipe)
        
        return unique_recipes
    
    def create_compound_order_model(self, order_df: pd.DataFrame, is_custom_starter_plate: bool = False):
        """
        Creates a compound order object and CSV file.
        
        Parameters
        ----------
        order_df: DataFrame
            The dataframe containing order information
        is_custom_starter_plate: bool
            Whether this is a custom starter materials plate
            
        Returns
        -------
        compound_order_obj: CompoundOrder
            The created compound order object
        """
        try:
            compound_order_obj = CompoundOrder()
            compound_order_obj.otsession_id = self.session.otsessionobj
            compound_order_obj.iscustomSMplate = is_custom_starter_plate
            
            # Convert the dataframe to CSV format
            csv_data = order_df.to_csv(encoding="utf-8", index=False)
            
            # Generate a filename based on session type and IDs
            filename = (
                f"{self.session.actionsessiontype}-session-orderplate-for-batch-"
                f"{self.session.batchobj.batchtag}-reactionstep-{self.session.reactionstep}-"
                f"sessionid-{str(self.session.otsessionobj.id)}.csv"
            )
            
            # Save the CSV file to storage
            order_csv = default_storage.save(
                f"compoundorders/{filename}",
                ContentFile(csv_data),
            )
            
            # Store the file path
            compound_order_obj.ordercsv = order_csv
            compound_order_obj.save()
            
            return compound_order_obj
            
        except Exception as e:
            logger.error(f"Error creating compound order model: {str(e)}")
            return None