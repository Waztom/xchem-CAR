"""
Manages pipette and tip rack selection, creation, and configuration.
"""

import logging
from statistics import median
from typing import List, Dict, Any

from ...models import Pipette, TipRack

logger = logging.getLogger(__name__)

class PipetteManager:
    """
    Manages pipette selection, tip rack creation, and related functionality.
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
        
        # Define available pipettes and their properties
        self.pipettes_available = [
            {
                "labware": "p10_single",
                "position": "right",
                "type": "single",
                "maxvolume": 10,
                "minvolume": 1,
            },
            {
                "labware": "p300_single",
                "position": "right",
                "type": "single",
                "maxvolume": 300,
                "minvolume": 30,
            },
            {
                "labware": "p1000_single_gen2",
                "position": "right",
                "type": "single",
                "maxvolume": 1000,
                "minvolume": 100,
            },
            {
                "labware": "p300_multi_gen2",
                "position": "left",
                "type": "multi",
                "maxvolume": 300,
                "minvolume": 30,
            },
            {
                "labware": "p20_multi_gen2",
                "position": "left",
                "type": "multi",
                "maxvolume": 20,
                "minvolume": 1,
            },
        ]
        
        # Define available tip racks
        self.tipracks_available = {
            "opentrons_96_tiprack_20ul": {"maxvolume": 20, "minvolume": 1},
            "opentrons_96_tiprack_300ul": {"maxvolume": 300, "minvolume": 20},
            "opentrons_96_tiprack_1000ul": {"maxvolume": 1000, "minvolume": 100},
        }
    
    def get_tip_rack_type(self, rounded_volumes: list) -> str:
        """
        Gets OT tiprack best suited for transferring volumes (ul)
        that minimises the number of transfers required.
        
        Parameters
        ----------
        rounded_volumes: list
            The list of rounded volumes needed for a set of actions
            
        Returns
        -------
        tipracktype: str
            The most suitable tiprack type
        """
        try:
            if not rounded_volumes:
                logger.warning("No volumes provided for tip rack selection")
                return "opentrons_96_tiprack_300ul"  # Default to medium tips
            
            # Find the median volume to determine most appropriate tip size
            median_volume = median(rounded_volumes)
            logger.debug(f"Median volume for tip selection: {median_volume}µL")
            
            # Choose appropriate tip rack based on volume needs
            if median_volume <= 20:
                return "opentrons_96_tiprack_20ul"
            elif median_volume <= 300:
                return "opentrons_96_tiprack_300ul"
            else:
                return "opentrons_96_tiprack_1000ul"
                
        except Exception as e:
            logger.error(f"Error selecting tip rack type: {str(e)}")
            return "opentrons_96_tiprack_300ul"  # Default to medium tips if there's an error
    
    def get_pipette_type(self, rounded_volumes: list, channel_type: str = "single") -> dict:
        """
        Gets the type of pipette that minimizes the number of transfers
        needed for transferring volumes (ul).
        
        Parameters
        ----------
        rounded_volumes: list
            The list of rounded volumes that need to be transferred
        channel_type: str
            The type of channel eg. single or multi. Default set to single
            
        Returns
        -------
        pipette_type: dict
            The pipette configuration needed
        """
        try:
            if not rounded_volumes:
                # Default to P300 if no volumes provided
                for pipette in self.pipettes_available:
                    if pipette["maxvolume"] == 300 and pipette["type"] == channel_type:
                        return pipette
                # If no matching pipette found, default to first compatible one
                for pipette in self.pipettes_available:
                    if pipette["type"] == channel_type:
                        return pipette
                return self.pipettes_available[0]
            
            # Get median volume for optimal pipette selection
            median_transfer_volume = median(rounded_volumes)
            logger.debug(f"Median volume for pipette selection: {median_transfer_volume}µL")
            
            # Filter for pipettes that match the channel type
            compatible_pipettes = [p for p in self.pipettes_available if p["type"] == channel_type]
            if not compatible_pipettes:
                logger.warning(f"No pipettes found for channel type: {channel_type}, using defaults")
                return self.pipettes_available[0]
            
            # Find pipettes with sufficient volume capacity
            suitable_pipettes = [p for p in compatible_pipettes if p["maxvolume"] >= median_transfer_volume]
            if not suitable_pipettes:
                # If no pipette has sufficient capacity, use the largest available
                largest_pipette = max(compatible_pipettes, key=lambda p: p["maxvolume"])
                logger.warning(
                    f"No pipettes can handle {median_transfer_volume}µL in a single transfer. "
                    f"Using largest available: {largest_pipette['labware']} ({largest_pipette['maxvolume']}µL)"
                )
                return largest_pipette
            
            # Compare pipettes based on transfer efficiency
            pipette_compare_dict = {}
            for pipette in suitable_pipettes:
                volume_difference = pipette["maxvolume"] - median_transfer_volume
                number_transfers = self.get_number_transfers(
                    pipette_volume=pipette["maxvolume"], rounded_volumes=rounded_volumes
                )
                pipette_compare_dict[pipette["labware"]] = {
                    "notransfers": number_transfers,
                    "volumedifference": volume_difference,
                }
            
            # Select pipette with minimum number of transfers
            if not pipette_compare_dict:
                return suitable_pipettes[0]
                
            minimum_transfers = min(
                [pipette_data["notransfers"] for pipette_data in pipette_compare_dict.values()]
            )
            
            optimal_pipettes = [
                pipette_name for pipette_name, pipette_data in pipette_compare_dict.items()
                if pipette_data["notransfers"] == minimum_transfers
            ]
            
            # If multiple pipettes have the same number of transfers, choose the one with
            # the smallest volume difference (most precise)
            if len(optimal_pipettes) > 1:
                minimum_volume_difference = min(
                    [
                        pipette_compare_dict[pipette_name]["volumedifference"]
                        for pipette_name in optimal_pipettes
                    ]
                )
                
                optimal_pipettes = [
                    pipette_name for pipette_name in optimal_pipettes
                    if pipette_compare_dict[pipette_name]["volumedifference"] == minimum_volume_difference
                ]
            
            # Get the full pipette info for the selected pipette
            selected_pipette = next(
                (pipette for pipette in self.pipettes_available if pipette["labware"] == optimal_pipettes[0]),
                suitable_pipettes[0]  # Fallback
            )
            
            logger.info(f"Selected pipette: {selected_pipette['labware']} for {channel_type} channel operations")
            return selected_pipette
            
        except Exception as e:
            logger.error(f"Error selecting pipette type: {str(e)}")
            # Default to P300 if there's an error
            for pipette in self.pipettes_available:
                if pipette["maxvolume"] == 300 and pipette["type"] == channel_type:
                    return pipette
            return self.pipettes_available[0]
    
    def get_number_transfers(self, pipette_volume: int, rounded_volumes: list) -> int:
        """
        Gets the number of transfers required for transferring
        a list of rounded volumes.
        
        Parameters
        ----------
        pipette_volume: int
            The pipette's maximum transfer volume (ul)
        rounded_volumes: list
            The list of rounded volumes that need to be transferred
            
        Returns
        -------
        number_transfers: int
            The number of transfers required for the pipette type used
        """
        try:
            number_transfers = 0
            
            for volume in rounded_volumes:
                if volume <= pipette_volume:
                    # Can transfer in a single operation
                    number_transfers += 1
                else:
                    # Need multiple transfers
                    transfers_needed = (volume + pipette_volume - 1) // pipette_volume  # Ceiling division
                    number_transfers += transfers_needed
            
            return number_transfers
            
        except Exception as e:
            logger.error(f"Error calculating number of transfers: {str(e)}")
            return len(rounded_volumes)  # Default to one transfer per volume
    
    def create_pipette_model(self):
        """
        Create a pipette object in the database.
        
        Returns
        -------
        pipette_obj: Pipette
            The created pipette object
        """
        try:
            if not hasattr(self.session, 'pipettetype'):
                logger.error("No pipette type selected for this session")
                return None
                
            pipette_obj = Pipette()
            pipette_obj.otsession_id = self.session.otsessionobj
            pipette_obj.position = self.session.pipettetype["position"]
            pipette_obj.maxvolume = self.session.pipettetype["maxvolume"]
            pipette_obj.type = self.session.pipettetype["type"]
            pipette_obj.name = f"{self.session.pipettetype['position']}_{self.session.pipettetype['labware']}"
            pipette_obj.labware = self.session.pipettetype["labware"]
            pipette_obj.save()
            
            logger.info(f"Created pipette model: {pipette_obj.name} (ID: {pipette_obj.id})")
            return pipette_obj
            
        except Exception as e:
            logger.error(f"Error creating pipette model: {str(e)}")
            return None
    
    def create_tiprack_model(self, name: str):
        """
        Creates TipRack object in the database.
        
        Parameters
        ----------
        name: str
            The name/type of the tiprack
            
        Returns
        -------
        tiprack_obj: TipRack or None
            The created tiprack object, or None if no deck slots available
        """
        try:
            index_slot = self.session.deck_manager.check_deck_slot_available()
            if not index_slot:
                logger.warning("No more deck slots available for tiprack")
                return None
                
            tiprack_obj = TipRack()
            tiprack_obj.otsession_id = self.session.otsessionobj
            tiprack_obj.deck_id = self.session.deckobj
            tiprack_obj.name = f"{name}_{index_slot}"
            tiprack_obj.index = index_slot
            tiprack_obj.labware = name
            
            # Set volume properties if available
            if name in self.tipracks_available:
                tiprack_obj.maxvolume = self.tipracks_available[name]["maxvolume"]
                tiprack_obj.minvolume = self.tipracks_available[name]["minvolume"]
                
            tiprack_obj.save()
            
            logger.info(f"Created tiprack model: {tiprack_obj.name} at position {index_slot}")
            return tiprack_obj
            
        except Exception as e:
            logger.error(f"Error creating tiprack model: {str(e)}")
            return None
    
    def create_tip_racks(self, tiprack_type: str, number_tipracks: int = 3):
        """
        Creates multiple tipracks of the specified type.
        
        Parameters
        ----------
        tiprack_type: str
            The type of tiprack needed
        number_tipracks: int
            The number of tip racks to create. Default is three.
            
        Returns
        -------
        tipracks: list
            List of created tiprack objects
        """
        try:
            tipracks = []
            
            logger.info(f"Creating {number_tipracks} tipracks of type {tiprack_type}")
            for i in range(number_tipracks):
                tiprack = self.create_tiprack_model(name=tiprack_type)
                if tiprack:
                    tipracks.append(tiprack)
                else:
                    logger.warning(f"Could only create {i} tipracks (deck full)")
                    break
            
            return tipracks
            
        except Exception as e:
            logger.error(f"Error creating tip racks: {str(e)}")
            return []
            
    def setup_pipettes_for_session(self, volumes: list, channel_type: str = "single"):
        """
        Set up all pipettes and tipracks needed for a session.
        
        Parameters
        ----------
        volumes: list
            The list of volumes to be transferred
        channel_type: str
            The type of channel to use (single or multi)
            
        Returns
        -------
        success: bool
            True if setup was successful, False otherwise
        """
        try:
            # Select optimal tiprack and pipette
            self.session.tipracktype = self.get_tip_rack_type(rounded_volumes=volumes)
            self.session.pipettetype = self.get_pipette_type(
                rounded_volumes=volumes, channel_type=channel_type
            )
            
            # Create tipracks
            tipracks = self.create_tip_racks(tiprack_type=self.session.tipracktype)
            if not tipracks:
                logger.error("Failed to create any tip racks")
                return False
                
            # Create pipette
            pipette = self.create_pipette_model()
            if not pipette:
                logger.error("Failed to create pipette")
                return False
                
            logger.info(f"Successfully set up pipettes for session with {len(tipracks)} tipracks")
            return True
            
        except Exception as e:
            logger.error(f"Error setting up pipettes for session: {str(e)}")
            return False