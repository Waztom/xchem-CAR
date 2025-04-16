"""Manages deck operations for OpenTrons sessions."""

import logging
from django.db.models import Q
from ...models import Deck, Plate, Well, Column

logger = logging.getLogger(__name__)

class DeckManager:
    """
    Manages deck slot allocation and resources.
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
    
    def create_deck_model(self):
        """
        Create a deck object.
        
        Returns
        -------
        deckobj: Deck
            The created deck object
        """
        deckobj = Deck()
        deckobj.otsession_id = self.session.otsessionobj
        deckobj.numberslots = 11
        deckobj.save()
        return deckobj
    
    def check_deck_slot_available(self):
        """
        Check if a deck slot is available and allocate it.
        
        Returns
        -------
        testslotavailable: int
            The index of the deck slot available
            
        Raises
        ------
        ValueError
            When no deck slots are available
        """
        testslotavailable = self.session.deckobj.indexslotavailable
        if testslotavailable <= self.session.deckobj.numberslots:
            self.session.deckobj.indexslotavailable = testslotavailable + 1
            self.session.deckobj.save()
            return testslotavailable
        else:
            self.session.deckobj.slotavailable = False
            self.session.deckobj.save()
            logger.error("No deck slots available")
            self.session.cleanup()
            raise ValueError("No deck slots available - cannot create more plates")
    
    