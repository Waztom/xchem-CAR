import logging
from django.db.models import QuerySet, Q

from .....models import (
    ActionSession,
    AddAction,
    ExtractAction,
    Plate,
)

logger = logging.getLogger(__name__)


class PlateQueryService:
    def __init__(self, session):
        """
        Initialize with a reference to the parent session.

        Parameters
        ----------
        session: BaseSession
            The parent session object
        """
        self.session = session

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
            role="startingmaterial",
            otsession_id__in=custom_session_ids,
        )

        return plates

    def get_action_session_by_plate_role(
        self, role: str, role_index: int = 1
    ) -> QuerySet[ActionSession]:
        """
        Get action sessions that target a specific plate role and index.

        Parameters
        ----------
        role: str
            The plate role to look for (e.g., "workup", "reaction")
        role_index: int, optional
            The plate role index (default 1)

        Returns
        -------
        action_session_queryset: QuerySet[ActionSession]
            Action sessions targeting the specified role/index
        """
        criterion1 = Q(id__in=self.session.actionsessionqueryset)
        criterion2 = Q(addaction__to_plate_role=role, addaction__to_plate_index=role_index)
        criterion3 = Q(extractaction__to_plate_role=role, extractaction__to_plate_index=role_index)

        action_session_queryset = ActionSession.objects.filter(
            criterion1 & (criterion2 | criterion3)
        ).distinct()

        return action_session_queryset

    def get_unique_to_plates(
        self, action_session_queryset: QuerySet[ActionSession], platetypes: list
    ) -> list:
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
        criterion2 = Q(role__in=["reaction", "workup", "spefilter"])

        otbatchprotocol_plate_queryset = Plate.objects.filter(criterion1 & criterion2)
        return otbatchprotocol_plate_queryset

    def get_input_plates_needed(
        self,
        searchsmiles: list,
        reaction_ids: list = None,
        groupreactionqueryset=None,
    ) -> list:
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
            otbatchprotocol_id=self.session.otbatchprotocolobj
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
        criterion4 = Q(role__in=["reaction", "workup", "spefilter"])

        # Find plates with wells matching our criteria
        for plate in otbatchprotocol_plates:
            well_matches = plate.well_set.filter(
                criterion1 & criterion2 & criterion3 & criterion4
            )
            if well_matches.exists():
                input_plates_needed.append(plate)

        logger.info(
            f"Found {len(input_plates_needed)} input plates with required materials"
        )
        return input_plates_needed
