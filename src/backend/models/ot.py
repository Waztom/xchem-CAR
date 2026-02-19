"""OpenTrons session models: OTProject, OTSession, Deck, Plate, Well, etc."""

from django.db import models

from .core import Project, Batch, Method, Reaction
from .actions import PlateRole, ActionSession


class OTProject(models.Model):
    """Django model to define an OTProject - an OT project will
       have one or more batch protocols for a project

    Parameters
    ----------
    project_id: ForeignKey
        Foreign key linking an OT protocol action to a reaction
    init_date: DateTimeField
        The date the OT project was created (autofield)
    name: CharField
        The name of the OT project
    """

    project_id = models.ForeignKey(Project, on_delete=models.CASCADE)
    init_date = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=150)


class OTBatchProtocol(models.Model):
    """Django model to define an OTBatchProtocol - OT protocols for a batch of
       targets

    Parameters
    ----------
    otproject_id: ForeignKey
        Foreign key linking a batch OT protocol with a OT project
    batch_id: ForeignKey
        Foreign key linking a batch OT protocol with batch of targets
    celery_taskid: CharField
        The Celery task id when a new OT batch protocol is created
    zipfile: FileField
        File link to a zip file of all the OT protocols required for executing the
        synthesis of a batch of targets
    """

    otproject_id = models.ForeignKey(OTProject, on_delete=models.CASCADE)
    batch_id = models.ForeignKey(
        Batch, related_name="otbatchprotocols", on_delete=models.CASCADE
    )
    celery_taskid = models.CharField(max_length=50)
    zipfile = models.FileField(upload_to="otbatchprotocols/", max_length=255, null=True)


class OTSession(models.Model):
    """Django model to define an OT Session - a OT session is a session (Reaction, Analysis)
       that needs to be executed on the OT

    Parameters
    ----------
    otbatchprotocol_id: ForeignKey
        Foreign key linking an OT session to a OT batch protocol
    reactionstep: IntegerField
        The reaction step that the session is being executed for
    sessiontype: CharField
        The type of session ebing excecuted
    """

    class SessionType(models.TextChoices):
        reaction = "reaction"
        workup = "workup"
        lcmsprep = "lcmsprep"

    otbatchprotocol_id = models.ForeignKey(
        OTBatchProtocol, related_name="otsessions", on_delete=models.CASCADE
    )
    reactionstep = models.IntegerField()
    sessiontype = models.CharField(
        choices=SessionType.choices, default=SessionType.reaction, max_length=10
    )


class Deck(models.Model):
    """Django model to define a Deck Session - an OT Deck

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a deck to an OT session
    numberslots: IntegerField
        The number of deck slots (default=11)
    slotavailable: BooleanField
        If a slot is still available on the deck (default=True)
    indexslotavailable: IntegerField
        The index (1-11) of the deck slot available
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="otdecks", on_delete=models.CASCADE
    )
    numberslots = models.IntegerField(default=11)
    slotavailable = models.BooleanField(default=True)
    indexslotavailable = models.IntegerField(default=1)


class Pipette(models.Model):
    """Django model to define a Pipette - an OT Pipette

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a deck to an OT session
    position: CharField
        The position (right or left) of the pipette for the OT session
    labware: CharField
        The name of the OT labware eg. p300_single
    type: CharField
        The type of pipette could be single or multi channel
    maxvolume: FloatField
        The maximum volume (uL) of the pipette
    name: CharField
        The name of the pipette - this is a combination of the OT labware name (p300_single)
        and the position of the pipette uses in the session (right) eg. right_p300_single
    """

    class Position(models.TextChoices):
        right = "Right"
        left = "Left"

    otsession_id = models.ForeignKey(
        OTSession, related_name="pipettes", on_delete=models.CASCADE
    )
    position = models.CharField(choices=Position.choices, max_length=10)
    labware = models.CharField(max_length=100)
    type = models.CharField(max_length=255)
    maxvolume = models.FloatField(default=300)
    name = models.CharField(max_length=255)


class TipRack(models.Model):
    """Django model to define a TipRack - an OT tiprack

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a tiprack to an OT session
    deck_id: ForeignKey
        Foreign key linking a tiprack to an OT deck
    labware: CharField
        The name of the OT labware eg. opentrons_96_tiprack_300ul
    index: IntegerField
        The deck index (1-11) of the tiprack
    name: CharField
        The name of the tiprack - this is a combination of the OT labware name (opentrons_96_tiprack_300ul)
        and the deck index of the tiprack eg. opentrons_96_tiprack_300ul_2
    """

    otsession_id = models.ForeignKey(OTSession, on_delete=models.CASCADE)
    deck_id = models.ForeignKey(Deck, on_delete=models.CASCADE)
    labware = models.CharField(max_length=255)
    index = models.IntegerField()
    name = models.CharField(max_length=255)


class Plate(models.Model):
    """Django model to define a Plate - an OT plate

    Parameters
    ----------
    otbatchprotocol_id: ForeignKey
        Foreign key linking a plate to an OT batch protocol
    deck_id: ForeignKey
        Foreign key linking a plate to an OT deck
    labware: CharField
        The name of the OT labware eg. labcyte_384_wellplate_100ul
    index: IntegerField
        The deck index (1-11) of the plate
    name: CharField
        The name of the plate
    role: CharField
        The functional role of the plate (reaction, workup, spefilter, etc.)
    role_index: PositiveIntegerField
        Instance number within the role (default 1). e.g. role=workup,
        role_index=2 is the second workup plate.
    maxwellvolume: FloatField
        The maximum plate well volume (uL)
    numberwells: IntegerField
        The number of plate wells
    wellavailable: BooleanField
        If a well is available on a plate
    numberwellsincolumn: IntegerField
        The number of wells in a column default is set to 8 for a 96 well plate
    indexswellavailable: IntegerField
        The index of the well available. Wells are occupied in increasing
        index starting from indices: A1, B1, C1 or 1, 2, 3 etc.
    columnavailable: BooleanField
        Wether a column of a plate is available. Example 4 rows of column 1 taken up
        by amidation reactions, column is then no longer avilable to any other reaction
        classes.
    indexscolumnavailable: IntegerField
        The index of the column available. Columns are occupied in increasing
        index starting from indices: A1, B1, C1 or 1, 2, 3 etc
    """

    otbatchprotocol_id = models.ForeignKey(
        OTBatchProtocol, related_name="otplates", on_delete=models.CASCADE
    )
    otsession_id = models.ForeignKey(
        OTSession,
        on_delete=models.CASCADE,
        null=True,
    )
    deck_id = models.ForeignKey(Deck, on_delete=models.CASCADE)
    labware = models.CharField(max_length=255)
    index = models.IntegerField()
    name = models.CharField(max_length=255, null=True)
    role = models.CharField(choices=PlateRole.choices, max_length=20, null=True)
    role_index = models.PositiveIntegerField(default=1)
    maxwellvolume = models.FloatField()
    numberwells = models.IntegerField()
    wellavailable = models.BooleanField(default=True)
    numberwellsincolumn = models.IntegerField(default=8)
    indexswellavailable = models.IntegerField(default=0)
    numbercolumns = models.IntegerField()
    columnavailable = models.BooleanField(default=True)
    indexcolumnavailable = models.IntegerField(default=0)


class Column(models.Model):
    """Django model to define a plate column

    Parameters
    ----------
    otbatchprotocol_id: ForeignKey
        Foreign key linking a plate to an OT batch protocol
    plate_id: ForeignKey
        Foreign key linking a well to a plate
    index: IntegerField
        The column index (0-11) on the plate
    role: CharField
        The functional role of the column (matches parent plate role)
    role_index: PositiveIntegerField
        Instance number within the role (matches parent plate role_index)
    reactionclass: CharField
        The reaction class eg. amidation. Each column can only contain
        one type of reaction class -> for multi-pipette handling
        and grouping reactions on plates
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="otcolumns", on_delete=models.CASCADE, null=True
    )
    plate_id = models.ForeignKey(Plate, on_delete=models.CASCADE)
    index = models.IntegerField()
    role = models.CharField(choices=PlateRole.choices, max_length=20, null=True)
    role_index = models.PositiveIntegerField(default=1)
    reactionclass = models.CharField(max_length=100)


class Well(models.Model):
    """Django model to define a Well - an OT plate well

    Parameters
    ----------
    otbatchprotocol_id: ForeignKey
        Foreign key linking a plate to an OT batch protocol
    plate_id: ForeignKey
        Foreign key linking a well to a plate
    method_id: ForeignKey
        Optional foreign key linking a well to a method
    reaction_id: ForeignKey
        Optional foreign key linking a well to a reaction
    index: IntegerField
        The well index (0-95) on the plate
    name: CharField
        The name of the well eg. A1, B1, C1
    type: CharField
        The type of well eg. analyse and reaction well.
        Uses PlateRole choices with a role_index for the plate
        instance number.
    role_index: PositiveIntegerField
        Instance number within the role (matches parent plate role_index)
    volume: FloatField
        The optional volume of the contents in the well (uL)
    smiles: CharField
        The optional SMILES of the well contents
    concentration: FloatField
        The optional concentration of the well contents
    solvent: CharField
        The optional solvent used for diluting the well contents
    reactantfornextstep: BooleanField
        Wether the contents are used for the next reaction step (default=False)
    available: BooleanField
        If the well is available w.r.t still containing it's contents
        vs. being empty (default=True)
    """

    otsession_id = models.ForeignKey(
        OTSession,
        related_name="otwells",
        on_delete=models.CASCADE,
        null=True,
    )
    plate_id = models.ForeignKey(Plate, on_delete=models.CASCADE)
    method_id = models.ForeignKey(Method, on_delete=models.CASCADE, null=True)
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE, null=True)
    column_id = models.ForeignKey(Column, on_delete=models.CASCADE, null=True)
    index = models.IntegerField()
    name = models.CharField(max_length=3, null=True)
    role = models.CharField(choices=PlateRole.choices, max_length=20, null=True)
    role_index = models.PositiveIntegerField(default=1)
    volume = models.FloatField(null=True)
    smiles = models.CharField(max_length=255, null=True)
    concentration = models.FloatField(null=True)
    solvent = models.CharField(max_length=255, null=True)
    reactantfornextstep = models.BooleanField(default=False)
    available = models.BooleanField(default=True)


class CompoundOrder(models.Model):
    """Django model to define a CompoundOrder - a csv
       file of ordering information that includes SMILES,
       plate name, well id, amount, solvent and concentration
       required for the synthesis of the reaction step

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    iscustomSMplate: BooleanField
        If the plate is a custom SM plate (default=False)
    ordercsv: FileField
        The csv file of the ordering information for executing
        a reaction step on the OpenTrons
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="compoundorders", on_delete=models.CASCADE
    )
    iscustomSMplate = models.BooleanField(default=False)
    ordercsv = models.FileField(upload_to="compoundorders/", max_length=255)


class OTScript(models.Model):
    """Django model to define a OTScript - a Python script
       for executing a reaction OpenTrons protocol

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    otscript: FileField
        The Python OpenTrons protocol file
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="otscripts", on_delete=models.CASCADE
    )
    otscript = models.FileField(upload_to="otscripts/", max_length=255)


class SolventPrep(models.Model):
    """Django model to define a SolventPrep - a csv
       file detailing the prepreation required for diluting previous
       reaction step products for use in the next reaction

    Parameters
    ----------
    otsession_id: ForeignKey
        Foreign key linking a plate to an OT session
    solventprepcsv: FileField
        The csv file with solvent amount (uL), plate name and well index
    """

    otsession_id = models.ForeignKey(
        OTSession, related_name="solventpreps", on_delete=models.CASCADE
    )
    solventprepcsv = models.FileField(upload_to="solventprep/", max_length=255)
