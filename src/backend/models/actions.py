"""Action-related models: ActionSession, AddAction, ExtractAction, etc."""

from django.db import models

from .core import Reaction


class PlateRole(models.TextChoices):
    """Functional role of a plate on the OT deck.

    The actual plate instance is identified by combining a role with a
    1-based ``role_index`` (default 1).  For example role='workup' +
    role_index=2 replaces the old hardcoded 'workup2'.
    """

    REACTION = "reaction"
    WORKUP = "workup"
    SPEFILTER = "spefilter"
    LCMS = "lcms"
    XCHEM = "xchem"
    NMR = "nmr"
    STARTINGMATERIAL = "startingmaterial"
    SOLVENT = "solvent"
    ANALYSE = "analyse"


class ActionSession(models.Model):
    class Type(models.TextChoices):
        reaction = "reaction"
        stir = "stir"
        workup = "workup"
        analyse = "analyse"

    class Driver(models.TextChoices):
        human = "human"
        robot = "robot"

    reaction_id = models.ForeignKey(
        Reaction, related_name="actionsessions", on_delete=models.CASCADE
    )
    sessionnumber = models.IntegerField()
    type = models.CharField(choices=Type.choices, max_length=10)
    driver = models.CharField(
        choices=Driver.choices, default=Driver.robot, max_length=10
    )
    continuation = models.BooleanField(default=False)


class AddAction(models.Model):
    """Django model to define a AddAction - the add action details

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    fromplatetype: CharField
        The plate the add action is moving material from
    toplatetype: CharField
        The plate the add action is moviong material to
    smiles: CharField
        Optional SMILES of the material being added
    calcunit: CharField
        The unit used for calculating the amount to add eg. molar eq. (moleq)
        and mass eq. (masseq)
    volume: FloatField
        The volume being added
    volumeunit: CharField
        The unit of the volume being added (default=uL)
    molecularweight: FloatField
        The molecular weight of the compound being added
    solvent: CharField
        Optional solvent used to dilute the material being added
    concentration: FloatField
        Optional concentration of the material solution prepared
    """

    class CalcUnit(models.TextChoices):
        moleq = "moleq"
        masseq = "masseq"
        uL = "uL"

    class VolumeUnit(models.TextChoices):
        uL = "uL"
        mL = "mL"

    class MassUnit(models.TextChoices):
        mg = "mg"
        g = "g"

    actionsession_id = models.ForeignKey(ActionSession, on_delete=models.CASCADE)
    reaction_id = models.ForeignKey(
        Reaction, related_name="addactions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    from_plate_role = models.CharField(
        choices=PlateRole.choices, max_length=20, default=PlateRole.STARTINGMATERIAL
    )
    from_plate_role_index = models.PositiveIntegerField(default=1)
    to_plate_role = models.CharField(
        choices=PlateRole.choices, max_length=20, default=PlateRole.REACTION
    )
    to_plate_role_index = models.PositiveIntegerField(default=1)
    smiles = models.CharField(max_length=255)
    calcunit = models.CharField(
        choices=CalcUnit.choices, default="moleq", max_length=10
    )
    volume = models.FloatField(null=True)
    volumeunit = models.CharField(
        choices=VolumeUnit.choices, default="uL", max_length=2
    )
    mass = models.FloatField(null=True)
    massunit = models.CharField(choices=MassUnit.choices, default="mg", max_length=2)

    molecularweight = models.FloatField()
    solvent = models.CharField(max_length=255, null=True)
    concentration = models.FloatField(null=True)


class ExtractAction(models.Model):
    """Django model to define an extract action

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    fromplatetype: CharField
        The plate the add action is moving material from
    toplatetype: CharField
        The plate the add action is moviong material to
    layer: CharField
        The layer to extract
    volume: FloatField
        The volume to extract
    volumeunit: CharField
        The unit of the volume being extracted (default=uL)
    molecularweight: FloatField
        The molecular weight of the compound being added
    bottomlayervolume:
        The volume of the bottom layer (uL)
    bottomlayervolumeunit:
        The otional unit used to calculate the volume of the bottom layer
    solvent: CharField
        Optional solvent used to dilute the material being added
    concentration: FloatField
        Optional concentration of the material solution prepared
    """

    class Layer(models.TextChoices):
        top = "top"
        bottom = "bottom"

    actionsession_id = models.ForeignKey(ActionSession, on_delete=models.CASCADE)
    reaction_id = models.ForeignKey(
        Reaction, related_name="extractactions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    from_plate_role = models.CharField(
        choices=PlateRole.choices, max_length=20, default=PlateRole.REACTION
    )
    from_plate_role_index = models.PositiveIntegerField(default=1)
    to_plate_role = models.CharField(
        choices=PlateRole.choices, max_length=20, default=PlateRole.WORKUP
    )
    to_plate_role_index = models.PositiveIntegerField(default=1)
    layer = models.CharField(choices=Layer.choices, default="bottom", max_length=10)
    smiles = models.CharField(max_length=255)
    volume = models.FloatField()
    volumeunit = models.CharField(default="uL", max_length=2)
    molecularweight = models.FloatField()
    bottomlayervolume = models.FloatField(null=True)
    bottomlayervolumeunit = models.CharField(default="uL", max_length=2)
    solvent = models.CharField(max_length=255, null=True)
    concentration = models.FloatField(null=True)


class MixAction(models.Model):
    """Django model to define a mix action

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    platetype: CharField
        The plate being mixed
    repetitions: IntField
        The number of mixes
    volume: FloatField
        The volume to mix
    """

    actionsession_id = models.ForeignKey(ActionSession, on_delete=models.CASCADE)
    reaction_id = models.ForeignKey(
        Reaction, related_name="mixactions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    plate_role = models.CharField(
        choices=PlateRole.choices, max_length=20, default=PlateRole.REACTION
    )
    plate_role_index = models.PositiveIntegerField(default=1)
    repetitions = models.IntegerField()


class StirAction(models.Model):
    """Django model to define a StirAction - the stir action details

    Parameters
    ----------
    actionsession_id: ForeignKey
        Foreign key linking an add action to an action session. An
        action session if a group of actions that represent a unit of
        operation executed by a robot or human eg. perfrom a reaction
        (liquid handling robot), stir (human)
        on a hot plate, analyse (human)
    reaction_id: ForeignKey
        Foreign key linking an add action to a reaction
    number: IntegerField
        The number of the action to be executed in a list of action numbers
    duration: FloatField
        The duration of the stir action
    durationunit: CharField
        The duration unit of the stir action (default=hours)
    temperature: IntegerField
        The temperature of the stir action (default=25)
    temperatureunit: CharField
        The temperature unit of the stir action (default=degC)
    stirringspeed: CharField
        The speed of the stir action (default=normal)
    """

    class TemperatureUnit(models.TextChoices):
        degcel = "degC"
        kelvin = "K"

    class Unit(models.TextChoices):
        seconds = "s"
        minutes = "m"
        hours = "h"

    class Speed(models.TextChoices):
        gentle = "gentle"
        normal = "normal"
        vigorous = "vigorous"

    actionsession_id = models.ForeignKey(
        ActionSession, related_name="stiractions", on_delete=models.CASCADE
    )
    reaction_id = models.ForeignKey(
        Reaction, related_name="stiractions", on_delete=models.CASCADE
    )
    number = models.IntegerField()
    plate_role = models.CharField(
        choices=PlateRole.choices, max_length=20, default=PlateRole.REACTION
    )
    plate_role_index = models.PositiveIntegerField(default=1)
    duration = models.FloatField()
    durationunit = models.CharField(
        choices=Unit.choices, default=Unit.hours, max_length=10
    )
    temperature = models.IntegerField()
    temperatureunit = models.CharField(
        choices=TemperatureUnit.choices, default=TemperatureUnit.degcel, max_length=10
    )
    stirringspeed = models.CharField(
        choices=Speed.choices, default=Speed.normal, max_length=10
    )
