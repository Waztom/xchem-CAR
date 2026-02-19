"""Versioned chemistry recipe models."""

from django.db import models

from .actions import PlateRole


class Recipe(models.Model):
    """A versioned chemistry recipe for a reaction class.

    Every save creates a new immutable row.  Edits produce a new row with
    ``parent`` pointing to the original and an auto-incremented ``version``.
    The unique key is ``(reaction_class, name, version)``.

    Parameters
    ----------
    reaction_class : CharField
        e.g. "Amidation", "Suzuki coupling".
    name : CharField
        Human-readable recipe name (e.g. "standard", "high-temp-v1").
    version : PositiveIntegerField
        Auto-incremented within (reaction_class, name).
    parent : ForeignKey (self, nullable)
        The recipe this version was derived from.  NULL for first version.
    created_at : DateTimeField
        Auto-populated creation timestamp.
    created_by : CharField
        Who created this version (optional).
    description : TextField
        Free-text notes.
    intramolecular : BooleanField
        Whether this reaction class supports an intramolecular pathway.
    estimated_yield : PositiveIntegerField
        Estimated yield percentage (optional).
    reaction_smarts : JSONField
        List of reaction SMARTS strings for product validation.
    references : TextField
        Literature or internal references (optional).
    """

    reaction_class = models.CharField(max_length=255, db_index=True)
    name = models.CharField(max_length=255)
    version = models.PositiveIntegerField(default=1)
    parent = models.ForeignKey(
        "self",
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="derived_versions",
    )
    created_at = models.DateTimeField(auto_now_add=True)
    created_by = models.CharField(max_length=255, blank=True, default="")
    description = models.TextField(blank=True, default="")
    intramolecular = models.BooleanField(default=False)
    estimated_yield = models.PositiveIntegerField(null=True, blank=True)
    reaction_smarts = models.JSONField(default=list, blank=True)
    references = models.TextField(null=True, blank=True)

    class Meta:
        unique_together = ["reaction_class", "name", "version"]
        ordering = ["reaction_class", "name", "-version"]

    def __str__(self):
        return f"{self.reaction_class} / {self.name} (v{self.version})"


class RecipeActionSession(models.Model):
    """A session (unit of execution) within a recipe.

    Sessions group actions into logical blocks — reaction, stir, workup,
    analyse.  Ordering is determined by position in the ``sessions`` list
    (stored as ``session_number``).

    Parameters
    ----------
    recipe : ForeignKey → Recipe
        Parent recipe.
    session_number : PositiveIntegerField
        1-indexed execution order within the recipe.
    session_type : CharField
        One of reaction | stir | workup | workup1 | workup2 | analyse.
    driver : CharField
        One of robot | human.
    continuation : BooleanField
        Whether this session continues from a previous session
        (e.g. add amine after pre-activation stir).
    """

    class SessionType(models.TextChoices):
        REACTION = "reaction"
        STIR = "stir"
        WORKUP = "workup"
        WORKUP1 = "workup1"
        WORKUP2 = "workup2"
        ANALYSE = "analyse"

    class SessionDriver(models.TextChoices):
        HUMAN = "human"
        ROBOT = "robot"

    recipe = models.ForeignKey(
        Recipe, related_name="action_sessions", on_delete=models.CASCADE
    )
    session_number = models.PositiveIntegerField(default=0)
    session_type = models.CharField(choices=SessionType.choices, max_length=20)
    driver = models.CharField(
        choices=SessionDriver.choices,
        default=SessionDriver.ROBOT,
        max_length=10,
    )
    continuation = models.BooleanField(default=False)

    class Meta:
        ordering = ["session_number"]
        unique_together = ["recipe", "session_number"]

    def save(self, *args, **kwargs):
        if not self.session_number:
            last = (
                RecipeActionSession.objects.filter(recipe=self.recipe)
                .order_by("-session_number")
                .values_list("session_number", flat=True)
                .first()
            )
            self.session_number = (last or 0) + 1
        super().save(*args, **kwargs)

    def __str__(self):
        return (
            f"Session {self.session_number} ({self.session_type}) "
            f"of {self.recipe}"
        )


class RecipeAddAction(models.Model):
    """An add-material step in a recipe session.

    Defines what material to add, how much, from/to which plate.
    The ``molecular_context`` field controls whether this action applies
    to the intermolecular pathway or intramolecular pathway.  Only relevant
    for reaction-type sessions where some recipes define different add
    actions depending on whether the specific reaction instance is
    intramolecular (single reactant) or intermolecular (two reactants).
    The consumer filters on ``Reaction.intramolecular`` at runtime.

    When ``molecular_context`` is NULL the action applies to **all**
    pathways (typical for workup / non-reaction sessions).

    Parameters
    ----------
    session : ForeignKey → RecipeActionSession
    action_number : PositiveIntegerField
        Execution order within the session (auto-assigned from list position).
    molecular_context : CharField or None
        intermolecular | intramolecular | NULL (applies to all pathways).
        Only relevant for reaction sessions.
    material_smarts : CharField
        SMARTS pattern to match a reactant by substructure.
    material_smiles : CharField
        Explicit SMILES for reagents, catalysts, solvents.
    equivalents : FloatField
        Amount value — meaning depends on quantity_unit.
    quantity_unit : CharField
        moleq | uL.
    solvent : CharField
        Solvent used to prepare the material solution.
    concentration : FloatField
        Concentration of the prepared solution (mol/L).
    density : FloatField
        Neat density (g/mL) for liquid reagents dispensed without solvent.
    from_plate_role : CharField
        Functional role of the source plate (default: startingmaterial).
    from_plate_role_index : PositiveIntegerField
        Instance number of the source plate (default: 1).
    to_plate_role : CharField
        Functional role of the destination plate (default: reaction).
    to_plate_role_index : PositiveIntegerField
        Instance number of the destination plate (default: 1).
    """

    class MolecularContext(models.TextChoices):
        INTERMOLECULAR = "intermolecular"
        INTRAMOLECULAR = "intramolecular"

    class QuantityUnit(models.TextChoices):
        MOLEQ = "moleq"
        MASSEQ = "masseq"
        UL = "uL"
        ML = "mL"
        MG = "mg"
        G = "g"
        MOLARITY = "M"
        MICROMOLARITY = "uM"

    session = models.ForeignKey(
        RecipeActionSession, related_name="add_actions", on_delete=models.CASCADE
    )
    action_number = models.PositiveIntegerField(default=0)
    molecular_context = models.CharField(
        choices=MolecularContext.choices,
        max_length=20,
        null=True,
        blank=True,
    )
    material_smarts = models.CharField(max_length=500, null=True, blank=True)
    material_smiles = models.CharField(max_length=500, null=True, blank=True)
    equivalents = models.FloatField()
    quantity_unit = models.CharField(
        choices=QuantityUnit.choices,
        default=QuantityUnit.MOLEQ,
        max_length=10,
    )
    solvent = models.CharField(max_length=255, null=True, blank=True)
    concentration = models.FloatField(null=True, blank=True)
    density = models.FloatField(null=True, blank=True)
    from_plate_role = models.CharField(
        choices=PlateRole.choices,
        default=PlateRole.STARTINGMATERIAL,
        max_length=20,
    )
    from_plate_role_index = models.PositiveIntegerField(default=1)
    to_plate_role = models.CharField(
        choices=PlateRole.choices,
        default=PlateRole.REACTION,
        max_length=20,
    )
    to_plate_role_index = models.PositiveIntegerField(default=1)

    class Meta:
        ordering = ["action_number"]

    def save(self, *args, **kwargs):
        if not self.action_number:
            from django.db.models import Max
            candidates = [
                self.session.add_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.stir_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.extract_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.mix_actions.aggregate(m=Max("action_number"))["m"] or 0,
            ]
            self.action_number = max(candidates) + 1
        super().save(*args, **kwargs)

    def __str__(self):
        material = self.material_smarts or self.material_smiles or "product"
        return f"Add {material} (#{self.action_number}) in {self.session}"


class RecipeStirAction(models.Model):
    """A stir/incubation step in a recipe session.

    Parameters
    ----------
    session : ForeignKey → RecipeActionSession
    action_number : PositiveIntegerField
    temperature : FloatField
        Target temperature.
    temperature_unit : CharField
        Default degC.
    duration : FloatField
        How long to stir.
    duration_unit : CharField
        Default hours.
    stirring_speed : CharField
        gentle | normal | vigorous.
    plate_role : CharField
        Functional role of the plate being stirred (default: reaction).
    plate_role_index : PositiveIntegerField
        Instance number of the plate (default: 1).
    """

    class StirSpeed(models.TextChoices):
        GENTLE = "gentle"
        NORMAL = "normal"
        VIGOROUS = "vigorous"

    class TemperatureUnit(models.TextChoices):
        DEGCEL = "degC"
        KELVIN = "K"

    class DurationUnit(models.TextChoices):
        SECONDS = "s"
        MINUTES = "m"
        HOURS = "h"

    session = models.ForeignKey(
        RecipeActionSession, related_name="stir_actions", on_delete=models.CASCADE
    )
    action_number = models.PositiveIntegerField(default=0)
    temperature = models.FloatField(default=25)
    temperature_unit = models.CharField(max_length=10, default="degC")
    duration = models.FloatField()
    duration_unit = models.CharField(max_length=10, default="h")
    stirring_speed = models.CharField(
        choices=StirSpeed.choices,
        default=StirSpeed.NORMAL,
        max_length=10,
    )
    plate_role = models.CharField(
        choices=PlateRole.choices,
        default=PlateRole.REACTION,
        max_length=20,
    )
    plate_role_index = models.PositiveIntegerField(default=1)

    class Meta:
        ordering = ["action_number"]

    def save(self, *args, **kwargs):
        if not self.action_number:
            from django.db.models import Max
            candidates = [
                self.session.add_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.stir_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.extract_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.mix_actions.aggregate(m=Max("action_number"))["m"] or 0,
            ]
            self.action_number = max(candidates) + 1
        super().save(*args, **kwargs)

    def __str__(self):
        return (
            f"Stir {self.temperature}{self.temperature_unit} "
            f"{self.duration}{self.duration_unit} on "
            f"{self.plate_role}{self.plate_role_index} "
            f"(#{self.action_number}) in {self.session}"
        )


class RecipeExtractAction(models.Model):
    """A liquid-liquid extraction step in a recipe session.

    Parameters
    ----------
    session : ForeignKey → RecipeActionSession
    action_number : PositiveIntegerField
    layer : CharField
        Which layer to extract (top | bottom).
    volume : FloatField
        Volume to extract (uL).
    bottom_layer_volume : FloatField
        Volume of the bottom layer (uL), used for volume calculations.
    smiles : CharField
        Optional SMILES of the extraction solvent.
    solvent : CharField
        Solvent name.
    concentration : FloatField
        Optional concentration of extracted material.
    from_plate_role : CharField
        Functional role of the source plate (default: reaction).
    from_plate_role_index : PositiveIntegerField
        Instance number of the source plate (default: 1).
    to_plate_role : CharField
        Functional role of the destination plate (default: workup).
    to_plate_role_index : PositiveIntegerField
        Instance number of the destination plate (default: 1).
    """

    class Layer(models.TextChoices):
        TOP = "top"
        BOTTOM = "bottom"

    session = models.ForeignKey(
        RecipeActionSession, related_name="extract_actions", on_delete=models.CASCADE
    )
    action_number = models.PositiveIntegerField(default=0)
    layer = models.CharField(
        choices=Layer.choices,
        default=Layer.BOTTOM,
        max_length=10,
    )
    volume = models.FloatField()
    bottom_layer_volume = models.FloatField(null=True, blank=True)
    smiles = models.CharField(max_length=500, null=True, blank=True)
    solvent = models.CharField(max_length=255, null=True, blank=True)
    concentration = models.FloatField(null=True, blank=True)
    from_plate_role = models.CharField(
        choices=PlateRole.choices,
        default=PlateRole.REACTION,
        max_length=20,
    )
    from_plate_role_index = models.PositiveIntegerField(default=1)
    to_plate_role = models.CharField(
        choices=PlateRole.choices,
        default=PlateRole.WORKUP,
        max_length=20,
    )
    to_plate_role_index = models.PositiveIntegerField(default=1)

    class Meta:
        ordering = ["action_number"]

    def save(self, *args, **kwargs):
        if not self.action_number:
            from django.db.models import Max
            candidates = [
                self.session.add_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.stir_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.extract_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.mix_actions.aggregate(m=Max("action_number"))["m"] or 0,
            ]
            self.action_number = max(candidates) + 1
        super().save(*args, **kwargs)

    def __str__(self):
        return (
            f"Extract {self.layer} layer {self.volume}uL "
            f"{self.from_plate_role}{self.from_plate_role_index}→"
            f"{self.to_plate_role}{self.to_plate_role_index} "
            f"(#{self.action_number}) in {self.session}"
        )


class RecipeMixAction(models.Model):
    """A mixing step in a recipe session.

    Parameters
    ----------
    session : ForeignKey → RecipeActionSession
    action_number : PositiveIntegerField
    plate_role : CharField
        Functional role of the plate to mix (default: reaction).
    plate_role_index : PositiveIntegerField
        Instance number of the plate (default: 1).
    repetitions : PositiveIntegerField
        Number of mix cycles.
    """

    session = models.ForeignKey(
        RecipeActionSession, related_name="mix_actions", on_delete=models.CASCADE
    )
    action_number = models.PositiveIntegerField(default=0)
    plate_role = models.CharField(
        choices=PlateRole.choices,
        default=PlateRole.REACTION,
        max_length=20,
    )
    plate_role_index = models.PositiveIntegerField(default=1)
    repetitions = models.PositiveIntegerField()

    class Meta:
        ordering = ["action_number"]

    def save(self, *args, **kwargs):
        if not self.action_number:
            from django.db.models import Max
            candidates = [
                self.session.add_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.stir_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.extract_actions.aggregate(m=Max("action_number"))["m"] or 0,
                self.session.mix_actions.aggregate(m=Max("action_number"))["m"] or 0,
            ]
            self.action_number = max(candidates) + 1
        super().save(*args, **kwargs)

    def __str__(self):
        return (
            f"Mix {self.repetitions}x on {self.plate_role}{self.plate_role_index} "
            f"(#{self.action_number}) in {self.session}"
        )
