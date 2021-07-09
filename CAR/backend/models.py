from django.db import models
from django.utils.text import slugify
import random
import string


def rand_slug():
    return "".join(random.choice(string.ascii_letters + string.digits) for _ in range(6))


class Project(models.Model):
    # The date it was made
    init_date = models.DateTimeField(auto_now_add=True)
    name = models.SlugField(max_length=100, db_index=True, unique=True)
    submitterorganisation = models.CharField(max_length=100)
    submittername = models.CharField(max_length=255)
    submitteremail = models.CharField(max_length=255)

    def save(self, *args, **kwargs):
        if not self.name:
            # Newly created object, so set slug
            self.name = slugify(
                self.submittername[0:3] + " " + self.submitterorganisation[0:3] + " " + rand_slug()
            )

        super(Project, self).save(*args, **kwargs)


class Target(models.Model):
    class Unit(models.TextChoices):
        mmol = "mmol"
        g = "g"
        mg = "mg"

    status = models.BooleanField(default=True)
    project_id = models.ForeignKey(Project, on_delete=models.CASCADE)
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to="targetimages/", max_length=255)
    name = models.CharField(max_length=255, db_index=True, unique=True)
    targetmass = models.FloatField()
    targetmols = models.FloatField()
    unit = models.CharField(choices=Unit.choices, default=Unit.mg, max_length=10)


class Method(models.Model):
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    status = models.BooleanField(default=True)
    nosteps = models.IntegerField(null=True)


class Reaction(models.Model):
    method_id = models.ForeignKey(Method, on_delete=models.CASCADE)
    reactionclass = models.CharField(max_length=255)
    reactionimage = models.FileField(
        upload_to="reactionimages/",
        max_length=255,
        null=True,
    )


class Product(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    name = models.CharField(max_length=255, unique=True)
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to="productimages/", max_length=255)


class AnalyseAction(models.Model):
    class QCMethod(models.TextChoices):
        LCMS = "LCMS"
        NMR = "NMR"
        XChem = "XChem"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    method = models.CharField(choices=QCMethod.choices, default=QCMethod.LCMS, max_length=10)


# Models to capture IBM actions
class IBMAddAction(models.Model):
    class Unit(models.TextChoices):
        mmol = "mmol"
        ml = "ml"
        ul = "ul"
        moleq = "moleq"

    class Atmosphere(models.TextChoices):
        nitrogen = "nitrogen"
        air = "air"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    additionorder = models.IntegerField(null=True)
    material = models.CharField(max_length=255)
    materialsmiles = models.CharField(max_length=255, null=True)
    materialquantity = models.FloatField()
    materialquantityunit = models.CharField(
        choices=Unit.choices,
        default=Unit.ul,
        max_length=10,
    )

    dropwise = models.BooleanField(default=False)
    atmosphere = models.CharField(choices=Atmosphere.choices, default=Atmosphere.air, max_length=10)

    # These are extras for helping robotic execution/calcs
    molecularweight = models.FloatField(null=True)
    materialimage = models.FileField(
        upload_to="addactionimages/",
        max_length=255,
        null=True,
    )


class IBMCollectLayerAction(models.Model):
    class Layer(models.TextChoices):
        organic = "organic"
        aqueous = "aqueous"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    layer = models.CharField(choices=Layer.choices, default=Layer.organic, max_length=15)


class IBMConcentrateAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()


class IBMDegasAction(models.Model):
    class Unit(models.TextChoices):
        seconds = "seconds"
        minutes = "minutes"
        hours = "hours"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    gas = models.CharField(max_length=100)
    duration = models.FloatField(null=True)
    durationunit = models.CharField(choices=Unit.choices, default=Unit.seconds, max_length=10)


class IBMDrySolidAction(models.Model):
    class Unit(models.TextChoices):
        seconds = "seconds"
        minutes = "minutes"
        hours = "hours"

    class Atmosphere(models.TextChoices):
        nitrogen = "nitrogen"
        vacuum = "vacuum"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    temperature = models.IntegerField(null=True)
    duration = models.FloatField(null=True)
    durationunit = models.CharField(choices=Unit.choices, default=Unit.seconds, max_length=10)
    atmosphere = models.CharField(
        choices=Atmosphere.choices, default=Atmosphere.nitrogen, max_length=10
    )


class IBMDrySolutionAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    dryingagent = models.CharField(max_length=100)


class IBMExtractAction(models.Model):
    class Unit(models.TextChoices):
        ml = "ml"
        mmol = "mmol"
        moleq = "moleq"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    solvent = models.CharField(max_length=100)
    solventquantity = models.FloatField(null=True)
    solventquantityunit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)
    numberofrepetitions = models.IntegerField(null=True)


class IBMFilterAction(models.Model):
    class PhaseToKeep(models.TextChoices):
        filtrate = "filtrate"
        precipitate = "precipitate"

    class Unit(models.TextChoices):
        ml = "ml"
        mmol = "mmol"
        moleq = "moleq"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    phasetokeep = models.CharField(
        choices=PhaseToKeep.choices, default=PhaseToKeep.filtrate, max_length=20
    )
    rinsingsolvent = models.CharField(max_length=255, null=True)
    rinsingsolventquantity = models.IntegerField(null=True)
    rinsingsolventquantityunit = models.CharField(
        choices=Unit.choices,
        default=Unit.ml,
        max_length=10,
    )

    extractionforprecipitatesolvent = models.CharField(max_length=255, null=True)
    extractionforprecipitatesolventquantity = models.IntegerField(null=True)
    extractionforprecipitatesolventquantityunit = models.CharField(
        choices=Unit.choices,
        default=Unit.ml,
        max_length=10,
    )


class IBMMakeSolutionAction(models.Model):
    class Unit(models.TextChoices):
        ml = "ml"
        mmol = "mmol"
        moleq = "moleq"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    solute = models.CharField(max_length=255)
    solutesmiles = models.CharField(max_length=255, null=True)
    soluteimage = models.FileField(upload_to="IBMmakesolnimages/", max_length=255, null=True)
    solvent = models.CharField(max_length=255)
    solventsmiles = models.CharField(max_length=255, null=True)
    solventimage = models.FileField(upload_to="IBMmakesolnimages/", max_length=255, null=True)
    solutequantity = models.FloatField()
    solutequantityunit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)
    solventquantity = models.FloatField()
    solventquantityunit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)


class IBMPartitionAction(models.Model):
    class Unit(models.TextChoices):
        ml = "ml"
        mmol = "mmol"
        moleq = "moleq"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()

    firstpartitionsolvent = models.CharField(max_length=255)
    secondpartitionsolvent = models.CharField(max_length=255)
    firstpartitionsolventquantity = models.IntegerField()
    secondpartitionsolventquantity = models.IntegerField()
    firstpartitionsolventquantityunit = models.CharField(
        choices=Unit.choices, default=Unit.ml, max_length=10
    )
    secondpartitionsolventquantityunit = models.CharField(
        choices=Unit.choices, default=Unit.ml, max_length=10
    )


class IBMpHAction(models.Model):
    class Unit(models.TextChoices):
        ml = "ml"
        mmol = "mmol"
        moleq = "moleq"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    material = models.CharField(max_length=100)
    materialquantity = models.FloatField(null=True)
    materialquantityunit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)
    pH = models.FloatField(null=True)
    dropwise = models.BooleanField(default=False)
    temperature = models.FloatField(null=True)


class IBMPhaseSeparationAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()


class IBMQuenchAction(models.Model):
    class Unit(models.TextChoices):
        ml = "ml"
        mmol = "mmol"
        moleq = "moleq"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()

    material = models.CharField(max_length=255)
    materialquantity = models.FloatField(null=True)
    materialquantityunit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)
    dropwise = models.BooleanField(default=False)
    temperature = models.IntegerField(null=True)


class IBMRefluxAction(models.Model):
    class Unit(models.TextChoices):
        seconds = "seconds"
        minutes = "minutes"
        hours = "hours"

    class Speed(models.TextChoices):
        gentle = "gentle"
        normal = "normal"
        vigorous = "vigorous"

    class Atmosphere(models.TextChoices):
        nitrogen = "nitrogen"
        air = "air"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    duration = models.FloatField(null=True)
    durationunit = models.CharField(choices=Unit.choices, default=Unit.seconds, max_length=10)
    stirringspeed = models.CharField(choices=Speed.choices, default=Speed.normal, max_length=10)
    deanstarkapparatus = models.BooleanField(default=False)
    atmosphere = models.CharField(choices=Atmosphere.choices, default=Atmosphere.air, max_length=10)


class IBMSetTemperatureAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    temperature = models.IntegerField()


class IBMStirAction(models.Model):
    class Unit(models.TextChoices):
        seconds = "seconds"
        minutes = "minutes"
        hours = "hours"

    class Speed(models.TextChoices):
        gentle = "gentle"
        normal = "normal"
        vigorous = "vigorous"

    class Atmosphere(models.TextChoices):
        nitrogen = "nitrogen"
        air = "air"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    duration = models.FloatField(null=True)
    durationunit = models.CharField(choices=Unit.choices, default=Unit.hours, max_length=10)
    temperature = models.IntegerField(null=True)
    stirringspeed = models.CharField(choices=Speed.choices, default=Speed.normal, max_length=10)
    atmosphere = models.CharField(choices=Atmosphere.choices, default=Atmosphere.air, max_length=10)


class IBMStoreAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    material = models.CharField(max_length=255)
    materialsmiles = models.CharField(max_length=255, null=True)
    materialimage = models.FileField(upload_to="storeimages/", max_length=255, null=True)


class IBMWaitAction(models.Model):
    class Unit(models.TextChoices):
        seconds = "seconds"
        minutes = "minutes"
        hours = "hours"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    duration = models.FloatField(null=True)
    durationunit = models.CharField(choices=Unit.choices, default=Unit.seconds, max_length=10)
    temperature = models.IntegerField(null=True)


class IBMWashAction(models.Model):
    class Unit(models.TextChoices):
        ml = "ml"

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    material = models.CharField(max_length=255)
    materialquantity = models.FloatField(null=True)
    materialquantityunit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)
    numberofrepetitions = models.IntegerField(null=True)
