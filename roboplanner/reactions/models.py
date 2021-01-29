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
                self.submittername[0:3]
                + " "
                + self.submitterorganisation[0:3]
                + " "
                + rand_slug()
            )

        super(Project, self).save(*args, **kwargs)


class Target(models.Model):
    class Unit(models.TextChoices):
        mmol = "mmol"
        g = "g"
        mg = "mg"

    project_id = models.ForeignKey(Project, on_delete=models.CASCADE)
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to="targetimages/", max_length=255)
    name = models.CharField(max_length=255, db_index=True, unique=True)
    targetmass = models.IntegerField()
    unit = models.CharField(choices=Unit.choices, default=Unit.mg, max_length=10)


class Method(models.Model):
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    nosteps = models.IntegerField(null=True)


class Reaction(models.Model):
    method_id = models.ForeignKey(Method, on_delete=models.CASCADE)
    reactionclass = models.CharField(max_length=255)


class Product(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    name = models.CharField(max_length=255, unique=True)
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to="productimages/", max_length=255)


class AddAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    name = models.CharField(max_length=255, unique=True)
    actionno = models.IntegerField()
    smiles = models.CharField(max_length=255)
    molequivalents = models.FloatField(default=1, null=True)
    molecularweight = models.FloatField()
    dropwise = models.BooleanField(default=False)
    image = models.FileField(upload_to="addactionimages/", max_length=255)


# Models to capture IBM actions
class MakeSolutionAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)

    class Unit(models.TextChoices):
        mmol = "mmol"
        ml = "ml"

    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    solute = models.CharField(max_length=100)
    solutequantity = models.IntegerField()
    soluteunit = models.CharField(choices=Unit.choices, default=Unit.mmol, max_length=10)

    solvent = models.CharField(max_length=100)
    solventquantity = models.IntegerField()
    solventunit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)


class StirAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)

    class Unit(models.TextChoices):
        seconds = "seconds"

    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    duration = models.IntegerField()
    unit = models.CharField(choices=Unit.choices, default=Unit.seconds, max_length=10)
    temperature = models.IntegerField()


class WashAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)

    class Unit(models.TextChoices):
        ml = "ml"

    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    material = models.CharField(max_length=100)
    norepetitions = models.IntegerField()
    unit = models.CharField(choices=Unit.choices, default=Unit.ml, max_length=10)


class DrySolutionAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    dryingagent = models.CharField(max_length=100)


class ConcentrateAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    concentrate = models.BooleanField(default=False)


class AnalyseAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)

    class QCMethod(models.TextChoices):
        LCMS = "LCMS"
        NMR = "NMR"
        XChem = "XChem"

    actiontype = models.CharField(max_length=100)
    actionno = models.IntegerField()
    method = models.CharField(
        choices=QCMethod.choices, default=QCMethod.LCMS, max_length=10
    )

