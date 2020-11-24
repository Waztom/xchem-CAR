from django.db import models

class Project(models.Model):
     # The date it was made
    init_date = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=255, db_index=True, unique=True)
    submittername = models.CharField(max_length=255)
    submitteremail = models.CharField(max_length=255)

class Target(models.Model):
    project_id = models.ForeignKey(Project, on_delete=models.CASCADE)
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to='targetimages/', max_length=255)
    name = models.CharField(max_length=255, db_index=True, unique=True)


class Method(models.Model):
    class Unit(models.TextChoices):
        mmol = 'mmol'
        g = 'g'
        mg = 'mg'

    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    nosteps = models.IntegerField(null=True)
    targetmass = models.IntegerField()
    unit = models.CharField(
        choices=Unit.choices,
        default=Unit.mg,
        max_length=10
    )


class Reaction(models.Model):
    method_id = models.ForeignKey(Method, on_delete=models.CASCADE)
    name = models.CharField(max_length=50) 
    productsmiles =  models.CharField(max_length=255, db_index=True, null=True)
    productimage = models.FileField(upload_to='productimages/', max_length=255)
    productname = models.CharField(max_length=100, unique=True)


class Reactant(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    smiles =  models.CharField(max_length=255, db_index=True, null=True)
    image= models.FileField(upload_to='reactantimages/', max_length=255)
    
# Models to capture actions and sequence

class AddAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    class Unit(models.TextChoices):
        mmol = 'mmol'
        ml = 'ml'
        mg = 'mg'

    material = models.CharField(max_length=255)
    quantity = models.IntegerField(null=True)
    unit = models.CharField(
        choices=Unit.choices,
        default=Unit.mmol,
        max_length=10
    )
    dropwise = models.BooleanField(default=False)


class MakeSolutionAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    class Unit(models.TextChoices):
        molar = 'mol/L'
        mass = 'g/L'

    solute = models.CharField(max_length=100)
    soluteammount = models.IntegerField()
    solvent = models.CharField(max_length=100)
    
    unit = models.CharField(
        choices=Unit.choices,
        default=Unit.molar,
        max_length=10
    )

    concentration = models.IntegerField()


class StirAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    class Unit(models.TextChoices):
        seconds = 's'
        minutes = 'min'
        hours = 'h'
    
    duration = models.IntegerField()
    
    unit = models.CharField(
        choices=Unit.choices,
        default=Unit.minutes,
        max_length=10
    )

class WashAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    class Unit(models.TextChoices):
        mmol = 'mmol'
        ml = 'ml'
        mg = 'mg'
    
    material = models.CharField(max_length=100)
    norepetitions = models.IntegerField() 
    unit = models.CharField(
        choices=Unit.choices,
        default=Unit.ml,
        max_length=10
    )

class DrySolutionAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)   
    dryingagent = models.CharField(max_length=100)


class ConcentrateAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    concentrate=models.BooleanField(default=False)


class AnalyseAction(models.Model):
    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    class QCMethod(models.TextChoices):
        LCMS = 'LCMS'
        NMR = 'NMR'
        XChem = 'XChem'

    method = models.CharField(
        choices=QCMethod.choices,
        default=QCMethod.LCMS,
        max_length=10)








    





#  moleq= models.IntegerField(null=True)
#     addorder=models.IntegerField(null=True)
#     concentration=models.IntegerField(null=True)

#     class Workup(models.TextChoices):
#         NA = 'None'    
#         HCl = 'HCl'
#         DCM = 'DCM' 

#     workup = models.CharField(
#         choices=Workup.choices,
#         default=Workup.NA,
#         max_length=10
#     )

# temperature = models.IntegerField(null=True)

#     class Solvent(models.TextChoices):
#         NA = 'None'
#         DMA = 'DMA'
#         DCM = 'DCM' 

#     solvent = models.CharField(
#         choices=Solvent.choices,
#         default=Solvent.NA,
#         max_length=10
#     )