from django.db import models

class Project(models.Model):
    name = models.CharField(max_length=255, db_index=True, unique=True)


class Target(models.Model):
    project_id = models.ForeignKey(Project, on_delete=models.CASCADE)
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to='targetimages/', max_length=255)
    name = models.CharField(max_length=255, db_index=True, unique=True)


class Method(models.Model):
    target_id = models.ForeignKey(Target, on_delete=models.CASCADE)
    nosteps = models.IntegerField(null=True)


class Reaction(models.Model):
    class Workup(models.TextChoices):
        NA = 'None'    
        HCl = 'HCl'
        DCM = 'DCM' 

    workup = models.CharField(
        choices=Workup.choices,
        default=Workup.NA,
        max_length=10
    )

    method_id = models.ForeignKey(Method, on_delete=models.CASCADE)
    name = models.CharField(max_length=50)
    temperature = models.IntegerField(null=True)
    productsmiles =  models.CharField(max_length=255, db_index=True, null=True)
    productimage = models.FileField(upload_to='productimages/', max_length=255)
    productname = models.CharField(max_length=100, unique=True)


class Reactant(models.Model):
    class Solvent(models.TextChoices):
        NA = 'None'
        DMA = 'DMA'
        DCM = 'DCM' 

    solvent = models.CharField(
        choices=Solvent.choices,
        default=Solvent.NA,
        max_length=10
    )

    reaction_id = models.ForeignKey(Reaction, on_delete=models.CASCADE)
    smiles =  models.CharField(max_length=255, db_index=True, null=True)
    image= models.FileField(upload_to='reactantimages/', max_length=255)
    moleq= models.IntegerField(null=True)
    addorder=models.IntegerField(null=True)
    concentration=models.IntegerField(null=True)
    





