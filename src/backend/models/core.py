"""Core domain models: Project, Batch, Target, Method, Reaction, and related."""

from django.db import models


class Project(models.Model):
    """Django model to define a Project - a compound synthesis project.

    Parameters
    ----------
    init_date: DateTimeField
        The date the target was initiated (autofield)
    name: SlugFieldField
        The name of the project created as combination of the first three letters of the
        submitter's name and submitter's organisation
    submittername: Charfield
        The name of the person creating the project
    submitterorganisation: Charfield
        The name of the person's organisation creating the project
    proteintarget: Charfield
        The target protein for the submitted compounds
    quotecost: FloatField
        Optional field for the total cost of the project using MCule as a supplier
    quoteurl: CharField
        Optional field for the MCule quote url
    """

    init_date = models.DateTimeField(auto_now_add=True)
    name = models.SlugField(max_length=100, db_index=True)
    submitterorganisation = models.CharField(max_length=100)
    submittername = models.CharField(max_length=255)
    proteintarget = models.CharField(max_length=100)
    quotedcost = models.FloatField(null=True)
    quoteurl = models.CharField(max_length=255, null=True)

    class Meta:
        ordering = ["-init_date"]


class Batch(models.Model):
    """Django model to define a Batch - a batch of compounds.

    Parameters
    ----------
    project_id: ForeignKey
        Foreign key linking a batch to it's project
    batch_id: ForeignKey
        Optional foreign key linking a sub-batch to it's parent-batch
    batchtag: Charfield
        The name of the batch
    """

    project_id = models.ForeignKey(
        Project, related_name="batches", on_delete=models.CASCADE
    )
    batch_id = models.ForeignKey("Batch", on_delete=models.CASCADE, null=True)
    batchtag = models.CharField(max_length=50)


class Target(models.Model):
    """Django model to define a Target - a target compound for synthesis.

    Parameters
    ----------
    batch_id: ForeignKey
        Foreign key linking a target to it's batch
    smiles: Charfield
        The SMILES of the target compound
    image: FileField
        File link to a stored version of the image file of the target compound
    name: Charfield
        The name of the compound (Must complete how this is made)
    concentration: FloatField
        The target synthesis concnetration (mM) for a target compound
    volume: FloatField
        The target synthesis volume (uL) for a target compound
    mols: FloatField
        The target mols of the compound to be synthesised
    """

    batch_id = models.ForeignKey(
        Batch, related_name="targets", on_delete=models.CASCADE
    )
    smiles = models.CharField(max_length=255, db_index=True)
    image = models.FileField(upload_to="targetimages/", max_length=255)
    name = models.CharField(max_length=255, db_index=True)
    concentration = models.FloatField()
    volume = models.FloatField()
    mass = models.FloatField()
    mols = models.FloatField()


class Method(models.Model):
    """Django model to define a Method - a retrosynthetic pathway of reactions
    for a target compound.

    Parameters
    ----------
    target_id: ForeignKey
        Foreign key linking a method to it's target compound
    nosteps: IntegerField
        The number of reaction steps for a method
    otchem: BooleanField
        Set to True if all the reactions in the method can be executed on the OpenTrons
    """

    target_id = models.ForeignKey(
        Target, related_name="methods", on_delete=models.CASCADE
    )
    nosteps = models.IntegerField()
    otchem = models.BooleanField(default=False)


class Reaction(models.Model):
    """Django model to define a Reaction - the reaction to make a compound.

    Parameters
    ----------
    method_id: ForeignKey
        Foreign key linking a reaction to it's method
    reactionclass: Charfield
        The name of the reaction
    recipe: CharField
        The encoded recipe type used for the reaction
    nnumber: IntegerField
        The number of the reaction eg. reaction 1 of 3
    intramolecular: BooleanField
        Set to True if the reaction is intermolecular
    image: FileField
        File link to a stored version of the image file of the reaction
    success: BooleanField
        The success of the reaction with default success set to True
    """

    method_id = models.ForeignKey(
        Method, related_name="reactions", on_delete=models.CASCADE
    )
    reactionclass = models.CharField(max_length=255)
    recipe = models.CharField(max_length=50, default="standard")
    recipe_id = models.ForeignKey(
        "Recipe",
        null=True,
        blank=True,
        on_delete=models.SET_NULL,
        related_name="reactions",
        help_text="Links to the exact Recipe version used for this reaction",
    )
    groupbycolumn = models.BooleanField(default=True)
    number = models.IntegerField()
    intramolecular = models.BooleanField(default=False)
    temperature = models.IntegerField(default=25)
    image = models.FileField(
        upload_to="reactionimages/",
        max_length=255,
        null=True,
    )
    success = models.BooleanField(default=True)


class PubChemInfo(models.Model):
    """Django model to define PubChemInfo - the PubChem info for a compound.

    Parameters
    ----------
    compoundid: IntegerField
        The PubChem compound id of a compound
    summaryurl: Charfield
        The PubChem compound url
    lcssurl: Charfield
        The PubChem Laboratory Chemical Safety Summary url
    smiles: Charfield
        The SMILES of the compound
    cas: Charfield
        Optional CAS number (if found) for a compound
    """

    compoundid = models.IntegerField()
    summaryurl = models.CharField(max_length=255)
    lcssurl = models.CharField(max_length=255)
    smiles = models.CharField(max_length=255)
    cas = models.CharField(max_length=50, null=True)


class Reactant(models.Model):
    """Django model to define a Reactant - the reactant compound of a reaction.

    Parameters
    ----------
    reaction_id: ForeignKey
        Optional foreign key linking a reactant to it's reaction
    pubcheminfo_id: ForeignKey
        Optional foreign key linking a reactant to it's pubchem info (if found)
    smiles: Charfield
        The SMILES of the reactant compound
    """

    reaction_id = models.ForeignKey(
        Reaction, related_name="reactants", on_delete=models.CASCADE, null=True
    )
    pubcheminfo_id = models.ForeignKey(
        PubChemInfo,
        related_name="reactantpubcheminfo",
        on_delete=models.PROTECT,
        null=True,
    )
    smiles = models.CharField(max_length=255)
    previousreactionproduct = models.BooleanField(default=False)


class Product(models.Model):
    """Django model to define a Product- the product compound of a reaction.

    Parameters
    ----------
    reaction_id: ForeignKey
        Foreign key linking a product to it's reaction
    pubcheminfo_id: ForeignKey
        Optional foreign key linking a product to it's pubchem info (if found)
    smiles: Charfield
        The SMILES of the product compound
    image: FileField
        File link to a stored version of the image file of the product compound
    """

    reaction_id = models.ForeignKey(
        Reaction, related_name="products", on_delete=models.CASCADE
    )
    pubcheminfo_id = models.ForeignKey(
        PubChemInfo,
        related_name="productpubcheminfo",
        on_delete=models.PROTECT,
        null=True,
    )
    smiles = models.CharField(max_length=255, db_index=True, null=True)
    image = models.FileField(upload_to="productimages/", max_length=255)


class CatalogEntry(models.Model):
    """Django model to define a CatalogEntry- the catalog information for a compound.

    Parameters
    ----------
    reactant_id: ForeignKey
        Optional Foreign key linking a catalog entry to a reactant
    target_id: ForeignKey
        Optional Foreign key linking a catalog entry to a product
    vendor: Charfield
        The vendor/supplier of the compound
    catalogid: Charfield
        The catalog/vendor id of the compound
    priceinfo: Charfield
        The catalog price info ($) for a compound. This can be a range -> "$100 < 1k/g"
    upperprice: IntegerField
        Optional upper price ($) of a compound. Highest price if range given.
    leadtime: IntegerField
        Optional lead time (weeks) for a compound to be delivered
    """

    reactant_id = models.ForeignKey(
        Reactant, related_name="catalogentries", on_delete=models.CASCADE, null=True
    )
    target_id = models.ForeignKey(
        Target, related_name="catalogentries", on_delete=models.CASCADE, null=True
    )

    vendor = models.CharField(max_length=100)
    catalogid = models.CharField(max_length=50)
    priceinfo = models.CharField(max_length=50)
    upperprice = models.IntegerField(null=True)
    leadtime = models.IntegerField(null=True)
