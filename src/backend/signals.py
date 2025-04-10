from django.db.models.signals import post_delete
from django.dispatch import receiver
import logging
from .models import Product, Reaction, Target
from .models import CompoundOrder, OTScript, SolventPrep, OTBatchProtocol
from .filemanager import delete_file_if_unused

logger = logging.getLogger(__name__)


@receiver(post_delete, sender=CompoundOrder)
def delete_compound_order_files(sender, instance, **kwargs):
    """Delete the file when CompoundOrder is deleted"""
    if instance.ordercsv:
        try:
            storage, path = instance.ordercsv.storage, instance.ordercsv.path
            storage.delete(path)
            logger.info(f"Deleted file {path}")
        except Exception as e:
            logger.error(f"Error deleting file: {str(e)}")


@receiver(post_delete, sender=OTScript)
def delete_otscript_files(sender, instance, **kwargs):
    """Delete the file when OTScript is deleted"""
    if instance.otscript:
        try:
            storage, path = instance.otscript.storage, instance.otscript.path
            storage.delete(path)
            logger.info(f"Deleted file {path}")
        except Exception as e:
            logger.error(f"Error deleting file: {str(e)}")


@receiver(post_delete, sender=SolventPrep)
def delete_solventprep_files(sender, instance, **kwargs):
    """Delete the file when SolventPrep is deleted"""
    if instance.solventprepcsv:
        try:
            storage, path = (
                instance.solventprepcsv.storage,
                instance.solventprepcsv.path,
            )
            storage.delete(path)
            logger.info(f"Deleted file {path}")
        except Exception as e:
            logger.error(f"Error deleting file: {str(e)}")


@receiver(post_delete, sender=OTBatchProtocol)
def delete_otbatchprotocol_files(sender, instance, **kwargs):
    """Delete the file when OTBatchProtocol is deleted"""
    if instance.zipfile:
        try:
            storage, path = instance.zipfile.storage, instance.zipfile.path
            storage.delete(path)
            logger.info(f"Deleted file {path}")
        except Exception as e:
            logger.error(f"Error deleting file: {str(e)}")


@receiver(post_delete, sender=Product)
def cleanup_product_image(sender, instance, **kwargs):
    """Check and delete product image if no longer used"""
    if hasattr(instance, "image") and instance.image:
        delete_file_if_unused(instance.image.name)


@receiver(post_delete, sender=Reaction)
def cleanup_reaction_image(sender, instance, **kwargs):
    """Check and delete reaction image if no longer used"""
    if hasattr(instance, "image") and instance.image:
        delete_file_if_unused(instance.image.name)


@receiver(post_delete, sender=Target)
def cleanup_target_image(sender, instance, **kwargs):
    """Check and delete target image if no longer used"""
    if hasattr(instance, "image") and instance.image:
        delete_file_if_unused(instance.image.name)
