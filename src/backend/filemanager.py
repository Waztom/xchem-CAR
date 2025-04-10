import logging
from django.db.models import Q
from django.core.files.storage import default_storage
from .models import Target, Product, Reaction

logger = logging.getLogger(__name__)

def delete_file_if_unused(file_path):
    """Delete a file only if no models are using it
    
    Parameters
    ----------
    file_path : str
        The path of the file to check and potentially delete
        
    Returns
    -------
    bool
        True if the file was deleted, False otherwise
    """
    if not file_path:
        return False
        
    # Check if the file is used by any model that might share images
    is_used = (
        Target.objects.filter(image=file_path).exists() or
        Reaction.objects.filter(image=file_path).exists() or
        Product.objects.filter(image=file_path).exists()
    )
    
    if is_used:
        logger.debug(f"File {file_path} is still in use by other models, not deleting")
        return False
    
    # If not used anywhere, delete the file
    try:
        default_storage.delete(file_path)
        logger.info(f"Deleted unused file: {file_path}")
        return True
    except Exception as e:
        logger.error(f"Error deleting unused file {file_path}: {str(e)}")
        return False