"""Tests for backend.signals module.

Covers all 7 ``post_delete`` signal handlers:
  - CompoundOrder  → deletes ordercsv file via storage
  - OTScript       → deletes otscript file via storage
  - SolventPrep    → deletes solventprepcsv file via storage
  - OTBatchProtocol→ deletes zipfile via storage
  - Product        → delegates to delete_file_if_unused
  - Reaction       → delegates to delete_file_if_unused
  - Target         → delegates to delete_file_if_unused

Tests verify:
  - Files are deleted on model deletion when file fields are populated
  - No deletion attempted when file fields are empty/None
  - Storage errors are caught and logged (not raised)
  - delete_file_if_unused is called with correct path for shared-image models
"""

from unittest import TestCase
from unittest.mock import patch, MagicMock, PropertyMock

from backend.signals import (
    delete_compound_order_files,
    delete_otscript_files,
    delete_solventprep_files,
    delete_otbatchprotocol_files,
    cleanup_product_image,
    cleanup_reaction_image,
    cleanup_target_image,
)


# ---------------------------------------------------------------------------
# CompoundOrder signal
# ---------------------------------------------------------------------------
class TestDeleteCompoundOrderFiles(TestCase):
    """Tests for delete_compound_order_files signal handler."""

    def test_deletes_file_when_ordercsv_present(self):
        """Should call storage.delete(path) when ordercsv is populated."""
        storage = MagicMock()
        instance = MagicMock()
        instance.ordercsv = MagicMock()
        instance.ordercsv.storage = storage
        instance.ordercsv.path = "compoundorders/test.csv"
        # bool(instance.ordercsv) must be True
        instance.ordercsv.__bool__ = lambda self: True

        delete_compound_order_files(sender=None, instance=instance)

        storage.delete.assert_called_once_with("compoundorders/test.csv")

    def test_no_delete_when_ordercsv_empty(self):
        """Should not attempt deletion when ordercsv is falsy."""
        instance = MagicMock()
        instance.ordercsv = None

        # Should not raise
        delete_compound_order_files(sender=None, instance=instance)

    def test_no_delete_when_ordercsv_blank(self):
        """Should not attempt deletion when ordercsv is empty string."""
        instance = MagicMock()
        instance.ordercsv = ""

        delete_compound_order_files(sender=None, instance=instance)

    def test_storage_error_is_caught(self):
        """Should log error but not raise when storage.delete fails."""
        storage = MagicMock()
        storage.delete.side_effect = OSError("disk full")
        instance = MagicMock()
        instance.ordercsv = MagicMock()
        instance.ordercsv.storage = storage
        instance.ordercsv.path = "compoundorders/fail.csv"
        instance.ordercsv.__bool__ = lambda self: True

        # Should not raise
        delete_compound_order_files(sender=None, instance=instance)
        storage.delete.assert_called_once()


# ---------------------------------------------------------------------------
# OTScript signal
# ---------------------------------------------------------------------------
class TestDeleteOTScriptFiles(TestCase):
    """Tests for delete_otscript_files signal handler."""

    def test_deletes_file_when_otscript_present(self):
        storage = MagicMock()
        instance = MagicMock()
        instance.otscript = MagicMock()
        instance.otscript.storage = storage
        instance.otscript.path = "otscripts/protocol.py"
        instance.otscript.__bool__ = lambda self: True

        delete_otscript_files(sender=None, instance=instance)

        storage.delete.assert_called_once_with("otscripts/protocol.py")

    def test_no_delete_when_otscript_empty(self):
        instance = MagicMock()
        instance.otscript = None

        delete_otscript_files(sender=None, instance=instance)

    def test_storage_error_is_caught(self):
        storage = MagicMock()
        storage.delete.side_effect = PermissionError("no access")
        instance = MagicMock()
        instance.otscript = MagicMock()
        instance.otscript.storage = storage
        instance.otscript.path = "otscripts/fail.py"
        instance.otscript.__bool__ = lambda self: True

        delete_otscript_files(sender=None, instance=instance)
        storage.delete.assert_called_once()


# ---------------------------------------------------------------------------
# SolventPrep signal
# ---------------------------------------------------------------------------
class TestDeleteSolventPrepFiles(TestCase):
    """Tests for delete_solventprep_files signal handler."""

    def test_deletes_file_when_solventprepcsv_present(self):
        storage = MagicMock()
        instance = MagicMock()
        instance.solventprepcsv = MagicMock()
        instance.solventprepcsv.storage = storage
        instance.solventprepcsv.path = "solventprep/prep.csv"
        instance.solventprepcsv.__bool__ = lambda self: True

        delete_solventprep_files(sender=None, instance=instance)

        storage.delete.assert_called_once_with("solventprep/prep.csv")

    def test_no_delete_when_solventprepcsv_empty(self):
        instance = MagicMock()
        instance.solventprepcsv = None

        delete_solventprep_files(sender=None, instance=instance)

    def test_storage_error_is_caught(self):
        storage = MagicMock()
        storage.delete.side_effect = OSError("IO error")
        instance = MagicMock()
        instance.solventprepcsv = MagicMock()
        instance.solventprepcsv.storage = storage
        instance.solventprepcsv.path = "solventprep/fail.csv"
        instance.solventprepcsv.__bool__ = lambda self: True

        delete_solventprep_files(sender=None, instance=instance)
        storage.delete.assert_called_once()


# ---------------------------------------------------------------------------
# OTBatchProtocol signal
# ---------------------------------------------------------------------------
class TestDeleteOTBatchProtocolFiles(TestCase):
    """Tests for delete_otbatchprotocol_files signal handler."""

    def test_deletes_file_when_zipfile_present(self):
        storage = MagicMock()
        instance = MagicMock()
        instance.zipfile = MagicMock()
        instance.zipfile.storage = storage
        instance.zipfile.path = "otbatchprotocols/batch.zip"
        instance.zipfile.__bool__ = lambda self: True

        delete_otbatchprotocol_files(sender=None, instance=instance)

        storage.delete.assert_called_once_with("otbatchprotocols/batch.zip")

    def test_no_delete_when_zipfile_empty(self):
        instance = MagicMock()
        instance.zipfile = None

        delete_otbatchprotocol_files(sender=None, instance=instance)

    def test_no_delete_when_zipfile_blank(self):
        instance = MagicMock()
        instance.zipfile = ""

        delete_otbatchprotocol_files(sender=None, instance=instance)

    def test_storage_error_is_caught(self):
        storage = MagicMock()
        storage.delete.side_effect = FileNotFoundError("already gone")
        instance = MagicMock()
        instance.zipfile = MagicMock()
        instance.zipfile.storage = storage
        instance.zipfile.path = "otbatchprotocols/gone.zip"
        instance.zipfile.__bool__ = lambda self: True

        delete_otbatchprotocol_files(sender=None, instance=instance)
        storage.delete.assert_called_once()


# ---------------------------------------------------------------------------
# Product image cleanup signal
# ---------------------------------------------------------------------------
class TestCleanupProductImage(TestCase):
    """Tests for cleanup_product_image signal handler."""

    @patch("backend.signals.delete_file_if_unused")
    def test_calls_delete_file_if_unused_with_image_name(self, mock_delete):
        instance = MagicMock()
        instance.image = MagicMock()
        instance.image.name = "productimages/product.svg"
        instance.image.__bool__ = lambda self: True

        cleanup_product_image(sender=None, instance=instance)

        mock_delete.assert_called_once_with("productimages/product.svg")

    @patch("backend.signals.delete_file_if_unused")
    def test_no_call_when_image_is_falsy(self, mock_delete):
        instance = MagicMock()
        instance.image = MagicMock()
        instance.image.__bool__ = lambda self: False

        cleanup_product_image(sender=None, instance=instance)

        mock_delete.assert_not_called()

    @patch("backend.signals.delete_file_if_unused")
    def test_no_call_when_no_image_attr(self, mock_delete):
        """Should not call delete when instance lacks image attribute."""
        instance = MagicMock(spec=[])  # no attributes at all

        cleanup_product_image(sender=None, instance=instance)

        mock_delete.assert_not_called()


# ---------------------------------------------------------------------------
# Reaction image cleanup signal
# ---------------------------------------------------------------------------
class TestCleanupReactionImage(TestCase):
    """Tests for cleanup_reaction_image signal handler."""

    @patch("backend.signals.delete_file_if_unused")
    def test_calls_delete_file_if_unused_with_image_name(self, mock_delete):
        instance = MagicMock()
        instance.image = MagicMock()
        instance.image.name = "reactionimages/rxn.svg"
        instance.image.__bool__ = lambda self: True

        cleanup_reaction_image(sender=None, instance=instance)

        mock_delete.assert_called_once_with("reactionimages/rxn.svg")

    @patch("backend.signals.delete_file_if_unused")
    def test_no_call_when_image_is_falsy(self, mock_delete):
        instance = MagicMock()
        instance.image = MagicMock()
        instance.image.__bool__ = lambda self: False

        cleanup_reaction_image(sender=None, instance=instance)

        mock_delete.assert_not_called()

    @patch("backend.signals.delete_file_if_unused")
    def test_no_call_when_no_image_attr(self, mock_delete):
        instance = MagicMock(spec=[])

        cleanup_reaction_image(sender=None, instance=instance)

        mock_delete.assert_not_called()


# ---------------------------------------------------------------------------
# Target image cleanup signal
# ---------------------------------------------------------------------------
class TestCleanupTargetImage(TestCase):
    """Tests for cleanup_target_image signal handler."""

    @patch("backend.signals.delete_file_if_unused")
    def test_calls_delete_file_if_unused_with_image_name(self, mock_delete):
        instance = MagicMock()
        instance.image = MagicMock()
        instance.image.name = "targetimages/mol.svg"
        instance.image.__bool__ = lambda self: True

        cleanup_target_image(sender=None, instance=instance)

        mock_delete.assert_called_once_with("targetimages/mol.svg")

    @patch("backend.signals.delete_file_if_unused")
    def test_no_call_when_image_is_falsy(self, mock_delete):
        instance = MagicMock()
        instance.image = MagicMock()
        instance.image.__bool__ = lambda self: False

        cleanup_target_image(sender=None, instance=instance)

        mock_delete.assert_not_called()

    @patch("backend.signals.delete_file_if_unused")
    def test_no_call_when_no_image_attr(self, mock_delete):
        instance = MagicMock(spec=[])

        cleanup_target_image(sender=None, instance=instance)

        mock_delete.assert_not_called()
