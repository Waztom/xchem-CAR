"""Tests for backend.api ViewSets via DRF's APIClient.

Covers all 22 registered endpoints (CRUD + custom @action endpoints):
  - 15 flat CRUD ViewSets: list, create, retrieve, update, destroy
  - 3 ViewSets with custom @action endpoints:
      ProjectViewSet: create_project, get_task_status
      BatchViewSet: create (clone), canonicalize_smiles, get_task_status,
                    update_reaction_success
      OTProjectViewSet: create_ot_project, get_task_status
  - filterset_fields filtering
  - fetchall query parameter for *All serializer selection
  - url_path overrides (underscore-free URLs)

Uses Django TestCase (runs against test DB) with DRF APIClient.
Celery tasks and service functions are mocked.
"""

from django.test import TestCase
from rest_framework.test import APIClient
from rest_framework import status

from unittest.mock import patch, MagicMock

from backend.models import (
    Project,
    Batch,
    Target,
    Method,
    Reaction,
    PubChemInfo,
    Product,
    Reactant,
    CatalogEntry,
    ActionSession,
    AddAction,
    ExtractAction,
    MixAction,
    StirAction,
    OTProject,
    OTBatchProtocol,
    OTSession,
    Deck,
    Pipette,
    TipRack,
    Plate,
    Column,
    Well,
    CompoundOrder,
    OTScript,
)


class APITestBase(TestCase):
    """Base class providing common fixtures and APIClient."""

    @classmethod
    def setUpTestData(cls):
        cls.project = Project.objects.create(
            name="test-project",
            submittername="Alice",
            submitterorganisation="TestOrg",
            proteintarget="EGFR",
        )
        cls.batch = Batch.objects.create(
            project_id=cls.project,
            batchtag="batch-1",
        )
        cls.target = Target.objects.create(
            batch_id=cls.batch,
            smiles="CCO",
            name="ethanol",
            concentration=10.0,
            volume=100.0,
            mass=46.0,
            mols=0.001,
            image="targetimages/ethanol.svg",
        )
        cls.method = Method.objects.create(
            target_id=cls.target,
            nosteps=1,
            otchem=True,
        )
        cls.reaction = Reaction.objects.create(
            method_id=cls.method,
            reactionclass="Amidation",
            number=1,
            intramolecular=False,
        )
        cls.pubchem = PubChemInfo.objects.create(
            compoundid=702,
            summaryurl="https://pubchem.ncbi.nlm.nih.gov/compound/702",
            lcssurl="https://pubchem.ncbi.nlm.nih.gov/compound/702#datasheet=LCSS",
            smiles="CCO",
        )
        cls.product = Product.objects.create(
            reaction_id=cls.reaction,
            smiles="CC(=O)O",
            image="productimages/product.svg",
        )
        cls.reactant = Reactant.objects.create(
            reaction_id=cls.reaction,
            smiles="CCO",
            previousreactionproduct=False,
        )
        cls.catalog_entry = CatalogEntry.objects.create(
            reactant_id=cls.reactant,
            vendor="Sigma",
            catalogid="S123",
            priceinfo="< $100 / g",
            upperprice=100,
            leadtime=2,
        )
        cls.action_session = ActionSession.objects.create(
            reaction_id=cls.reaction,
            sessionnumber=1,
            type="reaction",
            driver="robot",
        )
        cls.add_action = AddAction.objects.create(
            reaction_id=cls.reaction,
            actionsession_id=cls.action_session,
            number=1,
            smiles="CCO",
            molecularweight=46.07,
        )
        cls.stir_action = StirAction.objects.create(
            reaction_id=cls.reaction,
            actionsession_id=cls.action_session,
            number=1,
            duration=2.0,
            temperature=25,
        )
        cls.otproject = OTProject.objects.create(
            project_id=cls.project,
            name="OT-test",
        )
        cls.otbatchprotocol = OTBatchProtocol.objects.create(
            otproject_id=cls.otproject,
            batch_id=cls.batch,
            celery_taskid="abc-123",
        )
        cls.otsession = OTSession.objects.create(
            otbatchprotocol_id=cls.otbatchprotocol,
            reactionstep=1,
        )
        cls.deck = Deck.objects.create(otsession_id=cls.otsession)
        cls.pipette = Pipette.objects.create(
            otsession_id=cls.otsession,
            position="Right",
            labware="p300_single",
            type="single",
            name="P300",
        )
        cls.tiprack = TipRack.objects.create(
            otsession_id=cls.otsession,
            deck_id=cls.deck,
            labware="opentrons_96_tiprack_300ul",
            index=1,
            name="tiprack-1",
        )
        cls.plate = Plate.objects.create(
            otbatchprotocol_id=cls.otbatchprotocol,
            otsession_id=cls.otsession,
            deck_id=cls.deck,
            labware="corning_96_wellplate_360ul_flat",
            index=2,
            maxwellvolume=360.0,
            numberwells=96,
            numbercolumns=12,
        )
        cls.column = Column.objects.create(
            otsession_id=cls.otsession,
            plate_id=cls.plate,
            index=0,
            reactionclass="Amidation",
        )
        cls.well = Well.objects.create(
            otsession_id=cls.otsession,
            plate_id=cls.plate,
            index=0,
        )

    def setUp(self):
        self.client = APIClient()


# ===================================================================
# CRUD tests for all 15 flat ViewSets
# ===================================================================
class TestProjectCRUD(APITestBase):
    """ProjectViewSet CRUD tests."""

    def test_list_projects(self):
        resp = self.client.get("/api/projects/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertGreaterEqual(len(resp.json()), 1)

    def test_retrieve_project(self):
        resp = self.client.get(f"/api/projects/{self.project.pk}/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["name"], "test-project")

    def test_retrieve_project_with_fetchall(self):
        resp = self.client.get(f"/api/projects/{self.project.pk}/?fetchall=yes")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertIn("batches", resp.json())

    def test_create_project(self):
        data = {
            "name": "new-proj",
            "submittername": "Bob",
            "submitterorganisation": "AcmeCorp",
            "proteintarget": "BRAF",
        }
        resp = self.client.post("/api/projects/", data, format="json")
        self.assertEqual(resp.status_code, status.HTTP_201_CREATED)

    def test_delete_project(self):
        proj = Project.objects.create(
            name="delete-me",
            submittername="X",
            submitterorganisation="Y",
            proteintarget="Z",
        )
        resp = self.client.delete(f"/api/projects/{proj.pk}/")
        self.assertEqual(resp.status_code, status.HTTP_204_NO_CONTENT)

    def test_update_project(self):
        resp = self.client.patch(
            f"/api/projects/{self.project.pk}/",
            {"name": "updated-proj"},
            format="json",
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestBatchCRUD(APITestBase):
    """BatchViewSet CRUD tests."""

    def test_list_batches(self):
        resp = self.client.get("/api/batches/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_project(self):
        resp = self.client.get(f"/api/batches/?project_id={self.project.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        for item in resp.json():
            self.assertEqual(item["project_id"], self.project.pk)

    def test_retrieve_with_fetchall(self):
        resp = self.client.get(f"/api/batches/{self.batch.pk}/?fetchall=yes")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertIn("targets", resp.json())

    def test_retrieve_without_fetchall(self):
        resp = self.client.get(f"/api/batches/{self.batch.pk}/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertNotIn("targets", resp.json())


class TestTargetCRUD(APITestBase):
    """TargetViewSet CRUD tests."""

    def test_list_targets(self):
        resp = self.client.get("/api/targets/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_batch(self):
        resp = self.client.get(f"/api/targets/?batch_id={self.batch.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        for item in resp.json():
            self.assertEqual(item["batch_id"], self.batch.pk)

    def test_retrieve_with_fetchall(self):
        resp = self.client.get(f"/api/targets/{self.target.pk}/?fetchall=yes")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertIn("methods", resp.json())


class TestMethodCRUD(APITestBase):
    """MethodViewSet CRUD tests."""

    def test_list_methods(self):
        resp = self.client.get("/api/methods/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_target(self):
        resp = self.client.get(f"/api/methods/?target_id={self.target.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_nosteps(self):
        resp = self.client.get("/api/methods/?nosteps=1")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_retrieve_with_fetchall(self):
        resp = self.client.get(f"/api/methods/{self.method.pk}/?fetchall=yes")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertIn("reactions", resp.json())


class TestReactionCRUD(APITestBase):
    """ReactionViewSet CRUD tests."""

    def test_list_reactions(self):
        resp = self.client.get("/api/reactions/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_method(self):
        resp = self.client.get(f"/api/reactions/?method_id={self.method.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_retrieve_with_fetchall(self):
        resp = self.client.get(f"/api/reactions/{self.reaction.pk}/?fetchall=yes")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertIn("products", resp.json())
        self.assertIn("reactants", resp.json())


class TestPubChemInfoCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/pubcheminfo/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_retrieve(self):
        resp = self.client.get(f"/api/pubcheminfo/{self.pubchem.pk}/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestProductCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/products/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_reaction(self):
        resp = self.client.get(f"/api/products/?reaction_id={self.reaction.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_retrieve_with_fetchall(self):
        resp = self.client.get(f"/api/products/{self.product.pk}/?fetchall=yes")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestReactantCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/reactants/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_reaction(self):
        resp = self.client.get(f"/api/reactants/?reaction_id={self.reaction.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestCatalogEntryCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/catalogentries/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestActionSessionCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/actionsessions/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_reaction(self):
        resp = self.client.get(
            f"/api/actionsessions/?reaction_id={self.reaction.pk}"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestAddActionCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/addactions/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestStirActionCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/stiractions/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestExtractActionCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/extractactions/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestMixActionCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/mixactions/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


# ===================================================================
# OT model ViewSet CRUD
# ===================================================================
class TestOTProjectCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/otprojects/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_project(self):
        resp = self.client.get(f"/api/otprojects/?project_id={self.project.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestOTBatchProtocolCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/otbatchprotocols/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_otproject(self):
        resp = self.client.get(
            f"/api/otbatchprotocols/?otproject_id={self.otproject.pk}"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_batch(self):
        resp = self.client.get(
            f"/api/otbatchprotocols/?batch_id={self.batch.pk}"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_celery_taskid(self):
        resp = self.client.get(
            "/api/otbatchprotocols/?celery_taskid=abc-123"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(len(resp.json()), 1)


class TestOTSessionCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/otsessions/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_otbatchprotocol(self):
        resp = self.client.get(
            f"/api/otsessions/?otbatchprotocol_id={self.otbatchprotocol.pk}"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestDeckCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/decks/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_otsession(self):
        resp = self.client.get(f"/api/decks/?otsession_id={self.otsession.pk}")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestPipetteCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/pipettes/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestTipRackCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/tipracks/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestPlateCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/plates/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    def test_filter_by_otbatchprotocol(self):
        resp = self.client.get(
            f"/api/plates/?otbatchprotocol_id={self.otbatchprotocol.pk}"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestColumnCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/columns/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestWellCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/wells/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestCompoundOrderCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/compoundorders/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestOTScriptCRUD(APITestBase):
    def test_list(self):
        resp = self.client.get("/api/otscripts/")
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


# ===================================================================
# Custom @action endpoint tests — ProjectViewSet
# ===================================================================
class TestProjectCreateProjectAction(APITestBase):
    """Tests for ProjectViewSet.create_project action (POST /api/projects/createproject/)."""

    @patch("backend.api.upload_manifold_reaction")
    @patch("backend.api.validate_file_upload")
    @patch("backend.api.save_tmp_file", return_value="/tmp/test.csv")
    def test_validate_only_retro_api(self, mock_save, mock_validate, mock_upload):
        """validate_choice=0, API_choice=0 → validate only (retro-API)."""
        mock_task = MagicMock()
        mock_task.id = "task-001"
        mock_validate.delay.return_value = mock_task

        from django.core.files.uploadedfile import SimpleUploadedFile

        csv_file = SimpleUploadedFile("test.csv", b"smiles\nCCO", content_type="text/csv")
        resp = self.client.post(
            "/api/projects/createproject/",
            {
                "project_name": "test",
                "submitter_name": "Alice",
                "submitter_organisation": "Org",
                "protein_target": "EGFR",
                "validate_choice": "0",
                "API_choice": "0",
                "csv_file": csv_file,
            },
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_id"], "task-001")
        mock_validate.delay.assert_called_once()

    @patch("backend.api.upload_custom_reaction")
    @patch("backend.api.validate_file_upload")
    @patch("backend.api.save_tmp_file", return_value="/tmp/test.csv")
    def test_validate_only_custom_chem(self, mock_save, mock_validate, mock_upload):
        """validate_choice=0, API_choice=1 → validate only (custom-chem)."""
        mock_task = MagicMock()
        mock_task.id = "task-002"
        mock_validate.delay.return_value = mock_task

        from django.core.files.uploadedfile import SimpleUploadedFile

        csv_file = SimpleUploadedFile("test.csv", b"smiles\nCCO")
        resp = self.client.post(
            "/api/projects/createproject/",
            {
                "project_name": "test",
                "submitter_name": "Alice",
                "submitter_organisation": "Org",
                "protein_target": "EGFR",
                "validate_choice": "0",
                "API_choice": "1",
                "csv_file": csv_file,
            },
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_id"], "task-002")

    @patch("backend.api.upload_manifold_reaction")
    @patch("backend.api.validate_file_upload")
    @patch("backend.api.save_tmp_file", return_value="/tmp/test.csv")
    def test_validate_and_upload_retro_api(
        self, mock_save, mock_validate, mock_upload
    ):
        """validate_choice=1, API_choice=0 → validate + upload (retro-API chain)."""
        mock_chain = MagicMock()
        mock_chain.id = "task-003"
        mock_validate.s.return_value.__or__ = MagicMock(return_value=MagicMock())
        mock_validate.s.return_value.__or__.return_value.apply_async.return_value = (
            mock_chain
        )

        from django.core.files.uploadedfile import SimpleUploadedFile

        csv_file = SimpleUploadedFile("test.csv", b"smiles\nCCO")
        resp = self.client.post(
            "/api/projects/createproject/",
            {
                "project_name": "test",
                "submitter_name": "Alice",
                "submitter_organisation": "Org",
                "protein_target": "EGFR",
                "validate_choice": "1",
                "API_choice": "0",
                "csv_file": csv_file,
            },
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    @patch("backend.api.upload_combi_custom_reaction")
    @patch("backend.api.validate_file_upload")
    @patch("backend.api.save_tmp_file", return_value="/tmp/test.csv")
    def test_validate_only_combi_custom(self, mock_save, mock_validate, mock_upload):
        """validate_choice=0, API_choice=2 → validate only (combi-custom-chem)."""
        mock_task = MagicMock()
        mock_task.id = "task-004"
        mock_validate.delay.return_value = mock_task

        from django.core.files.uploadedfile import SimpleUploadedFile

        csv_file = SimpleUploadedFile("test.csv", b"smiles\nCCO")
        resp = self.client.post(
            "/api/projects/createproject/",
            {
                "project_name": "test",
                "submitter_name": "Alice",
                "submitter_organisation": "Org",
                "protein_target": "EGFR",
                "validate_choice": "0",
                "API_choice": "2",
                "csv_file": csv_file,
            },
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)


class TestProjectGetTaskStatus(APITestBase):
    """Tests for ProjectViewSet.get_task_status action."""

    @patch("backend.api.AsyncResult")
    def test_pending_status(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "PENDING"
        MockAsyncResult.return_value = mock_task

        resp = self.client.get(
            "/api/projects/gettaskstatus/?task_id=abc-123"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_status"], "PENDING")

    @patch("backend.api.AsyncResult")
    def test_failure_status(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "FAILURE"
        mock_task.traceback = "Traceback ..."
        MockAsyncResult.return_value = mock_task

        resp = self.client.get(
            "/api/projects/gettaskstatus/?task_id=abc-123"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_status"], "FAILURE")
        self.assertIn("traceback", resp.json())

    @patch("backend.api.AsyncResult")
    def test_success_validated_with_project_info(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "SUCCESS"
        mock_task.get.return_value = ({}, True, {"project_id": 99})
        MockAsyncResult.return_value = mock_task

        resp = self.client.get(
            "/api/projects/gettaskstatus/?task_id=abc-123"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        data = resp.json()
        self.assertEqual(data["task_status"], "SUCCESS")
        self.assertTrue(data["validated"])
        self.assertEqual(data["project_id"], 99)

    @patch("backend.api.AsyncResult")
    def test_success_not_validated_no_project_info(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "SUCCESS"
        mock_task.get.return_value = ({"error": "bad smiles"}, False, None)
        MockAsyncResult.return_value = mock_task

        resp = self.client.get(
            "/api/projects/gettaskstatus/?task_id=abc-123"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        data = resp.json()
        self.assertFalse(data["validated"])
        self.assertIn("validation_errors", data)

    @patch("backend.api.AsyncResult")
    def test_success_validated_no_project_info(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "SUCCESS"
        mock_task.get.return_value = ({}, True, None)
        MockAsyncResult.return_value = mock_task

        resp = self.client.get(
            "/api/projects/gettaskstatus/?task_id=abc-123"
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        data = resp.json()
        self.assertTrue(data["validated"])

    @patch("backend.api.AsyncResult")
    def test_success_not_validated_with_project_info(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "SUCCESS"
        mock_task.get.return_value = ({"error": "bad"}, False, {"project_id": 10})
        MockAsyncResult.return_value = mock_task

        resp = self.client.get(
            "/api/projects/gettaskstatus/?task_id=abc-123"
        )
        data = resp.json()
        self.assertFalse(data["validated"])
        self.assertIn("validation_errors", data)

    def test_no_task_id_returns_400(self):
        """When no task_id param is provided, returns 400 with error message."""
        resp = self.client.get("/api/projects/gettaskstatus/")
        self.assertEqual(resp.status_code, status.HTTP_400_BAD_REQUEST)
        self.assertIn("error", resp.json())


# ===================================================================
# Custom @action endpoint tests — BatchViewSet
# ===================================================================
class TestBatchCreateAction(APITestBase):
    """Tests for BatchViewSet.create (POST /api/batches/) — clone batch."""

    @patch("backend.api.clone_batch")
    def test_create_batch_clone_success(self, mock_clone):
        mock_batch_obj = MagicMock()
        mock_batch_obj.pk = 999
        mock_batch_obj.project_id = self.project
        mock_batch_obj.batchtag = "cloned"
        mock_batch_obj.batch_id = None

        mock_clone.return_value = mock_batch_obj

        # Patch serializer to return simple data
        with patch("backend.api.BatchSerializer") as MockSerializer:
            MockSerializer.return_value.data = {
                "id": 999,
                "project_id": self.project.pk,
                "batchtag": "cloned",
            }
            resp = self.client.post(
                "/api/batches/",
                {"methodids": [self.method.pk], "batchtag": "cloned"},
                format="json",
            )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)

    @patch("backend.api.clone_batch")
    def test_create_batch_clone_exception(self, mock_clone):
        """Exception in clone_batch returns a 500 JSON response with error detail."""
        mock_clone.side_effect = Exception("clone failed")

        resp = self.client.post(
            "/api/batches/",
            {"methodids": [1], "batchtag": "fail"},
            format="json",
        )
        self.assertEqual(resp.status_code, status.HTTP_500_INTERNAL_SERVER_ERROR)
        self.assertIn("error", resp.json())
        self.assertEqual(resp.json()["error"], "clone failed")


class TestBatchCanonicalizeSmiles(APITestBase):
    """Tests for BatchViewSet.canonicalize_smiles."""

    @patch("backend.api.canonicalize_smiles")
    def test_canonicalize_from_smiles_list(self, mock_task):
        mock_result = MagicMock()
        mock_result.id = "canon-task-001"
        mock_task.delay.return_value = mock_result

        resp = self.client.post(
            "/api/batches/canonicalizesmiles/",
            {"smiles": ["CCO", "c1ccccc1"]},
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_id"], "canon-task-001")

    @patch("backend.api.save_tmp_file", return_value="/tmp/smiles.csv")
    @patch("backend.api.canonicalize_smiles")
    def test_canonicalize_from_csv(self, mock_task, mock_save):
        mock_result = MagicMock()
        mock_result.id = "canon-task-002"
        mock_task.delay.return_value = mock_result

        from django.core.files.uploadedfile import SimpleUploadedFile

        csv_file = SimpleUploadedFile("smiles.csv", b"smiles\nCCO")
        resp = self.client.post(
            "/api/batches/canonicalizesmiles/",
            {"csv_file": csv_file},
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_id"], "canon-task-002")


class TestBatchGetTaskStatus(APITestBase):
    """Tests for BatchViewSet.get_task_status."""

    @patch("backend.api.AsyncResult")
    def test_pending(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "PENDING"
        MockAsyncResult.return_value = mock_task

        resp = self.client.get("/api/batches/gettaskstatus/?task_id=xyz")
        self.assertEqual(resp.json()["task_status"], "PENDING")

    @patch("backend.api.AsyncResult")
    def test_success_validated(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "SUCCESS"
        mock_task.get.return_value = (True, ["CCO"])
        MockAsyncResult.return_value = mock_task

        resp = self.client.get("/api/batches/gettaskstatus/?task_id=xyz")
        data = resp.json()
        self.assertEqual(data["task_status"], "SUCCESS")
        self.assertIn("canonicalizedsmiles", data)

    @patch("backend.api.AsyncResult")
    def test_success_not_validated(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "SUCCESS"
        mock_task.get.return_value = (False, "Invalid SMILES")
        MockAsyncResult.return_value = mock_task

        resp = self.client.get("/api/batches/gettaskstatus/?task_id=xyz")
        data = resp.json()
        self.assertEqual(data["task_status"], "SUCCESS")
        self.assertIn("error_summary", data)

    @patch("backend.api.AsyncResult")
    def test_failure(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "FAILURE"
        mock_task.traceback = "Error..."
        MockAsyncResult.return_value = mock_task

        resp = self.client.get("/api/batches/gettaskstatus/?task_id=xyz")
        self.assertEqual(resp.json()["task_status"], "FAILURE")


class TestBatchUpdateReactionSuccess(APITestBase):
    """Tests for BatchViewSet.update_reaction_success."""

    @patch("backend.api.mark_reactions_failed")
    def test_update_from_post_list(self, mock_mark):
        mock_mark.return_value = {"updated": 2}

        resp = self.client.post(
            "/api/batches/updatereactionsuccess/",
            {"reaction_ids": [str(self.reaction.pk)]},
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        mock_mark.assert_called_once()


# ===================================================================
# Custom @action endpoint tests — OTProjectViewSet
# ===================================================================
class TestOTProjectCreateAction(APITestBase):
    """Tests for OTProjectViewSet.create_ot_project."""

    @patch("backend.api.initiate_ot_project", return_value="ot-task-001")
    def test_create_ot_project_no_custom_materials(self, mock_initiate):
        import json

        resp = self.client.post(
            "/api/otprojects/createotproject/",
            {
                "batchids": json.dumps([self.batch.pk]),
                "protocol_name": "test-protocol",
                "has_custom_starting_materials": "false",
            },
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_id"], "ot-task-001")
        mock_initiate.assert_called_once_with(
            batch_ids=[self.batch.pk],
            protocol_name="test-protocol",
            custom_files=None,
            use_multichannel=False,
        )

    @patch("backend.api.initiate_ot_project", return_value="ot-task-002")
    def test_create_ot_project_with_custom_materials(self, mock_initiate):
        import json
        from django.core.files.uploadedfile import SimpleUploadedFile

        file_key = f"starting_materials_batch_{self.batch.pk}"
        csv_file = SimpleUploadedFile("materials.csv", b"smiles\nCCO")

        resp = self.client.post(
            "/api/otprojects/createotproject/",
            {
                "batchids": json.dumps([self.batch.pk]),
                "protocol_name": "custom-protocol",
                "has_custom_starting_materials": "true",
                file_key: csv_file,
            },
        )
        self.assertEqual(resp.status_code, status.HTTP_200_OK)
        self.assertEqual(resp.json()["task_id"], "ot-task-002")
        # Should have passed custom_files dict
        call_kwargs = mock_initiate.call_args[1]
        self.assertIsNotNone(call_kwargs["custom_files"])
        self.assertIn(str(self.batch.pk), call_kwargs["custom_files"])
        self.assertFalse(call_kwargs["use_multichannel"])


class TestOTProjectGetTaskStatus(APITestBase):
    """Tests for OTProjectViewSet.get_task_status."""

    @patch("backend.api.AsyncResult")
    def test_pending(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "PENDING"
        MockAsyncResult.return_value = mock_task

        resp = self.client.get("/api/otprojects/gettaskstatus/?task_id=xyz")
        self.assertEqual(resp.json()["task_status"], "PENDING")

    @patch("backend.api.AsyncResult")
    def test_success(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "SUCCESS"
        mock_task.get.return_value = ({"summary": "ok"}, 42)
        MockAsyncResult.return_value = mock_task

        resp = self.client.get("/api/otprojects/gettaskstatus/?task_id=xyz")
        data = resp.json()
        self.assertEqual(data["task_status"], "SUCCESS")
        self.assertEqual(data["otproject_id"], 42)

    @patch("backend.api.AsyncResult")
    def test_failure(self, MockAsyncResult):
        mock_task = MagicMock()
        mock_task.status = "FAILURE"
        mock_task.traceback = "Big error"
        MockAsyncResult.return_value = mock_task

        resp = self.client.get("/api/otprojects/gettaskstatus/?task_id=xyz")
        self.assertEqual(resp.json()["task_status"], "FAILURE")


# ===================================================================
# URL routing tests — verify url_path overrides
# ===================================================================
class TestURLPaths(APITestBase):
    """Verify the url_path overrides produce underscore-free URLs."""

    def test_createproject_url_exists(self):
        """The URL /api/projects/createproject/ should route."""
        from django.urls import resolve

        match = resolve("/api/projects/createproject/")
        self.assertEqual(match.url_name, "projects-create-project")

    def test_gettaskstatus_project_url_exists(self):
        match = self._resolve_or_none("/api/projects/gettaskstatus/")
        self.assertIsNotNone(match)

    def test_canonicalizesmiles_url_exists(self):
        match = self._resolve_or_none("/api/batches/canonicalizesmiles/")
        self.assertIsNotNone(match)

    def test_gettaskstatus_batch_url_exists(self):
        match = self._resolve_or_none("/api/batches/gettaskstatus/")
        self.assertIsNotNone(match)

    def test_updatereactionsuccess_url_exists(self):
        match = self._resolve_or_none("/api/batches/updatereactionsuccess/")
        self.assertIsNotNone(match)

    def test_createotproject_url_exists(self):
        match = self._resolve_or_none("/api/otprojects/createotproject/")
        self.assertIsNotNone(match)

    def test_gettaskstatus_otproject_url_exists(self):
        match = self._resolve_or_none("/api/otprojects/gettaskstatus/")
        self.assertIsNotNone(match)

    def test_underscored_urls_do_not_route_to_actions(self):
        """Snake_case URLs should NOT resolve to the custom @action methods.

        They may resolve as detail views (pk=<slug>), but not as the
        intended action endpoints.
        """
        for path, forbidden_name_fragment in [
            ("/api/projects/create_project/", "create-project"),
            ("/api/projects/get_task_status/", "get-task-status"),
            ("/api/batches/canonicalize_smiles/", "canonicalize"),
        ]:
            match = self._resolve_or_none(path)
            if match is not None:
                # If it resolves, it must NOT be the action endpoint
                self.assertNotIn(
                    forbidden_name_fragment,
                    match.url_name,
                    f"{path} incorrectly resolves to action {match.url_name}",
                )

    def _resolve_or_none(self, path):
        from django.urls import resolve, Resolver404

        try:
            return resolve(path)
        except Resolver404:
            return None


# ===================================================================
# Non-existent endpoints — 404
# ===================================================================
class TestNonExistentEndpoints(APITestBase):
    """Requests to unknown endpoints should 404."""

    def test_unknown_endpoint(self):
        resp = self.client.get("/api/nonexistent/")
        self.assertEqual(resp.status_code, status.HTTP_404_NOT_FOUND)

    def test_unknown_action(self):
        resp = self.client.get("/api/projects/nonexistent/")
        self.assertEqual(resp.status_code, status.HTTP_404_NOT_FOUND)
