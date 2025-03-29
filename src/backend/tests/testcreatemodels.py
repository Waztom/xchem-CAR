from django.test import TestCase
from .createmodels import (
    createProjectModel,
    createBatchModel,
    createReactionModel,
    createReactantModel,
    createProductModel,
    createCatalogEntryModel,
    createMethodModel,
    createTargetModel,
)

from .models import (
    Project,
    Batch,
    Reaction,
    Reactant,
    Product,
    CatalogEntry,
    Method,
    Target,
)


class CreateModelsTestCase(TestCase):
    def setUp(self) -> None:
        self.project_info = {
            "projectname": "Test",
            "submittername": "Tester",
            "submitterorganisation": "Rest Framework",
            "proteintarget": "MID2",
            "API_choice": 0,
        }

    def test_create_project_model(self):
        project_id = createProjectModel(self.project_info)
        project_obj = Project.objects.get(id=project_id)
        self.assertEqual(project_obj.name, "Test")
        self.assertEqual(project_obj.submittername, "Tester")
        self.assertEqual(project_obj.submitterorganisation, "Rest Framework")
        self.assertEqual(project_obj.proteintarget, "MID2")

    def test_create_batch_model(self):
        print(self.project_id)
        self.batch_id = createBatchModel(
            project_id=self.project_id, batch_tag="TestBatch"
        )
        batch_obj = Batch.objects.get(id=self.batch_id)
        self.assertEqual(batch_obj.project_id, self.project_id)
        self.assertEqual(batch_obj.batchtag, "TestBatch")
