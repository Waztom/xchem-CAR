# from django.test import TestCase
# from pytest import raises

# from celery.exceptions import Retry

# # for python 2: use mock.patch from `pip install mock`.
# from unittest.mock import patch

# from .models import *
# from .createmodels import *
# from .validate import ValidateFile

# class validateFileUploadTestCase(TestCase):
#     def setUp(self):
#         self.project_info = {
#             "projectname": "Test",
#             "submittername": "Tester",
#             "submitterorganisation": "Rest Framework",
#             "proteintarget": "MID2",
#             "API_choice": 0,
#         }
#         self.csvfile = open(
#             "/container/car/tests/testdata/indata/postera-ver1.5.csv", "rb"
#         )

#     @override_settings(CELERY_EAGER_PROPAGATES_EXCEPTIONS=True,
#                        CELERY_ALWAYS_EAGER=True,
#                        BROKER_BACKEND='memory')
#     def test_validation_celery_task(self):
#         ValidateFile(self.project_info, self.csvfile)
#         self.assertEqual(Project.objects.count(), 1)
#         self.assertEqual(Reaction.objects.count(), 1)
#         self.assertEqual(Reactant.objects.count(), 2)
#         self.assertEqual(Product.objects.count(), 1)
#         self.assertEqual(CatalogEntry.objects.count(), 1)
#         self.assertEqual(Method.objects.count(), 1)
#         self.assertEqual(Target.objects.count(), 1)
#         self.assertEqual(Batch.objects.count(), 1)

#     def test_create_project_model(self):
#         project = createProjectModel(self.project_info)
#         self.assertEqual(project.projectname, "Test")
#         self.assertEqual(project.submittername, "Tester")
#         self.assertEqual(project.submitterorganisation, "Rest Framework")
#         self.assertEqual(project.proteintarget, "MID2")
#         self.assertEqual(project.API_choice, 0)

#     def test_create_reaction_model(self):
#         project = createProjectModel(self.project_info)
#         reaction = createReactionModel(project)
#         self.assertEqual(reaction.project, project)

#     def test_create_reactant_model(self):
#         project = createProjectModel(self.project_info)
#         reaction = createReactionModel(project)
#         reactant = createReactantModel(reaction)
#         self.assertEqual(reactant.reaction, reaction)

#     def test_create_product_model(self):
#         project = createProjectModel(self.project_info)
#         reaction = createReactionModel(project)
#         product = createProductModel(reaction)
#         self.assertEqual(product.reaction, reaction)

#     def test_create_catalog_entry_model(self):
#         project = createProjectModel(self.project_info)
#         reaction = createReactionModel(project)
#         product = createProductModel(reaction)
#         catalog_entry = createCatalogEntryModel(product)
#         self.assertEqual(catalog_entry.product, product)

#     def test_create_method_model(self):
#         project = createProjectModel(self.project_info)
#         reaction = createReactionModel(project)
