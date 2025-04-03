# from rest_framework.test import APIClient, APITestCase
# from django.core.files.uploadedfile import SimpleUploadedFile
# from celery.contrib.testing.worker import start_worker
# from fragalysis.celery import app

# from django.core.files.base import ContentFile
# from django.core.files.storage import default_storage
# from django.conf import settings
# import os

# from rdkit.Chem import AllChem

# from .createmodels import (
#     CreateEncodedActionModels,
#     createBatchModel,
#     createCatalogEntryModel,
#     createMethodModel,
#     createProductModel,
#     createProjectModel,
#     createReactantModel,
#     createReactionModel,
#     createTargetModel,
# )
# from .tasks import validateFileUpload
# from .utils import getAddtionOrder

# # from .testdata.indata.testcreateproject import grouped_targets, routes


# def save_tmp_file(myfile):
#     path = default_storage.save(
#         "/container/car/tests/testdata/indata/tmpcsvfile.csv", ContentFile(myfile.read())
#     )
#     return path


# class APICreateProjectTestCase(APITestCase):
#     @classmethod
#     def setUpClass(cls):
#         super().setUpClass()
#         cls.celery_worker = start_worker(app, perform_ping_check=False)
#         cls.celery_worker.__enter__()

#     @classmethod
#     def tearDownClass(cls):
#         super().tearDownClass()
#         cls.celery_worker.__exit__(None, None, None)

#     def setUp(self) -> None:
#         """This is very close to the
#         uploadManifoldReaction Celery task. Ideal to
#         actually use that function for testing. Will need some
#         work to adapt the function.
#         """
#         super().setUp()
#         self.client = APIClient()
#         self.project_info = {
#             "projectname": "Test",
#             "submittername": "Tester",
#             "submitterorganisation": "Rest Framework",
#             "proteintarget": "MID2",
#             "API_choice": 0,
#         }
#         self.csvfile = open(
#             "/container/car/tests/testdata/indata/CAR-test-project-upload-1.2.csv", "rb"
#         )

#     def test_validation_celery_task(self):
#         self.project_info["validate_choice"] = 0
#         tmp_file = save_tmp_file(self.csvfile)
#         test_task = validateFileUpload.delay(csv_fp=tmp_file, validate_type="retor-API")
#         print("XXXXXXXXXXXXXXXXXXX")
#         print(test_task)
#         print("XXXXXXXXXXXXXXXXXXX")
#         self.assertEqual(test_task.status, "SUCCESS")
#         self.assertEqual(test_task.validated, True)


# #     def test_create_project(self):

# #         self.client.post()


# #         for batchtag, group in grouped_targets:
# #             batch_id = createBatchModel(project_id=self.project_id, batchtag=batchtag)
# #             for smiles, mass, route in zip(group["targets"], group["amount-required-mg"], routes):
# #                 target_id = createTargetModel(
# #                     batch_id=batch_id,
# #                     smiles=smiles,
# #                     mass=mass,
# #                 )

# #                 no_steps = len(route["reactions"])
# #                 reactions = route["reactions"]
# #                 method_id = createMethodModel(
# #                             target_id=target_id, nosteps=no_steps, otchem=True
# #                         )

# #                 for reaction in reversed(reactions):
# #                     reaction_name = reaction["name"]
# #                     recipes = encoded_recipes[reaction_name]["recipes"]
# #                     recipe_rxn_smarts = encoded_recipes[reaction_name][
# #                         "reactionSMARTS"
# #                     ]
# #                     reactant_smiles = reaction["reactantSmiles"]
# #                     product_smiles = reaction["productSmiles"]

# #                     if len(reactant_smiles) == 1:
# #                         actions = recipes["Intramolecular"]["actions"]
# #                         stir_action = [
# #                             action
# #                             for action in actions
# #                             if action["type"] == "stir"
# #                         ][0]
# #                         reaction_temperature = stir_action["content"][
# #                             "temperature"
# #                         ]["value"]
# #                         reactant_smiles_ordered = reactant_smiles
# #                     else:
# #                         actions = recipes["Standard"]["actions"]
# #                         stir_action = [
# #                             action
# #                             for action in actions
# #                             if action["type"] == "stir"
# #                         ][0]
# #                         reaction_temperature = stir_action["content"][
# #                             "temperature"
# #                         ]["value"]
# #                         reactant_smiles_ordered = getAddtionOrder(
# #                             product_smi=product_smiles,
# #                             reactant_SMILES=reactant_smiles,
# #                             reaction_SMARTS=recipe_rxn_smarts,
# #                         )
# #                         if not reactant_smiles_ordered:
# #                             continue

# #                     reaction_smarts = AllChem.ReactionFromSmarts(
# #                         "{}>>{}".format(
# #                             ".".join(reactant_smiles_ordered),
# #                             product_smiles,
# #                         ),
# #                         useSmiles=True,
# #                     )

# #                     reaction_id = createReactionModel(
# #                         method_id=method_id,
# #                         reaction_class=reaction_name,
# #                         reaction_temperature=reaction_temperature,
# #                         reaction_smarts=reaction_smarts,
# #                     )

# #                     createProductModel(
# #                         reaction_id=reaction_id,
# #                         product_smiles=product_smiles,
# #                     )

# #                     for reactant_smi in reactant_smiles_ordered:
# #                         reactant_id = createReactantModel(
# #                             reaction_id=reaction_id,
# #                             reactant_smiles=reactant_smi,
# #                         )

# #                         catalog_entries = [
# #                             molecule["catalogEntries"]
# #                             for molecule in route["molecules"]
# #                             if molecule["smiles"] == reactant_smi
# #                         ][0]
# #                         for catalog_entry in catalog_entries:
# #                             createCatalogEntryModel(
# #                                 catalog_entry=catalog_entry,
# #                                 reactant_id=reactant_id,
# #                             )

# #                     CreateEncodedActionModels(
# #                         actions=actions,
# #                         target_id=target_id,
# #                         reaction_id=reaction_id,
# #                         reactant_pair_smiles=reactant_smiles_ordered,
# #                         reaction_name=reaction_name,
# #                     )

# # # class APICreateOTProjectTestCase(APITestCase):
# # #     def setUp(self) -> None:
# # #         self.client = APIClient()
