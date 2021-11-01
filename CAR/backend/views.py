import os
from django.shortcuts import render
from django.views import View
from django.core.files.storage import default_storage
from django.core.files.base import ContentFile
from django.core.mail import send_mail
from django.http import JsonResponse
from django.conf import settings
from celery.result import AsyncResult
from .forms import API_CHOICES, UploadForm
from .tasks import (
    validateFileUpload,
    uploadIBMReaction,
    uploadManifoldReaction,
    uploadCustomReaction,
)
import pandas as pd


def save_tmp_file(myfile):
    name = myfile.name
    path = default_storage.save("tmp/" + name, ContentFile(myfile.read()))
    tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))
    return tmp_file


class UploadProject(View):
    def get(self, request):
        form = UploadForm()
        return render(request, "backend/upload.html", {"form": form})

    def post(self, request):
        form = UploadForm(request.POST, request.FILES)
        if form.is_valid():
            project_info = {}
            csvfile = request.FILES["csv_file"]
            project_info["submittername"] = request.POST["submitter_name"]
            project_info["submitterorganisation"] = request.POST["submitter_organisation"]
            project_info["submitteremail"] = request.POST["submitter_email"]
            validate_choice = request.POST["validate_choice"]
            API_choice = request.POST["API_choice"]
            tmp_file = save_tmp_file(csvfile)

            if str(validate_choice) == "0":

                if str(API_choice) == "2":
                    task_validate = validateFileUpload.delay(
                        csv_fp=tmp_file, validate_type="custom-chem"
                    )
                    context = {}
                    context["validate_task_id"] = task_validate.id
                    context["validate_task_status"] = task_validate.status

                    return render(request, "backend/upload.html", context)
                if str(API_choice) == "3":
                    task_validate = validateFileUpload.delay(
                        csv_fp=tmp_file, validate_type="combi-custom-chem"
                    )
                    context = {}
                    context["validate_task_id"] = task_validate.id
                    context["validate_task_status"] = task_validate.status

                    return render(request, "backend/upload.html", context)

                else:
                    task_validate = validateFileUpload.delay(
                        csv_fp=tmp_file, validate_type="retro-API"
                    )
                    context = {}
                    context["validate_task_id"] = task_validate.id
                    context["validate_task_status"] = task_validate.status

                    return render(request, "backend/upload.html", context)

            if str(validate_choice) == "1":
                if str(API_choice) == "0":
                    task_upload = (
                        validateFileUpload.s(
                            csv_fp=tmp_file,
                            validate_type="retro-API",
                            project_info=project_info,
                            validate_only=False,
                        )
                        | uploadManifoldReaction.s()
                    ).apply_async()

                    context = {}
                    context["upload_task_id"] = task_upload.id
                    context["upload_task_status"] = task_upload.status

                    return render(request, "backend/upload.html", context)
                    pass

                if str(API_choice) == "1":
                    task_upload = (
                        validateFileUpload.s(
                            csv_fp=tmp_file,
                            validate_type="retro-API",
                            project_info=project_info,
                            validate_only=False,
                        )
                        | uploadIBMReaction.s()
                    ).apply_async()

                    context = {}
                    context["upload_task_id"] = task_upload.id
                    context["upload_task_status"] = task_upload.status

                    return render(request, "backend/upload.html", context)

                if str(API_choice) == "2":
                    task_upload = (
                        validateFileUpload.s(
                            csv_fp=tmp_file,
                            validate_type="custom-chem",
                            project_info=project_info,
                            validate_only=False,
                        )
                        | uploadCustomReaction.s()
                    ).apply_async()

                    context = {}
                    context["upload_task_id"] = task_upload.id
                    context["upload_task_status"] = task_upload.status

                    return render(request, "backend/upload.html", context)

                if str(API_choice) == "3":
                    task_upload = (
                        validateFileUpload.s(
                            csv_fp=tmp_file,
                            validate_type="combi-custom-chem",
                            project_info=project_info,
                            validate_only=False,
                        )
                        | uploadCustomReaction.s()
                    ).apply_async()

                    context = {}
                    context["upload_task_id"] = task_upload.id
                    context["upload_task_status"] = task_upload.status

                    return render(request, "backend/upload.html", context)

        else:
            form = UploadForm()


class ValidateTaskView(View):
    """View to handle dynamic loading of validation results from `backend.tasks.validateFileUpload` - the validation of files
    uploaded to backend/upload
    Methods
    -------
    allowed requests:
        - GET: takes a task id, checks it's status and returns the status, and result if the task is complete
    url:
        validate_task/<validate_task_id>
    template:
        viewer/upload-cset.html or viewer/upload-tset.html
    """

    def get(self, request, validate_task_id):
        """Get method for `ValidateTaskView`. Takes a validate task id, checks it's status and returns the status,
        and result if the task is complete
        Parameters
        ----------
        request: request
            Context sent by `UploadCSet` or `UploadTset`
        validate_task_id: str
            task id provided by `UploadCSet` or `UploadTset`
        Returns
        -------
        response_data: JSON
            response data (dict) in JSON format:
                - if status = 'RUNNING':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                - if status = 'FAILURE':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                    - validate_traceback (str): task.traceback
                - if status = 'SUCCESS':
                    - validate_task_status (str): task.status
                    - validate_task_id (str): task.id
                    - html (str): html of task outcome - success message or html table of errors & fail message
        """

        task = AsyncResult(validate_task_id)
        response_data = {"validate_task_status": task.status, "validate_task_id": task.id}

        if task.status == "FAILURE":
            result = task.traceback
            response_data["validate_traceback"] = str(result)

            return JsonResponse(response_data)

        if task.status == "SUCCESS":
            results = task.get()
            validate_dict = results[0]
            validated = results[1]
            if validated:
                response_data[
                    "html"
                ] = "Your data was validated. \n It can now be uploaded using the upload option."
                response_data["validated"] = "Validated"

                return JsonResponse(response_data)

            if not validated:
                pd.set_option("display.max_colwidth", None)

                table = pd.DataFrame.from_dict(validate_dict)
                html_table = table.to_html()
                html_table += (
                    """<p> Your data was <b>not</b> validated. The table above shows errors</p>"""
                )

                response_data["html"] = html_table
                response_data["validated"] = "Not validated"

                return JsonResponse(response_data)

        return JsonResponse(response_data)


class UploadTaskView(View):
    """View to handle dynamic loading of upload results from `backend.tasks.UploadIBMReaction` - the upload of files
    for a computed set by a user at viewer/upload_cset or a target set by a user at viewer/upload_tset
    Methods
    -------
    allowed requests:
        - GET: takes a task id, checks it's status and returns the status, and result if the task is complete
    url:
        upload_task/<uploads_task_id>
    template:
        viewer/upload-cset.html or viewer/upload-tset.html
    """

    def get(self, request, upload_task_id):
        """Get method for `UploadTaskView`. Takes an upload task id, checks it's status and returns the status,
        and result if the task is complete
        Parameters
        ----------
        request: request
            Context sent by `UploadCSet` or `UploadTSet`
        upload_task_id: str
            task id provided by `UploadCSet` or `UploadTSet`
        Returns
        -------
        response_data: JSON
            response data (dict) in JSON format:
                - if status = 'RUNNING':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                - if status = 'FAILURE':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                    - upload_traceback (str): task.traceback
                - if status = 'SUCCESS':
                    - upload_task_status (str): task.status
                    - upload_task_id (str): task.id
                    - if results are a list (data was processed - validated or uploaded):
                        if this was a validation process
                        - validated (str): 'Not validated'
                        - html (str): html table of validation errors
                        if results are a validation/upload process:
                        - validated (str): 'Validated'
                        - results (dict): results
                        For compound sets ('cset')
                        - results['cset_download_url'] (str): download url for computed set sdf file
                        - results['pset_download_url'] (str): download url for computed set pdb files (zip)
                        For target sets ('tset')
                        - results['tset_download_url'] (str): download url for processed zip file
                    - if results are not string or list:
                        - processed (str): 'None'
                        - html (str): message to tell the user their data was not processed
        """
        task = AsyncResult(upload_task_id)
        response_data = {"upload_task_status": task.status, "upload_task_id": task.id}

        if task.status == "FAILURE":
            result = task.traceback
            response_data["upload_traceback"] = str(result)

            return JsonResponse(response_data)

        if task.status == "SUCCESS":

            results = task.get()
            validate_dict = results[0]
            validated = results[1]
            project_info = results[2]

            submitter_name = project_info["submittername"]
            submitter_email = project_info["submitteremail"]
            project_name = project_info["project_name"]

            if validated:
                send_mail(
                    "Project - {} - data uploaded to CAR".format(project_name),
                    "Hi {} - Good news, your data has been successfully processed and is available at https://car.xchem.diamond.ac.uk/".format(
                        submitter_name
                    ),
                    "waztom@gmail.com",
                    [submitter_email],
                    fail_silently=False,
                )

                response_data["validated"] = "Validated"
                return JsonResponse(response_data)

            if not validated:
                pd.set_option("display.max_colwidth", -1)

                table = pd.DataFrame.from_dict(validate_dict)
                html_table = table.to_html()
                html_table += (
                    """<p> Your data was <b>not</b> validated. The table above shows errors</p>"""
                )

                response_data["validated"] = "Not validated"
                response_data["html"] = html_table

                return JsonResponse(response_data)

            else:
                send_mail(
                    "Project - {} - data upload failed to CAR".format(project_name),
                    "Hi {} - Hmmmm, not so good news, your data failed to be processed and uploaded to CAR. Please try again or contact Warren by responding to this email".format(
                        submitter_name
                    ),
                    "waztom@gmail.com",
                    [submitter_email],
                    fail_silently=False,
                )
                html_table = """<p> Your data was <b>not</b> processed.</p>"""
                response_data["processed"] = "None"
                response_data["html"] = html_table
                return JsonResponse(response_data)

        return JsonResponse(response_data)
