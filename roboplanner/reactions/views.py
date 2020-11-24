from django.http import HttpResponseRedirect
from django.shortcuts import render
from django.views import View

from .forms import UploadForm


def save_tmp_file(myfile):
""" Save file in temporary location for validation/upload processing
"""

    name = myfile.name
    path = default_storage.save('tmp/' + name, ContentFile(myfile.read()))
    tmp_file = str(os.path.join(settings.MEDIA_ROOT, path))

class UploadProject(View):


    return tmp_file
    def get(self, request):
        form = UploadForm()
        return render(request, 'reactions/upload.html', {'form': form})

    def post(self, request):
        # check_services - from Fragalysis to check Celery stuff. May need it
        form = UploadForm(request.POST, request.FILES)
        context={}
        if form.is_valid():
            # Get info from form submitted
            csvfile = request.FILES['csv_file']
            projectname = request.POST['project_name']
            submittername =  request.POST['submitter_name']
            submitteremail =  request.POST['submitter_email']
            choice = request.POST['submit_choice']
            
            # Save csv to temp storage
            tmp_file = save_tmp_file(csvfile)

            # Settings for if validate option selected
            if str(choice) == '0':
                #### Got up to here - need to look at validate task and implement first
                # Start celery task
                task_validate = validate_csv.delay(tmp_file, target=target, zfile=zfile, update=update_set)

                context = {}
                context['validate_task_id'] = task_validate.id
                context['validate_task_status'] = task_validate.status

                # Update client side with task id and status
                return render(request, 'viewer/upload-cset.html', context)

            # if it's an upload, run the compound set task
            if str(choice) == '1':
                # Start chained celery tasks. NB first function passes tuple
                # to second function - see tasks.py
                task_upload = (
                            validate_compound_set.s(tmp_file, target=target, zfile=zfile, update=update_set) | process_compound_set.s()).apply_async()

                context = {}
                context['upload_task_id'] = task_upload.id
                context['upload_task_status'] = task_upload.status

                # Update client side with task id and status
                return render(request, 'viewer/upload-cset.html', context)
        
        else:
            form = UploadForm()


    context['form'] = form
    return render(request, 'reactions/upload.html', context)

# Create your views here.
