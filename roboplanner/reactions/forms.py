from django import forms

CHOICES = [
    (0, 'validate'),
    (1, 'upload'),
]

class UploadForm(forms.Form):
    submitter_name = forms.CharField(label='your_name', max_length=100)
    email = forms.EmailInput(label='your_email', max_length=100)
    project_name = forms.CharField(label='project_name', max_length=100)
    csv_file = forms.FileField(label='CSV file')
    submit_choice = forms.CharField(widget=forms.RadioSelect(choices=CHOICES))

