from django import forms

CHOICES = [
    (0, 'validate'),
    (1, 'upload'),
]

class UploadForm(forms.Form):
    submitter_name = forms.CharField(label='Your name', max_length=100)
    submitter_email = forms.EmailField(label='Your email')
    project_name = forms.CharField(label='Project name', max_length=100)
    csv_file = forms.FileField(label='CSV file')
    submit_choice = forms.CharField(widget=forms.RadioSelect(choices=CHOICES))

