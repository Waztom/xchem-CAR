from django import forms

validate_CHOICES = [
    (0, "validate"),
    (1, "upload"),
]

API_CHOICES = [(0, "Postera"), (1, "IBM"), (2, "Custom chemistry"), (3, "Combi custom chemistry")]


class UploadForm(forms.Form):
    submitter_name = forms.CharField(label="Your name", max_length=100)
    submitter_organisation = forms.CharField(label="Your organisation", max_length=100)
    submitter_email = forms.EmailField(label="Your email")
    csv_file = forms.FileField(label="CSV file")
    validate_choice = forms.CharField(widget=forms.RadioSelect(choices=validate_CHOICES))
    API_choice = forms.CharField(widget=forms.RadioSelect(choices=API_CHOICES))