cd CAR
python3 manage.py makemigrations backend
python3 manage.py migrate backend
python3 manage.py runserver
celery -A CAR worker -l info