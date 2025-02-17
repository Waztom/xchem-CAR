#!/bin/bash
echo "Running migrations..."
cd /code/CAR/
# initial migrations for existing stack
python3 manage.py makemigrations
python3 manage.py migrate

# echo "Running collectstatic..."
# python3 manage.py collectstatic --noinput -v 0 # collect static files
