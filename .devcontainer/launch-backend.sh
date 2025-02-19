#!/bin/bash

# Create docs directory
mkdir -p /code/CAR/docs

# Create log folder and file
mkdir -p /code/CAR/logs
touch /code/CAR/logs/logfile.log

# Run migrations
echo "Running migrations..."
cd /code/CAR/
python3 manage.py makemigrations
python3 manage.py migrate
