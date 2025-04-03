#!/bin/bash

# Create docs directory
mkdir -p /container/src/docs

# Create log folder and file
mkdir -p /container/src/logs
touch /container/src/logs/logfile.log

# Run migrations
echo "Running migrations..."
cd /container/src/
python3 manage.py makemigrations
python3 manage.py migrate
