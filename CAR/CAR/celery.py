from __future__ import absolute_import, unicode_literals
import os
from django.conf import settings  
from celery import Celery

# set default Django settings module for celery
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'CAR.settings')
 
app = Celery('CAR')
app.config_from_object('django.conf:settings', namespace='CELERY')  
app.autodiscover_tasks()

# Commented out below - worked with rabbit

# # set the default Django settings module for the 'celery' program.
# os.environ.setdefault("DJANGO_SETTINGS_MODULE", "CAR.settings")

# # Can use Redis - just need to install it
# # app = Celery('CAR', backend='redis://localhost:6379/0', broker='pyamqp://')
# # app = Celery("CAR", backend="db+sqlite:///results.db", broker="pyamqp://")

# # Using a string here means the worker doesn't have to serialize
# # the configuration object to child processes.
# # - namespace='CELERY' means all celery-related configuration keys
# #   should have a `CELERY_` prefix.
# app.config_from_object("django.conf:settings", namespace="CELERY")

# # Load task modules from all registered Django app configs.
# app.autodiscover_tasks()


@app.task(bind=True)
def debug_task(self):
    print(("Request: {0!r}".format(self.request)))
