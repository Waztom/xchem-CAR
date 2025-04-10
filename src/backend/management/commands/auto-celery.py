import shlex
import subprocess
import os
import signal
import time
from django.core.management.base import BaseCommand
from django.utils import autoreload


class Command(BaseCommand):
    help = "Starts Celery worker with auto-reload on code changes"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.celery_worker_process = None

    def handle(self, *args, **options):
        self.stdout.write("Starting celery worker with autoreload...")

        # For Django>=2.2
        autoreload.run_with_reloader(self._restart_celery)

    def _restart_celery(self):
        # Kill any existing worker processes
        self._kill_worker_processes()

        # Start a new worker process
        self._start_new_worker()

    def _kill_worker_processes(self):
        """Kill any running Celery workers more thoroughly"""
        self.stdout.write("Shutting down Celery workers...")

        # Kill the direct child process if we have it
        if self.celery_worker_process:
            try:
                os.killpg(os.getpgid(self.celery_worker_process.pid), signal.SIGTERM)
                # Give workers time to shut down gracefully
                time.sleep(1)
            except (ProcessLookupError, AttributeError):
                pass

        # Multiple approaches to ensure all workers get killed
        subprocess.call("pkill -9 -f 'celery worker'", shell=True)
        subprocess.call("pkill -9 -f 'celery -A CAR'", shell=True)

        # Wait for processes to terminate fully
        time.sleep(1)

    def _start_new_worker(self):
        """Start Celery worker as a detachable subprocess"""
        self.stdout.write("Starting new Celery workers...")

        # Use Popen instead of call to avoid blocking
        cmd = "celery -A CAR worker -l info"

        # Don't redirect stdout/stderr to see output in console
        self.celery_worker_process = subprocess.Popen(
            shlex.split(cmd),
            preexec_fn=os.setsid,  # Create a new process group
        )

        # Allow the worker time to initialize before returning
        time.sleep(1)
        self.stdout.write("Celery workers are now running...")
