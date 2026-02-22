"""Management command to delete orphaned OTProject records.

An OTProject becomes orphaned when all of its OTBatchProtocols have been
cascade-deleted (because their parent Batches were deleted).  This command
finds and removes those empty OTProjects so they no longer clutter the
OT Protocol History list.

Usage:
    python manage.py cleanup_orphan_otprojects
    python manage.py cleanup_orphan_otprojects --dry-run
"""

from django.core.management.base import BaseCommand

from backend.models import OTProject, OTBatchProtocol


class Command(BaseCommand):
    help = "Delete OTProject records that have no remaining OTBatchProtocols."

    def add_arguments(self, parser):
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="List orphans without deleting them.",
        )

    def handle(self, *args, **options):
        dry_run = options["dry_run"]

        orphans = OTProject.objects.exclude(
            pk__in=OTBatchProtocol.objects.values_list("otproject_id", flat=True)
        )
        count = orphans.count()

        if count == 0:
            self.stdout.write(self.style.SUCCESS("No orphaned OTProjects found."))
            return

        if dry_run:
            self.stdout.write(f"Found {count} orphaned OTProject(s) (dry run):")
            for otp in orphans:
                self.stdout.write(f"  id={otp.pk}  name={otp.name}  date={otp.init_date}")
        else:
            orphans.delete()
            self.stdout.write(
                self.style.SUCCESS(f"Deleted {count} orphaned OTProject(s).")
            )
