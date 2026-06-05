from django.core.management.base import BaseCommand

from backend.recipegenerator.generator import RecipeGenerator
from backend.recipegenerator.parsers import TemplateParser

"""

Usage::

    python3 manage.py generate_recipes <recipe_template.json> <design_matrix.csv>

"""

class Command(BaseCommand):
    help = "Generate recipe JSON files from a template and design matrix"

    def add_arguments(self, parser):
        parser.add_argument("template", help="Path to template JSON")
        parser.add_argument("design_matrix", help="Path to design matrix CSV")
        parser.add_argument(
            "-o",
            "--output",
            default="recipe_data",
            help="Output folder for generated JSON files",
        )
        parser.add_argument(
            "--recipe-name-column",
            default="recipe_name",
            help="Column to use for recipe names",
        )

    def handle(self, *args, **options):
        template_path = options["template"]
        design_matrix_path = options["design_matrix"]
        output_dir = options["output"]
        recipe_name_column = options["recipe_name_column"]

        template_parser = TemplateParser()
        template = template_parser.parse_json(template_path)

        generator = RecipeGenerator(template)
        recipes = generator.from_csv(
            design_matrix_path,
            recipe_name_column=recipe_name_column,
            output_dir=output_dir,
        )

        self.stdout.write(self.style.SUCCESS(f"Generated {len(recipes)} recipes"))