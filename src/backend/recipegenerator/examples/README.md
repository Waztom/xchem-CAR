# Recipe Generator Examples

This directory contains example templates and design matrices demonstrating the recipe generator workflow.

## Files

### `amidation_template.json`
A recipe template for screening amidation reaction conditions. Defines:
- Variable parameters (solvent, concentration, equivalents, temperature, duration)
- Action mappings for the standard Amidation recipe
- Named ordering presets for different addition sequences

### `amidation_design_matrix.csv`
A sample design matrix with 12 experiments exploring:
- Solvent screen (DMF, DMA, NMP, DMSO)
- Temperature screen (25°C, 40°C, 60°C)
- Coupling agent equivalents (1.2, 1.5, 2.0)
- Addition order variants
- Concentration effects

## Usage

```python
from backend.recipegenerator import RecipeTemplate, RecipeGenerator
from backend.recipegenerator.parsers import TemplateParser

# Load template
parser = TemplateParser()
template = parser.parse_json("amidation_template.json")

# Create generator
generator = RecipeGenerator(template)

# Generate recipes from design matrix
recipes = generator.from_csv(
    "amidation_design_matrix.csv",
    recipe_name_column="recipe_name"
)

# Each recipe is ready to use with CAR
for recipe_data in recipes:
    print(f"Generated: {recipe_data['name']}")
    # recipe_data['recipe'] contains the full encoded recipe dict
```

## Creating Your Own Templates

### Quick Start: Auto-generate from existing recipe

```python
from backend.recipegenerator.parsers import create_template_from_recipe

# Auto-generate a template with all possible variables mapped
template = create_template_from_recipe(
    reaction_class="Amidation",
    recipe_name="standard",
    template_name="my-amidation-screen"
)

# Save for customization
import json
with open("my_template.json", "w") as f:
    json.dump(template.to_dict(), f, indent=2)
```

### Manual Template Creation

1. Identify the base recipe (reaction class + recipe name)
2. Map each action you want to modify with a semantic ID
3. Define variables pointing to the parameters you want to change
4. Optionally define ordering presets

See `amidation_template.json` for a complete example.

## Design Matrix Format

The design matrix is a CSV or Excel file where:
- Each row is one experiment
- Columns match variable names from the template
- Special columns:
  - `experiment_id`: Optional identifier
  - `recipe_name`: Optional name for generated recipe
  - `order_preset`: Name of ordering preset to use
  - `{session_id}_action_order`: Comma-separated action IDs for custom ordering
  - `session_order`: Comma-separated session IDs for session reordering
  - `notes`: Ignored, for documentation

Columns not matching template variables are ignored (with a warning).
