# Recipe JSON Fixtures

This directory contains canonical recipe definitions in JSON format for use with the `load_recipe_json` management command and recipe model tests.

## Generating Recipes

To create a new recipe fixture:

1. **Start from a template:**
   - Use `amidation_standard.json` for a simple reaction with standard workup.
   - Use `mitsunobu_with_extraction.json` for a multi-stage workup example (liquid-liquid extraction, explicit plate routing).

2. **Edit the JSON:**
   - Only include fields that differ from the documented defaults. Omit `plate_index`, `from_plate_index`, `to_plate_index`, `quantity_unit`, `continuation`, etc. unless needed for multi-stage routing or non-default values.
   - Use `from_plate_role` and `to_plate_role` to specify plate types (e.g., `reaction`, `workup`, `spefilter`, `solvent`).
   - For multi-stage workup, set `to_plate_index` or `from_plate_index` explicitly (e.g., `workup1`, `workup2`).
   - Add comments (with `_comment` key) for clarity, but these are ignored by the loader.

3. **Validate the JSON:**
   - Run:
     ```bash
     python3 manage.py load_recipe_json --dry-run backend/tests/fixtures/your_recipe.json
     ```
   - Fix any errors reported.

4. **Load into the DB:**
   - Run:
     ```bash
     python3 manage.py load_recipe_json backend/tests/fixtures/your_recipe.json
     ```
   - Use `--update` to create a new version if the recipe already exists.

## Example Structure

```json
{
  "reaction_class": "Amidation",
  "name": "standard",
  "action_sessions": [
    {
      "session_type": "reaction",
      "driver": "robot",
      "actions": [
        {
          "type": "add",
          "material_smarts": "[#6:1](=[#8:2])-[#8;H1]",
          "equivalents": 1.0,
          "solvent": "DMA",
          "concentration": 0.5,
          "from_plate_role": "startingmaterial",
          "to_plate_role": "reaction"
        }
        // ... more actions ...
      ]
    }
    // ... more sessions ...
  ]
}
```

## Plate Routing

- For simple recipes, omit all `*_plate_index` fields (default is 1).
- For multi-stage workup, specify `to_plate_index` or `from_plate_index` as needed:
  - `extract: reaction → workup1`
  - `extract: reaction → workup2`
  - `add: workup2 → spefilter`

## Management Command Reference

- `load_recipe_json` supports:
  - `--dry-run`: Validate without writing
  - `--update`: Create new version if recipe exists
  - Multiple files at once

## Canonical Examples

- `amidation_standard.json`: Standard amidation, no explicit plate_index
- `mitsunobu_with_extraction.json`: Mitsunobu with multi-stage extraction, explicit plate_index routing

## See Also
- [../README.md](../README.md) for project overview
- [../management/commands/load_recipe_json.py](../management/commands/load_recipe_json.py) for loader implementation
