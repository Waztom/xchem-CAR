# Recipe JSON Fixtures

This directory contains canonical recipe definitions in JSON format for use with the `load_recipe_json` management command and recipe model tests.

## Generating Recipes

To create a new recipe fixture:

1. **Start from a template:**
   - Use `amidation_standard.json` for a simple reaction with standard workup.
   - Use `mitsunobu_with_extraction.json` for a multi-stage workup example (liquid-liquid extraction, explicit plate routing).

2. **Edit the JSON:**
   - Only include fields that differ from the documented defaults. Omit `plate_role_index`, `from_plate_role_index`, `to_plate_role_index`, `quantity_unit`, `continuation`, etc. unless needed for multi-stage routing or non-default values.
   - Use `from_plate_role` and `to_plate_role` to specify plate types (e.g., `reaction`, `workup`, `spefilter`, `solvent`).
   - For multi-stage workup, set `to_plate_role_index` or `from_plate_role_index` explicitly (e.g., index 1 for first workup plate, index 2 for second).
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
          "amount": 1.0,
          "solvent": "DMA",
          "concentration": 0.5,
          "from_plate_role": "startingmaterial",
          "to_plate_role": "reaction"
        }
      ]
    },
    {
      "session_type": "stir",
      "driver": "human",
      "actions": [
        {
          "type": "stir",
          "temperature": 25,
          "duration": 12.0,
          "plate_role": "reaction"
        }
      ]
    }
  ]
}
```

## Units Reference

### Quantity Units (`quantity_unit`)

Used with `amount` field in add actions to specify how a material quantity is measured:

| Value | Description | Example Use Case |
|-------|-------------|------------------|
| `moleq` | Molar equivalents relative to limiting reagent | `"amount": 1.5, "quantity_unit": "moleq"` — 1.5 eq |
| `masseq` | Mass equivalents relative to limiting reagent | `"amount": 2.0, "quantity_unit": "masseq"` — 2× mass |
| `uL` | Microliters (absolute volume) | `"amount": 50, "quantity_unit": "uL"` — 50 µL |
| `mL` | Milliliters (absolute volume) | `"amount": 0.5, "quantity_unit": "mL"` — 0.5 mL |
| `mg` | Milligrams (absolute mass) | `"amount": 10, "quantity_unit": "mg"` — 10 mg |
| `g` | Grams (absolute mass) | `"amount": 0.1, "quantity_unit": "g"` — 100 mg |
| `M` | Molarity (moles per liter) | `"amount": 0.1, "quantity_unit": "M"` — 0.1 M |
| `uM` | Micromolarity (µmol per liter) | `"amount": 100, "quantity_unit": "uM"` — 100 µM |

**Default:** `moleq`

### Duration Units (`duration_unit`)

Used with `duration` field in stir actions:

| Value | Description | Example |
|-------|-------------|---------|
| `s` | Seconds | `"duration": 30, "duration_unit": "s"` — 30 seconds |
| `m` | Minutes | `"duration": 15, "duration_unit": "m"` — 15 minutes |
| `h` | Hours | `"duration": 12, "duration_unit": "h"` — 12 hours |

**Default:** `h`

### Temperature Units (`temperature_unit`)

Used with `temperature` field in stir actions:

| Value | Description |
|-------|-------------|
| `degC` | Degrees Celsius (default) |
| `K` | Kelvin |

### Plate Roles

Valid values for `plate_role`, `from_plate_role`, `to_plate_role`:

| Value | Description |
|-------|-------------|
| `reaction` | Main reaction plate |
| `workup` | Workup plate (use with `plate_role_index` for multiple: workup 1, 2, etc.) |
| `spefilter` | SPE filter plate |
| `lcms` | LCMS sample plate |
| `xchem` | XChem crystallography plate |
| `nmr` | NMR sample plate |
| `startingmaterial` | Starting material source plate |
| `solvent` | Solvent reservoir plate |

## Field Reference

### Add Actions

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `amount` | Yes | — | Numeric amount — meaning depends on `quantity_unit` |
| `quantity_unit` | No | `moleq` | Unit for amount (see Units Reference above) |
| `solvent` | No | — | Solvent name (e.g., `DMA`, `THF`, `ACN`) |
| `concentration` | No | — | Solution concentration (mol/L) |
| `density` | No | — | Neat density (g/mL) for liquid reagents |
| `from_plate_role` | No | `startingmaterial` | Source plate role |
| `from_plate_role_index` | No | `1` | Source plate index (for multiple plates of same role) |
| `to_plate_role` | No | `reaction` | Destination plate role |
| `to_plate_role_index` | No | `1` | Destination plate index |

### Stir Actions

| Field | Required | Default | Description |
|-------|----------|---------|-------------|
| `temperature` | **Recommended** | `25` | Temperature value (always include for clarity) |
| `temperature_unit` | No | `degC` | Unit: `degC`, `K` |
| `duration` | Yes | — | Duration value |
| `duration_unit` | No | `h` | Unit: `s`, `m`, `h` |
| `stirring_speed` | No | `normal` | Speed: `gentle`, `normal`, `vigorous` |
| `plate_role` | No | `reaction` | Plate being stirred |
| `plate_role_index` | No | `1` | Plate role index |

## Plate Routing

- For simple recipes, omit all `*_plate_role_index` fields (default is 1).
- For multi-stage workup, specify `to_plate_role_index` or `from_plate_role_index` as needed:
  - `extract: reaction → workup (index 1)`
  - `extract: reaction → workup (index 2)`
  - `add: workup (index 2) → spefilter`

## Management Command Reference

- `load_recipe_json` supports:
  - `--dry-run`: Validate without writing
  - `--update`: Create new version if recipe exists
  - Multiple files at once

## Canonical Examples

- `amidation_standard.json`: Standard amidation, no explicit plate_role_index
- `mitsunobu_with_extraction.json`: Mitsunobu with multi-stage extraction, explicit plate_role_index routing

## See Also
- [../README.md](../README.md) for project overview
- [../management/commands/load_recipe_json.py](../management/commands/load_recipe_json.py) for loader implementation
