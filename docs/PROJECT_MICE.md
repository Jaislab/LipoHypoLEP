# Project mouse register

Recorded from the supplied `mice_genotypes.xlsx` workbook. Sex is recorded by worksheet (`female` or `male`). All 30 mouse IDs are unique across the two worksheets, and every record has a genotype and group.

| Group | Genotype as supplied | Female | Male | Total |
|---|---|---:|---:|---:|
| LD SALINE | `tg/LDLR` | 6 | 4 | 10 |
| LD LEPTIN | `tg/LDLR` | 6 | 4 | 10 |
| NON-LD control | `LDLR` | 6 | 4 | 10 |
| **Total** | | **18** | **12** | **30** |

Each group has the same sex composition: six females and four males.

## Individual mice

| Group | Female mouse IDs | Male mouse IDs |
|---|---|---|
| LD SALINE | 1021, 1076, 1122, 1129, 1131, 1180 | 1026, 1036, 1146, 1149 |
| LD LEPTIN | 1073, 1086, 1088, 1115, 1124, 1177 | 1029, 1030, 1071, 1104 |
| NON-LD control | 1022, 1097, 1130, 1156, 1159, 1160 | 1046, 1052, 1061, 1103 |

## Data and provenance

The [machine-readable register](../data/project/mice.csv) preserves workbook row order and contains:

- `mouse_id`: the supplied animal ID, represented as an exact decimal string.
- `sex`: the source worksheet name.
- `group`: the established project group label.
- `genotype_reported` and `group_reported`: the original cell values, unchanged.
- `source_sheet` and `source_row`: the location of the source record, including the header in row numbering.

Only group-label capitalization is standardized: `LD saline` → `LD SALINE`, `LD leptin` → `LD LEPTIN`, and `Non-LD control` → `NON-LD control`. Genotype abbreviations are retained exactly; no additional allele or zygosity information is inferred.

The [provenance record](../data/project/mice_provenance.json) includes the original workbook's filename and SHA-256, the CSV checksum, the label mapping and all group/sex counts. The original workbook remains local. All 30 animal records are retained, with no exclusions or invented values.

These are project animal counts. The workbook does not map mouse IDs to expression-file columns or identify which animals contributed to the historical DEG workbook. Sequencing metadata will need exact `sample_id`–`mouse_id` links when the expression files arrive; see the [input contract](INPUT_DATA_CONTRACT.md#sample-metadata-and-fpkm).
