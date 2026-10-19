# Experimental CSV input contract and read-only preflight

The planned handoff is **differential-expression results with sample-level FPKM columns included in the same CSV**. Raw counts are not required to draw the supplied DE estimates or run preranked GSEA from supplied signed statistics. A separate sample metadata CSV identifies the FPKM columns and study groups.

[validate_inputs.py](../analysis/validate_inputs.py) checks file structure and declared semantics. It does **not** run statistical analysis, read inputs into the demonstration renderer, normalize expression, compute new fold changes, change a sign, adjust a P value, strip gene versions, map identifiers, impute missing values, or modify any input. The existing simulation scripts remain demonstration code; this validator does not turn them into a production experimental-analysis pipeline.

## Three comparisons and direction

Supply a separate complete DE table for each available comparison. Positive `log2FoldChange` and positive signed `stat` must refer to higher expression in the numerator:

| Figure | Numerator | Denominator |
|---|---|---|
| 1 | `LD SALINE` | `NON-LD control` |
| 2 | `LD LEPTIN` | `LD SALINE` |
| 3 | `LD LEPTIN` | `NON-LD control` |

The contrast names are explicit CLI arguments; filenames and column contents are never used to guess direction. An opposite sign between nonzero `log2FoldChange` and `stat` fails preflight. An unsigned F statistic, unsigned likelihood-ratio statistic, absolute z score, or a statistic from a different contrast is not a valid signed ranking metric. The validator cannot prove the model or contrast merely from consistent signs; the bioinformatician must document their meaning upstream. It does not reconstruct a signed statistic from an unsigned one.

## DE CSV fields

The first column must be exactly `gene_id`. Required fields are:

| Field | Contract |
|---|---|
| `gene_id` | Unique, nonempty, exact gene identifier. Declare `--gene-id-type ensembl_mouse`, `symbol`, or `entrez_mouse`. |
| `log2FoldChange` | Supplied finite signed effect estimate for the declared contrast. |
| `stat` | Supplied finite signed test/ranking statistic with the same orientation. |
| `padj` | Supplied adjusted P value between 0 and 1, or an explicitly reported missing value as described below. |
| `pvalue` | Optional raw P-value column; if present, values must be between 0 and 1 except documented untestable rows. |
| `baseMean`, `lfcSE` | Optional recognized DE columns; if present, values must be finite and nonnegative except documented untestable rows. |
| Sample columns | With `--sample-expression fpkm`, one finite, nonnegative numeric column per metadata `sample_id`, named exactly by that identifier. |

All supplied DE statistics remain unchanged. The figure display rule is adjusted P < 0.05 **and** absolute log2 fold change ≥ 0.2. The absolute effect cutoff selects existing estimates; it neither shrinks effects nor changes the statistical null to a minimum effect of 0.2. The preflight does not apply that display filter or generate a DEG list.

For inline FPKM, headers other than required/recognized DE fields and metadata sample IDs must be declared using repeated `--annotation-column NAME` arguments. For example, declare `--annotation-column gene_symbol` if it is an annotation alongside Ensembl `gene_id`. This prevents a forgotten sample column from silently being treated as an annotation. Annotation values are not interpreted or transformed. Without the inline-FPKM declaration, extra columns are reported as unvalidated and heatmap-input readiness remains false.

Mouse Ensembl IDs must match `ENSMUSG` plus 11 digits, optionally followed by a numeric version suffix; versions are retained. Entrez IDs must be positive, unpadded integer identifiers. Symbol identity and species cannot be established from syntax alone. No namespace conversion, orthology conversion, alias repair, version stripping, or duplicate-gene aggregation occurs.

CSV must use UTF-8, optionally with a UTF-8 BOM, and commas. Quoting follows CSV rules. Duplicate/blank headers, duplicate gene IDs, wrong row widths, empty files, control characters in identifiers, and outer whitespace in identifiers or header names fail. Decimal and scientific notation are accepted for numeric cells, without surrounding whitespace or grouping separators. `NaN`, infinities, and nonnumeric values fail rather than being coerced. Finite decimals that would overflow float64 or turn a nonzero value into zero through float64 underflow also fail; source values are never replaced. No rows are silently discarded.

## Missing statistics and complete rankings

`NA` and an empty cell are the only recognized missing tokens; spelling and case are not repaired.

- **Missing `padj`:** allowed and reported with gene IDs and counts. Independent filtering can legitimately leave an adjusted P value missing. Such rows cannot enter DEG highlighting, but a finite signed statistic remains eligible for a GSEA ranking. The reason for missing adjusted P values must be documented upstream; the validator does not assume a reason or replace them with 1.
- **Missing `stat`:** fails by default. Only `--missing-stat-reason 'upstream explanation'` permits these rows as explicitly documented untestable entries. Their IDs and counts are reported as excluded from ranking. Under that same documented untestability condition, missing effects, raw P values, and optional numeric DE fields are retained and reported/allowed as appropriate. No missing statistic or effect is invented.
- **Finite zero P values:** accepted exactly and reported with a warning. A future plotting convention for a logarithmic axis must be disclosed separately; the validator does not replace zero values.

Preranked GSEA needs the **complete unfiltered tested-gene results**, not only significant DEGs. `--complete-tested-universe` records upstream confirmation of completeness. Without that flag, a valid DE file can pass structural preflight, but `readiness.gsea_rank_input` remains false. The validator cannot infer completeness from the number of rows or P-value distribution.

The ranking uses all available finite, correctly signed supplied statistics, including rows with missing `padj`; it must not use the figure's adjusted-P or fold-change filters. Documented missing-stat entries remain explicit exclusions. Different contrasts can legitimately have different tested universes because of upstream filtering or testability; this CLI validates one file at a time and never requires equal gene lists across the three files. Gene-set overlap, annotation coverage, species identity, and mapping quality still require later review against an appropriate mouse annotation release.

## Sample metadata and FPKM

Metadata must contain `sample_id` and `group` columns. Sample IDs must be unique and match inline FPKM column names exactly; there is no required prefix and no renaming. Include **every sample column in the file**, including samples from the third group if all groups are included in each DE CSV. Later question-specific heatmaps would use the two relevant groups.

Group labels must be exactly `LD SALINE`, `LD LEPTIN`, or `NON-LD control`; letter codes and spelling variants are not repaired. For a DE file, metadata can contain its two compared groups or all three groups. The numerator and denominator must each have at least two independent biological replicates; the third group is optional, but needs at least two when present. Standalone matrix preflight requires all three study groups. At least three independent biological replicates per represented group are recommended. The validator counts sample rows and warns that biological independence cannot be established from `sample_id` and `group`: one row must represent one independent biological mouse, rather than a technical replicate. Additional covariate columns are preserved but are not modeled by this preflight. Batch, sex, technical-replicate relationships, and other design factors must be described upstream.

FPKM columns must be explicitly declared with `--sample-expression fpkm` and paired metadata. All sample cells must be finite and nonnegative; missing or negative sample values are not imputed. A declaration with absent sample columns fails. Sample IDs cannot reuse recognized DE-statistic headers, because one column cannot supply both a statistic and sample expression.

FPKM is **display-only expression**, not input for refitting differential expression. A future heatmap can use a disclosed `log2(FPKM + c)` transformation followed by a disclosed row standardization, with `c` chosen only after inspecting the actual value distribution. This preflight performs neither transformation and selects no pseudocount. Neither supplied log2 fold changes/adjusted P values nor GSEA ranking statistics will be recomputed from FPKM.

## Example invocations

For Figure 1 with complete results and inline FPKM:

```bash
python analysis/validate_inputs.py LD_SALINE_vs_NON_LD_control.csv \
  --kind deg --gene-id-type ensembl_mouse \
  --numerator 'LD SALINE' --denominator 'NON-LD control' \
  --complete-tested-universe \
  --sample-expression fpkm --metadata samples.csv \
  --upstream-method 'Document the model and software/version' \
  --statistic-description 'Document the supplied signed statistic and contrast' \
  --adjustment-method 'Document adjustment method and tested family' \
  --annotation-release 'Document the mouse identifier release'
```

Use the corresponding numerator and denominator from the table for Figures 2 and 3. The provenance arguments record descriptions supplied by the bioinformatician; they are not automatically verified. If there is an additional annotation column, append `--annotation-column its_exact_header`. If untestable rows have missing statistics, append `--missing-stat-reason 'the documented upstream reason'`. Do not use that option to hide inconsistent or unsigned statistics.

For a DE-only file, omit both `--sample-expression` and `--metadata`. Such a file can supply DE/GSEA statistics, but cannot supply a sample-expression heatmap.

Secondary matrix preflights are also available:

```bash
python analysis/validate_inputs.py raw_counts.csv \
  --kind counts --gene-id-type ensembl_mouse --metadata samples.csv
python analysis/validate_inputs.py normalized_expression.csv \
  --kind expression --gene-id-type ensembl_mouse --metadata samples.csv
```

Matrices use `gene_id` followed only by sample columns, with an exact set match to metadata. `counts` means nonnegative integer raw counts, not TPM or FPKM; fractional values fail and are never rounded. Integer-valued numbers alone cannot prove units, so the raw-count declaration must be checked upstream. `expression` permits finite nonnegative fractional values for display; its units and any previous transformation must be supplied separately. Counts do not become heatmap-ready normalized expression by passing this check.

## JSON report and limits

The CLI writes JSON to **standard output only** and never creates a report file or overwrites inputs. Exit status is 0 for a structurally valid input, 1 for validation errors, and 2 for CLI usage errors. The report records:

- SHA-256 and byte length of the exact data/metadata bytes parsed, including any BOM.
- Kind, species/identifier declarations, contrast direction, completeness confirmation, expression unit, annotation declarations, and supplied upstream provenance.
- Row and unique-gene counts, sample/group counts where applicable, finite-statistic/plottable-row counts, and explicit missing-value exclusions by gene ID.
- Structured errors and warnings. Each list is limited to the first 100 issues; the corresponding total count records all issues.
- Separate `deg_plot_input`, `gsea_rank_input`, and `heatmap_expression_input` structural-readiness flags. A valid DE file without completeness confirmation is not GSEA-ready; a file without sample-expression validation is not heatmap-ready.

These readiness fields describe **input structure**, not a completed scientific analysis or verified identifier mapping. Hashes connect the report to exact submitted bytes but do not prove biological correctness. All rows remain in the source file, including explicitly excluded untestable rows. The CLI does not evaluate cross-file model consistency, power, sample quality, the truth of declared biological replicates, or the correctness of the upstream model/FDR method.

Run the boundary tests from the repository root:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s tests -p 'test_inputs.py' -v
```

The tests use temporary artificial CSV fixtures. They validate parsing, missingness rules, signed-statistic contradictions, exact metadata matching, inline-FPKM checks, numeric ranges, raw-count integrality, hashes, unchanged input bytes, and CLI JSON/exit behavior. They do not establish a production differential-expression or enrichment workflow.
