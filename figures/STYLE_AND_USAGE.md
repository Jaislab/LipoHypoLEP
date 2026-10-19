# Figure style and use

Two coordinated, data-free figures are supplied at **1536 × 1024 pixels (3:2)**.

1. `01_study_design`: experimental groups, tissue-to-RNA-seq workflow, three directional contrasts, and the intended joint interpretation.
2. `02_volcano_template`: empty, labeled volcano-plot panels for L vs S, S vs K, and L vs K. No observations, cutoffs, counts, or enrichment results have been invented.

## Visual specification

- White background; charcoal type; thin gray rules and axes.
- S, LD + saline: muted coral `#D3A297`, light tint `#F7EDE9`.
- L, LD + leptin: sage `#9CB8A4`, light tint `#EDF3EE`.
- K, non-LD control: blue-gray `#A4B5C3`, light tint `#EEF2F5`.
- **Nimbus Sans only**, using the actual URW regular and bold font files. The script checks every label's font and checks that every glyph is available. It fails if the required font is missing; it never substitutes another family.
- Main title 44 pixels; panel titles 29 pixels; major labels 26–28 pixels; body labels 21–24 pixels. Subscripts and superscripts are deliberately smaller.
- SVG and PDF text is outlined to preserve letterforms on other computers and in journal production. Source text remains editable in `make_figures.py`. The font audit records all labels.
- Native PNG is 1536 × 1024. The additional 300 dpi PNG is 4608 × 3072. Prefer vector PDF/SVG for layout; the anatomical illustration is an embedded raster asset, while typography and diagram elements are vector.

## Scientific conventions

- Positive log2 fold change means higher expression in the first named group.
- All cohorts are LDLR−/−; K is a non-lipodystrophic control, not wild-type.
- L vs K is labeled **Differences after leptin**. It can include persistent LD differences and treatment-associated differences.
- A nonsignificant L vs K comparison alone does not establish normalization. Joint categories require effect sizes, directions, uncertainty, and a stated statistical rule.
- The group-joining lines denote a common tissue and assay workflow, not pooling or paired sampling.
- Mouse-brain artwork is an AI-generated conceptual anatomical schematic and should receive anatomical review before publication. The green highlight marks the approximate hypothalamic region; it is not a measured spatial expression result or an anatomical reference atlas.
- The supplied transgene spelling should be confirmed before it is typeset in a publication caption. It has not been silently corrected or asserted in the figure.

## Regenerate

Run the following command from the repository root:

```bash
MPLCONFIGDIR=/tmp/lipohypolep-mpl python figures/make_figures.py
```

Install the Python dependencies with `python -m pip install -r figures/requirements.txt`. On Debian/Ubuntu, the `fonts-urw-base35` package supplies Nimbus Sans at the default font path. On another platform, set `NIMBUS_SANS_FONT_DIR` to the directory containing `NimbusSans-Regular.otf` and `NimbusSans-Bold.otf`; no substitute font is permitted. The `assets` directory must remain beside the script. Regeneration updates the exported figures and font audit in `figures/`.

When DEG tables arrive, replace placeholders with observations from those tables; select and disclose significance thresholds and adjust shared scales from the observed data. Do not treat this template as a completed data analysis.
