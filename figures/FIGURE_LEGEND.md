# Figure 1 | Synthetic demonstration of the hypothalamic bulk RNA-seq figure layout

**All observations in this figure are simulated. No experimental mouse measurements, biological findings, validated gene identities, or pathway results are shown.** The simulation demonstrates plot appearance and the proposed analysis workflow; it does not predict the outcome of the study. The group names are exactly **LD SALINE**, **LD LEPTIN**, and **NON-LD control**. All three study groups are on the LDLR−/− background; the control is non-lipodystrophic, not wild-type.

**a–c,** Volcano plots for LD LEPTIN vs LD SALINE, LD SALINE vs NON-LD control, and LD LEPTIN vs NON-LD control, respectively. Each point is one retained synthetic gene. Positive log2 fold change indicates higher expression in the first named group. Coral and blue identify higher and lower expression, respectively, at Benjamini–Hochberg-adjusted *P* < 0.05 and absolute estimated log2 fold change ≥ 1. Other genes are light gray. Dashed lines indicate these display thresholds. Scales are identical across the three panels; all retained genes are plotted without clipping the effect or significance range. The colors in these volcano panels encode effect direction, whereas the cohort colors in d and f encode group membership.

**d,** Principal-component analysis of the 1,500 retained genes with greatest sample variance after log2 normalization and gene centering, without variance scaling. Each symbol is one simulated mouse: LD SALINE, coral circles; LD LEPTIN, sage squares; NON-LD control, blue-gray triangles. The simulation includes six independent synthetic samples per group. Percentages on the axes indicate variance explained by the corresponding principal component. No confidence ellipse or inferential group-separation claim is made.

**e,** Estimated log2 fold changes for LD SALINE vs NON-LD control (horizontal) and LD LEPTIN vs NON-LD control (vertical), obtained from the same fitted model. Sage points identify genes meeting the panel-b display criteria; other genes are gray. The dashed diagonal indicates equal estimated effects in the two contrasts; horizontal and vertical reference lines indicate zero. The estimates share NON-LD control measurements and are statistically correlated. This is a descriptive comparison, not a test of rescue. Proximity to zero or nonsignificance after leptin does not establish equivalence or normalization, and genes are not assigned mechanistic response categories from this scatter plot.

**f,** Gene-centered expression heatmap for 60 retained genes ranked by the smallest adjusted *P* across the three contrasts. This ranking is a visualization-selection rule, not an additional cross-contrast significance test. Rows are ordered by average-linkage hierarchical clustering using Euclidean distance on gene-wise z-scores. Columns are individual simulated mice, grouped as LD SALINE, LD LEPTIN, and NON-LD control, with six samples in each block. Values are row-wise standardized log2-normalized expression (sample standard deviation, ddof = 1); the displayed color scale is limited to −2 to 2. No gene name implies a real molecular identity; the source tables use anonymous SIM identifiers. The selected gene list is supplied separately.

## Synthetic data and illustrative statistics

With random seed **20261019**, 8,000 anonymous genes were simulated for 18 samples using a negative-binomial distribution with gene-dependent abundance, gene-dependent overdispersion, and sample library-size variation. The generating model includes unchanged effects, near-complete and partial reversals, persistent effects, treatment-only changes, and amplified effects, with balanced positive and negative directions. Those labels describe simulation parameters only; they are stored separately from the inferred contrast tables and are not used to classify plotted observations.

Median-of-ratios size factors were estimated from genes with positive counts in every sample and rescaled to geometric mean 1. Genes with at least 10 raw counts in at least six samples were retained (**7,647 genes**). Expression was transformed as log2(count / size factor + 0.5). A common three-group ordinary least-squares model was fitted separately for each gene, using pooled residual variance and 15 residual degrees of freedom. Two-sided t-test *P* values were adjusted by the Benjamini–Hochberg procedure independently within each complete contrast family. All contrasts use the same fitted means and satisfy the exact effect identity: (LD LEPTIN vs NON-LD control) = (LD LEPTIN vs LD SALINE) + (LD SALINE vs NON-LD control), apart from floating-point rounding.

**This rapid log-scale analysis is an illustrative approximation, without empirical-Bayes moderation or precision weights. It is not DESeq2, edgeR, or limma-voom, and it does not establish the production analysis method for future experimental data.** The actual biological replicate counts, batch variables, covariates, sample quality, and experimental design will need to be supplied with the real dataset.

In this simulation only, the display criteria select 672 genes in LD LEPTIN vs LD SALINE, 1,690 in LD SALINE vs NON-LD control, and 1,094 in LD LEPTIN vs NON-LD control. These numbers are not results of the lipodystrophy experiment. The JSON summary gives numerical checks and separate generating truth.

## Artwork and reproducibility

The 1536 × 1024 pixel figure uses Nimbus Sans throughout. PDF fonts are embedded and searchable; the main SVG uses outlined lettering; the editable SVG contains embedded Nimbus Sans font data. The high-resolution PNG is 4608 × 3072 pixels. Synthetic-data provenance is also stored in PNG/PDF/SVG metadata, without adding explanatory prose to the figure canvas.

From the repository root:

```bash
python -m pip install -r figures/requirements.txt
python figures/synthetic_demo/simulate_and_analyze.py
python figures/synthetic_demo/render_figure.py
```

Install the URW Nimbus Sans fonts (Debian/Ubuntu package `fonts-urw-base35`). On another platform, set `NIMBUS_SANS_FONT_DIR` to the folder containing `NimbusSans-Regular.otf`, `NimbusSans-Bold.otf`, and `NimbusSans-Italic.otf`. Missing fonts cause an error; they are never replaced silently. Font conversion for robust PDF embedding is cached outside the repository. The data script is seed-reproducible, including deterministic gzip timestamps. Plot selections, group labels, dimensions, and PDF font identities are audited in the supplied files.
