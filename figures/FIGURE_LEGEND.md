# Shared methods for Figures 1–3

**All observations in these figures are simulated. No experimental mouse measurements, biological findings, validated gene identities, or pathway results are shown.** All three figures, including Figure 1, have been rebuilt with an assumed small-effect, low-variance scenario. The simulation demonstrates plot appearance and an illustrative analysis workflow; it is not fitted to brain or hypothalamic tissue and does not estimate this study's biological variability, expected effects, or statistical power.

The group names are exactly **LD SALINE**, **LD LEPTIN**, and **NON-LD control**. All three study groups are on the LDLR−/− background; the control is non-lipodystrophic, not wild-type. There are six independent synthetic samples per group. The questions and legends are separate from the artwork:

| Figure | Question | Contrast | Legend |
|---|---|---|---|
| **FIGURE 1** | How does untreated LD differ from non-LD controls? | LD SALINE vs NON-LD control | [Figure 1 legend](FIGURE_1_LEGEND.md) |
| **FIGURE 2** | What does leptin treatment change in LD mice? | LD LEPTIN vs LD SALINE | [Figure 2 legend](FIGURE_2_LEGEND.md) |
| **FIGURE 3** | What differences remain after leptin treatment? | LD LEPTIN vs NON-LD control | [Figure 3 legend](FIGURE_3_LEGEND.md) |

## Common panel definitions

**a, Volcano plot.** Each point is one of **8,000 retained synthetic genes**. Positive log2 fold change indicates higher expression in the first named group. Coral and blue identify higher and lower expression, respectively, at Benjamini–Hochberg-adjusted *P* < 0.05; other genes are light gray. There is **no fold-change threshold**. The dashed horizontal line marks adjusted *P* = 0.05, and the vertical line at zero is a reference for effect direction. All figures use the same axes and display every retained gene. An effect of ±0.2 log2 units is the intended scale of this demonstration, not an axis boundary, significance cutoff, or universal biological standard. Fitted estimates are neither clipped nor rescaled and can exceed ±0.2.

**b, Expression heatmap.** For the relevant contrast, the 40 genes with the smallest adjusted *P* values among **all tested genes** are shown, without requiring significance or a minimum effect size. Columns are the 12 individual simulated mice from the two named groups, with six samples per block. For each selected gene, log2-normalized expression is standardized across these 12 samples to a row-wise z-score using the sample standard deviation (ddof = 1). Rows are ordered by average-linkage hierarchical clustering with Euclidean distance on these z-scores. The color scale is limited to −2 to 2. Row standardization shows within-gene patterns; it does not show absolute effect magnitudes. Group annotation-strip colors identify LD SALINE (coral), LD LEPTIN (dark sage), and NON-LD control (light blue-gray); the heatmap colors encode standardized expression. These are separate color encodings. The selections are contrast-specific and are supplied as `figure_01_heatmap_genes.tsv`, `figure_02_heatmap_genes.tsv`, and `figure_03_heatmap_genes.tsv` in the [data directory](synthetic_demo/data).

**c, Selected gene effect estimates.** The six genes with the smallest adjusted *P* values in each direction among **all tested genes** are shown, ordered by estimated effect, without requiring significance or a minimum effect size. Points are estimated log2 fold changes. Horizontal intervals are the coefficient ± *t*<sub>0.975,15</sub> × its standard error; the vertical reference line is zero. These are **marginal, unadjusted 95% confidence intervals** from the illustrative fitted model. The genes were selected using the same data, and the intervals are not corrected for multiplicity or selection; 95% coverage is not asserted for this selected set. Selecting six genes per direction balances the illustration and does not show the overall distribution of effects. No SIM identifier implies a real molecular identity. The selections are supplied as `figure_01_forest_genes.tsv`, `figure_02_forest_genes.tsv`, and `figure_03_forest_genes.tsv` in the [data directory](synthetic_demo/data).

## Synthetic data and illustrative statistics

With random seed **20261019**, 8,000 anonymous genes were simulated for 18 samples using a negative-binomial distribution with sample library-size variation. Baseline abundance is lognormal with median 2,000 and log-scale standard deviation 0.7, clipped to 300–20,000. Gene-specific dispersion is log-uniform from 0.001 to 0.03, independent of baseline abundance and generating class. Base nonzero effect magnitudes are uniform from 0.08 to 0.22 log2 units, with balanced positive and negative directions. Derived amplified effects may exceed 0.22, and fitted estimates retain sampling noise. These parameters define an **assumed small-effect, low-variance visualization scenario**, rather than a scenario calibrated to hypothalamic bulk RNA-seq.

The generating model includes unchanged effects, near-complete and partial reversals, persistent effects, treatment-only changes, and amplified effects. Those labels describe simulation parameters only; they are stored separately from the inferred contrast tables and are not used to classify plotted observations or select plotted genes.

Median-of-ratios size factors were estimated from genes with positive counts in every sample and rescaled to geometric mean 1. Genes with at least 10 raw counts in at least six samples were retained (**all 8,000 genes**). Expression was transformed as log2(count / size factor + 0.5). A common three-group ordinary least-squares model was fitted separately for each gene, using pooled residual variance and 15 residual degrees of freedom. Two-sided t-test *P* values were adjusted by the Benjamini–Hochberg procedure independently within each complete contrast family. Every displayed effect estimate is a difference between fitted means on this transformed scale.

All contrasts use the same fitted means and satisfy the exact effect identity: (LD LEPTIN vs NON-LD control) = (LD LEPTIN vs LD SALINE) + (LD SALINE vs NON-LD control), apart from floating-point rounding. The contrasts share model estimates and are statistically correlated; they are not independent experiments.

**This rapid log-scale analysis is an illustrative approximation, without empirical-Bayes moderation or precision weights. It is not DESeq2, edgeR, or limma-voom, and it does not establish the production analysis method for future experimental data.** The actual biological replicate counts, batch variables, covariates, sample quality, and experimental design will need to be supplied with the real dataset.

In this simulation only, adjusted *P* < 0.05 identifies **871 genes** in LD SALINE vs NON-LD control, **210 genes** in LD LEPTIN vs LD SALINE, and **571 genes** in LD LEPTIN vs NON-LD control. These numbers are not results or power predictions for the lipodystrophy experiment. The [JSON summary](synthetic_demo/data/simulation_summary.json) gives numerical checks, fitted and true effect distributions, and separate generating truth.

## Interpretation of the post-treatment comparison

Figure 3 shows differences **present after leptin treatment**. These may include persistent LD-associated differences and treatment-induced differences. The comparison alone does not distinguish those origins or establish that a difference was uncorrected. Cross-contrast effect estimates and uncertainty would be needed for that assessment with experimental data. Proximity to zero or nonsignificance after treatment does not establish equivalence, normalization, or rescue.

## Color accessibility

Signed effects retain blue (`#779AAC`) and coral (`#C28A80`). Group annotations use LD SALINE coral (`#C28E85`), LD LEPTIN dark sage (`#527D6B`), and NON-LD control light blue-gray (`#A8BED2`). The group-color lightness was adjusted following simulated color-vision-deficiency checks. Full group names provide identification independent of color. The checks improve separation under the simulated conditions; they do not certify accessibility for every viewer or reproduction condition. See the [color-accessibility review and previews](accessibility/COLOR_ACCESSIBILITY.md) for the simulation method and evidence.

## Artwork and reproducibility

Each 1536 × 1024 pixel figure uses Nimbus Sans throughout. PDF fonts are embedded and searchable; main SVGs use outlined lettering; editable SVGs contain embedded Nimbus Sans font data. High-resolution PNGs are 4608 × 3072 pixels. The revised export names are `figure_01_small_effects_synthetic`, `figure_02_small_effects_synthetic`, and `figure_03_small_effects_synthetic`, with the appropriate file extensions. Synthetic-data provenance is also stored in PNG/PDF/SVG metadata, without adding explanatory prose to the figure canvas. PCA is omitted from the figures.

From the repository root:

```bash
python -m pip install -r figures/requirements.txt
python figures/synthetic_demo/simulate_and_analyze.py
python figures/synthetic_demo/render_figure.py
```

Install the URW Nimbus Sans fonts (Debian/Ubuntu package `fonts-urw-base35`). On another platform, set `NIMBUS_SANS_FONT_DIR` to the folder containing `NimbusSans-Regular.otf`, `NimbusSans-Bold.otf`, and `NimbusSans-Italic.otf`. Missing fonts cause an error; they are never replaced silently. Font conversion for robust PDF embedding is cached outside the repository. The data script is seed-reproducible, including deterministic gzip timestamps. Plot selections, group labels, dimensions, and PDF font identities are audited in the supplied files.
