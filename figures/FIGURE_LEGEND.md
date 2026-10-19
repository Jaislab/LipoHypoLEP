# Shared methods for Figures 1–3

**All observations in these figures are simulated. No experimental mouse measurements, biological findings, validated gene identities, or pathway results are shown.** The simulation demonstrates plot appearance and a proposed analysis workflow; it does not predict the outcome of the study.

The group names are exactly **LD SALINE**, **LD LEPTIN**, and **NON-LD control**. All three study groups are on the LDLR−/− background; the control is non-lipodystrophic, not wild-type. There are six independent synthetic samples per group. The questions and legends are separate from the artwork:

| Figure | Question | Contrast | Legend |
|---|---|---|---|
| **FIGURE 1** | How does untreated LD differ from non-LD controls? | LD SALINE vs NON-LD control | [Figure 1 legend](FIGURE_1_LEGEND.md) |
| **FIGURE 2** | What does leptin treatment change in LD mice? | LD LEPTIN vs LD SALINE | [Figure 2 legend](FIGURE_2_LEGEND.md) |
| **FIGURE 3** | What differences remain after leptin treatment? | LD LEPTIN vs NON-LD control | [Figure 3 legend](FIGURE_3_LEGEND.md) |

## Common panel definitions

**a, Volcano plot.** Each point is one of **7,647 retained synthetic genes**. Positive log2 fold change indicates higher expression in the first named group. Coral and blue identify higher and lower expression, respectively, at Benjamini–Hochberg-adjusted *P* < 0.05 and absolute estimated log2 fold change ≥ 1. Other genes are light gray. Dashed lines indicate these display criteria. The absolute fold-change rule is applied to the estimated coefficient; the *P* values test a zero effect, not a fold-change threshold. All figures use the same horizontal range, −5 to 5, and vertical range, −0.08 to 10; every retained gene is plotted within those limits.

**b, Expression heatmap.** For the relevant contrast, the 40 genes with the smallest adjusted *P* values among genes satisfying the panel-a display criteria are shown. Columns are the 12 individual simulated mice from the two named groups, with six samples per block. For each selected gene, log2-normalized expression is standardized across these 12 samples to a row-wise z-score using the sample standard deviation (ddof = 1). Rows are ordered by average-linkage hierarchical clustering with Euclidean distance on these z-scores. The color scale is limited to −2 to 2. Group annotation-strip colors identify LD SALINE (coral), LD LEPTIN (sage), and NON-LD control (blue-gray); the heatmap colors encode standardized expression. These are separate color encodings. The selections are contrast-specific and are supplied as `figure_01_heatmap_genes.tsv`, `figure_02_heatmap_genes.tsv`, and `figure_03_heatmap_genes.tsv` in the [data directory](synthetic_demo/data).

**c, Selected gene effect estimates.** The six genes with the smallest adjusted *P* values in each direction among genes satisfying the panel-a display criteria are shown, ordered by estimated effect. Points are estimated log2 fold changes. Horizontal intervals are the coefficient ± *t*<sub>0.975,15</sub> × its standard error; the vertical reference line is zero. These are **marginal, unadjusted 95% confidence intervals** from the illustrative fitted model. The genes were selected using the same data, and the intervals are not corrected for multiplicity or selection; 95% coverage is not asserted for this selected set. Selecting six genes per direction balances the illustration and does not show the overall distribution of effects. No SIM identifier implies a real molecular identity. The selections are supplied as `figure_01_forest_genes.tsv`, `figure_02_forest_genes.tsv`, and `figure_03_forest_genes.tsv` in the [data directory](synthetic_demo/data).

## Synthetic data and illustrative statistics

With random seed **20261019**, 8,000 anonymous genes were simulated for 18 samples using a negative-binomial distribution with gene-dependent abundance, gene-dependent overdispersion, and sample library-size variation. The generating model includes unchanged effects, near-complete and partial reversals, persistent effects, treatment-only changes, and amplified effects, with balanced positive and negative directions. Those labels describe simulation parameters only; they are stored separately from the inferred contrast tables and are not used to classify plotted observations or select plotted genes.

Median-of-ratios size factors were estimated from genes with positive counts in every sample and rescaled to geometric mean 1. Genes with at least 10 raw counts in at least six samples were retained (**7,647 genes**). Expression was transformed as log2(count / size factor + 0.5). A common three-group ordinary least-squares model was fitted separately for each gene, using pooled residual variance and 15 residual degrees of freedom. Two-sided t-test *P* values were adjusted by the Benjamini–Hochberg procedure independently within each complete contrast family. Every displayed effect estimate is a difference between fitted means on this transformed scale.

All contrasts use the same fitted means and satisfy the exact effect identity: (LD LEPTIN vs NON-LD control) = (LD LEPTIN vs LD SALINE) + (LD SALINE vs NON-LD control), apart from floating-point rounding. The contrasts share model estimates and are statistically correlated; they are not independent experiments.

**This rapid log-scale analysis is an illustrative approximation, without empirical-Bayes moderation or precision weights. It is not DESeq2, edgeR, or limma-voom, and it does not establish the production analysis method for future experimental data.** The actual biological replicate counts, batch variables, covariates, sample quality, and experimental design will need to be supplied with the real dataset.

In this simulation only, the display criteria select **1,690 genes** in LD SALINE vs NON-LD control, **672 genes** in LD LEPTIN vs LD SALINE, and **1,094 genes** in LD LEPTIN vs NON-LD control. These numbers are not results of the lipodystrophy experiment. The [JSON summary](synthetic_demo/data/simulation_summary.json) gives numerical checks and separate generating truth.

## Interpretation of the post-treatment comparison

Figure 3 shows differences **present after leptin treatment**. These may include persistent LD-associated differences and treatment-induced differences. The comparison alone does not distinguish those origins or establish that a difference was uncorrected. Cross-contrast effect estimates and uncertainty would be needed for that assessment with experimental data. Proximity to zero or nonsignificance after treatment does not establish equivalence, normalization, or rescue.

## Artwork and reproducibility

Each 1536 × 1024 pixel figure uses Nimbus Sans throughout. PDF fonts are embedded and searchable; main SVGs use outlined lettering; editable SVGs contain embedded Nimbus Sans font data. High-resolution PNGs are 4608 × 3072 pixels. Synthetic-data provenance is also stored in PNG/PDF/SVG metadata, without adding explanatory prose to the figure canvas. PCA is omitted from the figures.

From the repository root:

```bash
python -m pip install -r figures/requirements.txt
python figures/synthetic_demo/simulate_and_analyze.py
python figures/synthetic_demo/render_figure.py
```

Install the URW Nimbus Sans fonts (Debian/Ubuntu package `fonts-urw-base35`). On another platform, set `NIMBUS_SANS_FONT_DIR` to the folder containing `NimbusSans-Regular.otf`, `NimbusSans-Bold.otf`, and `NimbusSans-Italic.otf`. Missing fonts cause an error; they are never replaced silently. Font conversion for robust PDF embedding is cached outside the repository. The data script is seed-reproducible, including deterministic gzip timestamps. Plot selections, group labels, dimensions, and PDF font identities are audited in the supplied files.
