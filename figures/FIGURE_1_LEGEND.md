# FIGURE 1 | How does untreated LD differ from non-LD controls?

**Entirely synthetic demonstration. No experimental measurements or biological conclusions are shown.** The comparison is **LD SALINE vs NON-LD control**, with six simulated mice per group. All mice are on the LDLR−/− background. Positive log2 fold change indicates higher expression in LD SALINE.

**a,** Volcano plot of 7,647 retained synthetic genes. Coral and blue identify higher and lower expression in LD SALINE, respectively, at adjusted *P* < 0.05 and absolute estimated log2 fold change ≥ 1; other genes are gray. Dashed lines mark these display criteria. In this simulation, 1,690 genes meet both criteria: 868 higher and 822 lower. These counts are illustrative, not study results.

**b,** Expression heatmap of the 40 genes with the smallest adjusted *P* values among those meeting the panel-a criteria, restricted to the 12 LD SALINE and NON-LD control samples. Rows show gene-wise z-scores of log2-normalized expression across these samples (sample standard deviation, ddof = 1) and are ordered by average-linkage Euclidean clustering. The displayed scale is limited to −2 to 2. Columns are individual mice, grouped by the full study-group names.

**c,** Estimated effects for the six genes with the smallest adjusted *P* values in each direction among those meeting the panel-a criteria, sorted by effect. Intervals are marginal, unadjusted 95% confidence intervals from the common fitted model (coefficient ± *t*<sub>0.975,15</sub> × standard error). Gene selection uses the same data; these intervals do not provide multiplicity- or selection-adjusted coverage for the displayed set. SIM labels are anonymous synthetic identifiers.

See [shared methods and limitations](FIGURE_LEGEND.md) for simulation, normalization, the illustrative log-scale model, multiple-testing adjustment, and display-selection rules. Gene selections: [heatmap](synthetic_demo/data/figure_01_heatmap_genes.tsv) · [effect estimates](synthetic_demo/data/figure_01_forest_genes.tsv).
