# Shared methods for Figures 1–3

**All expression measurements and differential-expression statistics are synthetic. No experimental mouse measurements or tissue-derived biological findings are shown.** The figures retain the original synthetic simulation, with no automatic reduction or rescaling of effect magnitudes. It is not fitted to brain or hypothalamic tissue and does not estimate this study's biological variability, expected effects, or statistical power. Enrichment uses real native-mouse pathway memberships and an actual preranked GSEA calculation, with real gene symbols assigned randomly to the synthetic statistics solely for demonstration.

The group names are exactly **LD SALINE**, **LD LEPTIN**, and **NON-LD control**. All three study groups are on the LDLR−/− background; the control is non-lipodystrophic, not wild-type. There are six independent synthetic samples per group. The questions and legends are separate from the artwork:

| Figure | Question | Contrast | Legend |
|---|---|---|---|
| **FIGURE 1** | How does untreated LD differ from non-LD controls? | LD SALINE vs NON-LD control | [Figure 1 legend](FIGURE_1_LEGEND.md) |
| **FIGURE 2** | What does leptin treatment change in LD mice? | LD LEPTIN vs LD SALINE | [Figure 2 legend](FIGURE_2_LEGEND.md) |
| **FIGURE 3** | What differences remain after leptin treatment? | LD LEPTIN vs NON-LD control | [Figure 3 legend](FIGURE_3_LEGEND.md) |

## Common panel definitions

**a, Volcano plot.** Each point is one of **7,647 retained synthetic genes**. Positive log2 fold change indicates higher expression in the first named group. Warm and cool colors identify higher and lower expression, respectively, at **Benjamini–Hochberg-adjusted *P* < 0.05 and absolute estimated log2 fold change ≥ 0.2**. Other genes are light gray, including genes that meet only one criterion. Dashed lines mark adjusted *P* = 0.05 and the effect-size cutoffs −0.2 and 0.2; the line at zero is a direction reference. All figures use common axes and display every retained gene. The ±0.2 values are effect-size cutoffs, not axis boundaries: fitted estimates are unchanged and can exceed them. The *P* values test a zero effect; combining them with an estimated-effect cutoff does not turn them into tests against a 0.2 effect-size null.

**b, Expression heatmap.** For the relevant contrast, the 40 genes with the smallest adjusted *P* values among genes satisfying **both panel-a criteria** are shown. Columns are the 12 individual simulated mice from the two named groups, with six samples per block. For each selected gene, log2-normalized expression is standardized across these 12 samples to a row-wise z-score using the sample standard deviation (ddof = 1). Rows are ordered by average-linkage hierarchical clustering with Euclidean distance on these z-scores. The color scale is limited to −2 to 2, from the comparison's cool color through white to its warm color. Row standardization shows within-gene patterns; it does not show absolute effect magnitudes. The first group annotation strip uses the warm color and the second the cool color; full group labels identify each block. These colors are comparison-specific, rather than fixed colors for cohorts across the figure series. Contrast-specific selections are supplied as `figure_01_heatmap_genes.tsv`, `figure_02_heatmap_genes.tsv`, and `figure_03_heatmap_genes.tsv` in the [data directory](synthetic_demo/data).

**c, KEGG preranked GSEA; d, GO Biological Process preranked GSEA.** Both calculations use **all 7,647 retained genes ranked in descending order by the signed test statistic**, without adjusted-*P* or fold-change filtering. Horizontal position is normalized enrichment score (NES); positive NES indicates enrichment toward the first named group's higher-expression end of the ranking. Each panel shows the two largest positive and two most negative NES values in its collection, regardless of significance. This is a fixed display rule, not a selection of pathways to support a mechanism. Marker fill indicates empirical GSEA FDR q < 0.05: filled below that threshold, hollow otherwise. The adjacent q column reports GSEApy's empirical GSEA FDR, **not Benjamini–Hochberg adjusted *P***. Colors indicate NES direction. Displayed pathways need not be significant, and even a small q value in this random synthetic assignment is not biological evidence.

## Synthetic data and illustrative gene-level statistics

With random seed **20261019**, 8,000 anonymous genes were simulated for 18 samples using a negative-binomial distribution with sample library-size variation. Baseline abundance is lognormal with median 80 and log-scale standard deviation 1.3, clipped to 2–12,000. Gene-specific dispersion is (0.04 + 0.25 / (1 + baseline / 100)) multiplied by a lognormal draw with log mean 0 and log standard deviation 0.18. Base nonzero effect magnitudes are uniform from 0.8 to 2.2 log2 units, with balanced positive and negative directions; class-specific multipliers can reduce or amplify them. Fitted estimates retain sampling noise. These are the original demonstration parameters, not estimates calibrated to hypothalamic bulk RNA-seq. They are distinct from the absolute estimated log2 fold-change display cutoff of 0.2, and no effect magnitudes are rescaled to that cutoff.

The generating model includes unchanged effects, near-complete and partial reversals, persistent effects, treatment-only changes, and amplified effects. Those labels describe simulation parameters only; they are stored separately from the inferred contrast tables and are not used to classify plotted observations or select plotted genes or pathways.

Median-of-ratios size factors were estimated from genes with positive counts in every sample and rescaled to geometric mean 1. Genes with at least 10 raw counts in at least six samples were retained (**7,647 genes**). Expression was transformed as log2(count / size factor + 0.5). A common three-group ordinary least-squares model was fitted separately for each gene, using pooled residual variance and 15 residual degrees of freedom. Two-sided t-test *P* values were adjusted by the Benjamini–Hochberg procedure independently within each complete contrast family. Every displayed effect estimate is a difference between fitted means on this transformed scale.

All contrasts use the same fitted means and satisfy the exact effect identity: (LD LEPTIN vs NON-LD control) = (LD LEPTIN vs LD SALINE) + (LD SALINE vs NON-LD control), apart from floating-point rounding. The contrasts share model estimates and are statistically correlated; they are not independent experiments.

**These scripts generate and analyze synthetic data only. A reader and production analysis pipeline for experimental data have not been implemented or tested here.** The pooled ordinary least-squares analysis of log2(normalized counts + 0.5) is an illustrative approximation, without empirical-Bayes moderation or precision weights. Its estimates depend on the transformation and pseudocount and can compress effects at low expression; it is unsuitable as the production method for estimating experimental effects. It is not DESeq2, edgeR, or limma-voom. Future experimental plots should use validated effect estimates and signed ranking statistics supplied by an appropriate count-model analysis, such as DESeq2 or edgeR, with the actual replicate counts, batches, covariates, quality checks, and study design accounted for.

In this simulation, the joint gene-display criteria identify **1,826 genes** in LD SALINE vs NON-LD control (**936 higher, 890 lower**), **789 genes** in LD LEPTIN vs LD SALINE (**404 higher, 385 lower**), and **1,226 genes** in LD LEPTIN vs NON-LD control (**653 higher, 573 lower**). The [JSON summary](synthetic_demo/data/simulation_summary.json) records joint-cutoff and FDR-only counts separately. These are synthetic illustration results, not experimental results or power predictions for the lipodystrophy experiment.

Supporting gene-effect tables remain available but are **not panels in the current figures**. They contain six genes per direction ranked by adjusted *P* among genes passing both display cutoffs. Intervals are the estimated coefficient ± *t*<sub>0.975,15</sub> × standard error: marginal, unadjusted 95% intervals from the common model. They are not adjusted for multiplicity or selection from the same data, and 95% coverage is not asserted for the selected set. SIM labels are anonymous synthetic identifiers.

## Preranked GSEA computation and annotation provenance

The [GSEA source](synthetic_demo/run_gsea.py) performs weighted preranked GSEA with **GSEApy**, weight 1, **1,000 gene-set permutations**, seed **20261019**, two threads, and eligible gene-set sizes of **15–500 genes after intersection with the complete ranked universe**. Empirical FDR q values are calculated separately for each contrast and collection. These are gene-set permutations of the ranked statistic, not sample-label permutations. Finite-permutation nominal *P* values or empirical q values can equal zero; raw software values are retained. The two-decimal on-figure q values should be read together with the complete numerical tables.

The 7,647 retained SIM identifiers receive a seeded one-to-one random assignment to unique mouse gene symbols drawn from the union of the two annotation collections, independently of expression, effect size, test statistic, and simulation class. The assignment is fixed across contrasts. It supplies a demonstration identifier mapping, not an identification of measured mouse genes. Neither seeds nor memberships are selected to obtain significance, and no tissue mechanisms are inferred from the output. The [mapping](synthetic_demo/gsea_data/synthetic_to_mouse_symbol_map.tsv), [complete rankings and results](synthetic_demo/gsea_data), and [run provenance](synthetic_demo/gsea_data/gsea_provenance.json) are supplied, including leading-edge members, full result tables, checksums, software versions, and tested-set coverage.

- **KEGG:** 313 canonical native-mouse pathways from an archived NCBI BioSystems collection mirrored in the Stephens Lab `pathways` repository; 89 `mmu_M` modules are excluded. This is an archived collection associated with a 2020 retrieval, not a current KEGG release. The exact constituent KEGG release was not recorded by the mirror. See the [KEGG provenance](synthetic_demo/gsea_data/resources/kegg_mouse_provenance.json) and [GMT](synthetic_demo/gsea_data/resources/kegg_mouse_biosystems_archive_2020.gmt), including source commits and checksums. There is no human-to-mouse orthology conversion.
- **GO Biological Process:** 8,474 native-mouse sets derived from MGI annotations and the GO release recorded as **19 May 2026**, mirrored in `enrichlite`. Upstream processing excludes IEA, ND, and NOT-qualified records; propagates annotations through `is_a` and `part_of`; and retains sets containing 5–500 GAF objects. Local processing projects memberships to authoritative mouse protein-coding symbols. Terms excluded upstream for exceeding 500 objects remain absent even if an intersection with the ranked universe might have reduced their size. See the [GO BP provenance](synthetic_demo/gsea_data/resources/go_bp_mouse_provenance.json) and [GMT](synthetic_demo/gsea_data/resources/go_bp_mouse_2026_05_19_nonIEA_protein_coding.gmt). Gene Ontology and Mouse Genome Informatics supply the underlying annotations; the GO data are attributed under CC BY 4.0.

For experimental analysis, replace the arbitrary SIM-to-symbol mapping with validated measured-gene identifiers and review appropriate, versioned annotation resources. The archived KEGG version and restricted GO collection are documented limitations of this demonstration. Computational background: [Subramanian et al., GSEA, PNAS (2005)](https://doi.org/10.1073/pnas.0506580102); [Fang et al., GSEApy, Bioinformatics (2023)](https://doi.org/10.1093/bioinformatics/btac757).

## Interpretation of the post-treatment comparison

Figure 3 shows differences **present after leptin treatment**. These may include persistent LD-associated differences and treatment-induced differences. The comparison alone does not distinguish those origins or establish that a difference was uncorrected. Cross-contrast effect estimates and uncertainty would be needed for that assessment with experimental data. Proximity to zero or nonsignificance after treatment does not establish equivalence, normalization, or rescue.

## Comparison-specific color palettes

The warm/cool pairs follow the explicit RGB values printed in the supplied color-scheme guide:

| Figure | Higher expression / positive NES | Lower expression / negative NES |
|---|---|---|
| 1 | Yellow `#FDB338` | Blue `#025196` |
| 2 | Tan `#FFBE6A` | Turquoise `#40B0A6` |
| 3 | Orange `#EB6123` | Purple `#512888` |

Heatmaps run from the cool color through white to the warm color. Their annotation strips label the first group with the warm color and the second with the cool color. Thus a group's color can differ between figures; full labels supply its identity. Direction on the NES/log2-fold-change axes and GSEA marker fill provide additional encodings. See the [color-accessibility review and previews](accessibility/COLOR_ACCESSIBILITY.md) for modeled checks and limitations.

## Artwork, font evidence, and reproducibility

Each 1536 × 1024 pixel figure uses Nimbus Sans throughout. PDF fonts are embedded and searchable; main SVGs use outlined lettering; editable SVGs contain embedded Nimbus Sans font data. High-resolution PNGs are 4608 × 3072 pixels. The export names are `figure_01_gsea_synthetic`, `figure_02_gsea_synthetic`, and `figure_03_gsea_synthetic`, with the appropriate file extensions. Synthetic provenance is stored in the export metadata; questions and interpretive prose remain outside the artwork. PCA is omitted.

[Font verification](FONT_VERIFICATION.md) compares actual exported font programs, character metrics, and glyph outlines against genuine Nimbus Sans source fonts; the [verification report](synthetic_demo/export_font_verification.json) is tied to file hashes and must be regenerated after export. [Plot audits](synthetic_demo/font_and_plot_audit.json) record labels, selections, and panel definitions. A configured font-family string alone is not the verification evidence.

From the repository root:

```bash
python -m pip install -r figures/requirements.txt
python figures/synthetic_demo/simulate_and_analyze.py
python figures/synthetic_demo/run_gsea.py \
  --kegg-gmt figures/synthetic_demo/gsea_data/resources/kegg_mouse_biosystems_archive_2020.gmt \
  --go-bp-gmt figures/synthetic_demo/gsea_data/resources/go_bp_mouse_2026_05_19_nonIEA_protein_coding.gmt \
  --kegg-provenance figures/synthetic_demo/gsea_data/resources/kegg_mouse_provenance.json \
  --go-bp-provenance figures/synthetic_demo/gsea_data/resources/go_bp_mouse_provenance.json
python figures/synthetic_demo/render_figure.py
python figures/synthetic_demo/verify_export_fonts.py
```

Install the URW Nimbus Sans fonts (Debian/Ubuntu package `fonts-urw-base35`). On another platform, set `NIMBUS_SANS_FONT_DIR` to the folder containing `NimbusSans-Regular.otf`, `NimbusSans-Bold.otf`, and `NimbusSans-Italic.otf`. Missing fonts cause an error; they are never replaced silently. Font conversion for robust PDF embedding is cached outside the repository. The data script is seed-reproducible, including deterministic gzip timestamps.
