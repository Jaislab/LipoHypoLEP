# Male data dated 30 September 2026 — figure legends

These preliminary figures use the supplied, user-confirmed filtered male differential-expression lists. They display supplied estimates and probabilities, with no reconstructed expression values or new gene-level statistical tests. The complete tested-gene tables and individual sample expression were not supplied. Sample numbers and mouse-to-library assignments are not inferred from the project animal register.

## Figure 1

**How does untreated LD differ from non-LD controls?** LD SALINE versus NON-LD control, using `DEGmale_LDsaline_WT.csv` with its supplied direction. Positive effects denote higher expression in LD SALINE. **a**, All 456 supplied rows are plotted by supplied log2 fold change and −log10(adjusted P). Yellow and blue encode 84 higher and 257 lower genes meeting adjusted P < 0.05 and absolute log2 fold change ≥ 0.2. **b**, The ten eligible genes with the smallest adjusted P values, selected and ordered as described below. The large supplied Mir6236 and Lars2 effects remain visible in panel a.

## Figure 2

**What does leptin treatment change in LD mice?** LD LEPTIN versus LD SALINE, using `DEGmale_LDsaline_LDleptin.csv` with its supplied effects negated. Original estimates are retained separately; P values and adjusted P values are unchanged. **a**, All 11 supplied rows. Positive effects denote higher expression after leptin. Seven supplied genes meet both cutoffs, all lower after treatment and colored turquoise; tan would encode qualifying higher effects. **b**, All seven eligible genes: Lcn2, Tgm2, Apoe, Sgk3, Nbl1, Kcnj13 and Cyp26b1, ordered by displayed effect.

## Figure 3

**What differences remain after leptin treatment?** LD LEPTIN versus NON-LD control, using `DEGmale_LDleptin_WT.csv` with its supplied direction. Positive effects denote higher expression in LD LEPTIN. **a**, All 297 supplied rows; orange and purple encode 43 higher and 177 lower genes meeting both cutoffs. **b**, The ten eligible genes with the smallest adjusted P values. These differences can include persistent LD-associated and treatment-induced differences; this contrast alone cannot distinguish them.

## Shared methods

Dashed volcano guides mark log2 fold changes −0.2 and +0.2 and adjusted P = 0.05; the zero-effect guide is lighter. Colored points meet both cutoffs; gray points are supplied rows that fail the joint rule. All supplied rows already satisfy absolute log2 fold change > 0.2 and adjusted P < 0.1. Consequently, these are filtered-export volcano plots: the unexported central/nonsignificant population is unavailable. No unobserved genes or gray background cloud is invented.

Panel b selects up to ten genes by ascending supplied adjusted P, then Ensembl ID, among joint-cutoff hits. Selected genes are ordered by displayed fold change; point estimates have no added uncertainty intervals. Its dashed guides mark zero and ±0.2. All figures use explicit per-contrast symmetric effect axes; axis spans can differ between panels and figures. Estimates are neither compressed nor clipped. Panel a retains the full extent of every supplied estimate.

The joint rule is a display/selection criterion following the supplied upstream tests, not a test of a true effect exceeding 0.2. Original P values and adjusted P values are preserved, without recalculating multiple testing over these selected genes. An absent gene has unknown status and is not assigned an effect of zero or classified as rescued. The overlap of qualifying genes in the two LD-versus-control exports is descriptive, not an equivalence or treatment-rescue test.

KEGG and GO BP GSEA panels await complete tested-gene rankings, including nonsignificant genes and a correctly oriented signed test statistic. The prior workbook's pathways are not reused as current findings. Individual expression values are also needed for sample heatmaps; their absence does not justify expression reconstruction. These figures therefore contain two supported gene-level panels, without empty enrichment panels or methodological notes inside the artwork.

All visible text is 32 px Nimbus Sans on a 1536 × 1024 white canvas. PDF and SVG exports preserve typography; high-resolution PNG is 4608 × 3072. See the [analysis methods and findings](../../docs/MALE_2026_09_30_ANALYSIS.md), [input provenance](../../analysis/male_2026_09_30_results/analysis_manifest.json) and [export verification](figure_validation.json).
