# Legacy workbook: analysis scope and methods

This branch analyses `DGElist_AnkeData(1).xlsx` as historical experimental differential-expression results. It is separate from the synthetic figure template on `main`. The original workbook remains local; curated tables preserve the supplied numeric fields and record its SHA-256. Spreadsheet content is data, not analysis instructions.

The [initial results and weighting sensitivity](LEGACY_RESULTS.md) summarize the completed analyses and independent verification.

## Contrast mapping

| Figure | Workbook tab | Displayed contrast | Effect orientation |
|---|---|---|---|
| 1 | `saline_VS_wt` | LD SALINE vs NON-LD control | Supplied log2 fold change |
| 2 | `saline_VS_leptin` | LD LEPTIN vs LD SALINE | Negative of supplied log2 fold change |
| 3 | `leptin_VS_wt` | LD LEPTIN vs NON-LD control | Supplied log2 fold change |

`wt` denotes the non-lipodystrophic LDLR−/− controls described by the user, not a claim of a genetically wild-type LDLR background. Source effects remain in a separate column. Reversing the second contrast changes its direction, not its supplied P values, adjusted P values or SEs. Upstream statistical tests and effect-estimation methods have not been independently established.

The user confirmed that all tested genes were exported. We therefore use each tab's full tested universe: 15,518, 15,136 and 19,721 genes in figure order. The smaller universes are exact subsets of the largest at observed baseMean thresholds of approximately 10.39 and 12.84; the largest starts at approximately 0.68. This pattern documents differing analysis filters, without establishing which software procedure produced them. We do not force a shared universe or invent missing contrast results.

## What is preserved and what can be checked

The eight source columns are Ensembl ID, symbol, Entrez ID, baseMean, log2 fold change, SE, raw P and adjusted P. All 50,375 rows have finite, valid numeric fields and unique Ensembl IDs within their tab. The supplied adjusted P values agree with independent Benjamini–Hochberg calculations within approximately 4 × 10⁻¹⁵, separately within each complete tab. This is a verification, not a replacement of the supplied values.

Gene highlighting and selection use **adjusted P < 0.05 and absolute supplied log2 fold change ≥ 0.2**. In figure order:

| Figure | Joint-cutoff genes | Higher in first group | Lower in first group |
|---|---:|---:|---:|
| 1 | 574 | 142 | 432 |
| 2 | 750 | 512 | 238 |
| 3 | 1 | 1 | 0 |

These are estimate-based display criteria after zero-null testing; they are not a test that true effects exceed 0.2. No effect rescaling, truncation, refitting or new gene-level multiple-testing adjustment is performed.

The workbook has no sample-level expression, counts, replicate metadata or original signed test statistic. `baseMean` is a single per-gene summary and cannot supply individual mice or a sample-expression heatmap. The current legacy figures replace that panel with point estimates of selected genes. There are no invented replicates, expression values or uncertainty intervals.

## Effect estimates require upstream clarification

The supplied `log2FoldChange/lfcSE` ratios do not reproduce the source P values, so they are not relabelled as the original Wald statistic. The discrepancy is substantial for some genes; for example, the supplied Snx31 effect/SE implies a two-sided normal P near 0.606 while the source P is about 0.0000469. This is a compatibility diagnostic, not a replacement test.

The leptin-versus-control effects are unusually concentrated near zero: median absolute log2 fold change is approximately **3.49 × 10⁻⁶**, with **14,946 of 19,721** estimates below **10⁻⁵** in absolute value. The supplied contrasts also do not satisfy exact effect additivity on their shared genes. Contrast-specific shrinkage or other processing could explain these observations, but the workbook cannot identify the cause. Preserve the values and request upstream model, shrinkage and contrast definitions before drawing conclusions from their relative magnitudes.

The one gene passing the joint criteria after treatment is **Tcaf2**, with supplied log2 fold change approximately **3.459** and adjusted P approximately **3.17 × 10⁻¹³**. Few significant differences, compressed estimates or nonsignificance after treatment do not establish complete normalization, rescue or equivalence.

## Exploratory pathway enrichment

Because the original signed statistic is absent, the primary preranked GSEA score is explicitly derived as **sign(oriented supplied log2 fold change) × −log10(source raw P)**. It is a ranking proxy, not a recovered Wald statistic. A sensitivity run uses **sign(log2 fold change) × normal inverse-survival(P/2)**. These scores order genes identically here but give different weights to the enrichment walk; agreement is useful evidence of weighting stability, not a replacement for the missing upstream statistic.

The proxy uses effect direction and P value, not effect magnitude. Consequently, enrichment can be present even when very few supplied effects meet the 0.2 display cutoff, particularly in Figure 3. Such enrichment does not establish large residual effects, lack of correction or treatment rescue.

Every Ensembl row in each complete tab remains in its ranked background, including nonsignificant genes and genes outside a collection. Neither DEG cutoff filters GSEA. Weighted preranked GSEA uses weight 1, 1,000 gene-set permutations, seed 20261020, two threads, and set sizes 15–500 after intersection with the ranked universe. Empirical GSEA FDR q values are calculated within each contrast, collection and ranking metric; they are distinct from gene-level BH adjusted P values. Exact zero outputs are finite-permutation estimates, not proof of an exactly zero false-discovery probability; no unsupported lower-bound value is substituted. Gene-set permutations do not preserve within-pathway gene correlation or propagate biological-replicate/model uncertainty.

The existing pinned native-mouse resources provide archived KEGG BioSystems pathways associated with 2020 and GO BP annotations recorded as 19 May 2026, with the documented upstream size/non-IEA/protein-coding restrictions. These are real gene memberships. The synthetic demo's random gene assignment is never used.

Memberships are translated through the workbook's exact, case-preserved symbols to measured Ensembl IDs only where that mapping is globally unambiguous. Sixteen symbols map to more than one Ensembl ID. Their ambiguous memberships are excluded consistently; all their Ensembl rows remain in the ranked background and in DE tables. No gene is chosen because it has the strongest effect or smallest P value.

The 13 legacy mitochondrial labels match the archived KEGG collection but not the modern GO symbol names. This leaves a specific GO BP mitochondrial coverage gap. Blind alias conversion is unsafe: the cached alias table would map workbook `COX1`/`COX2` to `Ptgs1`/`Ptgs2`, which are different genes. Authoritative NCBI/Ensembl endpoints were unavailable in this environment. No unsupported alias substitution is made. This limitation matters for interpreting mitochondrial or respiratory pathways and should be resolved with a validated ID-based crosswalk before final analysis.

## Figure and table interpretation

Panels a and b preserve the source estimates and gene-level probabilities. Panel b shows up to three higher and three lower genes with the smallest adjusted P values among joint-cutoff hits, with stable Ensembl tie-breaking. Figure 3 correctly shows its single eligible gene. Point estimates have no added confidence intervals because the uncertainty-estimation method is unknown.

Panels c and d show up to the two largest positive and two most negative NES values, regardless of significance. If a direction has fewer than two sets, remaining positions use the strongest absolute NES among unselected sets. Filled markers encode GSEA FDR q < 0.05; open markers encode q ≥ 0.05. Full result tables and the weighting sensitivity must be consulted before interpretation. Pathways are selected by a stated numeric rule, not by a desired biological story.

All figure text uses **32 px Nimbus Sans** on the 1536 × 1024 canvas. Figure 1 uses yellow/blue; Figure 2 tan/turquoise; Figure 3 orange/purple. Positive effects and NES refer to the first named group. Methods and cautions remain in external legends rather than inside the artwork.

## Reproduce

Install `analysis/requirements.txt`, then run the workbook importer, legacy GSEA runner and legacy figure renderer from the repository root. The original `.xlsx` is a local input; its bytes are identified by the import provenance. The workflow writes only its separate legacy-results and legacy-figure directories. See the branch README for exact commands and result links.
