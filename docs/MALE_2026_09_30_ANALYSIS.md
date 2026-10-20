# Male differential-expression exports — 30 September 2026

This analysis uses the three latest CSVs supplied by the user. The user confirmed that these are **filtered lists** and that each filename denotes **first group minus second group**. The filenames identify male data. No individual expression values, original signed test statistics, model descriptions or mouse-to-library identifiers accompany these files.

## Scope and contrast mapping

| Figure | Source CSV | Displayed contrast | Source effect multiplier |
|---|---|---|---:|
| 1 | `DEGmale_LDsaline_WT.csv` | LD SALINE − NON-LD control | +1 |
| 2 | `DEGmale_LDsaline_LDleptin.csv` | LD LEPTIN − LD SALINE | −1 |
| 3 | `DEGmale_LDleptin_WT.csv` | LD LEPTIN − NON-LD control | +1 |

`WT` refers to the non-lipodystrophic LDLR−/− control group in the project context. The original [CSV files](../data/male_2026_09_30/source) are retained byte-for-byte. Curated tables retain all supplied rows, Ensembl IDs, symbols, raw P values, adjusted P values and original fold changes. Only the displayed Figure 2 effect is negated. There is no refit, new multiple-testing correction, effect compression or expression reconstruction.

All 764 source rows have valid, finite fields; there are no missing values or duplicate Ensembl IDs within a file. Every supplied row has absolute log2 fold change > 0.2 and adjusted P < 0.1. These observed bounds do not establish the full upstream export rule. Supplied adjusted P values are preserved, rather than recomputed over the selected rows.

## Gene-level results

The display criteria remain **adjusted P < 0.05 and absolute log2 fold change ≥ 0.2**.

| Figure | Supplied rows | Joint-cutoff hits among supplied rows | Higher in first group | Lower in first group |
|---|---:|---:|---:|---:|
| 1 | 456 | 341 | 84 | 257 |
| 2 | 11 | 7 | 0 | 7 |
| 3 | 297 | 220 | 43 | 177 |

These are counts within the supplied filtered exports. The number of all tested genes, completeness of qualifying-gene coverage and sample counts cannot be recovered from these files.

### What does leptin treatment change?

All seven supplied joint-cutoff treatment hits are lower in LD LEPTIN than LD SALINE:

| Gene | LD LEPTIN − LD SALINE log2 fold change | Supplied adjusted P |
|---|---:|---:|
| Lcn2 | −2.4284 | 0.006421 |
| Tgm2 | −0.8519 | 0.006421 |
| Apoe | −0.3476 | 0.018527 |
| Sgk3 | −0.7829 | 0.018527 |
| Nbl1 | −0.6422 | 0.018527 |
| Kcnj13 | −2.4201 | 0.018527 |
| Cyp26b1 | −1.4553 | 0.030317 |

This list does not establish absence of other treatment effects outside the export.

### What differences remain after leptin?

**134 genes** pass both cutoffs in both LD-versus-control exports. All have the same direction in those two comparisons: 11 higher and 123 lower. The [complete intersection table](../analysis/male_2026_09_30_results/baseline_posttreatment_shared_hits.tsv) retains both supplied effects and adjusted P values.

Examples lower in both comparisons include **Bcas1, Trf, Plp1, Ugt8a, Ermn and Mag**. Tcaf2 is higher in both, with supplied log2 fold changes approximately 3.282 in untreated LD and 3.186 after leptin. These are gene-level observations. They do not establish a pathway result, functional mechanism, statistical equivalence or complete rescue. An absent gene has unknown status in a filtered export; it is not assigned a zero effect or declared nonsignificant. There are no genes exported in all three files, so a three-contrast additivity check is unavailable.

The large untreated-LD estimates for **Mir6236** (−12.820) and **Lars2** (−6.805) are retained. Figure 1's volcano displays −4 to +4 log2 fold change, with labelled leftward arrows for these two genes at their actual −log10(adjusted P) heights. The other 454 supplied genes appear as points; all 456 remain in the analysis and counts. Sample expression, abundance, uncertainty estimates and upstream model details are needed to assess whether the extreme estimates reflect low-abundance instability or reproducible biological effects. No cause is inferred from the filtered results alone.

## Enrichment and expression panels

**KEGG/GO BP GSEA is pending complete tested-gene tables.** Filtering before ranking removes the necessary ranking/background information. Full results with a correctly oriented signed test statistic are needed, including nonsignificant genes. Individual expression values are a separate missing input: they would enable sample-expression displays but do not replace the complete DE rankings.

No GSEA or overrepresentation analysis is calculated from these selected lists. Neither the older workbook's background nor the union of the three current lists is substituted for the missing tested-gene universe. Older enrichment tables remain historical artifacts and are not results for these male CSVs.

Sample heatmaps, replicate-level uncertainty and PCA are not generated. The project register lists male mice, but this export does not identify which animals or libraries entered each contrast; no replicate count is assigned from that register.

## Figures and verification

Each figure has **a**, a volcano plot of the supplied rows, and **b**, up to ten joint-cutoff genes with the smallest supplied adjusted P values. Ensembl ID breaks adjusted-P ties; the displayed genes are ordered by effect. Figure 2 shows all seven eligible genes. External [figure legends](../figures/male_2026_09_30/FIGURE_LEGENDS.md) explain the filtered scope, selections and axes.

All text uses **32 px Nimbus Sans** on a 1536 × 1024 canvas, including off-scale gene labels. Figure 1 is yellow/blue, Figure 2 tan/turquoise and Figure 3 orange/purple. Explicit per-contrast axes and Figure 1's off-scale arrows represent every supplied gene. No estimates are rescaled. PDF embeds Nimbus; primary SVG outlines text; editable SVG embeds the fonts. Every exported PDF text span is checked, and actual embedded font metrics and glyph outlines are compared with installed URW Nimbus source fonts. Both PNG resolutions are compared with their companion PDFs. The [verification report](../figures/male_2026_09_30/figure_validation.json) binds inputs, code, source fonts and all 15 exports by SHA-256, and records the two arrows' true effects, P-value ordinates and annotation positions.

The [analysis manifest](../analysis/male_2026_09_30_results/analysis_manifest.json) records source hashes, directions, counts, row-preservation rules and the enrichment hold. Source CSV cells are treated as data, not analysis instructions. The older workbook analysis and original synthetic template are preserved separately.
