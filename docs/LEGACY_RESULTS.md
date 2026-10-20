# Historical workbook: initial results

These results use the supplied historical differential-expression tables, without refitting an expression model. All tested genes were exported, as confirmed by the user. The [methods](LEGACY_ANALYSIS.md) and [figure legends](../figures/legacy/FIGURE_LEGENDS.md) define the contrasts, selection rules and limitations.

## Gene-level results

| Figure | First group minus second group | Tested genes | Joint-cutoff hits | Higher | Lower |
|---|---|---:|---:|---:|---:|
| 1 | LD SALINE − NON-LD control | 15,518 | 574 | 142 | 432 |
| 2 | LD LEPTIN − LD SALINE | 15,136 | 750 | 512 | 238 |
| 3 | LD LEPTIN − NON-LD control | 19,721 | 1 | 1 | 0 |

The joint criteria are supplied adjusted P < 0.05 and absolute log2 fold change ≥ 0.2. Figure 2 reverses the source contrast's effect signs; original estimates are retained separately. No estimates are rescaled. Figure 3's single eligible gene is **Tcaf2** (log2 fold change 3.459; adjusted P approximately 3.17 × 10⁻¹³).

Most supplied Figure 3 effects are unusually compressed: the median absolute estimate is approximately 3.49 × 10⁻⁶. The original testing statistic is absent, and effect/SE ratios do not recreate source P values. These observations can arise when reported effects and SEs have undergone processing such as shrinkage while original testing probabilities are retained; the workbook cannot establish the procedure. Clarify the upstream model and effect estimation before interpreting relative effect magnitudes or claiming normalization after treatment.

## Pathway results and weighting sensitivity

The primary ranking is sign(log2 fold change) × −log10(raw P). The sensitivity ranking substitutes the two-sided normal-equivalent Z magnitude. Both are explicitly derived proxies using the same complete ordered genes; their weights differ. Neither is claimed to recover the missing original signed statistic.

| Figure | Collection | Tested sets | Primary q < 0.05 | Signed-Z q < 0.05 | Significant in both, same direction |
|---|---|---:|---:|---:|---:|
| 1 | KEGG | 275 | 63 | 76 | 62 |
| 1 | GO BP | 4,225 | 249 | 501 | 247 |
| 2 | KEGG | 274 | 46 | 60 | 43 |
| 2 | GO BP | 4,190 | 80 | 192 | 77 |
| 3 | KEGG | 287 | 1 | 11 | 1 |
| 3 | GO BP | 4,525 | 6 | 39 | 5 |

**Untreated LD versus control:** the primary Figure 1 panels show positive ECM-receptor interaction, complement/coagulation and collagen-related enrichment, and negative oxidative-phosphorylation and synaptic-vesicle enrichment. Each displayed Figure 1 pathway retains its direction and q < 0.05 under signed-Z weighting. These are associations in bulk tissue; they do not identify the responsible cell types or establish changes in physiological pathway activity.

**Leptin treatment in LD:** oxidative phosphorylation and synaptic-vesicle gene sets are enriched toward higher genes in LD LEPTIN versus LD SALINE, opposite their direction in untreated LD versus control. Ribosome and synaptic-translation sets are enriched toward lower genes after treatment. Each displayed Figure 2 pathway retains its direction and q < 0.05 under signed-Z weighting. Opposite enrichment directions are consistent with a treatment-associated reversal of those ranked signals, but do not demonstrate complete rescue or equivalence to controls.

**Post-treatment versus control:** neuroactive ligand-receptor interaction is the sole KEGG set with primary q < 0.05 (NES −1.974; q 0.00988); it remains significant with signed-Z weighting. Five GO BP sets are significant in both rankings with the same direction, including postsynaptic chemical transmission, excitatory postsynaptic potential, catecholamine metabolism and two cilium-related terms. Interpretation remains limited by the supplied effect estimates and weighting sensitivity. For example, NF-kappa B signaling has primary q 0.0829 versus signed-Z q 0.0159; the evidence crosses the significance threshold when weights change.

Across all six comparisons of rankings, NES Spearman correlations are 0.964–0.985, yet significance and selected terms can change substantially. Figure 3's significant-set Jaccard overlaps are only 0.091 for KEGG and 0.125 for GO BP; only one of its four displayed primary GO terms is also among the four displayed signed-Z terms. High rank correlation is therefore insufficient evidence that every pathway result is robust. The [complete sensitivity audit](../analysis/legacy_results/legacy_gsea_audit.json) records intersections, direction changes, display overlaps and the largest threshold disagreements.

KEGG disease names identify shared gene sets and are not diagnoses in these mice. GO terms overlap extensively and should not be counted as independent biological programs. Reported q = 0 values preserve finite-permutation estimates, not exactly zero false-discovery probabilities.

## Coverage and verification

Every tested Ensembl row remains in its ranking and background. Exact, globally unambiguous workbook symbols translate native mouse memberships to measured Ensembl IDs. Ambiguous memberships are excluded without discarding measured gene rows. The archived KEGG resource, modern GO projection, differing tested universes and legacy mitochondrial GO coverage gap are documented in the methods and manifests.

The independent reference audit recalculated **27,552 enrichment scores and leading-edge sets**, covering both ranking methods. Maximum absolute ES disagreement was approximately **1.25 × 10⁻¹⁵**, with **zero leading-edge or membership mismatches**. This verifies the observed weighted walks; it does not independently repeat null permutations, NES normalization or empirical FDR estimation.

The full test suite passes **128 tests** in the pinned Python 3.12 environment. All three final figures pass independent export verification: every actual PDF text span is **32 px Nimbus Sans**, embedded font metrics and glyph outlines match URW Nimbus Sans, and both PNG resolutions match their companion PDF content. See the [current export report](../figures/legacy/legacy_validation.json).

The original workbook stays local. The branch contains curated complete DE tables, complete rankings and enrichment results, leading-edge members, provenance, source code, figures and external legends. The original synthetic template remains unchanged on `main`.
