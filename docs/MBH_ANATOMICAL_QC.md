# Candidate panel for mediobasal hypothalamus coverage

Status: proposed anatomical QC panel, not a validated sample classifier or a normalization reference. Source: user-supplied results from a local HypoMap analysis (reported 2026-10-22). The underlying export, atlas release, species, cluster definitions, aggregation code and Lepr results have not been supplied or independently checked. The existing repository describes a mouse study; confirm that the local atlas analysis matches that species.

## Purpose and limits

Assess whether bulk samples contain evidence consistent with the intended mediobasal hypothalamic (MBH) compartments, including the median eminence (ME). Use multiple compartments and markers as descriptive evidence. Bulk expression alone cannot establish complete anatomical capture: it conflates cell abundance, per-cell expression, RNA quality and sampling depth. Histology, dissection records or spatial evidence are needed to substantiate completeness.

Assuming “Argp” referred to **Agrp**, Agrp and Pomc describe ARC neuronal populations; their transcripts do not by themselves establish retention of ME tissue. Keep these as treatment-sensitive context, outside any proposed treatment-independent QC score. Do not infer that all previously used genes are leptin-regulated without reviewing the actual list and evidence.

## Proposed candidates

All CPM values and cluster/region assignments below are supplied local-analysis summaries, not bulk thresholds. Ratios are approximate arithmetic from the displayed, rounded values. Their denominators and aggregation levels need confirmation before comparisons across rows.

| Intended compartment | Gene | Reported on-target CPM | Reported strongest off-target CPM | Ratio from displayed values | Qualification |
| --- | --- | --- | --- | --- | --- |
| ME fenestrated endothelium | Plvap | 2678, Plvap.Endothelial | 63, neural stem cells | 42.5× | Fenestrated-endothelial signal; ME localization requires spatial support and consideration of other fenestrated vascular sites. |
| ME floor / β2 tanycytes | Scn7a | 2907, Scn7a tanycytes | 162, ependymal | 17.9× | Verify β2/floor annotation and spatial distribution. |
| Tanycytes, all subtypes | Col23a1 | 667–1283, five tanycyte clusters | 151, meningeal fibroblasts | 4.4–8.5× | The weakest reported tanycyte cluster gives 4.4×, not 8×. |
| Tanycytes, all subtypes | Rax | 150–456, five tanycyte clusters | 54, one astrocyte cluster | 2.8–8.4× | The weakest reported tanycyte cluster gives 2.8×, not 8×. |
| ARC neurons | Crabp1 | 166, ARC; 138, posterior periventricular | ≤11, other reported regions | ≥15.1× for ARC; ≥12.5× for posterior periventricular | The supplied summary itself does not establish ARC exclusivity. |
| VMH + DMH neurons | Gpr149 | 97, VMH; 87, DMH | 32, zona incerta | 3.0× VMH; 2.7× DMH | Regional enrichment with off-target expression; MBH exclusivity is unproven. |
| VMH neurons | Fezf1 | 109, VMH | 30, anterior hypothalamus/PVH | 3.6× | Regional enrichment with off-target expression; MBH exclusivity is unproven. |

These are seven candidates to evaluate, not seven demonstrated treatment-stable markers. A hypothalamus-focused atlas cannot alone establish whole-brain exclusivity. Regional averages may also obscure expression in small off-target populations.

## Lepr and leptin-response validation

No Lepr coexpression or intervention results accompany the supplied panel. Do not label these genes or their expressing cells Lepr-negative yet. A gene does not itself express Lepr: the relevant question is whether marker-expressing cells express Lepr, and which receptor isoform the assay can resolve.

For each candidate, report the fraction of marker-positive cells with detectable Lepr, sample/donor-level pseudobulk Lepr abundance, cluster size, read depth and uncertainty. Define marker positivity and detection thresholds in advance. Examine joint expression within cells; separate cluster means cannot establish coexpression or its absence. Single-cell dropout and limited isoform resolution prevent equating undetected Lepr with absence of a functional receptor.

Even convincing receptor absence would not establish treatment independence: indirect neural, endocrine, inflammatory and vascular responses remain possible. Review relevant intervention evidence and assess changes across the actual treatment groups, accounting for biological replication. A nonsignificant treatment comparison does not demonstrate stability. Validate stability in independent data where possible; avoid selecting markers solely for flat expression in the same dataset used to estimate treatment effects.

## Preventing analytical bias

Keep the anatomical panel separate from size-factor estimation, reference-gene scaling and batch correction. Preserve supplied upstream DE estimates. Do not residualize treatment effects or automatically adjust DE models for marker scores that may be downstream of treatment. Likewise, do not automatically exclude samples because treatment changes a marker.

Use a prespecified, treatment-blinded anatomical review where feasible, corroborated by dissection/histology. Record compartment evidence and uncertainty rather than introducing an uncalibrated pass/fail cutoff. Missing evidence for one compartment should prompt review; high expression in another does not establish that it is present. Do not interpret atlas cell-cluster CPM as an expected bulk CPM or transfer it directly into a sample acceptance threshold.

Audit the upstream workflow before claiming bias occurred:

1. Recover the exact old marker list, inclusion/stratification rules, reference genes, scaling factors, batch-correction inputs, covariates and design matrix, together with raw counts and metadata.
2. Determine whether genes were used only for visualization, for selecting samples, for estimating normalization, or for fitting correction factors. These create different potential biases.
3. Check treatment/batch confounding and whether correction protected the biological contrast. Perfect confounding cannot be resolved by a different marker list.
4. In a separate analysis, compare defensible count-based models without marker-driven correction against the original workflow, retaining the treatment contrast and justified technical covariates. Examine effect estimates and sample inclusion, not only significant-gene counts.
5. Report sensitivity to uncertain dissections transparently; do not choose exclusions to obtain a desired treatment result.

The current checkout contains a synthetic median-of-ratios demonstration and a read-only experimental input preflight. It contains no identified experimental marker-based normalization, batch correction or anatomical stratification implementation. This inspection does not establish what happened in upstream processing.

## Evidence needed before adopting the panel

Retain the local HypoMap export and analysis code with atlas release/species, region and cell-type annotations, sample/donor counts, expression units, preprocessing and CPM aggregation method, and the full off-target search scope. Provide candidate expression and Lepr detection by cell type and region with biological-replicate support. Confirm anatomical localization using spatial evidence and evaluate treatment stability before defining any score or sample rule.

No samples are classified, excluded, normalized or batch-corrected by this document. No CPM acceptance cutoffs are defined.
