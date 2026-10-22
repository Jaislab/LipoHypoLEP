# LipoHypoLEP — anatomical QC of mediobasal hypothalamus samples

We are studying the effects of metreleptin on the brain using bulk RNA sequencing. To interpret those effects, we need to assess whether our dissections contain the intended mediobasal hypothalamic (MBH) compartments, including the median eminence (ME).

## Why we are revising the marker panel

Our previous approach may have introduced bias by stratifying samples using hypothalamic genes that respond to leptin. If these genes were also used for normalization or batch correction, part of the biological treatment response could have been treated as technical variation. Whether that happened still needs to be established from the upstream processing records.

An anatomical marker that changes with treatment can make an intact sample appear anatomically different. We therefore want a panel that covers ME vasculature, tanycytes and MBH neuronal populations, with particular attention to neuronal populations without detectable **Lepr** expression. The goal is to assess tissue coverage while preserving the metreleptin response we are trying to measure.

**Agrp** and **Pomc** remain useful context for ARC neuronal populations, but they are leptin-sensitive and their expression alone does not establish that the ME was retained.

## New candidate panel

The following expression summaries come from our **local analysis of HypoMap**. They describe enrichment in the reported atlas clusters or regions; they are not expression thresholds for bulk samples.

| Compartment | Gene | On-target CPM | Strongest reported off-target | Specificity from displayed CPM |
| --- | --- | --- | --- | --- |
| ME fenestrated endothelium | **Plvap** | 2678 (Plvap.Endothelial) | 63 (neural stem cells) | 42.5× |
| ME floor (β2 tanycytes) | **Scn7a** | 2907 (Scn7a tanycytes) | 162 (ependymal) | 17.9× |
| Tanycytes, all subtypes | **Col23a1** | 667–1283 across all five tanycyte clusters | 151 (meningeal fibroblasts) | 4.4–8.5× |
| Tanycytes, all subtypes | **Rax** | 150–456 across all five tanycyte clusters | 54 (one astrocyte cluster) | 2.8–8.4× |
| ARC neurons | **Crabp1** | 166 (ARC), 138 (posterior periventricular) | ≤11 in the other reported regions | ≥15.1× (ARC), ≥12.5× (posterior periventricular) |
| VMH + DMH neurons | **Gpr149** | 97 (VMH), 87 (DMH) | 32 (zona incerta) | 3.0× (VMH), 2.7× (DMH) |
| VMH neurons | **Fezf1** | 109 (VMH) | 30 (anterior hypothalamus/PVH) | 3.6× |

Ratios are calculated from the rounded CPM values shown. For Col23a1 and Rax, the ranges retain the variation across tanycyte clusters rather than reporting only the highest ratio. Specificity is relative to the reported comparison set, not established whole-brain exclusivity.

## What the panel tells us

**Plvap** provides a candidate readout of fenestrated endothelium associated with the ME. **Scn7a**, together with **Col23a1** and **Rax**, adds evidence for the tanycyte component. **Crabp1**, **Gpr149** and **Fezf1** extend coverage to ARC, VMH and DMH neuronal populations.

Together, these markers provide complementary evidence about dissection composition. Their combined presence would support representation of these compartments, but bulk expression cannot prove that the complete MBH was captured. Plvap is not intrinsically exclusive to ME vessels, and the neuronal candidates show expression outside their intended target regions. Anatomical interpretation should incorporate dissection records and spatial or histological evidence where available.

## How we will use it

Use this panel for descriptive anatomical QC, separate from normalization and batch correction. Review the compartments individually rather than allowing a strong signal in one compartment to compensate for missing evidence in another. Do not automatically exclude samples or adjust treatment effects using an unvalidated marker score.

The next validation step is to examine **Lepr coexpression in marker-positive cells**, spatial localization, and stability under leptin/metreleptin exposure. The supplied HypoMap summaries do not yet establish that these populations are Lepr-negative or that their markers are treatment-independent. Undetected Lepr can reflect limited single-cell sensitivity, and receptor-negative cells can still respond indirectly to treatment.

The panel is therefore our current anatomical-QC candidate set. Confirming its suitability requires retaining the HypoMap release, annotations and local analysis methods, checking the upstream use of the previous markers, and validating the new panel against the actual samples.
