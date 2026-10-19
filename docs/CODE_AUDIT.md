# Code audit before experimental data

Audited on 19 October 2026, starting from commit `052eca52cb9abc735c8e4f9c98cc5dd729244476`. Scope: every Python source file, numerical synthetic tables, six GSEA analyses and their annotation provenance, three figure export sets, font/accessibility evidence, dependencies, delivery bundle and upcoming CSV handoff. Confirmed issues were fixed and regression tests added; this is a source and numerical review, not a biological validation or an external security certification.

**The synthetic numerical calculations passed independent reconstruction. No effects, variability, fold changes or thresholds were rescaled.** The demo remains separate from experimental analysis. The incoming files will contain differential-expression results and per-sample FPKMs, so the supplied estimates will be preserved rather than refitted with the demo's model.

## Confirmed issues and resolutions

| Priority | Finding | Resolution |
|---|---|---|
| High | A failed GSEA rerun could retain an old `complete=true` marker alongside partially updated files. | Invalidate completion before processing; lock out concurrent writers; stage all outputs; publish completion last; verify all 30 output hashes and six expected contrast/collection combinations. |
| High | Source annotations could change between parsing and hashing; source provenance was optional despite a mouse-species claim. | Parse snapshotted source bytes, require declared mouse taxonomy and matching provenance checksums, and recheck source/input/code bytes before completion. This verifies declared provenance and byte integrity, not the truth of arbitrary metadata. |
| High | Simulation output could overwrite an unrelated nonempty folder; numerical or serialization failure could leave partial output. | Require an empty destination or an existing synthetic marker; reject symlink destinations; calculate/validate in a sibling stage; publish with handled-error rollback. |
| Medium | Invalid BH probabilities, unidentifiable median-ratio normalization and undefined OLS residual variance lacked explicit failure guards. | Validate ranges, dimensions, finite values, sample groups and identifiable variance before use. |
| Medium | Rendering required exactly 40 heatmap genes and 12 supporting genes and assumed four available pathways. | Select up to the stated maximum; handle empty or sparse results without inventing hits. |
| Medium | Zero adjusted P values could produce infinite volcano coordinates; constant heatmap rows could produce NaNs. | Use a documented positive floor for plotting zero probabilities only; preserve source probabilities; display constant expression rows at neutral zero z-score. |
| Medium | Heatmap group boundaries assumed six samples; fixed NES axes could clip larger scores. | Derive group widths and boundaries from metadata and expand NES axes when necessary. |
| Medium | Two-decimal q labels could turn a significant 0.049 into a displayed 0.05. | Keep threshold-aware labels and scientific notation for small positive q values; retain exact values in tables. |
| Medium | Default CSV float parsing could perturb the last binary digits of signed statistics and supporting estimates. | Use round-trip binary64 parsing for GSEA/rendering and test exact conservation into all six ranking calculations and saved ranks. Input DE tables remain unchanged; regenerated derived tables can differ in their final numerical digits. |
| Medium | Font checks used `assert`, which disappears under `python -O`; failed font/color reruns could leave old PASS reports. | Use explicit validation exceptions; atomic RUNNING/FAILED/PASS reports; hash-based freshness checks; inspect effective SVG fonts and both PNG resolutions against PDFs. |
| Medium | The prior cross-figure validation report had no saved executable counterpart; bundle creation was manual. | Add executable cross-file validation and a verified ZIP builder with CRC/SHA-256 checks, per-member manifest and a 40 MB limit. |
| Medium | There was no experimental CSV preflight or explicit handoff schema. | Add a read-only CSV validator and [input contract](INPUT_DATA_CONTRACT.md), covering supplied DE results, inline sample FPKMs, metadata, missingness, identifiers, contrast direction and full-universe declarations. |

## Independent numerical and export evidence

- Reconstructed size factors and log-expression from raw synthetic counts. Recomputed all contrasts using an independent design-matrix OLS fit, contrast covariance, Student t tests and SciPy's BH implementation; estimates, standard errors, statistics and probabilities agree. The current generator reproduces all ten quantitative output tables byte-for-byte.
- Independently recomputed weighted GSEA running-sum enrichment scores and leading-edge membership for **all 9,219 pathway results**, rather than relying on GSEApy's own checks. The initial maximum ES discrepancy was **1.09 × 10⁻¹³**, with no leading-edge membership discrepancies. Tests repeat this check on the current results. This does not independently reproduce empirical NES null distributions or q values; those remain GSEApy outputs.
- A clean Python 3.12.14 virtual environment installed the declared dependencies successfully; `pip check` found no broken requirements. The [tested dependency lock](../figures/requirements-lock.txt) records direct and transitive versions, including optional color checks. Nimbus source package: Debian `fonts-urw-base35` version `20200910-8`.
- Re-ran the full six-analysis, 1,000-permutation GSEA workflow, rendered all figures, independently verified embedded Nimbus glyph programs, compared both PNG resolutions with PDFs, and regenerated color-vision simulations. Current results are recorded in the machine-readable [cross-figure validation](../figures/synthetic_demo/question_figure_validation.json), [font verification](../figures/synthetic_demo/export_font_verification.json), and [color audit](../figures/accessibility/palette_simulation.json).
- Unit/regression checks cover independent numerical references, exact cutoff boundaries, failed reruns, changed inputs, corrupted output evidence, duplicate identifiers, invalid probabilities, absent discoveries, zero P values, constant expression, unequal sample counts, extreme NES, and genuinely exported edge-case figures.

**Final regression suite: 85 tests passed** in the clean environment. Categories: simulation/numerical integrity (16), GSEA runner/integrity (14), independent all-pathway ES/leading edges (3), renderer/precision/layout/fonts (19), export/accessibility evidence (14), and experimental CSV preflight (19). Python compilation and `git diff --check` also passed. Long labels that would collide or leave their panel now fail with an actionable layout error; missing Nimbus characters raise an error instead of silently disappearing.

The round-trip parser fix leaves the original synthetic counts, normalized expression, metadata, truth and three DE tables byte-identical to their pre-rescaling baseline (`2eb3599f9d4b1ffe2036dc2fa4d517edd0f647bc`). Derived ranking/supporting tables now preserve their source binary64 values exactly. Relative to the previous derived enrichment tables, ES differences are at floating-point scale (≤4.45 × 10⁻¹⁶); the largest NES/q differences are below 0.000876 in the final GO BP contrast, reflecting finite-permutation normalization. No displayed pathway, gene selection or significance classification changed.

## Experimental-data boundary

Use the supplied **log2** fold changes, adjusted P values and validated **signed** statistics. Apply `padj < 0.05` and `abs(log2FoldChange) >= 0.2` for gene highlighting/heatmap selection. This combines a zero-null significance test with an observed effect-size filter; it does not test that the true effect exceeds 0.2. Do not recompute BH silently or convert a linear fold-change column without an explicit recorded decision.

Use **all eligible tested genes**, including nonsignificant genes, for preranked GSEA. A significant-only DEG list is insufficient. Missing adjusted P values may reflect independent filtering and do not automatically invalidate an otherwise valid ranking statistic; missing or outlier-invalid tests need an explicit eligibility decision. An unsigned F/LR statistic must not be used directly as a signed ranking metric. Each experimental contrast can have a legitimately different eligible universe.

The FPKM sample columns can supply the heatmap's expression matrix. They must match sample metadata exactly. Any log transformation/pseudocount used for heatmap display will be chosen and documented after inspecting the values; it must never alter or replace the supplied DE estimates or ranking statistics. FPKM is not raw integer count data and will not be sent to a count model.

The experimental importer, identifier mapping and production rendering adapter cannot be verified against unseen CSVs. The new preflight checks structure and declared provenance; it cannot prove correct upstream modeling, mouse identity, biological independence, complete-universe claims or the true contrast direction. Once files arrive, inspect those declarations and column mappings before making figures. Do not route experimental CSVs through `simulate_and_analyze.py` or the synthetic random-symbol assignment in `run_gsea.py`.

For genuine mouse GSEA, use validated measured-gene identifiers and versioned mouse annotations; record unmapped/ambiguous identifiers and mapping coverage. Do not upper-case mouse symbols, silently collapse duplicates or randomly assign symbols. Review/refresh annotation resources before interpretation: the demo KEGG is an archived BioSystems collection associated with 2020, and the GO BP collection has documented upstream size and protein-coding restrictions. Preranked gene-set permutations do not preserve within-pathway gene correlation or propagate sample/model uncertainty; their q values describe that conditional enrichment procedure, rather than sample-label testing.

Figure 3 describes differences present after treatment. It cannot alone distinguish persistent LD effects from treatment-induced differences. Nonsignificance does not demonstrate normalization, equivalence or rescue. The three contrasts are related and must not be described as independent experiments; additivity of supplied shrunken effects is not a universal hard constraint.

## Repeat the audit

Use the reproduction commands in the [README](../README.md), then:

```bash
python -m unittest discover -s tests -v
python figures/synthetic_demo/verify_export_fonts.py
python figures/accessibility/audit_colors.py
python figures/synthetic_demo/validate_figures.py
python figures/synthetic_demo/build_bundle.py
```

After any modification, `validate_figures.py --check-report` rejects stale evidence. All reports concern the files identified by their hashes; a previous PASS is not evidence for later files. Simulation publication provides rollback for handled errors, but not crash-atomic directory exchange or concurrent simulation writers. GSEA's exclusive lock uses POSIX `fcntl`, appropriate for this Linux cloud environment. Render/export stages fail closed through their completion reports; they do not promise atomic replacement of an entire multi-file figure set. Long real-data labels and the final experimental figures still require visual inspection.
