#!/usr/bin/env python3
"""Recompute cross-file consistency checks for the published synthetic figures.

This verifies the demo's data/selection/provenance contract. It does not validate
experimental biology or replace visual inspection and exported-font checks.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
REPORT = ROOT / "question_figure_validation.json"
PALETTES = (("#FDB338", "#025196"), ("#FFBE6A", "#40B0A6"), ("#EB6123", "#512888"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_report(report):
    temporary = REPORT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, REPORT)


def check_report():
    report = json.loads(REPORT.read_text())
    require(report.get("status") == "PASS", "Figure validation is not a completed PASS")
    inputs = report.get("input_sha256", {})
    require(inputs and report.get("validator_sha256") == digest(Path(__file__)), "Stale validator evidence")
    for filename, expected in inputs.items():
        path = REPO / filename
        require(path.is_file() and digest(path) == expected, f"Stale figure validation: {filename}")
    for script in (ROOT / "verify_export_fonts.py", ROOT.parent / "accessibility/audit_colors.py"):
        subprocess.run([sys.executable, str(script), "--check-report"], check=True)
    return report


def verify():
    # Importing the renderer supplies its documented selection functions, while
    # the comparisons below independently verify those outputs against tables.
    from run_gsea import validate_completed_run

    gsea = ROOT / "gsea_data"
    data = ROOT / "data"
    provenance = validate_completed_run(gsea, data)
    plot = json.loads((ROOT / "font_and_plot_audit.json").read_text())
    require(plot.get("complete") is True, "Rendering did not complete")
    input_hashes = plot.get("input_sha256", {})
    require(input_hashes, "Renderer has no hashed input evidence")
    for filename, expected in input_hashes.items():
        require(digest(REPO / filename) == expected, f"Renderer input changed: {filename}")
    metadata = pd.read_csv(data / "sample_metadata.tsv", sep="\t")
    records = []
    require(len(plot["figures"]) == 3, "Exactly three dedicated figures required")
    for number, record in enumerate(plot["figures"], 1):
        table = pd.read_csv(data / record["contrast_table"], sep="\t", float_precision="round_trip")
        require(record["figure"] == number, "Figure sequence changed")
        require(record["absolute_log2fc_cutoff"] == 0.2 and record["adjusted_p_cutoff"] == 0.05,
                "The agreed DEG display criteria changed")
        require(record["pca_panel"] is False and record["panels"] ==
                ["volcano", "two-group heatmap", "KEGG preranked GSEA", "GO BP preranked GSEA"],
                "Unexpected panel composition")
        require(record["all_visible_typography"] == "Nimbus Sans", "Unexpected font selection")
        expected_exports={f"figure_{number:02d}_large_axis_synthetic{extension}" for extension in ('.png','.pdf')}
        require(set(record.get("versioned_export_sha256",{}))==expected_exports,
                "Missing distinctly named large-axis exports")
        for filename,expected_hash in record['versioned_export_sha256'].items():
            path=ROOT.parent/filename
            canonical=ROOT.parent/f"figure_{number:02d}_gsea_synthetic{path.suffix}"
            require(path.is_file() and digest(path)==expected_hash==digest(canonical),
                    f"Versioned export differs from verified artwork: {filename}")
        colors = record["comparison_colors"]
        require((colors["higher"], colors["lower"]) == PALETTES[number-1], "Requested palette changed")
        significant = table.loc[(table.padj < 0.05) & (table.log2FoldChange.abs() >= 0.2)]
        expected_heatmap = set(significant.sort_values("padj", kind="stable").head(40).gene_id)
        require(set(record["heatmap_genes"]) == expected_heatmap, "Heatmap genes disagree with joint criteria")
        expected_samples = [sample for group in record["heatmap_group_order"]
                            for sample in metadata.loc[metadata.group.eq(group), "sample_id"]]
        require(record["heatmap_samples"] == expected_samples, "Heatmap uses the wrong samples or order")
        ranking = pd.read_csv(gsea / f"figure_{number:02d}_complete_ranked_statistics.tsv.gz",
                              sep="\t", float_precision="round_trip")
        require(set(ranking.gene_id) == set(table.gene_id) and len(ranking) == len(table),
                "GSEA does not contain the full unfiltered tested universe")
        require(ranking.stat.is_monotonic_decreasing, "GSEA ranks are not signed-statistic descending")
        actual = ranking.set_index("gene_id").stat.reindex(table.gene_id).to_numpy()
        require(np.array_equal(actual, table.stat.to_numpy()), "GSEA statistics differ from supplied estimates")
        for library in ("kegg", "go_bp"):
            result = pd.read_csv(gsea / f"figure_{number:02d}_{library}_all_results.tsv.gz",
                                 sep="\t", float_precision="round_trip").set_index("term")
            shown = record["gsea"][library]
            require(len(shown["terms"]) <= 4, "Too many on-figure pathways")
            require(np.allclose(result.loc[shown["terms"], "NES"], shown["NES"], rtol=0, atol=1e-14),
                    f"{library}: plotted NES differs from complete results")
            require(np.allclose(result.loc[shown["terms"], "gsea_fdr_q"], shown["gsea_fdr_q"], rtol=0, atol=1e-14),
                    f"{library}: plotted q differs from complete results")
        records.append({"figure": number, "comparison": record["comparison"],
                        "tested_genes": len(table), "joint_DEG_genes": len(significant),
                        "gsea_full_unfiltered_rank_count": len(ranking),
                        "gsea_scores_and_q_match_complete_tables": True,
                        "requested_palette": colors, "panels": 4, "pca_absent": True})
    # These scripts read the actual fonts, outputs and their own freshness hashes.
    for script in (ROOT / "verify_export_fonts.py", ROOT.parent / "accessibility/audit_colors.py"):
        subprocess.run([sys.executable, str(script), "--check-report"], check=True)
    evidence = set(input_hashes)
    evidence.update(str((ROOT.parent/filename).relative_to(REPO)) for record in plot['figures']
                    for filename in record['versioned_export_sha256'])
    evidence.update(str(path.relative_to(REPO)) for path in (
        ROOT / "font_and_plot_audit.json", ROOT / "export_font_verification.json",
        ROOT.parent / "accessibility/palette_simulation.json", gsea / "gsea_provenance.json",
        ROOT / "run_gsea.py", ROOT / "simulate_and_analyze.py", ROOT.parent / "requirements.txt"))
    # Include every GSEA result/source hash rather than trusting just a summary.
    for path in gsea.rglob("*"):
        if path.is_file() and not any(part.startswith(".") for part in path.relative_to(gsea).parts):
            evidence.add(str(path.relative_to(REPO)))
    return {"status": "PASS", "verified_at_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "Synthetic figures only; numerical/table consistency plus current export evidence.",
            "validator_sha256": digest(Path(__file__)),
            "input_sha256": {filename: digest(REPO / filename) for filename in sorted(evidence)},
            "gsea_complete_analyses": len(provenance["analyses"]), "figures": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-report", action="store_true", help="Verify previously checked evidence is still current")
    args = parser.parse_args()
    if args.check_report:
        check_report()
        print("Figure validation: current PASS")
        return
    write_report({"status": "RUNNING"})
    try:
        report = verify()
    except Exception as exc:
        write_report({"status": "FAILED", "error": f"{type(exc).__name__}: {exc}"})
        raise
    write_report(report)
    print("PASS: three figures, joint 0.2/FDR criteria, complete GSEA ranks, palettes and current font/color evidence")


if __name__ == "__main__":
    main()
