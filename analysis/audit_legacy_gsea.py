#!/usr/bin/env python3
"""Independently audit legacy weighted enrichment scores and ranking sensitivity.

The defining running sum and leading edges are recalculated without GSEApy.
This does not independently reproduce null permutations, NES normalization,
or empirical FDR. The sensitivity comparison describes two weighting choices;
neither derived ranking recovers the unavailable original signed statistic.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from analysis.legacy_gsea import validate_completed_run, digest, write_manifest


def reference_walk(genes: np.ndarray, statistics: np.ndarray, members: set[str], positions: dict) -> tuple[float, set[str]]:
    """Calculate weight-1 ES and leading edge at the hit/miss boundaries."""
    hits = np.array(sorted(positions[gene] for gene in members if gene in positions), dtype=int)
    if not 0 < len(hits) < len(genes):
        raise ValueError("Reference enrichment requires both hits and misses")
    weights = np.abs(statistics[hits])
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("Reference enrichment requires positive finite hit weight")
    cumulative_weight = np.cumsum(weights / weights.sum())
    cumulative_misses = (hits - np.arange(len(hits))) / (len(genes) - len(hits))
    after_hit = cumulative_weight - cumulative_misses
    before_hit = np.r_[0., cumulative_weight[:-1]] - cumulative_misses
    peak, trough = int(after_hit.argmax()), int(before_hit.argmin())
    if after_hit[peak] > -before_hit[trough]:
        return float(after_hit[peak]), set(genes[hits[:peak + 1]])
    return float(before_hit[trough]), set(genes[hits[trough:]])


def audit_scores(output_dir: Path, data_dir: Path) -> dict:
    manifest = validate_completed_run(output_dir, data_dir)
    mapping = pd.read_csv(output_dir / "exact_symbol_mapping.tsv", sep="\t", dtype=str, keep_default_na=False)
    accepted = mapping.loc[mapping.mapping_status.eq("unambiguous_exact_supplied_symbol")]
    symbol_map = dict(zip(accepted.supplied_symbol, accepted.gene_id))
    inputs = {figure: pd.read_csv(data_dir / f"{figure}_deg.tsv.gz", sep="\t", float_precision="round_trip")
              for figure in ("figure_01", "figure_02", "figure_03")}
    source_annotations = pd.concat([table[["gene_id", "supplied_symbol"]] for table in inputs.values()]).drop_duplicates()
    if set(zip(mapping.gene_id, mapping.supplied_symbol)) != set(zip(source_annotations.gene_id, source_annotations.supplied_symbol)):
        raise ValueError("Reported mapping does not preserve every actual supplied identifier/symbol pair")
    symbol_counts = source_annotations.groupby("supplied_symbol").gene_id.nunique()
    if any(symbol_counts[symbol] != 1 for symbol in symbol_map):
        raise ValueError("Ambiguous source symbols were incorrectly assigned to pathway memberships")
    libraries = {}
    for library, metadata in manifest["sources"].items():
        path = Path(metadata["gmt_path"])
        path = path if path.is_absolute() else ROOT / path
        libraries[library] = {}
        for line in path.read_text(encoding="utf-8-sig").splitlines():
            if line.strip():
                fields = line.split("\t")
                libraries[library][fields[0]] = {symbol_map[symbol] for symbol in fields[2:] if symbol in symbol_map}
    analyses, cache = [], {}
    tie_audit = {}
    for analysis in manifest["analyses"]:
        rank_name = analysis["complete_rank_file"]
        if rank_name not in cache:
            ranking = pd.read_csv(output_dir / rank_name, sep="\t", float_precision="round_trip")
            source = inputs[analysis["figure_id"]].set_index("gene_id")
            column = {"proxy_signed_log10_p": "primaryproxy_signed_log10_p", "proxy_signed_z": "sensitivityproxy_signed_z"}[manifest["rank_metric"]]
            if (not ranking.gene_id.is_unique or set(ranking.gene_id) != set(source.index)
                    or not ranking.stat.is_monotonic_decreasing
                    or not np.array_equal(ranking.stat.to_numpy(), source.loc[ranking.gene_id, column].to_numpy())):
                raise ValueError("Ranking does not preserve every exact supplied gene/proxy statistic")
            genes = ranking.gene_id.to_numpy()
            cache[rank_name] = genes, ranking.stat.to_numpy(), {gene: index for index, gene in enumerate(genes)}
            counts = ranking.stat.value_counts()
            tied = ranking.loc[ranking.stat.duplicated(keep=False)]
            if any(group.gene_id.tolist() != sorted(group.gene_id.tolist()) for _, group in tied.groupby("stat")):
                raise ValueError("Exact ranking ties do not follow deterministic Ensembl ordering")
            tie_audit[rank_name] = {"duplicate_statistic_extra_rows": int(ranking.stat.duplicated().sum()),
                "tied_gene_rows": int((counts[counts > 1]).sum()), "tie_groups": int((counts > 1).sum()),
                "policy": "Descending signed proxy; lexicographic Ensembl identifier within exact ties; no jitter."}
        genes, statistics, positions = cache[rank_name]
        ranked_ids = set(positions)
        results = pd.read_csv(output_dir / analysis["all_results_file"], sep="\t", float_precision="round_trip", keep_default_na=False)
        errors, edge_errors, membership_errors = [], [], []
        for row in results.itertuples(index=False):
            members = libraries[analysis["library"]][row.term]
            score, leading = reference_walk(genes, statistics, members, positions)
            errors.append(abs(score - row.ES))
            reported = set(row.leading_edge_genes.split(";")) if row.leading_edge_genes else set()
            if reported != leading:
                edge_errors.append(row.term)
            reported_members = set(row.matched_genes.split(";")) if row.matched_genes else set()
            if reported_members != members & ranked_ids:
                membership_errors.append(row.term)
        maximum = max(errors, default=0.)
        if maximum > 2e-12 or edge_errors or membership_errors or len(results) != analysis["tested_gene_sets"]:
            raise ValueError(f"Independent legacy ES/leading-edge audit failed: {analysis['figure_id']}/{analysis['library']}")
        analyses.append({"figure_id": analysis["figure_id"], "library": analysis["library"],
            "rank_metric": manifest["rank_metric"], "tested_sets": len(results), "maximum_absolute_ES_error": maximum,
            "leading_edge_mismatches": len(edge_errors), "matched_membership_mismatches": len(membership_errors),
            "rank_ties": tie_audit[rank_name]})
    return {"rank_metric": manifest["rank_metric"], "manifest_sha256": digest(output_dir / "gsea_manifest.json"),
            "audited_sets": sum(item["tested_sets"] for item in analyses), "analyses": analyses}


def compare_metrics(primary_dir: Path, sensitivity_dir: Path, data_dir: Path) -> dict:
    primary = validate_completed_run(primary_dir, data_dir)
    sensitivity = validate_completed_run(sensitivity_dir, data_dir)
    if primary["input_table_sha256"] != sensitivity["input_table_sha256"]:
        raise ValueError("Sensitivity analyses use different input data")
    secondary = {(item["figure_id"], item["library"]): item for item in sensitivity["analyses"]}
    comparisons = []
    for analysis in primary["analyses"]:
        other = secondary[analysis["figure_id"], analysis["library"]]
        first = pd.read_csv(primary_dir / analysis["all_results_file"], sep="\t", float_precision="round_trip").set_index("term")
        second = pd.read_csv(sensitivity_dir / other["all_results_file"], sep="\t", float_precision="round_trip").set_index("term")
        if set(first.index) != set(second.index):
            raise ValueError("Sensitivity analyses test different pathway universes")
        second = second.loc[first.index]
        first_hits = set(first.index[first.gsea_fdr_q < .05])
        second_hits = set(second.index[second.gsea_fdr_q < .05])
        common_hits = first_hits & second_hits
        same_direction_hits = {term for term in common_hits if np.sign(first.loc[term, "NES"]) == np.sign(second.loc[term, "NES"])}
        disagreement = first.loc[(first.gsea_fdr_q < .05) != (second.gsea_fdr_q < .05), ["NES", "gsea_fdr_q"]].copy()
        disagreement.columns = ["primary_NES", "primary_q"]
        disagreement["sensitivity_NES"] = second.loc[disagreement.index, "NES"]
        disagreement["sensitivity_q"] = second.loc[disagreement.index, "gsea_fdr_q"]
        disagreement["minimum_q"] = disagreement[["primary_q", "sensitivity_q"]].min(axis=1)
        disagreement = disagreement.sort_values(["minimum_q"], kind="stable").drop(columns="minimum_q")
        top_first = set(pd.read_csv(primary_dir / analysis["display_file"], sep="\t").term)
        top_second = set(pd.read_csv(sensitivity_dir / other["display_file"], sep="\t").term)
        first_rank = pd.read_csv(primary_dir / analysis["complete_rank_file"], sep="\t")
        second_rank = pd.read_csv(sensitivity_dir / other["complete_rank_file"], sep="\t")
        comparisons.append({"figure_id": analysis["figure_id"], "library": analysis["library"], "tested_sets": len(first),
            "identical_gene_order": first_rank.gene_id.tolist() == second_rank.gene_id.tolist(),
            "NES_spearman_correlation": float(spearmanr(first.NES, second.NES).statistic),
            "NES_pearson_correlation": float(np.corrcoef(first.NES, second.NES)[0, 1]),
            "NES_direction_disagreements": int((np.sign(first.NES) != np.sign(second.NES)).sum()),
            "primary_q_below_0_05": len(first_hits), "sensitivity_q_below_0_05": len(second_hits),
            "significant_in_both": len(first_hits & second_hits), "primary_only": len(first_hits - second_hits),
            "significant_in_both_same_direction": len(same_direction_hits),
            "significant_in_both_opposite_direction": len(common_hits - same_direction_hits),
            "significant_in_both_same_direction_terms": sorted(same_direction_hits),
            "sensitivity_only": len(second_hits - first_hits),
            "significant_set_jaccard": len(first_hits & second_hits) / len(first_hits | second_hits) if first_hits | second_hits else 1.,
            "top_display_overlap": len(top_first & top_second), "top_primary_terms": sorted(top_first), "top_sensitivity_terms": sorted(top_second),
            "largest_threshold_disagreements_by_minimum_q": disagreement.head(10).reset_index().to_dict(orient="records")})
    return {"primary_metric": primary["rank_metric"], "sensitivity_metric": sensitivity["rank_metric"],
        "interpretation": "Same ordered genes, different monotone weighting of signed source P values; empirical q and NES can differ. Neither is the unavailable original test statistic.",
        "comparisons": comparisons}


def validate_audit_report(data_dir: Path, report_path: Path | None = None) -> dict:
    """Check report freshness against both complete runs without recalculating ES."""
    report_path = report_path or data_dir / "legacy_gsea_audit.json"
    report = json.loads(report_path.read_text())
    if report.get("complete") is not True or report.get("auditor_sha256") != digest(Path(__file__)):
        raise ValueError("Independent legacy GSEA audit is incomplete or its auditor source changed")
    directories = {"proxy_signed_log10_p": data_dir / "gsea", "proxy_signed_z": data_dir / "gsea/signed_z"}
    recorded = {item["rank_metric"]: item for item in report.get("score_audits", [])}
    if set(recorded) != set(directories) or len(report["score_audits"]) != 2:
        raise ValueError("Independent audit requires both primary and sensitivity ranking families")
    for metric, directory in directories.items():
        manifest = validate_completed_run(directory, data_dir)
        audit = recorded[metric]
        if (audit.get("manifest_sha256") != digest(directory / "gsea_manifest.json")
                or audit.get("audited_sets") != sum(item["tested_gene_sets"] for item in manifest["analyses"])
                or len(audit.get("analyses", [])) != 6):
            raise ValueError("Independent legacy GSEA audit is stale for a current enrichment run")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "analysis/legacy_results")
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    auditor_hash = digest(Path(__file__))
    primary = args.data_dir / "gsea"
    sensitivity = primary / "signed_z"
    report = {"complete": True, "auditor_sha256": auditor_hash,
        "independent_reference": "Weight-1 running sum at all hit/miss boundaries, independent of GSEApy.",
        "scope_limit": "Does not independently reproduce null permutations, NES normalization or empirical FDR calculation.",
        "score_audits": [audit_scores(primary, args.data_dir), audit_scores(sensitivity, args.data_dir)],
        "ranking_sensitivity": compare_metrics(primary, sensitivity, args.data_dir)}
    output = args.output or args.data_dir / "legacy_gsea_audit.json"
    if digest(Path(__file__)) != auditor_hash:
        raise ValueError("Independent auditor source changed during calculation")
    write_manifest(report, output)
    validate_audit_report(args.data_dir, output)
    print(json.dumps({"complete": True, "audited_sets": sum(item["audited_sets"] for item in report["score_audits"]), "output": output.name}))


if __name__ == "__main__":
    main()
