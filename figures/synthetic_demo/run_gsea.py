#!/usr/bin/env python3
"""Preranked mouse pathway demonstration using explicitly synthetic statistics.

SIM identifiers are randomly and reproducibly assigned to unique mouse symbols,
independently of every simulated effect or test statistic. Genuine mouse pathway
memberships are then used for an actual weighted preranked GSEA calculation.
This assignment does not turn the invented genes into measured mouse genes;
the results are a plotting demonstration and support no biological conclusions.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import io
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[1]
SEED = 20261019
CONTRASTS = (
    (1, "LD SALINE", "NON-LD control", "contrast_LD_SALINE_vs_NON_LD_control.tsv.gz"),
    (2, "LD LEPTIN", "LD SALINE", "contrast_LD_LEPTIN_vs_LD_SALINE.tsv.gz"),
    (3, "LD LEPTIN", "NON-LD control", "contrast_LD_LEPTIN_vs_NON_LD_control.tsv.gz"),
)
DISCLAIMER = (
    "Synthetic demonstration only. Invented SIM identifiers were assigned to real "
    "mouse symbols by a seeded one-to-one random mapping independent of expression, "
    "effect sizes, test statistics, simulation classes, and pathway membership. "
    "Neither gene assignments nor pathway results are findings from mouse tissue. "
    "No random seeds, pathways, or memberships were selected to obtain significance."
)


def digest(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def saved_path(path: Path) -> str:
    """Prefer checkout-relative paths while retaining external source locations."""
    absolute = path.resolve()
    try:
        return str(absolute.relative_to(REPOSITORY))
    except ValueError:
        return str(absolute)


def write_table(table: pd.DataFrame, path: Path) -> None:
    compression = {"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None
    table.to_csv(path, sep="\t", index=False, compression=compression, float_format="%.17g")


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets = {}
    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) < 3 or not columns[0].strip():
                raise ValueError(f"Invalid GMT row at {path}:{line_number}")
            term = columns[0].strip()
            if term in sets:
                raise ValueError(f"Duplicate GMT term {term!r} in {path}")
            # Preserve authoritative symbol case and remove duplicate memberships.
            genes = sorted({gene.strip() for gene in columns[2:] if gene.strip()})
            if not genes:
                raise ValueError(f"Empty gene set {term!r} in {path}")
            sets[term] = genes
    if not sets:
        raise ValueError(f"Empty GMT collection: {path}")
    return sets


def source_metadata(path: Path, provenance: Path | None) -> dict:
    result = {"gmt_path": saved_path(path), "gmt_sha256": digest(path)}
    if provenance is not None:
        result["provenance_path"] = saved_path(provenance)
        result["provenance_sha256"] = digest(provenance)
        supplied = json.loads(provenance.read_text())
        if supplied.get("sha256") != result["gmt_sha256"]:
            raise ValueError(f"GMT checksum disagrees with its source provenance: {path}")
        if supplied.get("ncbi_taxonomy_id") != 10090:
            raise ValueError(f"Source provenance does not identify native mouse annotation: {provenance}")
        result["provenance"] = supplied
    return result


def pick_display_rows(results: pd.DataFrame, count: int = 6) -> pd.DataFrame:
    """Choose up to three strongest NES in each direction without a q filter."""
    per_direction = count // 2
    positive = results.loc[results.NES > 0].sort_values(
        ["NES", "term"], ascending=[False, True], kind="stable"
    ).head(per_direction)
    negative = results.loc[results.NES < 0].sort_values(
        ["NES", "term"], ascending=[True, True], kind="stable"
    ).head(per_direction)
    selected = pd.concat([positive, negative])
    if len(selected) < min(count, len(results)):
        remaining = results.loc[~results.term.isin(selected.term)].copy()
        remaining["absolute_NES"] = remaining.NES.abs()
        remainder = remaining.sort_values(
            ["absolute_NES", "term"], ascending=[False, True], kind="stable"
        ).head(count - len(selected)).drop(columns="absolute_NES")
        selected = pd.concat([selected, remainder])
    return selected.sort_values(["NES", "term"], ascending=[False, True], kind="stable")


def normalize_result(prerank, sets: dict[str, list[str]], ranked_symbols: set[str]) -> pd.DataFrame:
    result = prerank.res2d.copy().rename(columns={
        "Name": "analysis", "Term": "term", "NOM p-val": "nominal_p",
        "FDR q-val": "gsea_fdr_q", "FWER p-val": "fwer_p",
        "Tag %": "tag_fraction", "Gene %": "rank_fraction", "Lead_genes": "leading_edge_genes",
    })
    required = {"term", "ES", "NES", "nominal_p", "gsea_fdr_q", "fwer_p", "leading_edge_genes"}
    missing = required.difference(result.columns)
    if missing:
        raise RuntimeError(f"Unrecognized GSEApy output; missing columns: {sorted(missing)}")
    for column in ("ES", "NES", "nominal_p", "gsea_fdr_q", "fwer_p"):
        result[column] = pd.to_numeric(result[column], errors="raise")
        if not np.isfinite(result[column]).all():
            raise RuntimeError(f"GSEA returned non-finite {column}")
    for column in ("nominal_p", "gsea_fdr_q", "fwer_p"):
        if not result[column].between(0, 1).all():
            raise RuntimeError(f"GSEA returned out-of-range {column}")
    result["collection_set_size"] = result.term.map(lambda term: len(sets[term]))
    result["ranked_set_size"] = result.term.map(lambda term: len(set(sets[term]) & ranked_symbols))
    result["matched_genes"] = result.term.map(
        lambda term: ";".join(sorted(set(sets[term]) & ranked_symbols))
    )
    result["leading_edge_genes"] = result.leading_edge_genes.fillna("")
    result["leading_edge_size"] = result.leading_edge_genes.map(lambda genes: len(genes.split(";")) if genes else 0)
    result["gsea_fdr_q_below_0_05"] = result.gsea_fdr_q < 0.05
    return result.sort_values(["gsea_fdr_q", "NES", "term"], ascending=[True, False, True], kind="stable")


def run(args: argparse.Namespace) -> dict:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/lipohypolep-gsea-mpl")
    import gseapy

    started = time.monotonic()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    libraries = {
        "kegg": (args.kegg_gmt, args.kegg_provenance),
        "go_bp": (args.go_bp_gmt, args.go_bp_provenance),
    }
    memberships = {key: read_gmt(paths[0]) for key, paths in libraries.items()}
    universe = sorted(set().union(*(set().union(*map(set, sets.values())) for sets in memberships.values())))

    contrasts = []
    tested_ids = None
    for number, numerator, denominator, filename in CONTRASTS:
        path = args.data_dir / filename
        payload = path.read_bytes()
        input_sha256 = hashlib.sha256(payload).hexdigest()
        table = pd.read_csv(
            io.BytesIO(payload), sep="\t", dtype={"gene_id": str},
            compression="gzip" if path.suffix == ".gz" else None,
        )
        if {"gene_id", "stat"}.difference(table.columns):
            raise ValueError(f"Missing gene_id/stat in {path}")
        if table.gene_id.isna().any() or not table.gene_id.str.fullmatch(r"SIM\d+").all():
            raise ValueError(
                "This demonstration accepts invented SIM identifiers only. Refusing "
                "random annotation assignment for real gene identifiers. Experimental "
                "GSEA must use measured mouse identifiers and validated identifier "
                "mapping; never apply this synthetic random assignment to real data."
            )
        if table.gene_id.duplicated().any() or not np.isfinite(table.stat).all():
            raise ValueError(f"Duplicate identifiers or invalid signed statistics in {path}")
        identifiers = set(table.gene_id)
        if tested_ids is None:
            tested_ids = identifiers
        elif identifiers != tested_ids:
            raise ValueError("All three contrasts must contain the same complete tested gene universe")
        contrasts.append((number, numerator, denominator, path, input_sha256, table))
    if len(universe) < len(tested_ids):
        raise ValueError(
            f"Only {len(universe)} unique symbols in the mouse annotation universe for "
            f"{len(tested_ids)} synthetic genes; unique assignment is impossible. "
            "Provide broader verified mouse collections; do not duplicate or reuse identifiers."
        )
    chosen_symbols = np.random.default_rng(args.seed).choice(universe, size=len(tested_ids), replace=False)
    mapping = pd.DataFrame({"synthetic_gene_id": sorted(tested_ids), "assigned_mouse_symbol": chosen_symbols})
    mapping["assignment_type"] = "random_synthetic_demonstration_only"
    mapping_path = output / "synthetic_to_mouse_symbol_map.tsv"
    write_table(mapping, mapping_path)
    symbol_map = mapping.set_index("synthetic_gene_id").assigned_mouse_symbol.to_dict()
    inverse_map = mapping.set_index("assigned_mouse_symbol").synthetic_gene_id.to_dict()
    ranked_symbols = set(chosen_symbols)

    summary = {
        "complete": False,
        "synthetic_disclaimer": DISCLAIMER,
        "species": "Mus musculus",
        "mapping": {
            "file": mapping_path.name, "sha256": digest(mapping_path), "seed": args.seed,
            "algorithm": "NumPy default_rng.choice, sorted library-union universe, no replacement; zipped with lexicographically sorted SIM IDs",
            "independent_of_statistics_and_simulation_classes": True,
            "available_unique_mouse_symbols": len(universe), "assigned_unique_symbols": len(mapping),
        },
        "parameters": {
            "ranking_metric": "unfiltered signed test statistic (stat), descending",
            "gene_level_p_or_fold_change_filter": None,
            "weight": 1.0, "min_size": args.min_size, "max_size": args.max_size,
            "permutation_num": args.permutations, "permutation_type": "gene_set/preranked",
            "seed": args.seed, "threads": args.threads,
            "multiple_testing": "GSEApy empirical GSEA FDR q, separately per contrast and library; not BH",
            "display_selection": "largest absolute NES, balanced up to 3 positive and 3 negative; no significance filter",
            "finite_permutation_note": "Software nominal P and empirical FDR q may equal zero at finite permutation resolution; raw values are preserved without replacement.",
        },
        "software": {
            "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
            "gseapy": importlib.metadata.version("gseapy"),
        },
        "sources": {}, "analyses": [],
    }
    for key, (path, provenance) in libraries.items():
        coverage = pd.DataFrame([
            {"term": term, "collection_set_size": len(genes),
             "ranked_set_size": len(set(genes) & ranked_symbols)}
            for term, genes in memberships[key].items()
        ])
        coverage["eligible"] = coverage.ranked_set_size.between(args.min_size, args.max_size)
        coverage_path = output / f"{key}_all_gene_set_coverage.tsv.gz"
        write_table(coverage, coverage_path)
        summary["sources"][key] = {
            **source_metadata(path, provenance),
            "total_gene_sets": len(coverage),
            "eligible_gene_sets": int(coverage.eligible.sum()),
            "below_min_size": int((coverage.ranked_set_size < args.min_size).sum()),
            "above_max_size": int((coverage.ranked_set_size > args.max_size).sum()),
            "coverage_file": coverage_path.name,
        }
        if not coverage.eligible.any():
            raise ValueError(f"No eligible gene sets for {key}")

    for number, numerator, denominator, path, input_sha256, table in contrasts:
        ranking = table[["gene_id", "stat"]].copy()
        ranking["mouse_symbol"] = ranking.gene_id.map(symbol_map)
        ranking = ranking.sort_values(["stat", "gene_id"], ascending=[False, True], kind="stable")
        ranking_path = output / f"figure_{number:02d}_complete_ranked_statistics.tsv.gz"
        write_table(ranking, ranking_path)
        rnk = ranking[["mouse_symbol", "stat"]]
        for key, sets in memberships.items():
            tick = time.monotonic()
            print(f"Figure {number}: {key}; {len(rnk)} complete ranked genes; {args.permutations} permutations", flush=True)
            pre = gseapy.prerank(
                rnk=rnk, gene_sets=sets, outdir=None, ascending=False,
                organism="mouse", pheno_pos=numerator, pheno_neg=denominator,
                min_size=args.min_size, max_size=args.max_size,
                permutation_num=args.permutations, weight=1.0,
                seed=args.seed, threads=args.threads, no_plot=True, verbose=False,
                method="permutation",
            )
            results = normalize_result(pre, sets, ranked_symbols)
            expected = summary["sources"][key]["eligible_gene_sets"]
            if len(results) != expected:
                raise RuntimeError(f"{key}: expected {expected} eligible sets but GSEA returned {len(results)}")
            results.insert(0, "library", key)
            results.insert(1, "contrast", f"{numerator} vs {denominator}")
            prefix = f"figure_{number:02d}_{key}"
            raw_path = output / f"{prefix}_gseapy_raw.tsv.gz"
            write_table(pre.res2d, raw_path)
            result_path = output / f"{prefix}_all_results.tsv.gz"
            write_table(results, result_path)
            display = pick_display_rows(results)
            display_path = output / f"{prefix}_display.tsv"
            write_table(display, display_path)
            leading = []
            for row in results.itertuples(index=False):
                for symbol in row.leading_edge_genes.split(";"):
                    if symbol:
                        if symbol not in inverse_map:
                            raise RuntimeError(f"Unmapped leading-edge symbol: {symbol}")
                        leading.append({"term": row.term, "assigned_mouse_symbol": symbol,
                                        "synthetic_gene_id": inverse_map[symbol]})
            leading_path = output / f"{prefix}_leading_edge.tsv.gz"
            write_table(pd.DataFrame(leading, columns=["term", "assigned_mouse_symbol", "synthetic_gene_id"]), leading_path)
            elapsed = round(time.monotonic() - tick, 3)
            summary["analyses"].append({
                "figure": number, "library": key, "numerator": numerator, "denominator": denominator,
                "source_contrast_file": saved_path(path), "source_contrast_sha256": input_sha256,
                "ranked_statistics_file": ranking_path.name, "tested_genes": len(ranking),
                "duplicate_rank_statistics": int(ranking.stat.duplicated().sum()),
                "tested_gene_sets": len(results),
                "gsea_fdr_q_below_0_05": int((results.gsea_fdr_q < 0.05).sum()),
                "gsea_fdr_q_below_0_25": int((results.gsea_fdr_q < 0.25).sum()),
                "all_results_file": result_path.name, "raw_results_file": raw_path.name,
                "display_file": display_path.name, "leading_edge_file": leading_path.name,
                "displayed_sets": display[["term", "NES", "nominal_p", "gsea_fdr_q", "ranked_set_size"]].to_dict(orient="records"),
                "runtime_seconds": elapsed,
            })
            (output / "gsea_provenance.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
            print(f"  {len(results)} tested sets; {int((results.gsea_fdr_q < .05).sum())} at GSEA FDR q < .05; {elapsed:.1f}s", flush=True)
    if any(digest(path) != input_sha256 for _, _, _, path, input_sha256, _ in contrasts):
        raise RuntimeError("Contrast inputs changed during GSEA; rerun all analyses on one consistent data version")
    summary["complete"] = True
    summary["runtime_seconds"] = round(time.monotonic() - started, 3)
    summary_path = output / "gsea_provenance.json"
    summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(f"Saved {summary_path}; total {summary['runtime_seconds']:.1f}s", flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kegg-gmt", type=Path, required=True, help="Verified mouse KEGG GMT")
    parser.add_argument("--go-bp-gmt", type=Path, required=True, help="Verified mouse GO biological-process GMT")
    parser.add_argument("--kegg-provenance", type=Path)
    parser.add_argument("--go-bp-provenance", type=Path)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "gsea_data")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    if args.permutations < 1 or args.threads < 1 or not 1 <= args.min_size <= args.max_size:
        parser.error("Require positive permutations/threads and 1 <= min-size <= max-size")
    run(args)


if __name__ == "__main__":
    main()
