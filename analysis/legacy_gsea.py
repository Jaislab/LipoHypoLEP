#!/usr/bin/env python3
"""Exploratory preranked GSEA for imported legacy mouse differential-expression results.

Every tested Ensembl row remains in its contrast's ranking. Only globally
unambiguous exact supplied-symbol mappings enter pathway memberships. The
signed-P rankings are derived proxies, never recovered original Wald statistics.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from analysis.legacy_workbook import FIGURES, validate_completed_import
from figures.synthetic_demo.run_gsea import read_gmt, source_metadata, normalize_result, pick_display_rows, write_manifest

HELPER_PATH = Path(read_gmt.__code__.co_filename).resolve()
RESOURCES = ROOT / "figures/synthetic_demo/gsea_data/resources"
LIBRARIES = {
    "kegg": (RESOURCES / "kegg_mouse_biosystems_archive_2020.gmt", RESOURCES / "kegg_mouse_provenance.json"),
    "go_bp": (RESOURCES / "go_bp_mouse_2026_05_19_nonIEA_protein_coding.gmt", RESOURCES / "go_bp_mouse_provenance.json"),
}
METRICS = {"proxy_signed_log10_p": "primaryproxy_signed_log10_p",
           "proxy_signed_z": "sensitivityproxy_signed_z"}
MT_LEGACY_SYMBOLS = {"COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8"}
MANIFEST = "gsea_manifest.json"


def digest(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_table(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False, float_format="%.17g",
                 compression={"method": "gzip", "mtime": 0} if path.suffix == ".gz" else None)


def load_tables(data_dir: Path) -> tuple[dict, dict[str, pd.DataFrame]]:
    provenance = validate_completed_import(data_dir)
    if provenance.get("tested_universe_confirmation", {}).get("all_tested_genes_exported") is not True:
        raise ValueError("Exploratory GSEA requires recorded confirmation that all tested genes were exported")
    tables = {}
    for figure, _, _, _, _ in FIGURES:
        filename = f"{figure}_deg.tsv.gz"
        payload = (data_dir / filename).read_bytes()
        if hashlib.sha256(payload).hexdigest() != provenance["data_artifact_sha256"][filename]:
            raise ValueError(f"Legacy input changed while loading its validated byte snapshot: {filename}")
        table = pd.read_csv(io.BytesIO(payload), compression="gzip", sep="\t", float_precision="round_trip",
                            dtype={"gene_id": str, "supplied_symbol": str, "entrez_id": str})
        required = {"gene_id", "supplied_symbol", "entrez_id", "log2FoldChange", "pvalue", *METRICS.values()}
        if required.difference(table.columns) or table.empty:
            raise ValueError(f"Missing legacy ranking columns or empty table: {figure}")
        if (table.gene_id.isna().any() or table.gene_id.duplicated().any()
                or not table.gene_id.str.fullmatch(r"ENSMUSG[0-9]+(?:\.[0-9]+)?").all()):
            raise ValueError("Legacy GSEA requires unique genuine mouse Ensembl gene identifiers")
        if (table.supplied_symbol.isna().any() or not table.supplied_symbol.eq(table.supplied_symbol.str.strip()).all()):
            raise ValueError("Supplied symbols must be nonmissing exact annotation strings without surrounding whitespace")
        if not table.pvalue.between(0, 1, inclusive="right").all() or not np.isfinite(table.log2FoldChange).all():
            raise ValueError("Finite effects and P values in (0,1] are required; no probability floor is invented")
        direction = np.sign(table.log2FoldChange.to_numpy())
        expected = {"primaryproxy_signed_log10_p": direction * -np.log10(table.pvalue.to_numpy()),
                    "sensitivityproxy_signed_z": direction * norm.isf(table.pvalue.to_numpy() / 2)}
        for column, values in expected.items():
            if not np.isfinite(values).all() or not np.array_equal(values, table[column].to_numpy()):
                raise ValueError(f"Imported exploratory ranking proxy is nonfinite or changed: {figure}/{column}")
        tables[figure] = table
    return provenance, tables


def exact_mapping(tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, dict[str, str]]:
    """Exclude ambiguous membership assignments without dropping any ranked row."""
    annotation = pd.concat([table[["gene_id", "supplied_symbol", "entrez_id"]]
                            for table in tables.values()], ignore_index=True).drop_duplicates()
    symbol_counts = annotation.groupby("supplied_symbol").gene_id.nunique()
    ambiguous = set(symbol_counts[symbol_counts > 1].index)
    gene_symbols = annotation.groupby("gene_id").supplied_symbol.nunique()
    gene_entrez = annotation.groupby("gene_id").entrez_id.nunique()
    inconsistent = set(gene_symbols[gene_symbols > 1].index) | set(gene_entrez[gene_entrez > 1].index)
    annotation["mapping_status"] = "unambiguous_exact_supplied_symbol"
    annotation.loc[annotation.supplied_symbol.isin(ambiguous), "mapping_status"] = "ambiguous_symbol_membership_excluded"
    annotation.loc[annotation.gene_id.isin(inconsistent), "mapping_status"] = "inconsistent_gene_annotation_membership_excluded"
    annotation = annotation.sort_values(["gene_id", "supplied_symbol", "entrez_id"], kind="stable").reset_index(drop=True)
    accepted = annotation.loc[annotation.mapping_status.eq("unambiguous_exact_supplied_symbol")]
    mapping = accepted.set_index("supplied_symbol").gene_id.to_dict()
    if len(mapping) != len(accepted) or len(set(mapping.values())) != len(mapping):
        raise RuntimeError("Exact mapping is not globally one-to-one")
    return annotation, mapping


def translate_library(sets: dict[str, list[str]], mapping: dict[str, str]) -> dict[str, list[str]]:
    return {term: sorted({mapping[symbol] for symbol in members if symbol in mapping})
            for term, members in sets.items()}


def validate_completed_run(output_dir: Path, data_dir: Path | None = None) -> dict:
    output_dir = Path(output_dir)
    path = output_dir / MANIFEST
    payload = path.read_bytes()
    summary = json.loads(payload)
    if summary.get("complete") is not True or summary.get("integrity_schema") != 1:
        raise ValueError("Legacy GSEA is incomplete or has an unsupported integrity schema")
    if summary.get("runner_sha256") != digest(Path(__file__)):
        raise ValueError("Legacy GSEA runner changed; regenerate results")
    if summary.get("pathway_helper_sha256") != digest(HELPER_PATH):
        raise ValueError("Legacy GSEA pathway helper source changed; regenerate results")
    if set(summary.get("sources", {})) != {"kegg", "go_bp"}:
        raise ValueError("Legacy GSEA requires exactly KEGG and GO BP annotation sources")
    metric = summary.get("rank_metric")
    expected = {(figure, library, metric) for figure, *_ in FIGURES for library in LIBRARIES}
    observed = [(item.get("figure_id"), item.get("library"), item.get("rank_metric")) for item in summary.get("analyses", [])]
    if metric not in METRICS or len(observed) != 6 or set(observed) != expected:
        raise ValueError("Legacy GSEA requires exactly six expected analyses for one ranking metric")
    hashes = summary.get("artifact_sha256", {})
    names = {"exact_symbol_mapping.tsv", "annotation_mapping_gaps.tsv", "kegg_coverage.tsv.gz", "go_bp_coverage.tsv.gz"}
    for item in summary["analyses"]:
        figure, library = item["figure_id"], item["library"]
        direction = next((numerator, denominator) for name, _, numerator, denominator, _ in FIGURES if name == figure)
        if (item.get("numerator"), item.get("denominator")) != direction:
            raise ValueError("Legacy GSEA contrast direction disagrees with the requested figure")
        expected_files = {"complete_rank_file": f"{figure}_complete_rank.tsv.gz",
                          "all_results_file": f"{figure}_{library}_all_results.tsv.gz",
                          "raw_results_file": f"{figure}_{library}_gseapy_raw.tsv.gz",
                          "display_file": f"{figure}_{library}_display.tsv",
                          "leading_edge_file": f"{figure}_{library}_leading_edge.tsv.gz"}
        if any(item.get(field) != value for field, value in expected_files.items()):
            raise ValueError("Unexpected legacy GSEA output filename")
        names.update(expected_files.values())
    if set(hashes) != names:
        raise ValueError("Legacy GSEA checksums do not cover every required artifact")
    for name, checksum in hashes.items():
        artifact = output_dir / name
        if artifact.is_symlink() or not artifact.is_file() or digest(artifact) != checksum:
            raise ValueError(f"Legacy GSEA output changed or is missing: {name}")
    for source in summary["sources"].values():
        for field, checksum in (("gmt_path", "gmt_sha256"), ("provenance_path", "provenance_sha256")):
            source_path = Path(source[field])
            source_path = source_path if source_path.is_absolute() else ROOT / source_path
            if not source_path.is_file() or digest(source_path) != source[checksum]:
                raise ValueError("Legacy GSEA annotation source changed")
    if data_dir is not None:
        imported = validate_completed_import(data_dir)
        if imported["data_artifact_sha256"] != summary["input_table_sha256"] or imported["input_sha256"] != summary["source_workbook_sha256"]:
            raise ValueError("Legacy GSEA inputs changed")
    if path.read_bytes() != payload:
        raise ValueError("Legacy GSEA completion changed during validation")
    return summary


def run_metric(args: argparse.Namespace, metric: str, output: Path) -> dict:
    """Publish one independently complete six-analysis ranking family."""
    import gseapy

    if metric not in METRICS or args.permutations < 1 or args.threads < 1 or not 1 <= args.min_size <= args.max_size:
        raise ValueError("Invalid ranking metric or GSEA parameters")
    output.mkdir(parents=True, exist_ok=True)
    with (output / ".gsea.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError("Another legacy GSEA process is writing this output") from error
        started = time.monotonic()
        runner_hash = digest(Path(__file__))
        helper_hash = digest(HELPER_PATH)
        write_manifest({"complete": False, "rank_metric": metric}, output / MANIFEST)
        imported, tables = load_tables(args.data_dir)
        annotation, symbol_map = exact_mapping(tables)
        libraries = {"kegg": (args.kegg_gmt, args.kegg_provenance), "go_bp": (args.go_bp_gmt, args.go_bp_provenance)}
        snapshots = {path.resolve(): path.read_bytes() for pair in libraries.values() for path in pair}
        native = {library: read_gmt(pair[0], snapshots[pair[0].resolve()]) for library, pair in libraries.items()}
        translated = {library: translate_library(sets, symbol_map) for library, sets in native.items()}
        sources = {library: source_metadata(pair[0], pair[1], snapshots[pair[0].resolve()], snapshots[pair[1].resolve()])
                   for library, pair in libraries.items()}
        summary = {
            "complete": False, "integrity_schema": 1, "rank_metric": metric,
            "ranking_status": "Exploratory signed-P proxy; original signed test statistic is unavailable, not recovered.",
            "all_tested_genes_exported_user_confirmed": True,
            "source_workbook_sha256": imported["input_sha256"],
            "source_workbook_basename": imported["input_file"],
            "input_table_sha256": imported["data_artifact_sha256"], "sources": sources, "analyses": [],
            "parameters": {"weight": 1.0, "permutation_num": args.permutations, "min_size": args.min_size,
                           "max_size": args.max_size, "threads": args.threads, "seed": args.seed,
                           "permutation_type": "gene-set/preranked", "gene_level_filter": None,
                           "multiple_testing": "GSEApy empirical GSEA FDR q separately per contrast/collection; not BH",
                           "display_selection": "up to two strongest NES in each direction, remainder strongest absolute NES; no significance filter"},
            "mapping_policy": "Exact supplied symbols translated to full Ensembl gene IDs through globally one-to-one mappings across all three contrasts; ambiguous symbol memberships excluded, every tested Ensembl row retained in ranking/background; no aliases, random assignment or strongest-statistic duplicate selection.",
            "limitations": ["Missing original signed test statistics require an explicitly exploratory ranking proxy.",
                            "Contrasts have complete but different tested universes and observed baseline-expression filters; compare pathway results cautiously.",
                            "KEGG is the archived BioSystems2020 collection, not a current KEGG release.",
                            "GO uses the pinned2026-05-19 non-IEA protein-coding projection with upstream size restrictions.",
                            "The 13 legacy mitochondrial symbols match archived KEGG but not current GO; no unvalidated alias conversion is performed.",
                            "No sample-level expression, counts, covariates or replicate metadata are available; this is not a refitted expression analysis."],
            "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
                         "gseapy": importlib.metadata.version("gseapy")},
            "started_utc": datetime.now(timezone.utc).isoformat(), "runner_sha256": runner_hash,
            "pathway_helper_sha256": helper_hash,
        }
        with tempfile.TemporaryDirectory(prefix=".legacy-gsea-stage-", dir=output) as temporary:
            stage = Path(temporary)
            write_table(annotation, stage / "exact_symbol_mapping.tsv")
            gap_rows, coverage_rows = [], {library: [] for library in libraries}
            for library, sets in native.items():
                symbols = set().union(*map(set, sets.values()))
                for row in annotation.itertuples(index=False):
                    if row.mapping_status != "unambiguous_exact_supplied_symbol" or row.supplied_symbol not in symbols:
                        gap_rows.append({"library": library, "gene_id": row.gene_id, "supplied_symbol": row.supplied_symbol,
                                         "reason": row.mapping_status if row.mapping_status != "unambiguous_exact_supplied_symbol" else "no_exact_symbol_in_collection",
                                         "legacy_mitochondrial_symbol": row.supplied_symbol in MT_LEGACY_SYMBOLS})
            write_table(pd.DataFrame(gap_rows, columns=["library", "gene_id", "supplied_symbol", "reason", "legacy_mitochondrial_symbol"]), stage / "annotation_mapping_gaps.tsv")
            symbol_lookup = annotation.drop_duplicates("gene_id").set_index("gene_id").supplied_symbol.to_dict()
            entrez_lookup = annotation.drop_duplicates("gene_id").set_index("gene_id").entrez_id.to_dict()
            for figure, _, numerator, denominator, _ in FIGURES:
                table = tables[figure]
                ranking = table[["gene_id", "supplied_symbol", METRICS[metric]]].rename(columns={METRICS[metric]: "stat"})
                ranking = ranking.sort_values(["stat", "gene_id"], ascending=[False, True], kind="stable")
                rank_file = f"{figure}_complete_rank.tsv.gz"
                write_table(ranking, stage / rank_file)
                ranked_ids = set(ranking.gene_id)
                for library, sets in translated.items():
                    eligible = {}
                    for term, members in sets.items():
                        size = len(set(members) & ranked_ids)
                        valid = args.min_size <= size <= args.max_size
                        coverage_rows[library].append({"figure_id": figure, "term": term,
                            "native_collection_set_size": len(native[library][term]), "mapped_collection_set_size": len(members),
                            "ranked_set_size": size, "eligible": valid})
                        if valid:
                            eligible[term] = members
                    if not eligible:
                        raise ValueError(f"No eligible pathways for {figure}/{library}")
                    tick = time.monotonic()
                    print(f"{metric} {figure} {library}: {len(ranking)} genes, {len(eligible)} pathways, {args.permutations} permutations", flush=True)
                    rnk = ranking[["gene_id", "stat"]].copy()
                    pre = gseapy.prerank(rnk=rnk, gene_sets=eligible, outdir=None, ascending=None,
                        organism="mouse", pheno_pos=numerator, pheno_neg=denominator, min_size=args.min_size,
                        max_size=args.max_size, permutation_num=args.permutations, weight=1.0, seed=args.seed,
                        threads=args.threads, no_plot=True, verbose=False, method="permutation")
                    if pre.ranking.index.tolist() != ranking.gene_id.tolist() or not np.array_equal(pre.ranking.to_numpy(), ranking.stat.to_numpy()):
                        raise RuntimeError("GSEA changed the complete deterministic ranking")
                    result = normalize_result(pre, eligible, ranked_ids)
                    if set(result.term) != set(eligible):
                        raise RuntimeError("GSEA did not test every eligible pathway exactly once")
                    result.insert(0, "library", library)
                    result.insert(1, "contrast", f"{numerator} vs {denominator}")
                    result["native_collection_set_size"] = result.term.map(lambda term: len(native[library][term]))
                    prefix = f"{figure}_{library}"
                    result_file, raw_file = f"{prefix}_all_results.tsv.gz", f"{prefix}_gseapy_raw.tsv.gz"
                    display_file, leading_file = f"{prefix}_display.tsv", f"{prefix}_leading_edge.tsv.gz"
                    write_table(result, stage / result_file)
                    write_table(pre.res2d, stage / raw_file)
                    display = pick_display_rows(result, count=4)
                    write_table(display, stage / display_file)
                    leading = [{"term": row.term, "gene_id": gene, "supplied_symbol": symbol_lookup[gene], "entrez_id": entrez_lookup[gene]}
                               for row in result.itertuples(index=False) for gene in row.leading_edge_genes.split(";") if gene]
                    write_table(pd.DataFrame(leading, columns=["term", "gene_id", "supplied_symbol", "entrez_id"]), stage / leading_file)
                    summary["analyses"].append({"figure_id": figure, "library": library, "rank_metric": metric,
                        "numerator": numerator, "denominator": denominator, "ranked_genes": len(ranking),
                        "all_tested_ensembl_rows_retained": len(ranking) == len(table),
                        "tested_gene_sets": len(result), "genes_with_pathway_membership": len(set().union(*map(set, sets.values())) & ranked_ids),
                        "gsea_fdr_q_below_0_05": int((result.gsea_fdr_q < .05).sum()),
                        "complete_rank_file": rank_file, "all_results_file": result_file, "raw_results_file": raw_file,
                        "display_file": display_file, "leading_edge_file": leading_file,
                        "displayed_sets": display[["term", "NES", "gsea_fdr_q"]].to_dict(orient="records"),
                        "runtime_seconds": round(time.monotonic() - tick, 3)})
                    print(f"  {int((result.gsea_fdr_q < .05).sum())} at empirical FDRq<.05; {time.monotonic()-tick:.1f}s", flush=True)
            for library, rows in coverage_rows.items():
                write_table(pd.DataFrame(rows), stage / f"{library}_coverage.tsv.gz")
            final_import = validate_completed_import(args.data_dir)
            if (final_import["data_artifact_sha256"] != summary["input_table_sha256"]
                    or final_import["input_sha256"] != summary["source_workbook_sha256"]):
                raise RuntimeError("Legacy input tables changed during enrichment")
            if (any(not path.is_file() or path.read_bytes() != payload for path, payload in snapshots.items())
                    or digest(Path(__file__)) != runner_hash or digest(HELPER_PATH) != helper_hash):
                raise RuntimeError("Legacy annotation or runner source changed during enrichment")
            summary["artifact_sha256"] = {path.name: digest(path) for path in sorted(stage.iterdir())}
            summary["runtime_seconds"] = round(time.monotonic() - started, 3)
            summary["completed_utc"] = datetime.now(timezone.utc).isoformat()
            for path in sorted(stage.iterdir()):
                path.replace(output / path.name)
            summary["complete"] = True
            write_manifest(summary, output / MANIFEST)
    return summary


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/lipohypolep-legacy-gsea-mpl")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "analysis/legacy_results")
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--rank-metric", choices=[*METRICS, "both"], default="both")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261020)
    for library, (gmt, provenance) in LIBRARIES.items():
        parser.add_argument(f"--{library.replace('_', '-')}-gmt", type=Path, default=gmt)
        parser.add_argument(f"--{library.replace('_', '-')}-provenance", type=Path, default=provenance)
    args = parser.parse_args()
    destination = args.output_dir or args.data_dir / "gsea"
    metrics = list(METRICS) if args.rank_metric == "both" else [args.rank_metric]
    for metric in metrics:
        output = destination / "signed_z" if metric == "proxy_signed_z" else destination
        run_metric(args, metric, output)
        validate_completed_run(output, args.data_dir)
        print(f"Completed {metric}: {output / MANIFEST}", flush=True)


if __name__ == "__main__":
    main()
