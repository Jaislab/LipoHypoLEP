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
import fcntl
import hashlib
import importlib.metadata
import io
import json
import os
import platform
import time
import tempfile
from datetime import datetime, timezone
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


def read_gmt(path: Path, payload: bytes | None = None) -> dict[str, list[str]]:
    """Parse the same byte snapshot that is recorded in source provenance."""
    sets = {}
    with io.StringIO((path.read_bytes() if payload is None else payload).decode("utf-8-sig")) as handle:
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


def source_metadata(path: Path, provenance: Path | None, payload: bytes | None = None,
                    provenance_payload: bytes | None = None) -> dict:
    result = {"gmt_path": saved_path(path), "gmt_sha256": hashlib.sha256(path.read_bytes() if payload is None else payload).hexdigest()}
    if provenance is not None:
        result["provenance_path"] = saved_path(provenance)
        provenance_payload = provenance.read_bytes() if provenance_payload is None else provenance_payload
        result["provenance_sha256"] = hashlib.sha256(provenance_payload).hexdigest()
        supplied = json.loads(provenance_payload)
        if supplied.get("sha256") != result["gmt_sha256"]:
            raise ValueError(f"GMT checksum disagrees with its source provenance: {path}")
        if supplied.get("ncbi_taxonomy_id") != 10090:
            raise ValueError(f"Source provenance does not identify native mouse annotation: {provenance}")
        result["provenance"] = supplied
    return result


def pick_display_rows(results: pd.DataFrame, count: int = 6) -> pd.DataFrame:
    """Choose up to three strongest NES in each direction without a q filter."""
    if not isinstance(count, int) or count < 0:
        raise ValueError("Display count must be a nonnegative integer")
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
    if result.term.isna().any() or result.term.duplicated().any() or not set(result.term).issubset(sets):
        raise RuntimeError("GSEA returned missing, duplicate, or unrecognized gene-set terms")
    for column in ("ES", "NES", "nominal_p", "gsea_fdr_q", "fwer_p"):
        result[column] = pd.to_numeric(result[column], errors="raise")
        if not np.isfinite(result[column]).all():
            raise RuntimeError(f"GSEA returned non-finite {column}")
    for column in ("nominal_p", "gsea_fdr_q", "fwer_p"):
        if not result[column].between(0, 1).all():
            raise RuntimeError(f"GSEA returned out-of-range {column}")
    if not result.ES.between(-1, 1).all() or ((result.ES * result.NES) < 0).any():
        raise RuntimeError("GSEA returned invalid enrichment-score magnitude or direction")
    for row in result.itertuples(index=False):
        leading = [] if pd.isna(row.leading_edge_genes) or row.leading_edge_genes == "" else str(row.leading_edge_genes).split(";")
        if len(set(leading)) != len(leading) or not set(leading).issubset(set(sets[row.term]) & ranked_symbols):
            raise RuntimeError(f"Invalid leading-edge members for {row.term}")
    result["collection_set_size"] = result.term.map(lambda term: len(sets[term]))
    result["ranked_set_size"] = result.term.map(lambda term: len(set(sets[term]) & ranked_symbols))
    result["matched_genes"] = result.term.map(
        lambda term: ";".join(sorted(set(sets[term]) & ranked_symbols))
    )
    result["leading_edge_genes"] = result.leading_edge_genes.fillna("")
    result["leading_edge_size"] = result.leading_edge_genes.map(lambda genes: len(genes.split(";")) if genes else 0)
    result["gsea_fdr_q_below_0_05"] = result.gsea_fdr_q < 0.05
    return result.sort_values(["gsea_fdr_q", "NES", "term"], ascending=[True, False, True], kind="stable")


def write_manifest(summary: dict, path: Path) -> None:
    """Replace a manifest atomically; incomplete runs must never look complete."""
    payload = json.dumps(summary, indent=2, allow_nan=False) + "\n"
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                     prefix=".gsea-manifest-", delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def artifact_names(summary: dict) -> set[str]:
    """Require exactly the six intended analyses and their expected filenames."""
    if set(summary.get("sources", {})) != {"kegg", "go_bp"}:
        raise ValueError("GSEA manifest requires KEGG and GO BP source collections")
    expected = {(number, library): (numerator, denominator, filename)
                for number, numerator, denominator, filename in CONTRASTS
                for library in ("kegg", "go_bp")}
    analyses = summary.get("analyses", [])
    observed = [(item.get("figure"), item.get("library")) for item in analyses]
    if len(observed) != len(expected) or set(observed) != set(expected):
        raise ValueError("GSEA manifest requires six unique expected contrast/library analyses")
    names = {"synthetic_to_mouse_symbol_map.tsv", "kegg_all_gene_set_coverage.tsv.gz",
             "go_bp_all_gene_set_coverage.tsv.gz"}
    if summary.get("mapping", {}).get("file") != "synthetic_to_mouse_symbol_map.tsv":
        raise ValueError("Unexpected GSEA mapping filename")
    for library, source in summary["sources"].items():
        if source.get("coverage_file") != f"{library}_all_gene_set_coverage.tsv.gz":
            raise ValueError("Unexpected GSEA coverage filename")
    for item in analyses:
        number, library = item["figure"], item["library"]
        numerator, denominator, filename = expected[number, library]
        if (item.get("numerator"), item.get("denominator")) != (numerator, denominator):
            raise ValueError("GSEA contrast direction disagrees with its figure")
        if Path(item.get("source_contrast_file", "")).name != filename:
            raise ValueError("Unexpected GSEA contrast source filename")
        fields = {"ranked_statistics_file": f"figure_{number:02d}_complete_ranked_statistics.tsv.gz",
                  "all_results_file": f"figure_{number:02d}_{library}_all_results.tsv.gz",
                  "raw_results_file": f"figure_{number:02d}_{library}_gseapy_raw.tsv.gz",
                  "display_file": f"figure_{number:02d}_{library}_display.tsv",
                  "leading_edge_file": f"figure_{number:02d}_{library}_leading_edge.tsv.gz"}
        if any(item.get(key) != value for key, value in fields.items()):
            raise ValueError("Unexpected GSEA output artifact filename")
        names.update(fields.values())
    return names


def validate_completed_run(output_dir: Path, data_dir: Path | None = None) -> dict:
    """Fail closed on incomplete, stale, missing, modified, or mixed GSEA outputs.

    Source paths recorded relative to the checkout are resolved against this
    checkout. ``data_dir`` optionally relocates the three contrast inputs.
    """
    output_dir = Path(output_dir)
    manifest_path = output_dir / "gsea_provenance.json"
    manifest_payload = manifest_path.read_bytes()
    summary = json.loads(manifest_payload)
    if summary.get("complete") is not True:
        raise ValueError("GSEA run is incomplete; rerun all six analyses before plotting")
    if summary.get("integrity_schema") != 1 or summary.get("runner_sha256") != digest(Path(__file__)):
        raise ValueError("GSEA integrity schema or runner source changed; regenerate the enrichment results")
    names = artifact_names(summary)
    hashes = summary.get("artifact_sha256", {})
    if set(hashes) != names:
        raise ValueError("GSEA artifact checksums are missing or do not cover every required output")
    for name in sorted(names):
        path = output_dir / name
        if not path.is_file() or digest(path) != hashes[name]:
            raise ValueError(f"GSEA output is missing or changed: {path}")
    if hashes[summary["mapping"]["file"]] != summary["mapping"].get("sha256"):
        raise ValueError("GSEA mapping checksum disagrees with its provenance")
    for source in summary["sources"].values():
        for path_field, hash_field in (("gmt_path", "gmt_sha256"),
                                       ("provenance_path", "provenance_sha256")):
            if path_field not in source:
                continue
            path = Path(source[path_field])
            path = path if path.is_absolute() else REPOSITORY / path
            if not path.is_file() or digest(path) != source.get(hash_field):
                raise ValueError(f"GSEA annotation source is missing or changed: {path}")
    checked = set()
    contrast_hashes = {}
    for analysis in summary["analyses"]:
        recorded = Path(analysis["source_contrast_file"])
        path = Path(data_dir) / recorded.name if data_dir is not None else (
            recorded if recorded.is_absolute() else REPOSITORY / recorded)
        if path in contrast_hashes and contrast_hashes[path] != analysis.get("source_contrast_sha256"):
            raise ValueError("GSEA libraries disagree about their contrast input checksum")
        contrast_hashes[path] = analysis.get("source_contrast_sha256")
        if path in checked:
            continue
        if not path.is_file() or digest(path) != analysis.get("source_contrast_sha256"):
            raise ValueError(f"GSEA contrast input is missing or changed: {path}")
        checked.add(path)
    if manifest_path.read_bytes() != manifest_payload:
        raise ValueError("GSEA completion manifest changed during validation")
    return summary


def _run(args: argparse.Namespace, output: Path) -> dict:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/lipohypolep-gsea-mpl")
    import gseapy

    started = time.monotonic()
    libraries = {
        "kegg": (args.kegg_gmt, args.kegg_provenance),
        "go_bp": (args.go_bp_gmt, args.go_bp_provenance),
    }
    if any(provenance is None for _, provenance in libraries.values()):
        raise ValueError("Both mouse GMT collections require source provenance with SHA256 and taxonomy 10090")
    snapshots = {path.resolve(): path.read_bytes()
                 for paths in libraries.values() for path in paths if path is not None}
    memberships = {key: read_gmt(paths[0], snapshots[paths[0].resolve()])
                   for key, paths in libraries.items()}
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
            float_precision="round_trip",
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
        table["stat"] = pd.to_numeric(table.stat, errors="raise")
        if table.empty or table.gene_id.duplicated().any() or not np.isfinite(table.stat).all():
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
            **source_metadata(path, provenance, snapshots[path.resolve()],
                              snapshots[provenance.resolve()] if provenance is not None else None),
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
                # The explicit stable stat/ID ordering above also resolves ties.
                rnk=rnk, gene_sets=sets, outdir=None, ascending=None,
                organism="mouse", pheno_pos=numerator, pheno_neg=denominator,
                min_size=args.min_size, max_size=args.max_size,
                permutation_num=args.permutations, weight=1.0,
                seed=args.seed, threads=args.threads, no_plot=True, verbose=False,
                method="permutation",
            )
            expected_index = pd.Index(rnk.mouse_symbol, name=pre.ranking.index.name)
            if (not pre.ranking.index.equals(expected_index)
                    or not np.array_equal(pre.ranking.to_numpy(), rnk.stat.to_numpy())):
                raise RuntimeError("GSEA changed the complete deterministic ranking")
            results = normalize_result(pre, sets, ranked_symbols)
            expected = summary["sources"][key]["eligible_gene_sets"]
            expected_terms = {term for term, genes in sets.items()
                              if args.min_size <= len(set(genes) & ranked_symbols) <= args.max_size}
            if len(results) != expected or set(results.term) != expected_terms:
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
            write_manifest(summary, output / "gsea_provenance.json")
            print(f"  {len(results)} tested sets; {int((results.gsea_fdr_q < .05).sum())} at GSEA FDR q < .05; {elapsed:.1f}s", flush=True)
    if any(digest(path) != input_sha256 for _, _, _, path, input_sha256, _ in contrasts):
        raise RuntimeError("Contrast inputs changed during GSEA; rerun all analyses on one consistent data version")
    if any(not path.is_file() or digest(path) != hashlib.sha256(payload).hexdigest()
           for path, payload in snapshots.items()):
        raise RuntimeError("Annotation sources changed during GSEA; rerun with consistent snapshots")
    summary["runtime_seconds"] = round(time.monotonic() - started, 3)
    return summary


def run(args: argparse.Namespace) -> dict:
    """Stage all artifacts and publish a complete manifest only after success."""
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # flock is released by the OS even after an interrupted process. Keeping its
    # inode in place avoids races between concurrent writers to this directory.
    with (output / ".gsea.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f"Another GSEA run is writing {output}") from error
        started = datetime.now(timezone.utc).isoformat()
        runner_sha256 = digest(Path(__file__))
        write_manifest({"complete": False, "status": "running_or_failed",
                        "started_utc": started}, output / "gsea_provenance.json")
        if args.permutations < 1 or args.threads < 1 or not 1 <= args.min_size <= args.max_size:
            raise ValueError("Require positive permutations/threads and 1 <= min-size <= max-size")
        with tempfile.TemporaryDirectory(prefix=".gsea-stage-", dir=output) as stage:
            staging = Path(stage)
            summary = _run(args, staging)
            if digest(Path(__file__)) != runner_sha256:
                raise RuntimeError("GSEA runner source changed during analysis; rerun with a consistent code version")
            names = artifact_names(summary)
            summary["artifact_sha256"] = {name: digest(staging / name) for name in sorted(names)}
            summary["started_utc"] = started
            summary["completed_utc"] = datetime.now(timezone.utc).isoformat()
            summary["runner_sha256"] = runner_sha256
            summary["integrity_schema"] = 1
            # If publication fails part way through, the completion marker remains
            # false, so mixed old/new files cannot be accepted by the renderer.
            for name in sorted(names):
                (staging / name).replace(output / name)
            summary["complete"] = True
            write_manifest(summary, output / "gsea_provenance.json")
        print(f"Saved {output / 'gsea_provenance.json'}; total {summary['runtime_seconds']:.1f}s", flush=True)
        return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kegg-gmt", type=Path, required=True, help="Verified mouse KEGG GMT")
    parser.add_argument("--go-bp-gmt", type=Path, required=True, help="Verified mouse GO biological-process GMT")
    parser.add_argument("--kegg-provenance", type=Path, required=True)
    parser.add_argument("--go-bp-provenance", type=Path, required=True)
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
