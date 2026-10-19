#!/usr/bin/env python3
"""Read-only structural preflight for mouse RNA-seq CSV handoffs.

This module validates inputs; it never fits a model, normalizes expression,
changes identifiers, changes statistics, writes files, or runs enrichment.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
import math
from pathlib import Path
import re
from typing import Iterator


GROUPS = ("LD SALINE", "LD LEPTIN", "NON-LD control")
KINDS = ("counts", "expression", "deg")
ID_TYPES = ("ensembl_mouse", "symbol", "entrez_mouse")
DE_REQUIRED = {"gene_id", "log2FoldChange", "stat", "padj"}
DE_OPTIONAL = {"pvalue", "baseMean", "lfcSE"}
MISSING = {"", "NA"}
NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z")


def issue(report: dict, level: str, code: str, message: str, **context) -> None:
    """Bound invalid-file reports without silently hiding their error count."""
    report[f"{level}_count"] += 1
    if len(report[level]) < 100:
        report[level].append({"code": code, "message": message, **context})


@dataclass
class Table:
    role: str
    header: list[str]
    reader: object
    report: dict

    def rows(self) -> Iterator[tuple[int, list[str]]]:
        try:
            for row in self.reader:
                line = self.reader.line_num
                if len(row) != len(self.header):
                    issue(self.report, "errors", "row_width", "CSV row has the wrong number of fields.",
                          file=self.role, line=line, expected=len(self.header), observed=len(row))
                    continue
                yield line, row
        except csv.Error as exc:
            issue(self.report, "errors", "malformed_csv", str(exc), file=self.role,
                  line=self.reader.line_num)


def identifier(value: str) -> bool:
    return bool(value) and value == value.strip() and not any(ord(c) < 32 or ord(c) == 127 for c in value)


def read_table(path: str | Path, role: str, report: dict) -> Table | None:
    path = Path(path)
    entry = {"path": str(path.absolute()), "sha256": None, "size_bytes": None}
    report["files"][role] = entry
    try:
        if not path.is_file():
            raise OSError("Input must be an existing regular file.")
        raw = path.read_bytes()
    except OSError as exc:
        issue(report, "errors", "file_read", str(exc), file=role)
        return None
    # Hash exactly the bytes that are parsed, including a possible UTF-8 BOM.
    entry.update(sha256=hashlib.sha256(raw).hexdigest(), size_bytes=len(raw))
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError:
        issue(report, "errors", "encoding", "Expected UTF-8 CSV, optionally with a UTF-8 BOM.", file=role)
        return None
    if "\x00" in text:
        issue(report, "errors", "nul_character", "CSV contains a NUL character.", file=role)
        return None
    reader = csv.reader(io.StringIO(text, newline=""), strict=True)
    try:
        header = next(reader)
    except StopIteration:
        issue(report, "errors", "empty_file", "CSV has no header.", file=role)
        return None
    except csv.Error as exc:
        issue(report, "errors", "malformed_csv", str(exc), file=role, line=reader.line_num)
        return None
    if not header or any(not identifier(name) for name in header):
        issue(report, "errors", "invalid_header", "Header names must be nonempty, without outer whitespace or control characters.", file=role)
        return None
    if len(set(header)) != len(header):
        issue(report, "errors", "duplicate_header", "Duplicate CSV column names are ambiguous.", file=role)
        return None
    return Table(role, header, reader, report)


def numeric(value: str, report: dict, *, line: int, column: str,
            file: str = "data", nonnegative: bool = False,
            probability: bool = False, integral: bool = False) -> Decimal | None:
    context = {"file": file, "line": line, "column": column}
    if value in MISSING:
        issue(report, "errors", "missing_numeric", "Required numeric value is missing.", **context)
        return None
    try:
        number = Decimal(value)
    except InvalidOperation:
        issue(report, "errors", "invalid_numeric", "Expected a decimal numeric value.", **context)
        return None
    if not number.is_finite():
        issue(report, "errors", "non_finite", "Numeric values must be finite.", **context)
        return None
    if not NUMBER.fullmatch(value):
        issue(report, "errors", "invalid_numeric", "Numeric cells cannot contain whitespace, grouping separators, or non-decimal syntax.", **context)
        return None
    floating = float(number)
    if not math.isfinite(floating):
        issue(report, "errors", "float_range", "Finite decimal value exceeds the downstream float64 range; no rescaling is performed.", **context)
        return None
    if number != 0 and floating == 0:
        issue(report, "errors", "numeric_underflow", "Nonzero decimal would underflow to zero in downstream float64 input; investigate upstream without silently replacing it.", **context)
        return None
    if (nonnegative or probability) and number < 0:
        issue(report, "errors", "numeric_range", "Expected a nonnegative value.", **context)
        return None
    if probability and number > 1:
        issue(report, "errors", "probability_range", "P values must lie between zero and one.", **context)
        return None
    if integral and number != number.to_integral_value():
        issue(report, "errors", "non_integer_count", "Raw counts must be integer-valued; no rounding is performed.", **context)
        return None
    return number


def gene_identifier(value: str, id_type: str, report: dict, line: int, seen: set[str]) -> None:
    if not identifier(value):
        issue(report, "errors", "invalid_gene_id", "Gene IDs must be nonempty, without outer whitespace or control characters.", line=line)
    elif id_type == "ensembl_mouse" and not re.fullmatch(r"ENSMUSG[0-9]{11}(?:\.[0-9]+)?", value):
        issue(report, "errors", "gene_id_type", "Expected a mouse Ensembl gene ID; version suffixes are retained.", line=line, gene_id=value)
    elif id_type == "entrez_mouse" and not re.fullmatch(r"[1-9][0-9]*", value):
        issue(report, "errors", "gene_id_type", "Expected an unpadded positive Entrez gene integer ID.", line=line, gene_id=value)
    if value in seen:
        issue(report, "errors", "duplicate_gene_id", "Duplicate gene IDs require upstream investigation; no aggregation is performed.", line=line, gene_id=value)
    seen.add(value)


def metadata_samples(path: str | Path, report: dict, required_groups: tuple[str, ...] = GROUPS) -> list[str]:
    table = read_table(path, "metadata", report)
    if table is None:
        return []
    if not {"sample_id", "group"}.issubset(table.header):
        issue(report, "errors", "metadata_columns", "Metadata requires sample_id and group columns.")
        return []
    positions = {name: table.header.index(name) for name in ("sample_id", "group")}
    samples, seen = [], set()
    counts = {group: 0 for group in GROUPS}
    for line, row in table.rows():
        sample, group = row[positions["sample_id"]], row[positions["group"]]
        if not identifier(sample):
            issue(report, "errors", "invalid_sample_id", "Sample IDs must be nonempty, without outer whitespace or control characters.", file="metadata", line=line)
        if sample in seen:
            issue(report, "errors", "duplicate_sample_id", "Metadata sample IDs must be unique.", file="metadata", line=line, sample_id=sample)
        seen.add(sample)
        samples.append(sample)
        if group not in GROUPS:
            issue(report, "errors", "invalid_group", "Group must exactly match one of the three study names.", file="metadata", line=line, observed=group)
        else:
            counts[group] += 1
    for group, count in counts.items():
        if count == 0 and group not in required_groups:
            continue
        if count < 2:
            issue(report, "errors", "insufficient_replicates", "At least two independent biological samples per group are required.", group=group, samples=count)
        elif count < 3:
            issue(report, "warnings", "few_replicates", "At least three independent biological replicates per group are recommended.", group=group, samples=count)
    issue(report, "warnings", "biological_independence_unverified", "Each metadata row must represent an independent biological mouse; independence and technical-replicate status cannot be verified from sample_id/group.")
    extras = sorted(set(table.header) - {"sample_id", "group"})
    if extras:
        issue(report, "warnings", "metadata_covariates_unvalidated", "Additional metadata fields are preserved but are not validated or fitted as covariates.", columns=extras)
    report["summary"]["samples"] = len(samples)
    report["summary"]["samples_per_group"] = counts
    return samples


def validate_inputs(path: str | Path, *, kind: str, gene_id_type: str | None = None,
                    metadata: str | Path | None = None, numerator: str | None = None,
                    denominator: str | None = None, complete_tested_universe: bool = False,
                    missing_stat_reason: str | None = None, sample_expression: str | None = None,
                    annotation_columns: list[str] | None = None, upstream_method: str | None = None,
                    statistic_description: str | None = None, annotation_release: str | None = None,
                    adjustment_method: str | None = None) -> dict:
    """Return a JSON-compatible report; input bytes and statistics are untouched."""
    report = {
        "schema_version": "1.0", "read_only": True, "kind": kind, "valid": False,
        "files": {}, "declarations": {"species": "Mus musculus", "gene_id_type": gene_id_type,
            "numerator": numerator, "denominator": denominator,
            "complete_tested_universe": complete_tested_universe,
            "missing_stat_reason": missing_stat_reason, "sample_expression": sample_expression,
            "annotation_columns": annotation_columns or [],
            "upstream_method": upstream_method, "statistic_description": statistic_description,
            "annotation_release": annotation_release, "adjustment_method": adjustment_method},
        "summary": {"data_rows": 0, "unique_gene_ids": 0},
        "exclusions": {"missing_effect": {"gene_ids": [], "count": 0},
                       "missing_padj": {"gene_ids": [], "count": 0},
                       "missing_stat": {"gene_ids": [], "count": 0, "reason": missing_stat_reason}},
        "readiness": {"deg_plot_input": False, "gsea_rank_input": False, "heatmap_expression_input": False},
        "errors": [], "warnings": [], "errors_count": 0, "warnings_count": 0,
    }
    if kind not in KINDS:
        issue(report, "errors", "input_kind", "Explicit kind must be counts, expression, or deg.")
    if gene_id_type not in ID_TYPES:
        issue(report, "errors", "gene_id_declaration", "Declare mouse gene IDs as ensembl_mouse, symbol, or entrez_mouse.")
    if missing_stat_reason is not None and not missing_stat_reason.strip():
        issue(report, "errors", "missing_stat_reason", "Missing-stat exclusions require a nonempty upstream reason.")
    if kind == "deg":
        if numerator not in GROUPS or denominator not in GROUPS or numerator == denominator:
            issue(report, "errors", "contrast_declaration", "DEG input requires distinct exact study-group numerator and denominator names.")
        if not complete_tested_universe:
            issue(report, "warnings", "universe_not_confirmed", "The complete unfiltered tested-gene universe has not been confirmed; this file is not eligible as GSEA ranking input.")
        if sample_expression not in (None, "fpkm"):
            issue(report, "errors", "expression_unit", "Inline sample expression must be explicitly declared as fpkm.")
        if bool(metadata) != bool(sample_expression):
            issue(report, "errors", "inline_expression_declaration", "Inline FPKM validation requires both --metadata and --sample-expression fpkm.")
    elif kind in ("counts", "expression") and metadata is None:
        issue(report, "errors", "metadata_required", "Matrix inputs require a metadata CSV.")
    if kind != "deg" and (numerator is not None or denominator is not None or sample_expression is not None or complete_tested_universe or missing_stat_reason is not None):
        issue(report, "errors", "inapplicable_option", "Contrast, universe, missing-stat, and inline-FPKM options apply only to DEG inputs.")
    if metadata is not None and Path(metadata).absolute() == Path(path).absolute():
        issue(report, "errors", "file_role_collision", "Data and sample metadata must be separate files.")

    table = read_table(path, "data", report)
    required_groups = tuple(g for g in (numerator, denominator) if g in GROUPS) if kind == "deg" else GROUPS
    samples = metadata_samples(metadata, report, required_groups) if metadata is not None else []
    if table is None:
        return report
    if table.header[0] != "gene_id":
        issue(report, "errors", "first_column", "The first data column must be exactly gene_id.")
        return report
    if kind == "deg" and not DE_REQUIRED.issubset(table.header):
        issue(report, "errors", "deg_columns", "DEG input requires gene_id, log2FoldChange, stat, and padj; pvalue is optional.", missing=sorted(DE_REQUIRED - set(table.header)))
        return report
    columns = {name: i for i, name in enumerate(table.header)}
    matrix_samples = table.header[1:] if kind in ("counts", "expression") else samples if sample_expression == "fpkm" else []
    if kind in ("counts", "expression"):
        if not matrix_samples:
            issue(report, "errors", "no_sample_columns", "A matrix needs sample columns after gene_id.")
        if set(matrix_samples) != set(samples):
            issue(report, "errors", "sample_mismatch", "Matrix headers and metadata sample IDs must match exactly.", missing_from_data=sorted(set(samples)-set(matrix_samples)), missing_from_metadata=sorted(set(matrix_samples)-set(samples)))
    elif sample_expression == "fpkm":
        annotations = set(annotation_columns or [])
        if set(samples) & (DE_REQUIRED | DE_OPTIONAL):
            issue(report, "errors", "sample_column_collision", "A metadata sample ID cannot also name a required or recognized DE-statistic column.", columns=sorted(set(samples) & (DE_REQUIRED | DE_OPTIONAL)))
        if annotations & (DE_REQUIRED | DE_OPTIONAL | set(samples)):
            issue(report, "errors", "annotation_collision", "An annotation column cannot also be a DE or sample-expression column.")
        missing = set(samples) - set(table.header)
        unexpected = set(table.header) - DE_REQUIRED - DE_OPTIONAL - set(samples) - annotations
        if missing or unexpected:
            issue(report, "errors", "sample_mismatch", "Inline FPKM sample columns must match metadata exactly; declare extra non-expression annotations explicitly.", missing_from_data=sorted(missing), undeclared_columns=sorted(unexpected))
        absent_annotations = annotations - set(table.header)
        if absent_annotations:
            issue(report, "errors", "annotation_missing", "Declared annotation columns are absent.", columns=sorted(absent_annotations))
        if not samples:
            issue(report, "errors", "no_sample_columns", "Inline FPKM declaration requires sample columns and valid sample metadata.")
    elif kind == "deg":
        extras = sorted(set(table.header) - DE_REQUIRED - DE_OPTIONAL)
        if extras:
            issue(report, "warnings", "extra_columns_unvalidated", "Additional columns are preserved but unvalidated. Declare inline FPKM and metadata to validate sample-expression columns.", columns=extras)
    seen: set[str] = set()
    ranked, plottable, zero_pvalues = 0, 0, 0
    for line, row in table.rows():
        report["summary"]["data_rows"] += 1
        gene = row[0]
        gene_identifier(gene, gene_id_type, report, line, seen)
        if kind == "deg":
            raw_stat, raw_padj = row[columns["stat"]], row[columns["padj"]]
            documented_untestable = raw_stat in MISSING and bool(missing_stat_reason)
            if row[columns["log2FoldChange"]] in MISSING and documented_untestable:
                report["exclusions"]["missing_effect"]["gene_ids"].append(gene)
                effect = None
            else:
                effect = numeric(row[columns["log2FoldChange"]], report, line=line, column="log2FoldChange")
            if raw_stat in MISSING:
                report["exclusions"]["missing_stat"]["gene_ids"].append(gene)
                if not missing_stat_reason:
                    issue(report, "errors", "missing_stat", "A missing signed statistic requires --missing-stat-reason; the row is never silently dropped.", line=line, gene_id=gene)
                statistic = None
            else:
                statistic = numeric(raw_stat, report, line=line, column="stat")
                ranked += statistic is not None
            if raw_padj in MISSING:
                report["exclusions"]["missing_padj"]["gene_ids"].append(gene)
            else:
                adjusted = numeric(raw_padj, report, line=line, column="padj", probability=True)
                plottable += adjusted is not None and effect is not None
                zero_pvalues += adjusted == 0
            if effect is not None and statistic is not None and ((effect < 0 < statistic) or (statistic < 0 < effect)):
                issue(report, "errors", "statistic_direction", "Signed statistic and log2 fold change disagree in direction; investigate upstream, never auto-flip or take an absolute value.", line=line, gene_id=gene)
            if "pvalue" in columns:
                raw_pvalue = row[columns["pvalue"]]
                if not (raw_pvalue in MISSING and raw_stat in MISSING and missing_stat_reason):
                    pvalue = numeric(raw_pvalue, report, line=line, column="pvalue", probability=True)
                    zero_pvalues += pvalue == 0
            for optional in ("baseMean", "lfcSE"):
                if optional in columns:
                    if not (row[columns[optional]] in MISSING and documented_untestable):
                        numeric(row[columns[optional]], report, line=line, column=optional, nonnegative=True)
        for sample in matrix_samples:
            if sample in columns:
                numeric(row[columns[sample]], report, line=line, column=sample,
                        nonnegative=True, integral=kind == "counts")
    report["summary"]["unique_gene_ids"] = len(seen)
    if not report["summary"]["data_rows"]:
        issue(report, "errors", "no_data_rows", "CSV must contain at least one well-formed data row.")
    for label, exclusion in report["exclusions"].items():
        exclusion["count"] = len(exclusion["gene_ids"])
        if exclusion["count"]:
            message = {
                "missing_padj": "Missing adjusted P values are excluded from DEG highlighting, not from GSEA when a finite signed statistic is present. Document upstream filtering/untestability reasons.",
                "missing_stat": "Explicitly documented missing-stat rows are ineligible for GSEA; no rows or input files were modified.",
                "missing_effect": "Documented untestable rows with missing effects cannot be plotted as DE genes; they are retained and reported, never silently dropped.",
            }[label]
            issue(report, "warnings", label, message, count=exclusion["count"])
    if zero_pvalues:
        issue(report, "warnings", "zero_pvalues", "Zero P/adjusted-P values are retained exactly; any later plotting convention must be disclosed rather than changing source statistics.", cells=zero_pvalues)
    issue(report, "warnings", "identity_unverified", "Mouse species, identifier identity/release, upstream statistical model, contrast truth, and biological independence cannot be established from numeric CSV contents. Declarations require upstream confirmation.")
    if kind == "deg":
        report["summary"].update(rows_with_finite_signed_stat=ranked, rows_with_effect_and_padj=plottable)
        if not sample_expression:
            issue(report, "warnings", "heatmap_expression_missing", "DE summary statistics alone do not supply sample-level heatmap values. Inline FPKM columns with sample metadata can supply display-only expression.")
        if not all((upstream_method, statistic_description, annotation_release, adjustment_method)):
            issue(report, "warnings", "upstream_provenance_incomplete", "Confirm the upstream model, signed-statistic definition, adjusted-P method/family, and mouse identifier annotation release before analysis.")
    elif kind == "counts":
        issue(report, "warnings", "raw_count_units_unverified", "--kind counts declares raw integer counts, not TPM/FPKM. Integer values alone cannot prove units; no differential-expression model is fitted here.")
    report["valid"] = report["errors_count"] == 0
    report["readiness"].update(
        deg_plot_input=report["valid"] and kind == "deg" and plottable > 0,
        gsea_rank_input=report["valid"] and kind == "deg" and complete_tested_universe and ranked >= 2,
        heatmap_expression_input=report["valid"] and (kind == "expression" or (kind == "deg" and sample_expression == "fpkm")),
    )
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Unmodified UTF-8 comma-separated data file")
    parser.add_argument("--kind", choices=KINDS, required=True)
    parser.add_argument("--gene-id-type", choices=ID_TYPES, required=True)
    parser.add_argument("--metadata", type=Path, help="CSV with sample_id and exact study group")
    parser.add_argument("--numerator", choices=GROUPS)
    parser.add_argument("--denominator", choices=GROUPS)
    parser.add_argument("--complete-tested-universe", action="store_true", help="Explicit upstream confirmation; cannot be inferred from the CSV")
    parser.add_argument("--missing-stat-reason", help="Document upstream reason before allowing empty/NA statistics")
    parser.add_argument("--sample-expression", choices=("fpkm",), help="Validate inline FPKM columns named by --metadata")
    parser.add_argument("--annotation-column", action="append", default=[], dest="annotation_columns", help="Explicit extra non-DE, non-FPKM annotation column; repeat as needed")
    parser.add_argument("--upstream-method")
    parser.add_argument("--statistic-description")
    parser.add_argument("--annotation-release")
    parser.add_argument("--adjustment-method")
    args = parser.parse_args(argv)
    report = validate_inputs(args.input, kind=args.kind, gene_id_type=args.gene_id_type,
        metadata=args.metadata, numerator=args.numerator, denominator=args.denominator,
        complete_tested_universe=args.complete_tested_universe,
        missing_stat_reason=args.missing_stat_reason, sample_expression=args.sample_expression,
        annotation_columns=args.annotation_columns, upstream_method=args.upstream_method,
        statistic_description=args.statistic_description, annotation_release=args.annotation_release,
        adjustment_method=args.adjustment_method)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
