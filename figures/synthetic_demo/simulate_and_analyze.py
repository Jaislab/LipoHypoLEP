#!/usr/bin/env python3
"""Generate an explicitly synthetic bulk RNA-seq plotting demonstration.

This is a reproducible illustration, not an experiment or a production RNA-seq
pipeline. Counts follow a negative-binomial model; inference uses an approximate
linear model on log2-normalized counts, without empirical-Bayes moderation or
precision weights. It must not be described as DESeq2, edgeR, or limma-voom.
All effects are simulation parameters and all gene identifiers are invented.
This generator never reads experimental count files. Its log-count OLS and
pseudocount are illustrative only; real differential-expression statistics must
come from a separate validated count-model analysis, not this module.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd
import scipy
from scipy import stats


SEED = 20261019
N_GENES = 8000
N_REPLICATES = 6
GROUPS = ("NON-LD control", "LD SALINE", "LD LEPTIN")
CONTRASTS = (
    ("LD LEPTIN", "LD SALINE", "contrast_LD_LEPTIN_vs_LD_SALINE.tsv.gz"),
    ("LD SALINE", "NON-LD control", "contrast_LD_SALINE_vs_NON_LD_control.tsv.gz"),
    ("LD LEPTIN", "NON-LD control", "contrast_LD_LEPTIN_vs_NON_LD_control.tsv.gz"),
)
CLASS_COUNTS = {
    "unchanged": 5200,
    "near_complete_reversal": 640,
    "partial_reversal": 800,
    "persistent": 800,
    "treatment_only": 320,
    "amplified": 240,
}


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Adjust one contrast's complete tested family, preserving row order."""
    pvalues = np.asarray(pvalues, dtype=float)
    if pvalues.ndim != 1 or not np.all(np.isfinite(pvalues)):
        raise ValueError("P values must be a finite one-dimensional array.")
    if np.any((pvalues < 0) | (pvalues > 1)):
        raise ValueError("P values must be between zero and one.")
    if not pvalues.size:
        return np.empty(0, dtype=float)
    order = np.argsort(pvalues, kind="stable")
    adjusted = pvalues[order] * len(pvalues) / np.arange(1, len(pvalues) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def median_ratio_size_factors(counts: np.ndarray) -> np.ndarray:
    """Estimate demo size factors, failing when median-of-ratios is undefined."""
    counts = np.asarray(counts)
    if counts.ndim != 2 or not all(counts.shape):
        raise ValueError("Counts must be a nonempty genes-by-samples matrix.")
    if not np.all(np.isfinite(counts)) or np.any(counts < 0):
        raise ValueError("Counts must be finite and nonnegative.")
    if np.any(counts != np.floor(counts)):
        raise ValueError("The synthetic count matrix must contain integer counts.")
    positive = counts[np.all(counts > 0, axis=1)].astype(float)
    if not len(positive):
        raise ValueError("Median-of-ratios requires a gene positive in every sample.")
    geometric_means = np.exp(np.mean(np.log(positive), axis=1))
    factors = np.median(positive / geometric_means[:, None], axis=0)
    if not np.all(np.isfinite(factors)) or np.any(factors <= 0):
        raise ValueError("Median-of-ratios produced invalid size factors.")
    factors /= np.exp(np.mean(np.log(factors)))
    return factors


def fit_group_model(
    expression: np.ndarray, sample_groups: np.ndarray, groups: tuple[str, ...]
) -> tuple[dict[str, np.ndarray], np.ndarray, int]:
    """Fit the illustrative group-only OLS model, rejecting undefined tests."""
    expression = np.asarray(expression, dtype=float)
    sample_groups = np.asarray(sample_groups)
    if expression.ndim != 2 or not all(expression.shape):
        raise ValueError("Expression must be a nonempty genes-by-samples matrix.")
    if not np.all(np.isfinite(expression)):
        raise ValueError("Expression values must be finite.")
    if sample_groups.ndim != 1 or len(sample_groups) != expression.shape[1]:
        raise ValueError("Sample groups must match the expression columns.")
    if not groups or len(set(groups)) != len(groups):
        raise ValueError("Model groups must be nonempty and unique.")
    if set(sample_groups) != set(groups):
        raise ValueError("Each declared group must have samples, with no unknown groups.")
    residual_df = expression.shape[1] - len(groups)
    if residual_df <= 0:
        raise ValueError("The group model requires positive residual degrees of freedom.")
    means = {
        group: np.mean(expression[:, sample_groups == group], axis=1)
        for group in groups
    }
    fitted = np.column_stack([means[group] for group in sample_groups])
    residual_variance = np.sum((expression - fitted) ** 2, axis=1) / residual_df
    if not np.all(np.isfinite(residual_variance)) or np.any(residual_variance <= 0):
        raise ValueError("The illustrative t test requires positive residual variance for every gene.")
    return means, residual_variance, residual_df


def write_table(frame: pd.DataFrame, destination: Path) -> None:
    """Write reproducible tables; compressed files have a fixed gzip timestamp."""
    compression = {"method": "gzip", "mtime": 0} if destination.suffix == ".gz" else None
    frame.to_csv(
        destination,
        sep="\t",
        index=False,
        compression=compression,
        float_format="%.17g",
    )


def balanced_effects(rng: np.random.Generator, size: int) -> np.ndarray:
    signs = np.r_[np.ones(size // 2), -np.ones(size - size // 2)]
    rng.shuffle(signs)
    return signs * rng.uniform(0.8, 2.2, size)


def _generate_demo_into(output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    gene_ids = np.array([f"SIM{i:05d}" for i in range(1, N_GENES + 1)])
    classes = np.concatenate(
        [np.repeat(name, number) for name, number in CLASS_COUNTS.items()]
    )
    rng.shuffle(classes)

    # Baseline expression and biological overdispersion are independently drawn
    # for each invented gene; every class contains balanced effect directions.
    baseline = np.clip(rng.lognormal(np.log(80.0), 1.3, N_GENES), 2.0, 12000.0)
    dispersion = (0.04 + 0.25 / (1.0 + baseline / 100.0)) * rng.lognormal(
        0.0, 0.18, N_GENES
    )
    effect_saline = np.zeros(N_GENES)
    effect_leptin = np.zeros(N_GENES)
    for name in CLASS_COUNTS:
        if name == "unchanged":
            continue
        locations = np.flatnonzero(classes == name)
        effect = balanced_effects(rng, len(locations))
        if name == "treatment_only":
            effect_leptin[locations] = effect
            continue
        effect_saline[locations] = effect
        if name == "near_complete_reversal":
            effect_leptin[locations] = effect * rng.uniform(0.0, 0.2, len(locations))
        elif name == "partial_reversal":
            effect_leptin[locations] = effect * rng.uniform(0.35, 0.7, len(locations))
        elif name == "persistent":
            effect_leptin[locations] = effect
        elif name == "amplified":
            effect_leptin[locations] = effect * rng.uniform(1.25, 1.65, len(locations))

    sample_groups = np.repeat(np.array(GROUPS), N_REPLICATES)
    sample_ids = [
        f"{group} {replicate:02d}"
        for group in GROUPS
        for replicate in range(1, N_REPLICATES + 1)
    ]
    library_factors = rng.lognormal(0.0, 0.15, len(sample_ids))
    library_factors /= np.exp(np.mean(np.log(library_factors)))
    group_effects = np.column_stack(
        [np.zeros(N_GENES), effect_saline, effect_leptin]
    )
    sample_effects = np.repeat(group_effects, N_REPLICATES, axis=1)
    mean_counts = baseline[:, None] * library_factors[None, :] * np.exp2(sample_effects)
    inverse_dispersion = 1.0 / dispersion[:, None]
    probability = inverse_dispersion / (inverse_dispersion + mean_counts)
    counts = rng.negative_binomial(inverse_dispersion, probability)

    # Standard median-of-ratios normalization uses genes positive in every sample.
    # Counts themselves, rather than the generating factors, determine normalization.
    positive_in_every_sample = np.all(counts > 0, axis=1)
    estimated_factors = median_ratio_size_factors(counts)
    normalized = counts / estimated_factors[None, :]
    retained = np.sum(counts >= 10, axis=1) >= N_REPLICATES
    if not np.any(retained):
        raise ValueError("No synthetic genes passed the count filter.")
    tested_ids = gene_ids[retained]
    log_expression = np.log2(normalized[retained] + 0.5)
    base_mean = np.mean(normalized[retained], axis=1)

    metadata = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "group": sample_groups,
            "group_label": sample_groups,
            "replicate": np.tile(np.arange(1, N_REPLICATES + 1), len(GROUPS)),
            "simulated_library_size_factor": library_factors,
            "estimated_size_factor": estimated_factors,
        }
    )
    count_frame = pd.DataFrame(counts, columns=sample_ids)
    count_frame.insert(0, "gene_id", gene_ids)
    expression_frame = pd.DataFrame(log_expression, columns=sample_ids)
    expression_frame.insert(0, "gene_id", tested_ids)
    write_table(count_frame, output_dir / "counts.tsv.gz")
    write_table(metadata, output_dir / "sample_metadata.tsv")
    write_table(expression_frame, output_dir / "normalized_log2_expression.tsv.gz")

    # One group-only OLS model is used for all contrasts. Its residual degrees of
    # freedom are 18 - 3 = 15. Contrasts share model estimates and covariance.
    group_means, residual_variance, residual_df = fit_group_model(
        log_expression, sample_groups, GROUPS
    )
    contrast_se = np.sqrt(residual_variance * (2.0 / N_REPLICATES))
    estimates = {}
    contrast_summary = {}
    true_effects_by_group = {
        "NON-LD control": np.zeros(N_GENES),
        "LD SALINE": effect_saline,
        "LD LEPTIN": effect_leptin,
    }
    for numerator, denominator, filename in CONTRASTS:
        name = f"{numerator} vs {denominator}"
        estimate = group_means[numerator] - group_means[denominator]
        statistic = estimate / contrast_se
        pvalue = 2.0 * stats.t.sf(np.abs(statistic), df=residual_df)
        padj = benjamini_hochberg(pvalue)
        frame = pd.DataFrame(
            {
                "gene_id": tested_ids,
                "baseMean": base_mean,
                "log2FoldChange": estimate,
                "lfcSE": contrast_se,
                "stat": statistic,
                "pvalue": pvalue,
                "padj": padj,
            }
        )
        write_table(frame, output_dir / filename)
        estimates[name] = estimate
        true_difference = (
            true_effects_by_group[numerator] - true_effects_by_group[denominator]
        )[retained]
        true_null = np.abs(true_difference) < 1e-12
        selected = (padj < 0.05) & (np.abs(estimate) >= 0.2)
        contrast_summary[name] = {
            "file": filename,
            "tested_genes": int(len(tested_ids)),
            "adjusted_p_below_0_05": int(np.sum(padj < 0.05)),
            "adjusted_p_below_0_05_and_absolute_log2fc_at_least_0_2": int(np.sum(selected)),
            "higher_expression_at_joint_cutoffs": int(np.sum(selected & (estimate > 0))),
            "lower_expression_at_joint_cutoffs": int(np.sum(selected & (estimate < 0))),
            "true_null_genes": int(np.sum(true_null)),
            "true_null_unadjusted_p_below_0_05_fraction": float(
                np.mean(pvalue[true_null] < 0.05)
            ),
            "estimated_vs_true_effect_correlation": float(
                np.corrcoef(estimate, true_difference)[0, 1]
            ),
            "estimated_vs_true_effect_mean_absolute_error": float(
                np.mean(np.abs(estimate - true_difference))
            ),
            "pvalues_finite_and_bounded": bool(
                np.all(np.isfinite(pvalue)) & np.all((pvalue >= 0) & (pvalue <= 1))
            ),
            "adjusted_pvalues_finite_and_bounded": bool(
                np.all(np.isfinite(padj)) & np.all((padj >= 0) & (padj <= 1))
            ),
        }

    # Truth stays separate so invented mechanistic classes cannot accidentally be
    # presented as classifications inferred from differential-expression P values.
    truth = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "simulation_class": classes,
            "baseline_mean_NON_LD_control": baseline,
            "negative_binomial_dispersion": dispersion,
            "true_log2FoldChange_LD_SALINE_vs_NON_LD_control": effect_saline,
            "true_log2FoldChange_LD_LEPTIN_vs_NON_LD_control": effect_leptin,
            "true_log2FoldChange_LD_LEPTIN_vs_LD_SALINE": effect_leptin - effect_saline,
            "retained_after_count_filter": retained,
        }
    )
    write_table(truth, output_dir / "simulation_truth.tsv.gz")

    # PCA is unscaled after gene centering, using 1,500 most variable tested genes.
    gene_variance = np.var(log_expression, axis=1, ddof=1)
    pca_gene_indices = np.argsort(-gene_variance, kind="stable")[:1500]
    pca_expression = log_expression[pca_gene_indices]
    centered = (pca_expression - pca_expression.mean(axis=1, keepdims=True)).T
    left, singular_values, right = np.linalg.svd(centered, full_matrices=False)
    # Resolve arbitrary SVD signs with a deterministic loading-based convention.
    for component in range(len(singular_values)):
        largest_loading = np.argmax(np.abs(right[component]))
        if right[component, largest_loading] < 0:
            left[:, component] *= -1
            right[component] *= -1
    pca_scores = left * singular_values[None, :]
    explained_variance = singular_values**2 / (len(sample_ids) - 1)
    explained_ratio = explained_variance / np.sum(explained_variance)
    scores_frame = metadata[["sample_id", "group"]].copy()
    for component in range(3):
        scores_frame[f"PC{component + 1}"] = pca_scores[:, component]
    variance_frame = pd.DataFrame(
        {
            "component": [f"PC{i + 1}" for i in range(len(singular_values))],
            "explained_variance": explained_variance,
            "explained_variance_ratio": explained_ratio,
        }
    )
    write_table(scores_frame, output_dir / "pca_scores.tsv")
    write_table(variance_frame, output_dir / "pca_variance.tsv")
    write_table(pd.DataFrame({"gene_id": tested_ids[pca_gene_indices]}), output_dir / "pca_genes.tsv")

    contrast_identity_error = float(
        np.max(
            np.abs(
                estimates["LD LEPTIN vs NON-LD control"]
                - estimates["LD LEPTIN vs LD SALINE"]
                - estimates["LD SALINE vs NON-LD control"]
            )
        )
    )
    validations = {
        "counts_shape_is_8000_by_18": bool(counts.shape == (N_GENES, 18)),
        "counts_are_nonnegative_integers": bool(
            np.issubdtype(counts.dtype, np.integer) & np.all(counts >= 0)
        ),
        "six_replicates_per_group": bool(
            all(np.sum(sample_groups == group) == N_REPLICATES for group in GROUPS)
        ),
        "all_gene_ids_are_unique_and_synthetic": bool(
            len(np.unique(gene_ids)) == N_GENES
            and all(identifier.startswith("SIM") for identifier in gene_ids)
        ),
        "normalization_factors_are_finite_and_positive": bool(
            np.all(np.isfinite(estimated_factors)) & np.all(estimated_factors > 0)
        ),
        "residual_degrees_of_freedom_equal_15": bool(residual_df == 15),
        "contrast_identity_max_absolute_error": contrast_identity_error,
        "contrast_identity_within_1e_minus_12": bool(contrast_identity_error < 1e-12),
        "all_pvalues_finite_and_bounded": bool(
            all(item["pvalues_finite_and_bounded"] for item in contrast_summary.values())
        ),
        "all_adjusted_pvalues_finite_and_bounded": bool(
            all(item["adjusted_pvalues_finite_and_bounded"] for item in contrast_summary.values())
        ),
        "pca_variance_ratios_sum_to_one": bool(np.isclose(explained_ratio.sum(), 1.0)),
    }
    summary = {
        "data_status": "Entirely synthetic illustration; no experimental measurements or biological findings.",
        "seed": SEED,
        "software": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "genes_simulated": N_GENES,
        "genes_retained": int(np.sum(retained)),
        "genes_filtered": int(np.sum(~retained)),
        "samples": len(sample_ids),
        "replicates_per_group": {group: N_REPLICATES for group in GROUPS},
        "simulation": {
            "distribution": "Negative binomial; variance = mean + dispersion * mean^2.",
            "baseline": "Lognormal(log(80), 1.3), clipped to [2, 12000].",
            "dispersion": "(0.04 + 0.25 / (1 + baseline / 100)) * Lognormal(0, 0.18).",
            "direction": "Balanced positive and negative effects within each non-null truth class.",
            "truth_class_counts": CLASS_COUNTS,
            "truth_class_fractions": {name: number / N_GENES for name, number in CLASS_COUNTS.items()},
            "class_interpretation": "Generating parameters only; never classifications inferred from statistical significance.",
        },
        "analysis": {
            "normalization": "Median of ratios using genes positive in every sample; factors rescaled to geometric mean 1.",
            "genes_used_for_size_factor_estimation": int(np.sum(positive_in_every_sample)),
            "filter": "Raw count >= 10 in at least 6 samples, irrespective of group.",
            "transformation": "log2(count / estimated_size_factor + 0.5)",
            "model": "One pooled group-only OLS model on transformed expression, with common per-gene residual variance.",
            "residual_degrees_of_freedom": residual_df,
            "contrasts": "First group minus second group; all contrasts derived from the same fitted group means.",
            "tests": "Two-sided Student t tests; Benjamini-Hochberg adjustment independently per contrast across all retained genes.",
            "limitation": "Illustrative approximate log-scale inference, without empirical-Bayes moderation or precision weights; not DESeq2, edgeR, or limma-voom.",
            "contrast_covariance": "All three contrasts share model estimates. Per-gene contrast covariance follows the common design matrix; contrasts are not independent.",
            "plot_selection_example": "DEG display uses adjusted P < 0.05 and absolute log2 fold change >= 0.2. GSEA ranks all tested genes without either DEG filter. This does not define truth classes or establish normalization.",
        },
        "pca": {
            "gene_selection": "1,500 tested genes with highest sample variance of transformed expression.",
            "centering_and_scaling": "Genes centered; no variance scaling.",
            "genes_used": len(pca_gene_indices),
            "explained_variance_ratio": {f"PC{i + 1}": float(explained_ratio[i]) for i in range(3)},
        },
        "contrasts": contrast_summary,
        "validations": validations,
    }
    failed = [name for name, value in validations.items() if isinstance(value, bool) and not value]
    if failed:
        raise RuntimeError(f"Synthetic demonstration validation failed: {failed}")
    with (output_dir / "simulation_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")
    return summary


def generate_demo(output_dir: Path) -> dict:
    """Stage a validated demo before publication; refuse unmarked data folders.

    Existing synthetic auxiliary files are retained. A calculation or serialization
    failure leaves the previous output untouched; the two final directory renames
    are rolled back if publication fails. Concurrent writers are not supported.
    """
    output_dir = Path(output_dir).absolute()
    if output_dir.is_symlink():
        raise ValueError("The synthetic output directory must not be a symbolic link.")
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError("The synthetic output destination must be a directory.")
    if output_dir.exists() and any(output_dir.iterdir()):
        marker = output_dir / "simulation_summary.json"
        try:
            previous = json.loads(marker.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise ValueError(
                "Refusing to overwrite a nonempty directory without a synthetic demo marker."
            ) from error
        if not isinstance(previous, dict) or not str(previous.get("data_status", "")).startswith(
            "Entirely synthetic illustration;"
        ):
            raise ValueError("Refusing to overwrite a directory not marked as synthetic.")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=f".{output_dir.name}-staging-", dir=output_dir.parent
    ) as temporary:
        staging_root = Path(temporary)
        staged = staging_root / "data"
        if output_dir.exists():
            shutil.copytree(output_dir, staged, symlinks=True)
            # Never follow an existing output-file symlink while regenerating.
            for path in staged.iterdir():
                if path.is_symlink():
                    raise ValueError("Synthetic output directories must not contain symbolic links.")
        else:
            staged.mkdir()
        summary = _generate_demo_into(staged)
        backup = staging_root / "previous"
        if output_dir.exists():
            os.replace(output_dir, backup)
        try:
            os.replace(staged, output_dir)
        except OSError:
            if backup.exists():
                os.replace(backup, output_dir)
            raise
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
        help="Destination for reproducible synthetic tables.",
    )
    arguments = parser.parse_args()
    summary = generate_demo(arguments.output_dir)
    print(
        json.dumps(
            {
                "output_directory": str(arguments.output_dir.resolve()),
                "genes_retained": summary["genes_retained"],
                "contrasts": summary["contrasts"],
                "validations": summary["validations"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
