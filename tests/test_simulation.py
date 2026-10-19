"""Independent numerical and failure-path checks for the synthetic demo only."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "figures" / "synthetic_demo" / "simulate_and_analyze.py"
SPEC = importlib.util.spec_from_file_location("simulation_under_test", SOURCE)
simulation = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(simulation)
DATA = SOURCE.parent / "data"


class AdjustmentTests(unittest.TestCase):
    def test_bh_matches_independent_scipy_reference_with_ties_and_boundaries(self):
        pvalues = np.r_[0.0, 1.0, 0.01, 0.01, 0.2, 0.9,
                        np.random.default_rng(417).uniform(size=1000)]
        np.testing.assert_allclose(
            simulation.benjamini_hochberg(pvalues),
            stats.false_discovery_control(pvalues, method="bh"),
            rtol=2e-15, atol=0,
        )

    def test_empty_bh_family_is_supported(self):
        self.assertEqual(simulation.benjamini_hochberg(np.array([])).shape, (0,))

    def test_invalid_pvalues_are_rejected(self):
        for values in ([np.nan], [np.inf], [-0.01], [1.01], [[0.1, 0.2]]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                simulation.benjamini_hochberg(values)


class ModelGuardTests(unittest.TestCase):
    def test_median_ratio_factors_recover_known_multipliers(self):
        counts = np.array([[10, 20, 40], [20, 40, 80], [0, 10, 0]])
        np.testing.assert_allclose(
            simulation.median_ratio_size_factors(counts), [0.5, 1.0, 2.0],
            rtol=1e-14, atol=0,
        )

    def test_invalid_or_unidentifiable_normalization_is_rejected(self):
        bad = (
            np.empty((0, 3)), np.empty((3, 0)), [[1, 0], [0, 1]],
            [[1, -1]], [[1, np.nan]], [[1, np.inf]], [[1, 1.5]],
        )
        for counts in bad:
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                simulation.median_ratio_size_factors(counts)

    def test_group_model_matches_independent_matrix_ols(self):
        groups = np.array(["A", "A", "A", "B", "B", "C", "C", "C", "C"])
        expression = np.random.default_rng(9).normal(size=(25, len(groups)))
        means, variance, degrees = simulation.fit_group_model(expression, groups, ("A", "B", "C"))
        design = np.column_stack([groups == group for group in ("A", "B", "C")]).astype(float)
        coefficients, _, _, _ = np.linalg.lstsq(design, expression.T, rcond=None)
        residuals = expression.T - design @ coefficients
        np.testing.assert_allclose(np.stack(list(means.values())), coefficients, rtol=1e-13, atol=1e-14)
        np.testing.assert_allclose(variance, np.sum(residuals**2, axis=0) / 6, rtol=1e-13)
        self.assertEqual(degrees, 6)

    def test_undefined_group_models_fail_explicitly(self):
        cases = (
            (np.ones((1, 4)), ["A", "A", "B", "B"], ("A", "B")),
            (np.array([[1.0, np.nan, 2, 3]]), ["A", "A", "B", "B"], ("A", "B")),
            (np.array([[1.0, 2]]), ["A", "B"], ("A", "B")),
            (np.array([[1.0, 2, 3]]), ["A", "A"], ("A",)),
            (np.array([[1.0, 2, 3]]), ["A", "A", "A"], ("A", "B")),
            (np.array([[1.0, 2, 3]]), ["A", "A", "A"], ("A", "A")),
        )
        for expression, groups, model_groups in cases:
            with self.subTest(groups=groups), self.assertRaises(ValueError):
                simulation.fit_group_model(expression, np.array(groups), model_groups)


class PublicationSafetyTests(unittest.TestCase):
    def test_unmarked_nonempty_directory_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "experimental_data"
            destination.mkdir()
            original = destination / "counts.csv"
            original.write_text("gene,sample\nRealGene,23\n")
            with self.assertRaisesRegex(ValueError, "Refusing to overwrite"):
                simulation.generate_demo(destination)
            self.assertEqual(original.read_text(), "gene,sample\nRealGene,23\n")

    def test_failed_generation_keeps_entire_previous_dataset(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "demo"
            destination.mkdir()
            (destination / "simulation_summary.json").write_text(json.dumps({
                "data_status": "Entirely synthetic illustration; fixture."
            }))
            (destination / "counts.tsv.gz").write_bytes(b"previous-valid-data")
            before = {p.name: p.read_bytes() for p in destination.iterdir()}

            def fail_mid_generation(staged):
                (staged / "counts.tsv.gz").write_bytes(b"partial-invalid-data")
                raise ValueError("injected numerical or serialization failure")

            with mock.patch.object(simulation, "_generate_demo_into", side_effect=fail_mid_generation):
                with self.assertRaisesRegex(ValueError, "injected"):
                    simulation.generate_demo(destination)
            self.assertEqual({p.name: p.read_bytes() for p in destination.iterdir()}, before)

    def test_failed_publication_rolls_back_previous_dataset(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "demo"
            destination.mkdir()
            (destination / "simulation_summary.json").write_text(json.dumps({
                "data_status": "Entirely synthetic illustration; fixture."
            }))
            (destination / "sentinel.txt").write_text("previous")
            original_replace = simulation.os.replace

            def replace_with_injected_failure(source, target):
                if Path(source).name == "data":
                    raise OSError("injected publication failure")
                return original_replace(source, target)

            with mock.patch.object(simulation, "_generate_demo_into", return_value={}), \
                 mock.patch.object(simulation.os, "replace", side_effect=replace_with_injected_failure):
                with self.assertRaisesRegex(OSError, "injected"):
                    simulation.generate_demo(destination)
            self.assertEqual((destination / "sentinel.txt").read_text(), "previous")

    def test_successful_publication_preserves_auxiliary_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "demo"
            destination.mkdir()
            (destination / "simulation_summary.json").write_text(json.dumps({
                "data_status": "Entirely synthetic illustration; fixture."
            }))
            (destination / "selected_plot_genes.tsv").write_text("gene_id\nSIM00001\n")
            (destination / "counts.tsv.gz").write_bytes(b"previous")

            def publish_fixture(staged):
                (staged / "counts.tsv.gz").write_bytes(b"reproduced")
                return {"complete": True}

            with mock.patch.object(simulation, "_generate_demo_into", side_effect=publish_fixture):
                self.assertEqual(simulation.generate_demo(destination), {"complete": True})
            self.assertEqual((destination / "counts.tsv.gz").read_bytes(), b"reproduced")
            self.assertEqual((destination / "selected_plot_genes.tsv").read_text(), "gene_id\nSIM00001\n")

    def test_output_file_symlinks_cannot_redirect_writes(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "demo"
            destination.mkdir()
            (destination / "simulation_summary.json").write_text(json.dumps({
                "data_status": "Entirely synthetic illustration; fixture."
            }))
            original = Path(temporary) / "real_counts.csv"
            original.write_bytes(b"experimental")
            (destination / "counts.tsv.gz").symlink_to(original)
            with self.assertRaisesRegex(ValueError, "symbolic links"):
                simulation.generate_demo(destination)
            self.assertEqual(original.read_bytes(), b"experimental")

    def test_symlink_destination_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "target"
            target.mkdir()
            link = Path(temporary) / "link"
            link.symlink_to(target, target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "symbolic link"):
                simulation.generate_demo(link)


class ExistingDataNumericalAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.counts = pd.read_csv(DATA / "counts.tsv.gz", sep="\t", index_col="gene_id")
        cls.metadata = pd.read_csv(DATA / "sample_metadata.tsv", sep="\t", float_precision="round_trip")
        cls.expression = pd.read_csv(DATA / "normalized_log2_expression.tsv.gz", sep="\t", index_col="gene_id", float_precision="round_trip")

    def test_stored_normalization_and_log_expression_match_raw_counts(self):
        counts = self.counts.to_numpy(dtype=float)
        positive = counts[np.all(counts > 0, axis=1)]
        ratios = positive / np.exp(np.mean(np.log(positive), axis=1))[:, None]
        factors = np.median(ratios, axis=0)
        factors /= np.exp(np.log(factors).mean())
        np.testing.assert_allclose(factors, self.metadata["estimated_size_factor"], rtol=1e-14)
        retained = np.sum(counts >= 10, axis=1) >= 6
        self.assertTrue(self.counts.index[retained].equals(self.expression.index))
        np.testing.assert_allclose(
            self.expression, np.log2(counts[retained] / factors[None, :] + 0.5),
            rtol=1e-14, atol=0,
        )

    def test_all_contrasts_match_independent_design_matrix_inference(self):
        groups = self.metadata["group"].to_numpy()
        design = np.column_stack([groups == g for g in simulation.GROUPS]).astype(float)
        coefficients, _, rank, _ = np.linalg.lstsq(design, self.expression.to_numpy().T, rcond=None)
        residual_df = len(groups) - rank
        residuals = self.expression.to_numpy().T - design @ coefficients
        variances = (residuals**2).sum(axis=0) / residual_df
        covariance = np.linalg.inv(design.T @ design)
        estimates = {}
        for numerator, denominator, filename in simulation.CONTRASTS:
            contrast = np.array([int(g == numerator) - int(g == denominator) for g in simulation.GROUPS])
            estimate = contrast @ coefficients
            se = np.sqrt(variances * (contrast @ covariance @ contrast))
            statistic = estimate / se
            pvalue = 2 * stats.t.sf(np.abs(statistic), residual_df)
            padj = stats.false_discovery_control(pvalue, method="bh")
            stored = pd.read_csv(DATA / filename, sep="\t", index_col="gene_id", float_precision="round_trip")
            self.assertTrue(stored.index.equals(self.expression.index))
            for column, independent in (("log2FoldChange", estimate), ("lfcSE", se), ("stat", statistic), ("pvalue", pvalue), ("padj", padj)):
                with self.subTest(contrast=filename, column=column):
                    np.testing.assert_allclose(stored[column], independent, rtol=5e-12, atol=5e-14)
            estimates[(numerator, denominator)] = stored["log2FoldChange"].to_numpy()
        np.testing.assert_allclose(
            estimates[("LD LEPTIN", "NON-LD control")],
            estimates[("LD LEPTIN", "LD SALINE")] + estimates[("LD SALINE", "NON-LD control")],
            rtol=0, atol=1e-12,
        )

    def test_current_generator_reproduces_every_existing_numerical_output_byte(self):
        names = ["counts.tsv.gz", "sample_metadata.tsv", "normalized_log2_expression.tsv.gz",
                 "simulation_truth.tsv.gz", "pca_scores.tsv", "pca_variance.tsv", "pca_genes.tsv"
                 ] + [item[2] for item in simulation.CONTRASTS]
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "reproduced"
            simulation.generate_demo(destination)
            for name in names:
                with self.subTest(file=name):
                    self.assertEqual(
                        hashlib.sha256((destination / name).read_bytes()).hexdigest(),
                        hashlib.sha256((DATA / name).read_bytes()).hexdigest(),
                    )
            # JSON key insertion order has no numerical meaning.
            self.assertEqual(
                json.loads((destination / "simulation_summary.json").read_text()),
                json.loads((DATA / "simulation_summary.json").read_text()),
            )


if __name__ == "__main__":
    unittest.main()
