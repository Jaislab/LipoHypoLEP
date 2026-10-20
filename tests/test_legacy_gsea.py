"""Real-identifier, mapping, provenance and failure-safety legacy GSEA tests."""

import argparse
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from analysis import legacy_gsea as gsea
from analysis import audit_legacy_gsea as independent_audit
from analysis.audit_legacy_gsea import audit_scores, compare_metrics, reference_walk, validate_audit_report
from analysis.legacy_workbook import FIGURES, STATUS


class LegacyGseaTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.data = self.root / "data"
        self.data.mkdir()
        self.output = self.root / "gsea"
        self.calls = []
        self.callback = None
        self.tables = {}
        symbols = ["DUP", "DUP", "A", "B", "C", "COX1", "UNANNOTATED", "D"]
        pvalues = np.array([.5, .25, .01, .02, .03, .04, .9, .05])
        effects = np.array([.1, -.1, .02, -.02, .3, -.3, .001, -.001])
        for offset, (figure, *_) in enumerate(FIGURES):
            count = 8 - offset
            table = pd.DataFrame({"gene_id": [f"ENSMUSG{i:011d}" for i in range(1, count + 1)],
                "supplied_symbol": symbols[:count], "entrez_id": [str(i) for i in range(1, count + 1)],
                "log2FoldChange": effects[:count], "pvalue": pvalues[:count],
                "primaryproxy_signed_log10_p": np.sign(effects[:count]) * -np.log10(pvalues[:count]),
                "sensitivityproxy_signed_z": np.sign(effects[:count]) * norm.isf(pvalues[:count] / 2)})
            self.tables[figure] = table
            gsea.write_table(table, self.data / f"{figure}_deg.tsv.gz")
        self.publish_import()
        self.kegg = self.root / "kegg.gmt"
        self.go = self.root / "go.gmt"
        self.kegg.write_text("path_one\tna\tA\tDUP\tCOX1\npath_two\tna\tB\tC\tD\n")
        # The legacy mitochondrial symbol is intentionally absent from GO; no
        # alias mapping may turn COX1 into either mt-Co1 or the unrelated Ptgs1.
        self.go.write_text("go_one\tna\tA\tDUP\tmt-Co1\tPtgs1\ngo_two\tna\tB\tC\tD\n")
        self.kegg_provenance = self.root / "kegg.json"
        self.go_provenance = self.root / "go.json"
        for gmt, provenance in ((self.kegg, self.kegg_provenance), (self.go, self.go_provenance)):
            provenance.write_text(json.dumps({"sha256": gsea.digest(gmt), "ncbi_taxonomy_id": 10090}))
        self.args = argparse.Namespace(data_dir=self.data, kegg_gmt=self.kegg, go_bp_gmt=self.go,
            kegg_provenance=self.kegg_provenance, go_bp_provenance=self.go_provenance,
            permutations=10, threads=1, seed=20261020, min_size=1, max_size=7)

    def tearDown(self):
        self.temporary.cleanup()

    def publish_import(self, confirmed=True):
        manifest = {"complete": True, "data_status": STATUS, "input_file": "fixture.xlsx",
            "input_sha256": "a" * 64, "tested_universe_confirmation": {"all_tested_genes_exported": confirmed},
            "data_artifact_sha256": {f"{figure}_deg.tsv.gz": gsea.digest(self.data / f"{figure}_deg.tsv.gz")
                                     for figure, *_ in FIGURES},
            "figures": [{"figure_id": figure, "data_file": f"{figure}_deg.tsv.gz", "source_sheet": sheet,
                         "numerator": numerator, "denominator": denominator, "fold_change_multiplier": multiplier}
                        for figure, sheet, numerator, denominator, multiplier in FIGURES]}
        (self.data / "legacy_provenance.json").write_text(json.dumps(manifest))

    def fake_prerank(self, **kwargs):
        self.calls.append(kwargs)
        if self.callback is not None:
            self.callback()
        rank = kwargs["rnk"].copy()
        ids = set(rank.gene_id)
        rows = []
        for index, (term, genes) in enumerate(kwargs["gene_sets"].items()):
            matched = sorted(set(genes) & ids)
            direction = 1 if index % 2 == 0 else -1
            rows.append({"Term": term, "ES": direction * .4, "NES": direction * 1.5,
                         "NOM p-val": .01, "FDR q-val": .125, "FWER p-val": .3,
                         "Lead_genes": matched[0]})
        return types.SimpleNamespace(res2d=pd.DataFrame(rows),
            ranking=pd.Series(rank.stat.to_numpy(), index=rank.gene_id.to_numpy()))

    def run_fixture(self, metric="proxy_signed_log10_p", output=None):
        with patch.dict(sys.modules, {"gseapy": types.SimpleNamespace(prerank=self.fake_prerank)}), \
             patch.object(gsea.importlib.metadata, "version", return_value="test-engine"), \
             contextlib.redirect_stdout(io.StringIO()):
            return gsea.run_metric(self.args, metric, output or self.output)

    def test_exact_mapping_excludes_every_ambiguous_membership_keeps_ranked_rows(self):
        annotation, mapping = gsea.exact_mapping(self.tables)
        self.assertNotIn("DUP", mapping)
        self.assertEqual(len(annotation), 8)
        self.assertEqual((annotation.mapping_status == "ambiguous_symbol_membership_excluded").sum(), 2)
        translated = gsea.translate_library({"term": ["DUP", "A", "UNKNOWN"]}, mapping)
        self.assertEqual(translated, {"term": ["ENSMUSG00000000003"]})
        result = self.run_fixture()
        self.assertEqual([len(call["rnk"]) for call in self.calls], [8, 8, 7, 7, 6, 6])
        for call in self.calls:
            self.assertIn("ENSMUSG00000000001", set(call["rnk"].gene_id))
            self.assertIn("ENSMUSG00000000002", set(call["rnk"].gene_id))
            self.assertFalse(any({"ENSMUSG00000000001", "ENSMUSG00000000002"} & set(genes)
                                 for genes in call["gene_sets"].values()))
        self.assertEqual(len(result["artifact_sha256"]), 31)

    def test_inconsistent_gene_annotations_are_excluded_without_picking_a_symbol(self):
        tables = {key: value.copy() for key, value in self.tables.items()}
        tables["figure_02"].loc[2, "supplied_symbol"] = "OTHER"
        annotation, mapping = gsea.exact_mapping(tables)
        self.assertNotIn("A", mapping)
        self.assertNotIn("OTHER", mapping)
        rows = annotation.loc[annotation.gene_id.eq("ENSMUSG00000000003")]
        self.assertTrue(rows.mapping_status.eq("inconsistent_gene_annotation_membership_excluded").all())

    def test_no_alias_substitution_for_legacy_mitochondrial_symbols(self):
        self.run_fixture()
        mitochondrial = "ENSMUSG00000000006"
        self.assertTrue(any(mitochondrial in genes for genes in self.calls[0]["gene_sets"].values()))
        self.assertFalse(any(mitochondrial in genes for genes in self.calls[1]["gene_sets"].values()))
        gaps = pd.read_csv(self.output / "annotation_mapping_gaps.tsv", sep="\t")
        gap = gaps.loc[gaps.library.eq("go_bp") & gaps.supplied_symbol.eq("COX1")]
        self.assertEqual(len(gap), 1)
        self.assertTrue(gap.legacy_mitochondrial_symbol.iloc[0])

    def test_complete_runs_keep_no_deg_filter_and_real_identifiers_in_both_metrics(self):
        for metric in gsea.METRICS:
            self.calls.clear()
            out = self.output / metric
            summary = self.run_fixture(metric, out)
            gsea.validate_completed_run(out, self.data)
            self.assertTrue(summary["all_tested_genes_exported_user_confirmed"])
            for call, (figure, *_expected) in zip(self.calls[::2], FIGURES):
                expected = self.tables[figure].sort_values([gsea.METRICS[metric], "gene_id"], ascending=[False, True], kind="stable")
                np.testing.assert_array_equal(call["rnk"].stat.to_numpy(), expected[gsea.METRICS[metric]].to_numpy())
                self.assertEqual(set(call["rnk"].gene_id), set(self.tables[figure].gene_id))
                self.assertIsNone(call["ascending"])
            display = pd.read_csv(out / "figure_01_kegg_display.tsv", sep="\t")
            self.assertTrue(display.gsea_fdr_q.eq(.125).all())
            self.assertGreater(len(display), 0)

    def test_absent_tested_universe_confirmation_and_changed_proxy_are_rejected(self):
        self.publish_import(confirmed=False)
        with self.assertRaisesRegex(ValueError, "all tested genes"):
            gsea.load_tables(self.data)
        changed = self.tables["figure_01"].copy()
        changed.loc[0, "primaryproxy_signed_log10_p"] += .1
        gsea.write_table(changed, self.data / "figure_01_deg.tsv.gz")
        self.publish_import()
        with self.assertRaisesRegex(ValueError, "proxy.*changed"):
            gsea.load_tables(self.data)

    def test_failure_invalidates_prior_completion_and_does_not_mix_results(self):
        self.run_fixture()
        path = self.output / "figure_01_kegg_all_results.tsv.gz"
        original = path.read_bytes()
        self.callback = lambda: (_ for _ in ()).throw(RuntimeError("fit failed"))
        with self.assertRaisesRegex(RuntimeError, "fit failed"):
            self.run_fixture()
        self.assertFalse(json.loads((self.output / gsea.MANIFEST).read_text())["complete"])
        self.assertEqual(path.read_bytes(), original)
        self.assertFalse(list(self.output.glob(".legacy-gsea-stage-*")))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            gsea.validate_completed_run(self.output, self.data)

    def test_modified_results_input_tables_and_sources_are_rejected(self):
        self.run_fixture()
        for path in (self.output / "figure_01_kegg_display.tsv", self.output / "exact_symbol_mapping.tsv",
                     self.data / "figure_01_deg.tsv.gz", self.go):
            original = path.read_bytes()
            try:
                path.write_bytes(original + b"\n")
                with self.subTest(path=path.name), self.assertRaises(ValueError):
                    gsea.validate_completed_run(self.output, self.data)
            finally:
                path.write_bytes(original)

    def test_primary_and_sensitivity_completion_are_independent(self):
        self.run_fixture()
        primary_bytes = (self.output / gsea.MANIFEST).read_bytes()
        self.callback = lambda: (_ for _ in ()).throw(RuntimeError("sensitivity failed"))
        with self.assertRaisesRegex(RuntimeError, "sensitivity failed"):
            self.run_fixture("proxy_signed_z", self.output / "signed_z")
        self.assertEqual((self.output / gsea.MANIFEST).read_bytes(), primary_bytes)
        gsea.validate_completed_run(self.output, self.data)

    def test_missing_library_wrong_direction_and_changed_helper_are_rejected(self):
        original = self.run_fixture()
        path = self.output / gsea.MANIFEST
        changes = (
            lambda value: value["sources"].pop("go_bp"),
            lambda value: value["analyses"][0].__setitem__("numerator", "NON-LD control"),
            lambda value: value.__setitem__("pathway_helper_sha256", "0" * 64),
        )
        for change in changes:
            value = json.loads(json.dumps(original))
            change(value)
            gsea.write_manifest(value, path)
            with self.assertRaises(ValueError):
                gsea.validate_completed_run(self.output, self.data)
        gsea.write_manifest(original, path)

    @unittest.skipUnless(importlib.util.find_spec("gseapy"), "Install the pinned GSEA dependency for integration coverage")
    def test_actual_engine_accepts_real_ensembl_background_and_translated_memberships(self):
        with contextlib.redirect_stdout(io.StringIO()):
            summary = gsea.run_metric(self.args, "proxy_signed_log10_p", self.output)
        self.assertEqual(len(summary["analyses"]), 6)
        gsea.validate_completed_run(self.output, self.data)
        ranks = pd.read_csv(self.output / "figure_01_complete_rank.tsv.gz", sep="\t")
        self.assertEqual(set(ranks.gene_id), set(self.tables["figure_01"].gene_id))
        independent = audit_scores(self.output, self.data)
        self.assertEqual(independent["audited_sets"], 12)
        self.assertTrue(all(item["leading_edge_mismatches"] == 0 for item in independent["analyses"]))

    def test_independent_boundary_reference_agrees_with_full_gene_running_sum(self):
        generator = np.random.default_rng(20261020)
        genes = np.array([f"gene{i:03d}" for i in range(200)])
        statistics = np.sort(generator.normal(size=200))[::-1]
        positions = {gene: index for index, gene in enumerate(genes)}
        for _ in range(20):
            hit_indices = np.sort(generator.choice(200, size=15, replace=False))
            members = set(genes[hit_indices])
            increments = np.full(200, -1 / 185)
            increments[hit_indices] = np.abs(statistics[hit_indices]) / np.abs(statistics[hit_indices]).sum()
            running = np.cumsum(increments)
            peak, trough = int(running.argmax()), int(running.argmin())
            expected = running[peak] if running[peak] > -running[trough] else running[trough]
            edge = set(genes[hit_indices[hit_indices <= peak]]) if expected > 0 else set(genes[hit_indices[hit_indices > trough]])
            actual, actual_edge = reference_walk(genes, statistics, members, positions)
            self.assertAlmostEqual(actual, expected, places=12)
            self.assertEqual(actual_edge, edge)

    @unittest.skipUnless(importlib.util.find_spec("gseapy"), "Install the pinned GSEA dependency for integration coverage")
    def test_two_actual_rankings_have_an_independent_freshness_bound_sensitivity_audit(self):
        primary = self.data / "gsea"
        sensitivity = primary / "signed_z"
        with contextlib.redirect_stdout(io.StringIO()):
            gsea.run_metric(self.args, "proxy_signed_log10_p", primary)
            gsea.run_metric(self.args, "proxy_signed_z", sensitivity)
        comparison = compare_metrics(primary, sensitivity, self.data)
        self.assertEqual(len(comparison["comparisons"]), 6)
        self.assertTrue(all(row["identical_gene_order"] for row in comparison["comparisons"]))
        self.assertTrue(all(0 <= row["significant_set_jaccard"] <= 1 for row in comparison["comparisons"]))
        report = {"complete": True, "auditor_sha256": gsea.digest(Path(independent_audit.__file__)),
                  "score_audits": [audit_scores(primary, self.data), audit_scores(sensitivity, self.data)],
                  "ranking_sensitivity": comparison}
        path = self.data / "legacy_gsea_audit.json"
        gsea.write_manifest(report, path)
        validate_audit_report(self.data)
        report["score_audits"][0]["manifest_sha256"] = "0" * 64
        gsea.write_manifest(report, path)
        with self.assertRaisesRegex(ValueError, "stale"):
            validate_audit_report(self.data)


if __name__ == "__main__":
    unittest.main()
