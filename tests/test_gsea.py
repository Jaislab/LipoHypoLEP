"""Numerical contracts and failure-safety tests for the synthetic GSEA runner."""

import argparse
import contextlib
import fcntl
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


SOURCE = Path(__file__).resolve().parents[1] / "figures/synthetic_demo/run_gsea.py"
SPEC = importlib.util.spec_from_file_location("audit_run_gsea", SOURCE)
gsea = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gsea)


class GseaTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.data = self.root / "data"
        self.output = self.root / "results"
        self.data.mkdir()
        self.calls = []
        self.callback = None
        self.corrupt_ranking = False
        for _, _, _, filename in gsea.CONTRASTS:
            # Ties, both signs, and a zero statistic must all remain in the rank.
            table = pd.DataFrame({"gene_id": [f"SIM{i:03d}" for i in range(6)],
                                  "stat": [2.0, 2.0, 0.0, -0.5, -1.0, -2.0]})
            gsea.write_table(table, self.data / filename)
        gmt = self.root / "mouse.gmt"
        gmt.write_text("set_a\tna\tA\tB\tC\tD\nset_b\tna\tE\tF\tG\tH\n")
        provenance = self.root / "mouse.json"
        provenance.write_text(json.dumps({"sha256": gsea.digest(gmt), "ncbi_taxonomy_id": 10090}))
        self.args = argparse.Namespace(data_dir=self.data, output_dir=self.output,
            kegg_gmt=gmt, go_bp_gmt=gmt, kegg_provenance=provenance,
            go_bp_provenance=provenance, seed=20261019, min_size=1, max_size=5,
            threads=1, permutations=10)

    def tearDown(self):
        self.temporary.cleanup()

    def fake_prerank(self, **kwargs):
        self.calls.append(kwargs)
        if self.callback is not None:
            self.callback()
        rank = kwargs["rnk"].copy()
        symbols = set(rank.iloc[:, 0])
        rows = []
        for index, (term, genes) in enumerate(kwargs["gene_sets"].items()):
            matched = sorted(set(genes) & symbols)
            if kwargs["min_size"] <= len(matched) <= kwargs["max_size"]:
                direction = 1 if index % 2 == 0 else -1
                rows.append({"Name": "prerank", "Term": term, "ES": direction * .4,
                    "NES": direction * 1.5, "NOM p-val": .01, "FDR q-val": .125,
                    "FWER p-val": .3, "Lead_genes": matched[0]})
        ranking = pd.Series(rank.iloc[:, 1].to_numpy(), index=rank.iloc[:, 0].to_numpy())
        if self.corrupt_ranking:
            ranking = ranking.iloc[::-1]
        return types.SimpleNamespace(res2d=pd.DataFrame(rows), ranking=ranking)

    def run_fixture(self):
        fake = types.SimpleNamespace(prerank=self.fake_prerank)
        with patch.dict(sys.modules, {"gseapy": fake}), \
             patch.object(gsea.importlib.metadata, "version", return_value="test-engine"), \
             contextlib.redirect_stdout(io.StringIO()):
            return gsea.run(self.args)

    def read_manifest(self):
        return json.loads((self.output / "gsea_provenance.json").read_text())

    def test_gmt_parser_deduplicates_membership_and_rejects_duplicate_terms(self):
        self.assertEqual(gsea.read_gmt(Path("snapshot"), b"term\tna\tB\tA\tA\n"), {"term": ["A", "B"]})
        for payload in (b"term\tna\tA\nterm\tna\tB\n", b"term\tna\t\n", b"\n"):
            with self.subTest(payload=payload), self.assertRaises(ValueError):
                gsea.read_gmt(Path("snapshot"), payload)

    def test_mouse_source_provenance_rejects_wrong_species_and_hash(self):
        for data in ({"ncbi_taxonomy_id": 9606, "sha256": gsea.digest(self.args.kegg_gmt)},
                     {"ncbi_taxonomy_id": 10090, "sha256": "bad"}):
            with self.subTest(data=data), self.assertRaises(ValueError):
                gsea.source_metadata(self.args.kegg_gmt, self.args.kegg_provenance,
                    provenance_payload=json.dumps(data).encode())

    def test_result_normalization_preserves_empirical_q_and_rejects_bad_members(self):
        result = types.SimpleNamespace(res2d=pd.DataFrame([{"Term": "term", "ES": .5,
            "NES": 1.4, "NOM p-val": .001, "FDR q-val": 0.0,
            "FWER p-val": .1, "Lead_genes": "A"}]))
        normalized = gsea.normalize_result(result, {"term": ["A", "B"]}, {"A", "C"})
        self.assertEqual(normalized.gsea_fdr_q.iloc[0], 0.0)
        self.assertEqual(normalized.ranked_set_size.iloc[0], 1)
        for column, value in (("Lead_genes", "B"), ("Lead_genes", "A;A"),
                              ("Term", "unknown"), ("NES", float("inf")),
                              ("FDR q-val", 1.1), ("ES", -.5)):
            bad = result.res2d.copy()
            bad.loc[0, column] = value
            with self.subTest(column=column, value=value), self.assertRaises(RuntimeError):
                gsea.normalize_result(types.SimpleNamespace(res2d=bad), {"term": ["A", "B"]}, {"A"})

    def test_display_selection_accepts_one_direction_few_sets_and_no_significance(self):
        rows = pd.DataFrame({"term": ["a", "b", "c"], "NES": [1., 3., 2.], "gsea_fdr_q": [1.] * 3})
        self.assertEqual(gsea.pick_display_rows(rows).term.tolist(), ["b", "c", "a"])
        self.assertEqual(len(gsea.pick_display_rows(rows, count=0)), 0)
        with self.assertRaises(ValueError):
            gsea.pick_display_rows(rows, count=-1)

    def test_complete_run_hashes_every_artifact_and_keeps_complete_stable_rank(self):
        result = self.run_fixture()
        self.assertTrue(result["complete"])
        self.assertEqual(len(result["artifact_sha256"]), 30)
        self.assertEqual(len(self.calls), 6)
        for call in self.calls:
            self.assertIsNone(call["ascending"])
            self.assertEqual(call["rnk"].stat.tolist(), [2., 2., 0., -.5, -1., -2.])
            self.assertEqual(len(call["rnk"]), 6)
        verified = gsea.validate_completed_run(self.output, self.data)
        self.assertEqual(verified["artifact_sha256"], result["artifact_sha256"])

    def test_binary64_statistics_survive_input_and_rank_output_without_precision_loss(self):
        # The default pandas parser changes the third/fourth values by one or
        # more final binary digits despite the sufficient 17-digit text export.
        expected = np.array([2.0000000000000004, 1.2345678901234567,
                             .12345678901234566, -.12345678901234566,
                             -1.2345678901234567, -2.0000000000000004])
        table = pd.DataFrame({"gene_id": [f"SIM{i:03d}" for i in range(6)], "stat": expected})
        for _, _, _, filename in gsea.CONTRASTS:
            gsea.write_table(table, self.data / filename)
        summary = self.run_fixture()
        for call in self.calls:
            np.testing.assert_array_equal(call["rnk"].stat.to_numpy(), expected)
        for name in {item["ranked_statistics_file"] for item in summary["analyses"]}:
            rank = pd.read_csv(self.output / name, sep="\t", float_precision="round_trip")
            np.testing.assert_array_equal(rank.stat.to_numpy(), expected)

    def test_failed_rerun_invalidates_old_completion_and_preserves_old_artifacts(self):
        result = self.run_fixture()
        display = self.output / result["analyses"][0]["display_file"]
        old_bytes = display.read_bytes()
        self.callback = lambda: (_ for _ in ()).throw(RuntimeError("engine failed"))
        with self.assertRaisesRegex(RuntimeError, "engine failed"):
            self.run_fixture()
        self.assertFalse(self.read_manifest()["complete"])
        self.assertEqual(display.read_bytes(), old_bytes)
        self.assertFalse(list(self.output.glob(".gsea-stage-*")))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            gsea.validate_completed_run(self.output, self.data)

    def test_preprocessing_failure_also_invalidates_old_completion(self):
        self.run_fixture()
        self.args.go_bp_provenance = None
        with self.assertRaisesRegex(ValueError, "require source provenance"):
            self.run_fixture()
        self.assertFalse(self.read_manifest()["complete"])

    def test_no_random_annotation_assignment_for_real_gene_identifiers(self):
        for _, _, _, filename in gsea.CONTRASTS:
            gsea.write_table(pd.DataFrame({"gene_id": ["ENSMUSG00000000001"], "stat": [1.]}), self.data / filename)
        with patch.object(gsea.np.random, "default_rng", side_effect=AssertionError("must not assign real IDs")), \
             self.assertRaisesRegex(ValueError, "invented SIM identifiers only"):
            self.run_fixture()

    def test_concurrent_writer_does_not_invalidate_a_completed_run(self):
        self.run_fixture()
        before = (self.output / "gsea_provenance.json").read_bytes()
        with (self.output / ".gsea.lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError, "Another GSEA run"):
                self.run_fixture()
        self.assertEqual((self.output / "gsea_provenance.json").read_bytes(), before)

    def test_changed_outputs_sources_and_manifest_are_rejected(self):
        result = self.run_fixture()
        paths = [self.output / result["analyses"][0]["display_file"],
                 self.output / result["mapping"]["file"], self.args.kegg_gmt,
                 self.args.kegg_provenance, self.data / gsea.CONTRASTS[0][3]]
        for path in paths:
            original = path.read_bytes()
            try:
                path.write_bytes(original + b"\n")
                with self.subTest(path=path.name), self.assertRaises(ValueError):
                    gsea.validate_completed_run(self.output, self.data)
            finally:
                path.write_bytes(original)
        manifest = self.output / "gsea_provenance.json"
        for change in (lambda data: data["artifact_sha256"].pop(next(iter(data["artifact_sha256"]))),
                       lambda data: data["analyses"].__setitem__(0, data["analyses"][1]),
                       lambda data: data["analyses"][1].__setitem__("source_contrast_sha256", "bad")):
            bad = json.loads(json.dumps(result))
            change(bad)
            gsea.write_manifest(bad, manifest)
            with self.assertRaises(ValueError):
                gsea.validate_completed_run(self.output, self.data)
        gsea.write_manifest(result, manifest)

    def test_inputs_changing_during_analysis_never_receive_completion(self):
        for changed in (self.args.kegg_gmt, self.data / gsea.CONTRASTS[0][3]):
            original = changed.read_bytes()
            self.callback = lambda: changed.write_bytes(original + b"\n")
            try:
                with self.subTest(path=changed.name), self.assertRaisesRegex(RuntimeError, "changed during GSEA"):
                    self.run_fixture()
                self.assertFalse(self.read_manifest()["complete"])
            finally:
                changed.write_bytes(original)

    def test_engine_cannot_drop_or_reorder_ranked_genes(self):
        self.corrupt_ranking = True
        with self.assertRaisesRegex(RuntimeError, "complete deterministic ranking"):
            self.run_fixture()
        self.assertFalse(self.read_manifest()["complete"])

    @unittest.skipUnless(importlib.util.find_spec("gseapy"), "Install figures/requirements.txt for the GSEApy integration check")
    def test_actual_weighted_enrichment_scores_match_independent_running_sum(self):
        import gseapy

        symbols = list("ABCDEFGHIJ")
        statistics = np.array([5., 4., 3., 2., 1., -1., -2., -3., -4., -5.])
        sets = {"positive": list("ABC"), "negative": list("HIJ"), "mixed": list("AI")}
        rank = pd.DataFrame({"mouse_symbol": symbols, "stat": statistics})
        result = gseapy.prerank(rnk=rank, gene_sets=sets, outdir=None, ascending=None,
            min_size=2, max_size=8, permutation_num=100, weight=1., seed=20261019,
            threads=1, no_plot=True, verbose=False, method="permutation")
        tidy = gsea.normalize_result(result, sets, set(symbols)).set_index("term")
        for term, members in sets.items():
            hits = np.isin(symbols, members)
            weighted_hits = np.where(hits, np.abs(statistics) / np.abs(statistics[hits]).sum(), 0.)
            misses = (~hits) / (~hits).sum()
            running_sum = np.cumsum(weighted_hits - misses)
            expected = max(running_sum) if max(running_sum) > abs(min(running_sum)) else min(running_sum)
            self.assertAlmostEqual(tidy.loc[term, "ES"], expected, places=12)
        self.assertEqual(result.ranking.index.tolist(), symbols)
        self.assertEqual(result.ranking.to_numpy().tolist(), statistics.tolist())


if __name__ == "__main__":
    unittest.main()
