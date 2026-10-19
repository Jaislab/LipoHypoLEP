"""Independently verify exported weighted GSEA scores and leading edges.

The reference uses the defining running-sum calculation, without importing
GSEApy or the runner under test. It does not independently reproduce null
permutations, normalized enrichment scores, or empirical FDR estimates.
"""

from __future__ import annotations

import json
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "figures" / "synthetic_demo" / "gsea_data"


def reference_weighted_walk(symbols, statistics, pathway, positions=None):
    """Return weight-1 ES and its leading edge directly from a ranked list."""
    if positions is None:
        positions = {symbol: index for index, symbol in enumerate(symbols)}
    hits = np.array(sorted(positions[gene] for gene in pathway if gene in positions))
    if not 0 < len(hits) < len(symbols):
        raise ValueError("The reference needs both pathway hits and misses")
    weights = np.abs(statistics[hits])
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        raise ValueError("Pathway hit weights must have a positive finite sum")
    increments = np.full(len(symbols), -1.0 / (len(symbols) - len(hits)))
    increments[hits] = weights / weights.sum()
    running = np.cumsum(increments)
    peak, trough = int(running.argmax()), int(running.argmin())
    score = running[peak] if running[peak] > -running[trough] else running[trough]
    leading_hits = hits[hits <= peak] if score > 0 else hits[hits > trough]
    return float(score), frozenset(symbols[leading_hits])


class ReferenceWalkTests(unittest.TestCase):
    def test_positive_and_negative_extremes_have_the_expected_leading_edge(self):
        symbols = np.array(list("ABCDEFGH"))
        statistics = np.array([4.0, 3.0, 2.0, 1.0, -1.0, -2.0, -3.0, -4.0])
        positive, positive_edge = reference_weighted_walk(symbols, statistics, {"A", "B"})
        negative, negative_edge = reference_weighted_walk(symbols, statistics, {"G", "H"})
        self.assertAlmostEqual(positive, 1.0, places=14)
        self.assertAlmostEqual(negative, -1.0, places=14)
        self.assertEqual(positive_edge, frozenset({"A", "B"}))
        self.assertEqual(negative_edge, frozenset({"G", "H"}))


class ExportedScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        provenance = json.loads((DATA / "gsea_provenance.json").read_text())
        if not provenance.get("complete") or len(provenance.get("analyses", [])) != 6:
            raise AssertionError("The reference audit requires six completed GSEA analyses")
        if provenance["parameters"]["weight"] != 1.0:
            raise AssertionError("This reference implements the declared weight-1 calculation")

        collections = {}
        for library, source in provenance["sources"].items():
            path = Path(source["gmt_path"])
            if not path.is_absolute():
                path = ROOT / path
            collection = {}
            for line in path.read_text(encoding="utf-8-sig").splitlines():
                if line.strip():
                    fields = line.split("\t")
                    collection[fields[0]] = set(fields[2:])
            collections[library] = collection

        ranks = {}
        cls.maximum_error = 0.0
        cls.edge_mismatches = []
        cls.audited_sets = 0
        cls.expected_sets = sum(analysis["tested_gene_sets"] for analysis in provenance["analyses"])
        for analysis in provenance["analyses"]:
            rank_file = analysis["ranked_statistics_file"]
            if rank_file not in ranks:
                table = pd.read_csv(DATA / rank_file, sep="\t", usecols=["mouse_symbol", "stat"],
                                    float_precision="round_trip")
                if not table.mouse_symbol.is_unique or not table.stat.is_monotonic_decreasing:
                    raise AssertionError("Ranked symbols must be unique and signed statistics descending")
                symbols = table.mouse_symbol.to_numpy()
                positions = {symbol: index for index, symbol in enumerate(symbols)}
                ranks[rank_file] = (symbols, table.stat.to_numpy(), positions)
            symbols, statistics, positions = ranks[rank_file]
            results = pd.read_csv(
                DATA / analysis["all_results_file"], sep="\t",
                usecols=["term", "ES", "leading_edge_genes"], keep_default_na=False,
                float_precision="round_trip",
            )
            if len(results) != analysis["tested_gene_sets"] or not results.term.is_unique:
                raise AssertionError("Result rows must cover every declared tested set exactly once")
            for row in results.itertuples(index=False):
                score, edge = reference_weighted_walk(
                    symbols, statistics, collections[analysis["library"]][row.term], positions,
                )
                cls.maximum_error = max(cls.maximum_error, abs(score - row.ES))
                exported_edge = frozenset(row.leading_edge_genes.split(";")) if row.leading_edge_genes else frozenset()
                if edge != exported_edge:
                    cls.edge_mismatches.append((analysis["figure"], analysis["library"], row.term))
                cls.audited_sets += 1

    def test_all_enrichment_scores_match_the_independent_running_sum(self):
        self.assertEqual(self.audited_sets, self.expected_sets)
        self.assertEqual(self.audited_sets, 9219, "Pinned demonstration collection changed; review the reference fixture")
        self.assertLessEqual(self.maximum_error, 2e-12)

    def test_all_positive_and_negative_leading_edges_match_the_reference(self):
        self.assertEqual(self.edge_mismatches, [])


if __name__ == "__main__":
    unittest.main()
