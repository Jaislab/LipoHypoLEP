import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from analysis import analyze_male_2026_09_30 as a


def source(rows):
    return (','.join(a.FIELDS) + '\n' + '\n'.join(','.join(map(str, row)) for row in rows) + '\n').encode()


def gene(number):
    return f'ENSMUSG{number:011d}'


class MaleAnalysisTests(unittest.TestCase):
    def test_treatment_reversal_preserves_original_effect_and_probabilities(self):
        table = a.read_source(source([(gene(1), 'Lcn2', 2.4284321933091206, 8.020379548114328e-7, .006420714847242925)]), 'input.csv', -1)
        row = table.iloc[0]
        self.assertEqual(row.log2FoldChange, -2.4284321933091206)
        self.assertEqual(row.source_log2FoldChange, 2.4284321933091206)
        self.assertEqual(row.pvalue, 8.020379548114328e-7)
        self.assertEqual(row.padj, .006420714847242925)

    def test_joint_cutoff_boundaries(self):
        table = a.read_source(source([(gene(1), 'A', .2, .001, .049), (gene(2), 'B', -.2, .001, .049),
                                      (gene(3), 'C', .2, .001, .05), (gene(4), 'D', .199, .001, .01)]), 'input.csv', 1)
        self.assertEqual(table.loc[a.joint_hits(table), 'gene_id'].tolist(), [gene(1), gene(2)])

    def test_duplicate_ensembl_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate Ensembl'):
            a.read_source(source([(gene(1), 'A', .2, .001, .01)] * 2), 'input.csv', 1)

    def test_invalid_numeric_and_identifier_fields_are_rejected(self):
        for row in [('human', 'A', .2, .01, .02), (gene(1), 'A', 'nan', .01, .02),
                    (gene(1), 'A', .2, 0, .02), (gene(1), 'A', .2, .01, 1.1),
                    (gene(1), 'A', '1e-999', .01, .02)]:
            with self.subTest(row=row), self.assertRaises(ValueError):
                a.read_source(source([row]), 'input.csv', 1)

    def test_gene_selection_uses_adjusted_p_and_stable_id_ties(self):
        table = a.read_source(source([(gene(i), f'Gene{i}', .5, .001, .01) for i in range(12, 0, -1)]), 'input.csv', 1)
        picked = a.selected_genes(table)
        self.assertEqual(len(picked), 10)
        self.assertEqual(set(picked.gene_id), {gene(i) for i in range(1, 11)})

    def test_zero_hits_are_retained_without_inventing_selections(self):
        table = a.read_source(source([(gene(1), 'A', .1, .1, .2)]), 'input.csv', 1)
        self.assertEqual(len(a.selected_genes(table)), 0)

    def run_fixture(self, directory):
        source_dir, output_dir = Path(directory) / 'sources', Path(directory) / 'outputs'
        source_dir.mkdir()
        for _, _, _, filename, _ in a.SPECS:
            (source_dir / filename).write_bytes(source([(gene(1), 'A', .5, .001, .04), (gene(2), 'B', -.3, .002, .06)]))
        with contextlib.redirect_stdout(io.StringIO()):
            manifest = a.run(source_dir, output_dir)
        return source_dir, output_dir, manifest

    def test_filtered_exports_preserve_adjusted_p_and_gate_enrichment(self):
        with tempfile.TemporaryDirectory() as temp:
            _, output, manifest = self.run_fixture(temp)
            self.assertEqual([r['joint_cutoff_hits'] for r in manifest['figures']], [1, 1, 1])
            self.assertEqual(manifest['figures'][1]['lower'], 1)
            self.assertFalse(manifest['enrichment']['performed'])
            import pandas as pd
            table = pd.read_csv(output / 'figure_02_supplied.tsv', sep='\t', float_precision='round_trip')
            self.assertEqual(table.padj.tolist(), [.04, .06])
            self.assertEqual(table.log2FoldChange.tolist(), [-.5, .3])
            a.validate_results(output)

    def test_output_and_source_changes_invalidate_saved_results(self):
        for changed_source in (True, False):
            with self.subTest(source=changed_source), tempfile.TemporaryDirectory() as temp:
                inputs, output, _ = self.run_fixture(temp)
                path = inputs / a.SPECS[0][3] if changed_source else output / 'figure_01_supplied.tsv'
                path.write_bytes(path.read_bytes() + b'\n')
                with self.assertRaisesRegex(ValueError, 'changed'):
                    a.validate_results(output)

    def test_failed_rerun_invalidates_previous_completion(self):
        with tempfile.TemporaryDirectory() as temp:
            inputs, output, _ = self.run_fixture(temp)
            (inputs / a.SPECS[0][3]).write_text('bad,input\n')
            with self.assertRaises(ValueError):
                a.run(inputs, output)
            self.assertFalse(json.loads((output / 'analysis_manifest.json').read_text())['complete'])


if __name__ == '__main__':
    unittest.main()
