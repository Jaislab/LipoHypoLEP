"""Legacy figures preserve supplied gene estimates and verified export typography."""
from pathlib import Path
import tempfile
import unittest

import fitz
import numpy as np
import pandas as pd

from analysis import render_legacy_figures as render


def contrast(effects=(.2, -.2, .199, -.3), probabilities=(.049, .049, .001, .05)):
    return pd.DataFrame({
        'gene_id': [f'ENSMUSG{i:011d}' for i in range(len(effects))],
        'supplied_symbol': [f'Gene{i}' for i in range(len(effects))],
        'log2FoldChange': effects, 'pvalue': np.asarray(probabilities) / 2,
        'padj': probabilities,
        # Deliberately meaningless legacy SE values: no intervals may be inferred.
        'lfcSE': [999] * len(effects),
    })


class SelectionTests(unittest.TestCase):
    def test_cutoff_boundaries_preserve_original_estimates(self):
        table = contrast()
        original = table.copy(deep=True)
        selected = render.select_gene_effects(table)
        self.assertEqual(selected.gene_id.tolist(), table.gene_id[:2].tolist())
        self.assertEqual(selected.log2FoldChange.tolist(), [.2, -.2])
        pd.testing.assert_frame_equal(table, original)

    def test_top_hits_are_selected_per_direction_without_effect_rescaling(self):
        table = contrast(effects=(10, 2, 1, .3, -8, -2, -1, -.3),
                         probabilities=(.04, .03, .02, .01, .04, .03, .02, .01))
        selected = render.select_gene_effects(table)
        self.assertEqual(selected.log2FoldChange.tolist(), [2, 1, .3, -.3, -1, -2])
        np.testing.assert_array_equal(selected.lfcSE, np.full(6, 999))
        self.assertGreater(render.effect_limit(selected), 2)

    def test_duplicate_symbols_retain_distinct_ensembl_rows(self):
        selected = contrast(effects=(1, -1), probabilities=(.01, .01))
        selected['supplied_symbol'] = ['Duplicate', 'Duplicate']
        labels = render.gene_effect_labels(selected)
        self.assertEqual(len(set(labels)), 2)
        for identifier, label in zip(selected.gene_id, labels):
            self.assertIn(identifier[-6:], label)
            self.assertTrue(label.startswith('Duplicate'))

    def test_sparse_and_no_hit_contrasts_do_not_invent_gene_rows(self):
        self.assertTrue(render.select_gene_effects(contrast(probabilities=(1,) * 4)).empty)
        selected = render.select_gene_effects(contrast(effects=(3.45903607, 0, 0, 0),
                                                       probabilities=(3.1687e-13, 1, 1, 1)))
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected.log2FoldChange.iloc[0], 3.45903607)


class InputTests(unittest.TestCase):
    def test_failed_rerun_invalidates_a_previous_completed_audit(self):
        import json

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'legacy_plot_audit.json').write_text('{"complete": true}')
            with self.assertRaises(ValueError):
                render.run(root / 'missing_data', root / 'missing_gsea', root)
            self.assertFalse(json.loads((root / 'legacy_plot_audit.json').read_text())['complete'])

    def test_all_audited_pathways_are_kept_when_one_nes_direction_is_scarce(self):
        display = pd.DataFrame({'term': ['One', 'Two', 'Three', 'Four'],
                                'NES': [-1, -2, -3, -4], 'gsea_fdr_q': [.1, .2, .3, .4]})
        shown = render.gsea_display(display)
        self.assertEqual(shown.term.tolist(), ['One', 'Two', 'Three', 'Four'])
        with self.assertRaises(ValueError):
            render.gsea_display(pd.concat([display, pd.DataFrame({'term': ['Five'], 'NES': [1],
                                                                 'gsea_fdr_q': [.01]})]))

    def test_reader_rejects_duplicate_genes_invalid_probabilities_and_missing_labels(self):
        for mutation in ('duplicate_gene', 'invalid_p', 'missing_symbol', 'nonfinite_effect'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as directory:
                table = contrast()
                if mutation == 'duplicate_gene':
                    table.loc[1, 'gene_id'] = table.loc[0, 'gene_id']
                elif mutation == 'invalid_p':
                    table.loc[0, 'padj'] = 1.01
                elif mutation == 'missing_symbol':
                    table.loc[0, 'supplied_symbol'] = ''
                else:
                    table.loc[0, 'log2FoldChange'] = np.inf
                path = Path(directory) / 'contrast.tsv'
                table.to_csv(path, sep='\t', index=False)
                with self.assertRaises(ValueError):
                    render.read_contrast(path)

    def test_reader_keeps_binary64_cutoff_boundaries(self):
        table = contrast(effects=(np.nextafter(.2, 0), np.nextafter(.2, 1),
                                  1.1554300887446989, .30661085401655512))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'contrast.tsv'
            table.to_csv(path, sep='\t', index=False, float_format='%.17g')
            read = render.read_contrast(path)
            np.testing.assert_array_equal(read.log2FoldChange, table.log2FoldChange)


class ExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = Path('/usr/share/fonts/opentype/urw-base35/NimbusSans-Regular.otf')
        if not source.is_file():
            raise unittest.SkipTest('Nimbus Sans source fonts are unavailable')
        cls.fonts = render.style.prepare_nimbus()
        render.style.configure_theme()

    def test_sparse_realistic_panel_exports_only_estimates_with_uniform_nimbus(self):
        table = contrast(effects=(3.45903607, .001, -.001, 0),
                         probabilities=(3.1687e-13, .9, 1, 1))
        table.loc[0, 'supplied_symbol'] = 'Tcaf2'
        enrichment = pd.DataFrame({
            'term': ['Synaptic signaling', 'Regulation of hormone secretion',
                     'Inflammatory response', 'Protein processing'],
            'NES': [2.1, 1.6, -1.8, -2.2], 'gsea_fdr_q': [.01, .049, .3, .9],
        })
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            details = render.render_question(render.SPECS[2], table,
                                              {'kegg': enrichment, 'go_bp': enrichment},
                                              self.fonts, output, 4, 14, 2.6)
            self.assertEqual(details['gene_effects']['display_labels'], ['Tcaf2'])
            self.assertEqual(details['gene_effects']['log2FoldChange'], [3.45903607])
            self.assertIsNone(details['gene_effects']['confidence_intervals'])
            self.assertFalse(details['heatmap_panel'])
            self.assertFalse(details['sample_expression_available'])
            self.assertEqual(len(details['export_sha256']), 5)
            self.assertEqual({record['size_px'] for record in details['labels']}, {32})
            self.assertEqual(details['gsea']['kegg']['gsea_fdr_q'], [.01, .049, .3, .9])
            self.assertNotIn('<text', (output / 'figure_03_legacy.svg').read_text())
            editable = (output / 'figure_03_legacy_editable.svg').read_text()
            self.assertIn('<text', editable)
            self.assertEqual(editable.count('@font-face'), 3)
            with fitz.open(output / 'figure_03_legacy.pdf') as pdf:
                spans = [span for page in pdf for block in page.get_text('dict')['blocks']
                         for line in block.get('lines', []) for span in line['spans'] if span['text'].strip()]
            self.assertTrue(spans)
            self.assertTrue(all(span['font'].startswith('NimbusSans-') for span in spans))
            self.assertTrue(all(abs(span['size'] * 100 / 72 - 32) < .01 for span in spans))


if __name__ == '__main__':
    unittest.main()
