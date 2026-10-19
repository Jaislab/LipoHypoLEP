"""Scientific display edge cases, independent of the synthetic effect generator."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "demo_renderer", REPO / "figures/synthetic_demo/render_figure.py")
render = importlib.util.module_from_spec(spec)
spec.loader.exec_module(render)


def contrast(effects=(.2, -.2, .199, -.3), probabilities=(.049, .049, .001, .05)):
    return pd.DataFrame({"log2FoldChange": effects, "padj": probabilities,
                         "lfcSE": [.1] * len(effects)},
                        index=pd.Index([f"SIM{i}" for i in range(len(effects))], name="gene_id"))


def inputs():
    table = contrast()
    metadata = pd.DataFrame({"group": ["LD SALINE", "LD SALINE", "LD LEPTIN", "NON-LD control"]},
                            index=pd.Index(["s1", "s2", "l1", "k1"], name="sample_id"))
    expression = pd.DataFrame(np.arange(16).reshape(4, 4), index=table.index, columns=metadata.index)
    return [table.copy() for _ in range(3)], expression, metadata


class CutoffTests(unittest.TestCase):
    def test_cutoff_is_inclusive_for_effect_exclusive_for_probability(self):
        table = contrast()
        original = table.copy(deep=True)
        heatmap, forest = render.selected_rows(table)
        self.assertEqual(heatmap.index.tolist(), ["SIM0", "SIM1"])
        self.assertEqual(set(forest.index), {"SIM0", "SIM1"})
        pd.testing.assert_frame_equal(table, original)

    def test_empty_and_one_sided_selections_are_valid(self):
        empty, forest = render.selected_rows(contrast(probabilities=(1, 1, 1, 1)))
        self.assertTrue(empty.empty)
        self.assertTrue(forest.empty)
        heatmap, forest = render.selected_rows(contrast(effects=(.3, .4, .5, .6)))
        self.assertEqual(len(heatmap), 3)
        self.assertEqual(len(forest), 3)

    def test_effects_are_not_rescaled_or_clipped(self):
        table = contrast(effects=(5.4, -3.1, .001, .3), probabilities=(.001,) * 4)
        heatmap, _ = render.selected_rows(table)
        self.assertEqual(heatmap.log2FoldChange.tolist(), [5.4, -3.1, .3])
        self.assertGreaterEqual(render.volcano_limits([table])[0], 5.4)


class ValidationTests(unittest.TestCase):
    def test_valid_unequal_group_sizes(self):
        render.validate_inputs(*inputs())

    def test_misaligned_or_duplicate_identifiers_fail(self):
        for mutation in ("contrast_order", "duplicate_gene", "expression_order", "sample_set", "group_name"):
            with self.subTest(mutation=mutation):
                tables, expression, metadata = inputs()
                if mutation == "contrast_order":
                    tables[1] = tables[1].iloc[::-1]
                elif mutation == "duplicate_gene":
                    tables[0].index = ["same"] * 4
                elif mutation == "expression_order":
                    expression = expression.iloc[::-1]
                elif mutation == "sample_set":
                    expression = expression.rename(columns={"s1": "unknown"})
                else:
                    metadata.loc["s1", "group"] = "L"
                with self.assertRaises(ValueError):
                    render.validate_inputs(tables, expression, metadata)

    def test_invalid_probabilities_errors_and_nonfinite_data_fail(self):
        for column, value in (("padj", -1), ("padj", 1.1), ("padj", np.nan),
                              ("log2FoldChange", np.inf), ("lfcSE", -.1)):
            with self.subTest(column=column, value=value):
                tables, expression, metadata = inputs()
                tables[0].loc["SIM0", column] = value
                with self.assertRaises(ValueError):
                    render.validate_inputs(tables, expression, metadata)
        tables, expression, metadata = inputs()
        expression.iloc[0, 0] = np.nan
        with self.assertRaises(ValueError):
            render.validate_inputs(tables, expression, metadata)

    def test_failed_rerun_invalidates_previous_success_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "demo"
            root.mkdir()
            audit = root / "font_and_plot_audit.json"
            audit.write_text('{"complete": true}')
            with patch.object(render, "ROOT", root), patch.object(render, "DATA", root / "missing_data"):
                with self.assertRaises((ValueError, FileNotFoundError)):
                    render.main()
            self.assertFalse(audit.exists())


class NumericalDisplayTests(unittest.TestCase):
    def test_table_reader_preserves_full_binary64_values(self):
        # These estimates include a cutoff boundary and values requiring 17 digits.
        values=np.array([np.nextafter(.2,0),np.nextafter(.2,1),1.1554300887446989,
                         .30661085401655512,np.nextafter(.05,0)])
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'precision.tsv'
            pd.DataFrame({'estimate':values}).to_csv(path,sep='\t',index=False,float_format='%.17g')
            np.testing.assert_array_equal(render.read_table(path).estimate.to_numpy(),values)

    def test_zero_probability_is_finite_without_mutating_source(self):
        probabilities = np.array([0, .01, 1])
        y = render.volcano_ordinates(probabilities)
        self.assertTrue(np.isfinite(y).all())
        self.assertGreater(y[0], 300)
        np.testing.assert_allclose(y[1:], [2, 0])
        self.assertEqual(probabilities[0], 0)

    def test_zero_effect_all_nonsignificant_limits_retain_cutoff_guides(self):
        self.assertEqual(render.volcano_limits([contrast(effects=(0,) * 4, probabilities=(1,) * 4)]), (1, 2))

    def test_constant_and_single_heatmap_rows_do_not_make_nan(self):
        expression = pd.DataFrame([[2, 2, 2], [1, 2, 3]], index=["constant", "variable"])
        z, genes, constant = render.heatmap_values(expression, expression.index, expression.columns)
        self.assertEqual(constant, ["constant"])
        self.assertTrue(np.isfinite(z).all())
        np.testing.assert_array_equal(z[genes.index("constant")], [0, 0, 0])
        np.testing.assert_allclose(z[genes.index("variable")], [-1, 0, 1])
        z, genes, constant = render.heatmap_values(expression, ["constant"], expression.columns)
        self.assertEqual(z.shape, (1, 3))

    def test_empty_heatmap_is_supported(self):
        expression = pd.DataFrame([[1, 2]], index=["g"])
        z, genes, constant = render.heatmap_values(expression, [], expression.columns)
        self.assertEqual(z.shape, (0, 2))
        self.assertEqual((genes, constant), ([], []))

    def test_q_labels_preserve_threshold_and_positive_values(self):
        self.assertEqual(render.q_label(.049), "<0.05")
        self.assertEqual(render.q_label(.049999999), "<0.05")
        self.assertEqual(render.q_label(.05), "0.05")
        self.assertEqual(render.q_label(.001), "0.001")
        self.assertEqual(render.q_label(0), "0")
        self.assertNotEqual(float(render.q_label(1e-120)), 0)


class GseaDisplayTests(unittest.TestCase):
    def test_long_terms_fail_layout_check_without_reducing_font_or_truncating_text(self):
        font_dir=Path('/usr/share/fonts/opentype/urw-base35')
        if not (font_dir/'NimbusSans-Regular.otf').is_file():
            self.skipTest('Nimbus Sans source fonts are unavailable')
        fig=render.plt.figure(figsize=(render.W/render.DPI,render.H/render.DPI),dpi=render.DPI)
        ax=fig.add_axes([538/render.W,116/render.H,150/render.W,340/render.H])
        long_term='regulation of hypothalamic neuronal cellular metabolic response '*6
        full_label=render.textwrap.fill(long_term,width=42,break_long_words=False,break_on_hyphens=False)
        ax.set_yticks(np.arange(4),labels=[full_label,'Second pathway','Third pathway','Fourth pathway'])
        ax.set(ylim=(3.5,-.5))
        for label in ax.get_yticklabels():
            label.set_fontproperties(render.fm.FontProperties(fname=str(font_dir/'NimbusSans-Regular.otf'),size=render.px(21)))
        fig._gsea_layout_regions=[{'axes':ax,'collection':'GO BP','q_labels':[],'left':40,'right':768}]
        try:
            fig.canvas.draw()
            with self.assertRaisesRegex(RuntimeError,'GO BP GSEA layout.*allocate more panel space'):
                render.validate_gsea_layout(fig)
            self.assertEqual(ax.get_yticklabels()[0].get_text(),full_label)
            self.assertEqual(ax.get_yticklabels()[0].get_fontsize(),render.px(21))
        finally:
            render.plt.close(fig)

    def test_one_sided_and_empty_gsea_are_supported_without_inventing_hits(self):
        table = pd.DataFrame({"term": ["a", "b", "c"], "NES": [1, 3, 2], "gsea_fdr_q": [.2, .4, .8]})
        shown = render.displayed_gsea(table)
        self.assertEqual(shown.term.tolist(), ["b", "c"])
        self.assertEqual(render.displayed_gsea(table.iloc[:0]).shape, (0, 3))

    def test_invalid_gsea_probabilities_or_duplicate_terms_fail(self):
        for table in (pd.DataFrame({"term": ["a"], "NES": [np.inf], "gsea_fdr_q": [.1]}),
                      pd.DataFrame({"term": ["a"], "NES": [1], "gsea_fdr_q": [1.1]}),
                      pd.DataFrame({"term": ["a", "a"], "NES": [1, -1], "gsea_fdr_q": [.1, .2]})):
            with self.assertRaises(ValueError):
                render.displayed_gsea(table)

    def test_actual_layout_handles_no_degs_unequal_replicates_and_large_nes(self):
        tables, expression, metadata = inputs()
        table = tables[0].assign(padj=1)
        font_dir = Path("/usr/share/fonts/opentype/urw-base35")
        if not (font_dir / "NimbusSans-Regular.otf").is_file():
            self.skipTest("Nimbus Sans source fonts are unavailable")
        fonts = {style: font_dir / f"NimbusSans-{style}.otf" for style in ("Regular", "Bold", "Italic")}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "synthetic_demo"
            (root / "data").mkdir(parents=True)
            (root / "gsea_data").mkdir()
            for collection in ("kegg", "go_bp"):
                pd.DataFrame({"term": ["GO:000001 short pathway"], "NES": [4.7], "gsea_fdr_q": [.049]}).to_csv(
                    root / "gsea_data" / f"figure_01_{collection}_display.tsv", sep="\t", index=False)
            captured = {}

            def capture(fig, number, question, supplied_fonts, details):
                captured["figure"] = fig
                return details

            with patch.object(render, "ROOT", root), patch.object(render, "OUT", root.parent), \
                 patch.object(render, "DATA", root / "data"), patch.object(render, "export_figure", capture):
                details = render.render_question(render.SPECS[0], table, expression, metadata, fonts, 1, 2, 1, 9)
            fig = captured["figure"]
            try:
                self.assertEqual(details["heatmap_genes"], [])
                self.assertEqual(details["heatmap_samples"], ["s1", "s2", "k1"])
                self.assertGreater(details["gsea"]["kegg"]["absolute_nes_limit"], 4.7)
                strips = [ax for ax in fig.axes if len(ax.patches)]
                self.assertEqual(len(strips), 2)
                self.assertAlmostEqual(strips[0].get_position().width * render.W, 298)
                self.assertAlmostEqual(strips[1].get_position().width * render.W, 148)
            finally:
                render.plt.close(fig)

    def test_missing_nimbus_character_is_an_error(self):
        font_dir = Path("/usr/share/fonts/opentype/urw-base35")
        font = font_dir / "NimbusSans-Regular.otf"
        if not font.is_file():
            self.skipTest("Nimbus Sans source fonts are unavailable")
        render.configure_theme()
        fig = render.plt.figure()
        fig.text(.5, .5, "\U0001f9ec", fontproperties=render.fm.FontProperties(fname=str(font)))
        try:
            with self.assertRaisesRegex(UserWarning, "Glyph .* missing from font"):
                fig.canvas.draw()
        finally:
            render.plt.close(fig)

    def test_five_format_export_keeps_nimbus_with_zero_p_and_no_degs(self):
        font_dir = Path("/usr/share/fonts/opentype/urw-base35")
        if not (font_dir / "NimbusSans-Regular.otf").is_file():
            self.skipTest("Nimbus Sans source fonts are unavailable")
        _, expression, metadata = inputs()
        table = contrast(effects=(0,) * 4, probabilities=(0, 1, 1, 1))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "demo"
            (root / "data").mkdir(parents=True)
            (root / "gsea_data").mkdir()
            for collection in ("kegg", "go_bp"):
                pd.DataFrame({"term": ["GO:000001 Short pathway"], "NES": [4.7], "gsea_fdr_q": [.049]}).to_csv(
                    root / "gsea_data" / f"figure_01_{collection}_display.tsv", sep="\t", index=False)
            with patch.object(render, "ROOT", root), patch.object(render, "OUT", root.parent), \
                 patch.object(render, "DATA", root / "data"):
                fonts = render.prepare_nimbus()
                render.configure_theme()
                xlim, ylim = render.volcano_limits([table])
                details = render.render_question(render.SPECS[0], table, expression, metadata, fonts,
                                                 xlim, ylim, 1, 9)
            self.assertEqual(len(list(root.parent.glob("figure_01*"))), 7)
            for extension in ('.png','.pdf'):
                self.assertEqual((root.parent/f'figure_01_gsea_synthetic{extension}').read_bytes(),
                                 (root.parent/f'figure_01_large_axis_synthetic{extension}').read_bytes())
            self.assertTrue(all("NimbusSans" in font for font in details["pdf_fonts"]))
            self.assertEqual(details["heatmap_genes"], [])
            self.assertEqual(ylim, 324)


if __name__ == "__main__":
    unittest.main()
