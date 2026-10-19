"""Failure-path tests for independent export evidence; no figure regeneration."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
VERIFIER_PATH = REPOSITORY / 'figures' / 'synthetic_demo' / 'verify_export_fonts.py'
spec = importlib.util.spec_from_file_location('export_verifier', VERIFIER_PATH)
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)
color_spec = importlib.util.spec_from_file_location(
    'color_auditor', REPOSITORY / 'figures' / 'accessibility' / 'audit_colors.py')
color_auditor = importlib.util.module_from_spec(color_spec)
try:
    color_spec.loader.exec_module(color_auditor)
except ModuleNotFoundError as error:
    if error.name != 'daltonlens':
        raise
    color_auditor = None


class ExportVerificationTests(unittest.TestCase):
    def test_checks_remain_active_with_python_optimized(self):
        with tempfile.TemporaryDirectory() as temporary:
            svg = Path(temporary) / 'unexpected_live_text.svg'
            svg.write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>Arial</text></svg>')
            code = ('import sys; from pathlib import Path; sys.path.insert(0, sys.argv[1]); '
                    'import verify_export_fonts as v; v.audit_outlined_svg(Path(sys.argv[2]))')
            result = subprocess.run([sys.executable, '-O', '-c', code,
                                     str(VERIFIER_PATH.parent), str(svg)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('FontVerificationError', result.stderr)
            self.assertIn('still contains live text', result.stderr)

    def test_outlined_svg_rejects_unresolved_nimbus_glyph(self):
        with tempfile.TemporaryDirectory() as temporary:
            svg = Path(temporary) / 'missing_glyph.svg'
            svg.write_text('<svg xmlns="http://www.w3.org/2000/svg" '
                           'xmlns:xlink="http://www.w3.org/1999/xlink">'
                           '<g id="text_1"><use xlink:href="#NimbusSans-Regular-41"/></g></svg>')
            with self.assertRaisesRegex(verifier.FontVerificationError, 'Missing outlined glyph'):
                verifier.audit_outlined_svg(svg)

    def test_embedded_program_rejects_false_family_before_outline_checks(self):
        with mock.patch.object(verifier, 'TTFont', return_value={}), \
             mock.patch.object(verifier, 'names', return_value={'1': ['Arial'], '2': ['Regular']}):
            with self.assertRaisesRegex(verifier.FontVerificationError, 'Unexpected embedded family'):
                verifier.audit_program(b'fixture', {})

    def test_mathtext_svg_checks_tspan_family_without_requiring_parent_family(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'mathtext.svg'
            fonts = ''.join('@font-face {font-family:"Nimbus Sans"; '
                            'src:url(data:font/ttf;base64,AAAA);}' for _ in range(3))
            path.write_text('<svg xmlns="http://www.w3.org/2000/svg"><style>' + fonts +
                            '</style><text>\n<tspan style="font-family: Nimbus Sans">x</tspan>'
                            '\n</text></svg>')
            programs = [{'style': style} for style in verifier.STYLES]
            with mock.patch.object(verifier, 'audit_program', side_effect=programs):
                self.assertEqual(verifier.audit_editable_svg(path, {})['live_text_elements'], 1)

    def test_svg_rejects_unstyled_glyphs_despite_other_nimbus_declarations(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'fallback.svg'
            fonts = ''.join('@font-face {font-family:"Nimbus Sans"; '
                            'src:url(data:font/ttf;base64,AAAA);}' for _ in range(3))
            path.write_text('<svg xmlns="http://www.w3.org/2000/svg"><style>' + fonts +
                            '</style><text>fallback</text>'
                            '<text style="font-family: Nimbus Sans">verified</text></svg>')
            programs = [{'style': style} for style in verifier.STYLES]
            with mock.patch.object(verifier, 'audit_program', side_effect=programs):
                with self.assertRaisesRegex(verifier.FontVerificationError, 'no verified Nimbus'):
                    verifier.audit_editable_svg(path, {})

    def test_png_rejects_wrong_dimensions_before_pdf_access(self):
        from PIL import Image
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'wrong_size.png'
            Image.new('RGB', (10, 10)).save(path)
            with self.assertRaisesRegex(verifier.FontVerificationError, 'Unexpected image size'):
                verifier.audit_png(path, Path(temporary) / 'missing.pdf')

    def test_failed_run_invalidates_prior_pass(self):
        reports = []
        with mock.patch.object(sys, 'argv', ['verify_export_fonts.py']), \
             mock.patch.object(verifier, 'write_report', side_effect=lambda record: reports.append(record)), \
             mock.patch.object(verifier, 'artifact_snapshot', side_effect=RuntimeError('missing export')):
            with self.assertRaisesRegex(RuntimeError, 'missing export'):
                verifier.main()
        self.assertEqual([record['status'] for record in reports], ['RUNNING', 'FAILED'])
        self.assertIn('missing export', reports[-1]['error'])

    def test_atomic_report_rejects_nonfinite_values(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'report.json'
            path.write_text('{"status":"old"}\n')
            with self.assertRaises(ValueError):
                verifier.write_report({'value': float('nan')}, path)
            self.assertEqual(path.read_text(), '{"status":"old"}\n')
            self.assertEqual(list(Path(temporary).iterdir()), [path])

    def make_fresh_report(self, root, out, source):
        (root / 'render_figure.py').write_text('# test renderer\n')
        fonts = []
        for style in verifier.STYLES:
            path = source / f'NimbusSans-{style}.otf'
            path.write_bytes(style.encode())
            fonts.append({'style': style, 'sha256': verifier.sha256(path.read_bytes())})
        figures = []
        for number in (1, 2, 3):
            stem = f'figure_{number:02d}_gsea_synthetic'
            hashes = {}
            for extension in ('.pdf', '.svg', '_editable.svg', '.png', '_300dpi.png'):
                path = out / f'{stem}{extension}'
                path.write_bytes(path.name.encode())
                hashes[path.name] = verifier.sha256(path.read_bytes())
            figures.append({'figure': number, 'file_stem': stem, 'export_sha256': hashes})
        return {'status': 'PASS', 'source_fonts': fonts, 'figures': figures,
                'renderer_sha256': verifier.sha256((root / 'render_figure.py').read_bytes()),
                'verifier_sha256': verifier.sha256(VERIFIER_PATH.read_bytes())}

    def test_saved_pass_is_rejected_after_export_edit(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root, out, source = (base / name for name in ('code', 'exports', 'fonts'))
            for directory in (root, out, source):
                directory.mkdir()
            report = self.make_fresh_report(root, out, source)
            path = base / 'report.json'
            path.write_text(json.dumps(report))
            with mock.patch.multiple(verifier, ROOT=root, OUT=out, SOURCE=source):
                self.assertEqual(verifier.verify_report_current(path)['status'], 'PASS')
                (out / 'figure_02_gsea_synthetic.pdf').write_bytes(b'replaced')
                with self.assertRaisesRegex(verifier.FontVerificationError, 'Export changed or missing'):
                    verifier.verify_report_current(path)

    def test_report_requires_three_distinct_figures(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            root, out, source = (base / name for name in ('code', 'exports', 'fonts'))
            for directory in (root, out, source):
                directory.mkdir()
            report = self.make_fresh_report(root, out, source)
            report['figures'] = [deepcopy(report['figures'][0]) for _ in range(3)]
            path = base / 'report.json'
            path.write_text(json.dumps(report))
            with mock.patch.multiple(verifier, ROOT=root, OUT=out, SOURCE=source):
                with self.assertRaisesRegex(verifier.FontVerificationError, 'Incomplete export evidence'):
                    verifier.verify_report_current(path)


@unittest.skipIf(color_auditor is None, 'Install figures/accessibility/requirements.txt for color-audit tests')
class ColorAuditTests(unittest.TestCase):
    def test_grayscale_uses_linear_light_luminance(self):
        import numpy as np
        # Rec.709 red luminance is 0.2126729; its sRGB encoding rounds to 127.
        result = color_auditor.transform(np.array([[[255, 0, 0]]], dtype=np.uint8), 'gray', None)
        np.testing.assert_array_equal(result, np.array([[[127, 127, 127]]], dtype=np.uint8))

    def test_lab_maps_black_and_white_to_reference_endpoints(self):
        import numpy as np
        values = color_auditor.lab(np.array([[0, 0, 0], [255, 255, 255]], dtype=np.uint8))
        np.testing.assert_allclose(values[0], [0, 0, 0], atol=1e-10)
        np.testing.assert_allclose(values[1], [100, 0, 0], atol=2e-5)

    def test_failed_color_run_invalidates_prior_pass(self):
        reports = []
        with mock.patch.object(sys, 'argv', ['audit_colors.py']), \
             mock.patch.object(color_auditor, 'write_report', side_effect=lambda record: reports.append(record)), \
             mock.patch.object(color_auditor, 'source_snapshot', side_effect=RuntimeError('missing figure')):
            with self.assertRaisesRegex(RuntimeError, 'missing figure'):
                color_auditor.main()
        self.assertEqual([record['status'] for record in reports], ['RUNNING', 'FAILED'])

    def test_color_report_requires_exact_source_figure_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            fonts = base / 'fonts'
            fonts.mkdir()
            regular = fonts / 'NimbusSans-Regular.otf'
            regular.write_bytes(b'test font evidence')
            report = {'status': 'PASS', 'auditor_sha256': color_auditor.sha256(Path(color_auditor.__file__)),
                      'font_sha256': color_auditor.sha256(regular), 'input_export_sha256': {}}
            path = base / 'report.json'
            path.write_text(json.dumps(report))
            with mock.patch.object(color_auditor, 'FONT_DIR', fonts):
                with self.assertRaisesRegex(RuntimeError, 'Incomplete color-audit input evidence'):
                    color_auditor.verify_report_current(path)


if __name__ == '__main__':
    unittest.main()
