"""Independent real-PDF typography and stale-evidence failure paths."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import fitz

from analysis import verify_legacy_outputs as verifier


class LegacyVerificationTests(unittest.TestCase):
    def make_pdf(self, path, sizes=(32,), nimbus=True, canvas=(1536, 1024)):
        with fitz.open() as document:
            page = document.new_page(width=canvas[0] * .72, height=canvas[1] * .72)
            name = 'NimbusFixture' if nimbus else 'helv'
            if nimbus:
                page.insert_font(fontname=name, fontfile=str(verifier.font_audit.SOURCE / 'NimbusSans-Regular.otf'))
            for number, size in enumerate(sizes):
                page.insert_text((72, 72 + 40 * number), 'KEGG' if number == 0 else 'q',
                                 fontname=name, fontsize=size * .72)
            document.save(path)

    def test_every_actual_pdf_span_is_checked(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'actual_nimbus.pdf'
            self.make_pdf(path, sizes=(32, 32))
            report = verifier.audit_typography(path)
            self.assertEqual(report['nonempty_text_spans_checked'], 2)
            self.assertTrue(report['every_span_is_32_native_px'])
            self.assertEqual(report['actual_font_names'], ['NimbusSans-Regular'])

    def test_small_annotation_cannot_hide_behind_a_correct_heading(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'mixed_size.pdf'
            self.make_pdf(path, sizes=(32, 24))
            with self.assertRaisesRegex(ValueError, "Actual PDF text 'q'.*expected 32px"):
                verifier.audit_typography(path)

    def test_actual_fallback_font_is_rejected_even_at_correct_size(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'fallback.pdf'
            self.make_pdf(path, nimbus=False)
            with self.assertRaisesRegex(ValueError, 'Actual PDF font fallback'):
                verifier.audit_typography(path)

    def test_wrong_pdf_canvas_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'wrong_canvas.pdf'
            self.make_pdf(path, canvas=(1000, 1000))
            with self.assertRaisesRegex(ValueError, 'PDF canvas'):
                verifier.audit_typography(path)

    def test_declared_checksum_cannot_be_missing(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'plot_data.tsv'
            path.write_text('gene\tvalue\nA\t1\n')
            with self.assertRaisesRegex(ValueError, 'Missing or malformed declared checksum'):
                verifier.checked_file(path, None)

    def test_changed_declared_artifact_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'plot_data.tsv'
            path.write_text('original\n')
            expected = verifier.digest(path)
            path.write_text('changed\n')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                verifier.checked_file(path, expected)

    def test_failed_rerun_replaces_old_pass(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            report = output / verifier.REPORT_NAME
            report.write_text('{"status":"PASS"}\n')
            with mock.patch.object(verifier, 'collect_state', side_effect=RuntimeError('incomplete GSEA')), \
                 mock.patch.object(verifier.font_audit, 'main', side_effect=AssertionError('synthetic entry point called')):
                with self.assertRaisesRegex(RuntimeError, 'incomplete GSEA'):
                    verifier.run(output / 'data', output / 'gsea', output)
            failed = json.loads(report.read_text())
            self.assertEqual(failed['status'], 'FAILED')
            self.assertIn('incomplete GSEA', failed['error'])

    def test_saved_pass_becomes_stale_when_completion_manifest_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            manifest = output / 'gsea_manifest.json'
            manifest.write_text('{"complete":true,"version":1}\n')

            def collect(*args):
                return {'input_sha256': {'gsea_manifest': verifier.digest(manifest)}}, {}, {}, {}

            record = {'status': 'PASS', 'integrity_schema': 1, 'font_pixels': 32,
                      'paths': {'data_dir': str(output / 'data'), 'gsea_dir': str(output / 'gsea'),
                                'output_dir': str(output)},
                      'state_sha256': collect()[0],
                      'figures': [{'figure': number, 'typography': {'every_span_is_32_native_px': True}}
                                  for number in (1, 2, 3)]}
            report = output / verifier.REPORT_NAME
            report.write_text(json.dumps(record))
            with mock.patch.object(verifier, 'collect_state', side_effect=collect):
                self.assertEqual(verifier.verify_report_current(report)['status'], 'PASS')
                manifest.write_text('{"complete":true,"version":2}\n')
                with self.assertRaisesRegex(ValueError, 'validation is stale'):
                    verifier.verify_report_current(report)


if __name__ == '__main__':
    unittest.main()
