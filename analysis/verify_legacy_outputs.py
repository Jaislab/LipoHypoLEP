#!/usr/bin/env python3
"""Independently verify legacy figure files, fonts, typography and provenance.

This only reads the shared font-audit functions. It never calls the synthetic
renderer or its fixed-path verification entry point, and never edits demo files.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys
import tempfile

import fitz
from fontTools.ttLib import TTFont

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))
from analysis.legacy_workbook import FIGURES, validate_completed_import
from analysis.legacy_gsea import validate_completed_run
from analysis.audit_legacy_gsea import validate_audit_report
from figures.synthetic_demo import verify_export_fonts as font_audit

DEFAULT_DATA = REPOSITORY / 'analysis/legacy_results'
DEFAULT_OUTPUT = REPOSITORY / 'figures/legacy'
REPORT_NAME = 'legacy_validation.json'
SUFFIXES = ('.png', '_300dpi.png', '.pdf', '.svg', '_editable.svg')
FONT_PIXELS = 32
CHECKSUM_UNSET = object()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def saved_path(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(REPOSITORY))
    except ValueError:
        return str(path)


def resolve_saved(path):
    path = Path(path)
    return path if path.is_absolute() else REPOSITORY / path


def write_report(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', encoding='utf-8', dir=path.parent,
                                     suffix='.json', delete=False) as handle:
        pending = Path(handle.name)
        try:
            json.dump(record, handle, indent=2, allow_nan=False)
            handle.write('\n')
        except BaseException:
            pending.unlink(missing_ok=True)
            raise
    pending.replace(path)


def checked_file(path, expected=CHECKSUM_UNSET):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), f'Missing or symbolic-link validation input: {path}')
    actual = digest(path)
    if expected is not CHECKSUM_UNSET:
        require(isinstance(expected, str) and len(expected) == 64
                and all(character in '0123456789abcdef' for character in expected),
                f'Missing or malformed declared checksum: {path}')
        require(actual == expected, f'File changed or checksum mismatch: {path}')
    return actual


def collect_state(data_dir, gsea_dir, output_dir):
    """Validate producer manifests and gather exact current dependency hashes."""
    data_dir, gsea_dir, output_dir = map(Path, (data_dir, gsea_dir, output_dir))
    provenance = validate_completed_import(data_dir)
    gsea = validate_completed_run(gsea_dir, data_dir)
    independent = validate_audit_report(data_dir)
    score_audits = {item['rank_metric']: item for item in independent['score_audits']}
    require(score_audits['proxy_signed_log10_p']['manifest_sha256'] == digest(gsea_dir / 'gsea_manifest.json'),
            'Independent ES audit belongs to a different primary GSEA run')
    sensitivity_dir = data_dir / 'gsea/signed_z'
    sensitivity = validate_completed_run(sensitivity_dir, data_dir)
    require(provenance.get('importer_sha256') == digest(REPOSITORY / 'analysis/legacy_workbook.py'),
            'Legacy importer changed after import; regenerate imported results')
    require(provenance.get('tested_universe_confirmation', {}).get('all_tested_genes_exported') is True,
            'The tested-gene universe must be confirmed in import provenance')
    require(gsea.get('rank_metric') == 'proxy_signed_log10_p',
            'Legacy figures must use the primary signed-log-P proxy results')
    render_path = output_dir / 'legacy_plot_audit.json'
    render_bytes = render_path.read_bytes()
    rendered = json.loads(render_bytes)
    require(rendered.get('complete') is True and rendered.get('legacy_exploratory') is True,
            'Legacy rendering is incomplete or is not marked exploratory')
    require(rendered.get('font_pixels') == FONT_PIXELS, 'Legacy renderer declared a different font-size policy')
    require(rendered.get('gsea_primary_analyses') == 6, 'Legacy rendering must use all six primary analyses')
    require(rendered.get('gsea_manifest_sha256') == checked_file(gsea_dir / 'gsea_manifest.json'),
            'Legacy renderer used a different GSEA manifest')
    required_inputs = [REPOSITORY / 'analysis/render_legacy_figures.py',
                       REPOSITORY / 'figures/synthetic_demo/render_figure.py',
                       data_dir / 'legacy_provenance.json', gsea_dir / 'gsea_manifest.json',
                       REPOSITORY / 'analysis/legacy_workbook.py', REPOSITORY / 'analysis/legacy_gsea.py']
    required_inputs += [data_dir / f'{figure}_deg.tsv.gz' for figure, *_ in FIGURES]
    required_inputs += [gsea_dir / f'figure_{number:02d}_{library}_display.tsv'
                        for number in (1, 2, 3) for library in ('kegg', 'go_bp')]
    declared = rendered.get('input_sha256', {})
    require(set(declared) == {saved_path(path) for path in required_inputs},
            'Legacy renderer input manifest is incomplete or unexpected')
    inputs = {saved_path(path): checked_file(path, declared[saved_path(path)]) for path in required_inputs}
    inputs[saved_path(render_path)] = hashlib.sha256(render_bytes).hexdigest()
    for family_dir, family in ((gsea_dir, gsea), (sensitivity_dir, sensitivity)):
        manifest_path = family_dir / 'gsea_manifest.json'
        inputs[saved_path(manifest_path)] = checked_file(manifest_path)
        for filename, checksum in family['artifact_sha256'].items():
            path = family_dir / filename
            inputs[saved_path(path)] = checked_file(path, checksum)
        for source in family['sources'].values():
            for filename, checksum in (('gmt_path', 'gmt_sha256'), ('provenance_path', 'provenance_sha256')):
                path = resolve_saved(source[filename])
                inputs[saved_path(path)] = checked_file(path, source[checksum])
    independent_path = data_dir / 'legacy_gsea_audit.json'
    inputs[saved_path(independent_path)] = checked_file(independent_path)
    code_paths = (Path(__file__), Path(font_audit.__file__), REPOSITORY / 'figures/synthetic_demo/run_gsea.py',
                  REPOSITORY / 'analysis/audit_legacy_gsea.py')
    code = {saved_path(path): checked_file(path) for path in code_paths}
    artifacts = {}
    figures = rendered.get('figures', [])
    require(len(figures) == 3, 'Legacy rendering must describe exactly three ordered figures')
    for number, (figure, expected) in enumerate(zip(figures, FIGURES), 1):
        stem = f'figure_{number:02d}_legacy'
        require(figure.get('figure') == number and figure.get('file_stem') == stem,
                'Legacy figure identifiers are incomplete, duplicated or out of order')
        require(figure.get('comparison') == f'{expected[2]} vs {expected[3]}',
                f'Legacy Figure {number} comparison direction is incorrect')
        require(figure.get('canvas_pixels') == [1536, 1024] and figure.get('font_pixels') == FONT_PIXELS,
                f'Legacy Figure {number} declared dimensions or typography differ')
        require(figure.get('absolute_log2fc_cutoff') == .2 and figure.get('adjusted_p_cutoff') == .05,
                f'Legacy Figure {number} display cutoffs changed')
        require(figure.get('heatmap_panel') is False and figure.get('pca_panel') is False
                and figure.get('sample_expression_available') is False,
                'Unavailable sample-level expression must not be represented as a heatmap or PCA')
        expected_exports = {f'{stem}{suffix}' for suffix in SUFFIXES}
        exported = figure.get('export_sha256', {})
        require(set(exported) == expected_exports, f'Incomplete exports for legacy Figure {number}')
        for filename, checksum in exported.items():
            path = output_dir / filename
            artifacts[saved_path(path)] = checked_file(path, checksum)
        plot_data = output_dir / f'figure_{number:02d}_gene_effects.tsv'
        require(figure.get('plot_data_file') == saved_path(plot_data), 'Unexpected legacy plot-data filename')
        artifacts[saved_path(plot_data)] = checked_file(plot_data, figure.get('plot_data_sha256'))
    require(render_path.read_bytes() == render_bytes, 'Legacy renderer completion changed during validation')
    sources = {f'NimbusSans-{style}.otf': checked_file(font_audit.SOURCE / f'NimbusSans-{style}.otf')
               for style in font_audit.STYLES}
    return {'input_sha256': inputs, 'verification_code_sha256': code,
            'artifact_sha256': artifacts, 'source_font_sha256': sources}, provenance, gsea, rendered


def audit_typography(path):
    """Check every actual PDF span, with no exceptions for pathway labels or math."""
    with fitz.open(path) as document:
        require(len(document) == 1, 'Legacy figure PDF must have exactly one page')
        page = document[0]
        require(abs(page.rect.width / .72 - 1536) < .01 and abs(page.rect.height / .72 - 1024) < .01,
                'Legacy PDF canvas is not 1536 × 1024 native pixels')
        spans = [span for block in page.get_text('dict')['blocks']
                 for line in block.get('lines', []) for span in line['spans'] if span['text'].strip()]
        require(spans, 'Legacy PDF has no searchable text spans')
        for span in spans:
            require(span['font'].startswith('NimbusSans-'), f'Actual PDF font fallback in {span["text"]!r}')
            require(abs(span['size'] / .72 - FONT_PIXELS) < .01,
                    f'Actual PDF text {span["text"]!r} is {span["size"] / .72:.5f}px; expected 32px')
        return {'nonempty_text_spans_checked': len(spans), 'every_span_is_32_native_px': True,
                'native_dpi': 100, 'pdf_points_to_native_pixels': 1 / .72,
                'actual_font_names': sorted({span['font'] for span in spans}),
                'actual_native_pixel_size_counts': dict(Counter(f'{span["size"] / .72:.4f}' for span in spans))}


def verify_report_current(path=DEFAULT_OUTPUT / REPORT_NAME):
    """Cheap freshness check; rerun independent export inspection after any change."""
    path = Path(path)
    report = json.loads(path.read_text())
    require(report.get('status') == 'PASS', 'Legacy validation is not a completed PASS')
    require(report.get('integrity_schema') == 1 and report.get('font_pixels') == FONT_PIXELS,
            'Unsupported legacy verification schema or typography policy')
    paths = report['paths']
    output_dir = resolve_saved(paths['output_dir']).resolve()
    require(output_dir == path.parent.resolve(), 'Legacy validation belongs to a different export directory')
    current, *_ = collect_state(resolve_saved(paths['data_dir']), resolve_saved(paths['gsea_dir']), output_dir)
    require(current == report.get('state_sha256'),
            'Legacy validation is stale: a manifest, source, font, dependency or output changed')
    figures = report.get('figures', [])
    require(len(figures) == 3 and {figure.get('figure') for figure in figures} == {1, 2, 3},
            'Legacy validation does not cover all three figures')
    require(all(figure.get('typography', {}).get('every_span_is_32_native_px') is True for figure in figures),
            'Legacy validation lacks complete uniform-size evidence')
    return report


def run(data_dir=DEFAULT_DATA, gsea_dir=None, output_dir=DEFAULT_OUTPUT):
    data_dir, output_dir = Path(data_dir).resolve(), Path(output_dir).resolve()
    gsea_dir = Path(gsea_dir or data_dir / 'gsea').resolve()
    report_path = output_dir / REPORT_NAME
    write_report(report_path, {'status': 'RUNNING', 'started_at_utc': datetime.now(timezone.utc).isoformat()})
    expected, opened = {}, []
    try:
        before, imported, gsea, rendered = collect_state(data_dir, gsea_dir, output_dir)
        for style in font_audit.STYLES:
            font = TTFont(font_audit.SOURCE / f'NimbusSans-{style}.otf')
            opened.append(font)
            require(font_audit.names(font)['1'] == ['Nimbus Sans'], f'Trusted {style} source is not Nimbus Sans')
            expected[style] = (font, font_audit.converted_glyphs(font))
        figures = []
        for number in (1, 2, 3):
            stem = f'figure_{number:02d}_legacy'
            paths = {suffix: output_dir / f'{stem}{suffix}' for suffix in SUFFIXES}
            figures.append({'figure': number, 'file_stem': stem,
                            'export_sha256': {path.name: digest(path) for path in paths.values()},
                            'typography': audit_typography(paths['.pdf']),
                            'pdf': font_audit.audit_pdf(paths['.pdf'], expected),
                            'outlined_svg': font_audit.audit_outlined_svg(paths['.svg']),
                            'editable_svg': font_audit.audit_editable_svg(paths['_editable.svg'], expected),
                            'png': font_audit.audit_png(paths['.png'], paths['.pdf']),
                            'high_resolution_png': font_audit.audit_png(paths['_300dpi.png'], paths['.pdf'],
                                                                       expected_size=(4608, 3072))})
        after, *_ = collect_state(data_dir, gsea_dir, output_dir)
        require(before == after, 'A legacy source, manifest, font or export changed during verification')
        report = {'status': 'PASS', 'integrity_schema': 1,
                  'verified_at_utc': datetime.now(timezone.utc).isoformat(),
                  'scope': 'Actual legacy exports, producer completion manifests and exact dependency hashes; no synthetic reports are read or written.',
                  'font_pixels': FONT_PIXELS, 'state_sha256': after,
                  'paths': {'data_dir': saved_path(data_dir), 'gsea_dir': saved_path(gsea_dir), 'output_dir': saved_path(output_dir)},
                  'source_workbook_sha256': imported['input_sha256'],
                  'gsea_rank_metric': gsea['rank_metric'], 'gsea_primary_analyses': len(gsea['analyses']),
                  'independent_gsea_audit_complete': True,
                  'independent_gsea_audit_sha256': after['input_sha256'][saved_path(data_dir / 'legacy_gsea_audit.json')],
                  'independent_gsea_audit_scope': 'Weighted ES, leading edges and membership independently recalculated for both rankings; sensitivity comparison completed. Null permutations, NES normalization and empirical FDR are not independently reproduced.',
                  'renderer_complete': rendered['complete'], 'figures': figures,
                  'font_identity_method': 'Actual embedded PDF and editable-SVG font metrics and glyph outlines independently compared with installed URW Nimbus Sans source OTFs; outlined SVG glyph references checked.',
                  'software': {package: version(package) for package in ('PyMuPDF', 'fonttools', 'Pillow', 'numpy')}}
        write_report(report_path, report)
    except BaseException as error:
        write_report(report_path, {'status': 'FAILED', 'verified_at_utc': datetime.now(timezone.utc).isoformat(),
                                   'error': f'{type(error).__name__}: {error}'})
        raise
    finally:
        for font in opened:
            font.close()
    for figure in figures:
        print(f'Figure {figure["figure"]}: PASS; every actual PDF span 32px Nimbus Sans; '
              f'embedded PDF/SVG fonts and both PNG resolutions verified.')
    print(f'Report: {report_path}')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--gsea-dir', type=Path, help='Defaults to DATA_DIR/gsea')
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--check-report', action='store_true')
    args = parser.parse_args()
    if args.check_report:
        verify_report_current(args.output_dir / REPORT_NAME)
        print('PASS: saved legacy export verification matches current manifests, sources, fonts and outputs.')
    else:
        run(args.data_dir, args.gsea_dir, args.output_dir)


if __name__ == '__main__':
    main()
