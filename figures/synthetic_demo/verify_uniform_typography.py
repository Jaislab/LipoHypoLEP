#!/usr/bin/env python3
"""Independently enforce one Nimbus Sans size on every exported PDF text span.

The review montage uses unchanged native-resolution crops from actual PDFs.
It is separate from the manuscript artwork. No renderer audit claims are used.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import io
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

os.environ.setdefault('MPLCONFIGDIR', '/tmp/lipohypolep-mpl-v2')
import fitz
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

REPOSITORY = Path(__file__).resolve().parents[2]
OUTPUT = REPOSITORY / 'figures' / 'typography'
REPORT = OUTPUT / 'uniform_type_measurements.json'
DPI = 100
CANVAS = (1536, 1024)
DATA_ROOTS = ('figures/synthetic_demo/data', 'figures/synthetic_demo/gsea_data')
FONT_DIR = Path(os.environ.get('NIMBUS_SANS_FONT_DIR', '/usr/share/fonts/opentype/urw-base35'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def git_bytes(reference, path):
    return subprocess.check_output(['git', 'show', f'{reference}:{path}'], cwd=REPOSITORY)


def atomic_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('wb', dir=path.parent, delete=False) as handle:
        pending = Path(handle.name)
        try:
            handle.write(data)
        except BaseException:
            pending.unlink(missing_ok=True)
            raise
    pending.replace(path)


def save_report(report):
    atomic_write(REPORT, (json.dumps(report, indent=2, allow_nan=False) + '\n').encode())


def current_data_files():
    return {str(path.relative_to(REPOSITORY)) for root in DATA_ROOTS
            for path in (REPOSITORY / root).rglob('*')
            if path.is_file() and path.name != '.gsea.lock'}


def verify_unchanged_data(reference):
    files = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', reference, *DATA_ROOTS],
                                    cwd=REPOSITORY, text=True).splitlines()
    require(files and set(files) == current_data_files(), 'Analysis-data file inventory changed')
    evidence = {}
    for filename in files:
        current = sha256((REPOSITORY / filename).read_bytes())
        require(sha256(git_bytes(reference, filename)) == current, f'Analysis data changed: {filename}')
        evidence[filename] = current
    return evidence


def measure_pdf(data, expected_px=None):
    """Visit every nonempty actual span, including pathway labels and q values."""
    with fitz.open(stream=data, filetype='pdf') as document:
        require(len(document) == 1, 'Expected a one-page figure PDF')
        page = document[0]
        scale = DPI / 72
        require(abs(page.rect.width * scale - CANVAS[0]) < .01 and
                abs(page.rect.height * scale - CANVAS[1]) < .01,
                'PDF canvas differs from 1536 × 1024 native pixels')
        spans, lines = [], []
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                nonempty = [span for span in line['spans'] if span['text'].strip()]
                if nonempty:
                    lines.append(''.join(span['text'] for span in line['spans']).strip())
                    spans.extend(nonempty)
        require(spans, 'PDF contains no searchable text spans')
        for span in spans:
            require(span['font'].startswith('NimbusSans-'),
                    f'Unexpected PDF font {span["font"]!r} in {span["text"]!r}')
            native_size = span['size'] * scale
            if expected_px is not None:
                require(abs(native_size - expected_px) < .01,
                        f'Actual PDF span {span["text"]!r} is {native_size:.5f}px; '
                        f'every span must be {expected_px:g}px')
        normalized = [re.sub(r'\s+', '', line).casefold() for line in lines]
        critical = ('log2foldchange', '−log10(adjustedp)', 'genes', 'z-score',
                    'nes', 'kegg', 'gobp', 'q', 'a', 'b', 'c', 'd')
        require(all(value in normalized for value in critical), 'An essential plot label is missing')
        require(normalized.count('nes') == 2 and normalized.count('q') == 2,
                'Expected both KEGG and GO BP text labels')
        sizes = Counter(f'{span["size"] * scale:.4f}'.rstrip('0').rstrip('.') for span in spans)
        fonts = Counter(span['font'] for span in spans)
        regions = Counter(('upper_panels' if (span['bbox'][1] + span['bbox'][3]) / 2 < page.rect.height / 2
                           else 'lower_KEGG' if (span['bbox'][0] + span['bbox'][2]) / 2 < page.rect.width / 2
                           else 'lower_GO_BP') for span in spans)
        measurements = {'pdf_sha256': sha256(data), 'nonempty_text_spans_checked': len(spans),
                        'all_span_sizes_native_px_and_counts': dict(sorted(sizes.items(), key=lambda pair: float(pair[0]))),
                        'all_span_sizes_pdf_pt_and_counts': dict(sorted(Counter(f'{span["size"]:.6f}' for span in spans).items())),
                        'actual_font_name_counts': dict(fonts), 'span_counts_by_page_region': dict(regions),
                        'essential_text_labels_present': list(critical),
                        'text_without_whitespace_sha256': sha256(''.join(normalized).encode())}
        if expected_px is not None:
            measurements['every_nonempty_span_matches_expected_size'] = True
        pixmap = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False, colorspace=fitz.csRGB)
        image = Image.frombytes('RGB', (pixmap.width, pixmap.height), pixmap.samples)
    return measurements, image


def make_review(before, after, before_measurement, expected_px):
    # Stack both complete GSEA rows at native scale. No crop is resampled.
    box = (0, 496, CANVAS[0], CANVAS[1])
    crops = [image.crop(box) for image in (before, after)]
    width, height = crops[0].size
    margin, header = 24, 90
    canvas_width = width + 2 * margin
    canvas_height = 2 * (height + header) + 3 * margin
    bold, regular = FONT_DIR / 'NimbusSans-Bold.otf', FONT_DIR / 'NimbusSans-Regular.otf'
    require(bold.is_file() and regular.is_file(), 'Nimbus Sans source fonts are required')
    for path in (bold, regular):
        require(fm.FontProperties(fname=str(path)).get_name() == 'Nimbus Sans', 'Unexpected review font family')
    fig = plt.figure(figsize=(canvas_width / DPI, canvas_height / DPI), dpi=DPI, facecolor='white')
    locations = []
    for index, crop in enumerate(crops):
        header_top = margin + index * (header + height + margin)
        crop_top = header_top + header
        locations.append((margin, crop_top, margin + width, crop_top + height))
        fig.figimage(np.asarray(crop), xo=margin, yo=canvas_height - crop_top - height,
                     origin='upper', resize=False)
        heading = 'Previous PDF · mixed sizes' if index == 0 else 'Revised PDF · one font size'
        previous_sizes = ', '.join(before_measurement['all_span_sizes_native_px_and_counts'])
        detail = (f'Measured sizes across the figure: {previous_sizes} px' if index == 0
                  else f'Every text span: {expected_px:g} px · Nimbus Sans')
        for text, offset, font in ((heading, 0, bold), (detail, 43, regular)):
            fig.text(margin / canvas_width, 1 - (header_top + offset) / canvas_height, text, va='top',
                     fontproperties=fm.FontProperties(fname=str(font), size=32 * 72 / DPI), color='#262A2B')
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=DPI, facecolor='white')
    plt.close(fig)
    data = buffer.getvalue()
    with Image.open(io.BytesIO(data)) as montage:
        require(montage.size == (canvas_width, canvas_height), 'Unexpected review-image dimensions')
        for crop, location in zip(crops, locations):
            require(np.array_equal(np.asarray(montage.crop(location).convert('RGB')), np.asarray(crop)),
                    'Review crop changed or was resampled')
    return data, {'source_figure': 1, 'source_crop_native_pixels': list(box),
                  'source_pixel_to_review_pixel': '1:1', 'source_crops_pixel_identical': True,
                  'review_canvas_pixels': [canvas_width, canvas_height]}


def verify_report_current(path=REPORT, expected_px=32):
    """Fail closed if any current export, analysis input, font or proof changed."""
    report = json.loads(path.read_text())
    require(report.get('status') == 'PASS', 'Uniform-typography report is not a completed PASS')
    require(report.get('expected_native_pixel_size') == expected_px, 'Saved uniform-size policy differs from the required size')
    require(report.get('verifier_sha256') == sha256(Path(__file__).read_bytes()),
            'Uniform-typography verifier changed; rerun it')
    figures = report.get('figures', [])
    require(len(figures) == 3 and {entry.get('figure') for entry in figures} == {1, 2, 3},
            'Incomplete current PDF evidence')
    for entry in figures:
        filename = f'figures/figure_{entry["figure"]:02d}_gsea_synthetic.pdf'
        require(sha256((REPOSITORY / filename).read_bytes()) == entry['current']['pdf_sha256'],
                f'Current PDF changed: {filename}')
        require(entry['current'].get('every_nonempty_span_matches_expected_size') is True,
                'Missing completed whole-PDF size verification')
        histogram = entry['current'].get('all_span_sizes_native_px_and_counts', {})
        require(histogram and all(abs(float(size)-expected_px) < .01 for size in histogram),
                'Saved span-size evidence is not uniform')
        require(sum(histogram.values()) == entry['current']['nonempty_text_spans_checked'],
                'Saved span-size evidence is incomplete')
    data = report.get('unchanged_analysis_sha256', {})
    require(data and set(data) == current_data_files(), 'Analysis-data evidence is incomplete or stale')
    for filename, checksum in data.items():
        require(sha256((REPOSITORY / filename).read_bytes()) == checksum, f'Analysis data changed: {filename}')
    for filename, checksum in report.get('review_font_sha256', {}).items():
        require(sha256((FONT_DIR / filename).read_bytes()) == checksum, f'Review font changed: {filename}')
    require(set(report.get('review_font_sha256', {})) == {'NimbusSans-Regular.otf', 'NimbusSans-Bold.otf'},
            'Incomplete review-font evidence')
    require(sha256((OUTPUT / 'uniform_type_before_after.png').read_bytes()) == report['review_image']['sha256'],
            'Uniform-typography review image changed')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', default='46faeea', help='Previous mixed-size figure Git revision')
    parser.add_argument('--expected-px', type=float, default=32, help='Required size of every nonempty PDF span')
    parser.add_argument('--check-report', action='store_true', help='Check saved evidence against current files')
    args = parser.parse_args()
    if args.check_report:
        verify_report_current()
        print('PASS: uniform-typography evidence matches current PDFs, analysis data, fonts and proof.')
        return
    save_report({'status': 'RUNNING'})
    try:
        reference = subprocess.check_output(['git', 'rev-parse', args.reference], cwd=REPOSITORY, text=True).strip()
        figures, images = [], []
        for number in (1, 2, 3):
            filename = f'figures/figure_{number:02d}_gsea_synthetic.pdf'
            previous, before = measure_pdf(git_bytes(reference, filename))
            current, after = measure_pdf((REPOSITORY / filename).read_bytes(), args.expected_px)
            figures.append({'figure': number, 'previous': previous, 'current': current})
            if number == 1:
                images = [before, after]
        unchanged = verify_unchanged_data(reference)
        png, crop = make_review(*images, figures[0]['previous'], args.expected_px)
        atomic_write(OUTPUT / 'uniform_type_before_after.png', png)
        report = {'status': 'PASS', 'verified_at_utc': datetime.now(timezone.utc).isoformat(),
                  'reference_commit': reference, 'current_version': 'Working-tree PDF exports identified by SHA-256',
                  'measurement_method': 'Every nonempty actual PDF text span read independently with PyMuPDF; renderer claims are not used.',
                  'scope': 'All figure text, including pathway terms, q values, headings, panel letters, group names, axes and numeric ticks; no font-size exceptions.',
                  'expected_native_pixel_size': args.expected_px, 'native_dpi': DPI,
                  'pdf_points_to_native_pixels': DPI / 72, 'canvas_native_pixels': list(CANVAS),
                  'figures': figures, 'analysis_data_unchanged': True,
                  'unchanged_analysis_file_count': len(unchanged), 'unchanged_analysis_sha256': unchanged,
                  'review_image': {**crop, 'file': 'uniform_type_before_after.png', 'sha256': sha256(png)},
                  'review_font_sha256': {f'NimbusSans-{style}.otf': sha256((FONT_DIR / f'NimbusSans-{style}.otf').read_bytes())
                                         for style in ('Regular', 'Bold')},
                  'software': {package: version(package) for package in ('PyMuPDF', 'matplotlib', 'numpy', 'Pillow')},
                  'verifier_sha256': sha256(Path(__file__).read_bytes())}
        # Recheck mutable sources before exposing completed evidence.
        for entry in figures:
            filename = f'figures/figure_{entry["figure"]:02d}_gsea_synthetic.pdf'
            require(sha256((REPOSITORY / filename).read_bytes()) == entry['current']['pdf_sha256'],
                    f'PDF changed during verification: {filename}')
        require(all(sha256((REPOSITORY / filename).read_bytes()) == checksum for filename, checksum in unchanged.items()),
                'Analysis data changed during verification')
        save_report(report)
        print(f'PASS: every nonempty actual PDF text span is Nimbus Sans at {args.expected_px:g}px in all three figures; '
              f'{len(unchanged)} analysis files unchanged.')
        print(f'Review image: {OUTPUT / "uniform_type_before_after.png"}')
    except BaseException as error:
        save_report({'status': 'FAILED', 'error': f'{type(error).__name__}: {error}'})
        raise


if __name__ == '__main__':
    main()
