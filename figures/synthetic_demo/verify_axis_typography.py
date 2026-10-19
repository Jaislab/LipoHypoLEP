#!/usr/bin/env python3
"""Measure exported PDF typography and make a native-scale review comparison.

This inspects actual PDF text spans, independently of the renderer's settings.
The comparison image is a review artifact, separate from manuscript artwork.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
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
DPI = 100
WIDTH, HEIGHT = 1536, 1024
NUMBER = re.compile(r'^[−-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+−-]?\d+)?$')
LABEL_ROLES = {'log2foldchange': 'volcano_x', '−log10(adjustedp)': 'volcano_y',
               'genes': 'heatmap_y', 'z-score': 'heatmap_colorbar'}
EXPECTED_LABEL_ROLES = set(LABEL_ROLES.values()) | {'kegg_x', 'go_bp_x'}
EXPECTED_TICK_ROLES = {'volcano_x', 'volcano_y', 'heatmap_colorbar', 'kegg_x', 'go_bp_x'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
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
    atomic_write(OUTPUT / 'axis_type_measurements.json',
                 (json.dumps(report, indent=2, allow_nan=False) + '\n').encode())


def measure_pdf(data):
    with fitz.open(stream=data, filetype='pdf') as document:
        require(len(document) == 1, 'Expected a one-page figure PDF')
        page = document[0]
        scale = DPI / 72
        require(abs(page.rect.width * scale - WIDTH) < .01 and
                abs(page.rect.height * scale - HEIGHT) < .01,
                'PDF canvas differs from the requested 1536 × 1024 native pixels')
        lines = []
        for block in page.get_text('dict')['blocks']:
            for line in block.get('lines', []):
                spans = [span for span in line['spans'] if span['text'].strip()]
                text = ''.join(span['text'] for span in line['spans']).strip()
                if not text or not spans:
                    continue
                require(all(span['font'].startswith('NimbusSans-') for span in spans),
                        f'Unexpected actual PDF font in {text!r}')
                box = list(line['bbox'])
                lines.append({'text': text, 'bbox_pdf_pt': box,
                              'center_x': (box[0] + box[2]) / 2,
                              'center_y': (box[1] + box[3]) / 2,
                              'span_sizes_pdf_pt': [span['size'] for span in spans],
                              'font_names': sorted({span['font'] for span in spans})})
        labels = {}
        for line in lines:
            normalized = re.sub(r'\s+', '', line['text']).casefold()
            role = LABEL_ROLES.get(normalized)
            if normalized == 'nes':
                role = 'kegg_x' if line['center_x'] < page.rect.width / 2 else 'go_bp_x'
            if role:
                require(role not in labels, f'Ambiguous axis title: {role}')
                labels[role] = line
        require(set(labels) == EXPECTED_LABEL_ROLES, f'Missing PDF axis titles: {EXPECTED_LABEL_ROLES-set(labels)}')
        numeric = [line for line in lines if NUMBER.fullmatch(line['text'])]
        top_left = [line for line in numeric if line['center_x'] < page.rect.width / 2
                    and line['center_y'] < page.rect.height / 2]
        require(top_left, 'No upper-left numeric tick labels')
        lowest = max(line['center_y'] for line in top_left)
        volcano_x = [line for line in top_left if abs(line['center_y'] - lowest) < 1]
        ticks = {'volcano_x': volcano_x,
                 'volcano_y': [line for line in top_left if line not in volcano_x],
                 'heatmap_colorbar': [line for line in numeric
                                     if line['center_x'] > page.rect.width / 2
                                     and line['center_y'] < page.rect.height / 2]}
        for role, left_half in [('kegg_x', True), ('go_bp_x', False)]:
            candidates = [line for line in numeric if line['center_y'] > page.rect.height / 2
                          and (line['center_x'] < page.rect.width / 2) == left_half]
            require(candidates, f'No lower-panel numeric labels for {role}')
            lowest = max(line['center_y'] for line in candidates)
            # q values sit beside pathway rows above the bottom NES tick row.
            ticks[role] = [line for line in candidates if abs(line['center_y'] - lowest) < 1]
        for role, entries in ticks.items():
            require(len(entries) >= 3, f'Too few spatially identified numeric ticks: {role}')
            for entry in entries:
                require(len(entry['span_sizes_pdf_pt']) == 1,
                        f'Unexpected multi-size numeric tick: {entry["text"]!r}')
        def public(line):
            result = {key: value for key, value in line.items() if not key.startswith('center_')}
            # Math subscripts retain their intended smaller size; ordinary title
            # glyphs determine the full axis-label size.
            result['full_size_pdf_pt'] = max(line['span_sizes_pdf_pt'])
            result['full_size_native_px'] = result['full_size_pdf_pt'] * scale
            result['span_sizes_native_px'] = [value * scale for value in line['span_sizes_pdf_pt']]
            return result
        pixmap = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=False, colorspace=fitz.csRGB)
        image = Image.frombytes('RGB', (pixmap.width, pixmap.height), pixmap.samples)
    return {'pdf_sha256': digest(data), 'canvas_native_pixels': list(image.size),
            'axis_titles': {role: public(line) for role, line in labels.items()},
            'numeric_ticks': {role: [public(line) for line in entries] for role, entries in ticks.items()}}, image


def verify_unchanged_data(reference):
    roots = ('figures/synthetic_demo/data', 'figures/synthetic_demo/gsea_data')
    files = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', reference, *roots],
                                    cwd=REPOSITORY, text=True).splitlines()
    require(files, 'No reference analysis data found')
    current_files = {str(path.relative_to(REPOSITORY)) for root in roots
                     for path in (REPOSITORY / root).rglob('*')
                     if path.is_file() and path.name != '.gsea.lock'}
    require(set(files) == current_files, 'Analysis-data file inventory changed')
    evidence = {}
    for filename in files:
        previous = digest(git_bytes(reference, filename))
        current = digest((REPOSITORY / filename).read_bytes())
        require(previous == current, f'Analysis data changed: {filename}')
        evidence[filename] = current
    return evidence


def comparison_image(before, after, before_measurement, after_measurement):
    # One identical crop rectangle at 100 DPI: each source pixel remains one
    # montage pixel. Include complete volcano axis titles and numeric ticks.
    y0 = max(0, int(np.floor(min(before_measurement['axis_titles']['volcano_y']['bbox_pdf_pt'][1],
                                after_measurement['axis_titles']['volcano_y']['bbox_pdf_pt'][1]) * DPI / 72)) - 6)
    y1 = int(np.ceil(max(before_measurement['axis_titles']['volcano_x']['bbox_pdf_pt'][3],
                         after_measurement['axis_titles']['volcano_x']['bbox_pdf_pt'][3]) * DPI / 72)) + 2
    box = (0, y0, 755, y1)
    crops = [image.crop(box) for image in (before, after)]
    width, height = crops[0].size
    gap, header = 24, 108
    canvas_width, canvas_height = 2 * width + 3 * gap, height + header + gap
    font_dir = Path(os.environ.get('NIMBUS_SANS_FONT_DIR', '/usr/share/fonts/opentype/urw-base35'))
    regular, bold = font_dir / 'NimbusSans-Regular.otf', font_dir / 'NimbusSans-Bold.otf'
    require(regular.is_file() and bold.is_file(), 'Nimbus Sans is required for review headings')
    fig = plt.figure(figsize=(canvas_width / DPI, canvas_height / DPI), dpi=DPI, facecolor='white')
    for index, (crop, measurement, heading) in enumerate(zip(
            crops, (before_measurement, after_measurement), ('Previous PDF', 'Revised PDF'))):
        x = gap + index * (width + gap)
        fig.figimage(np.asarray(crop), xo=x, yo=gap, origin='upper', resize=False)
        label = measurement['axis_titles']['volcano_x']['full_size_native_px']
        tick = measurement['numeric_ticks']['volcano_x'][0]['full_size_native_px']
        fig.text(x / canvas_width, 1 - 22 / canvas_height, heading, va='top',
                 fontproperties=fm.FontProperties(fname=str(bold), size=28 * 72 / DPI), color='#262A2B')
        fig.text(x / canvas_width, 1 - 63 / canvas_height,
                 f'Axis labels {label:g} px  ·  numeric ticks {tick:g} px', va='top',
                 fontproperties=fm.FontProperties(fname=str(regular), size=24 * 72 / DPI), color='#525858')
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=DPI, facecolor='white')
    plt.close(fig)
    data = buffer.getvalue()
    with Image.open(io.BytesIO(data)) as montage:
        require(montage.size == (canvas_width, canvas_height), 'Unexpected review canvas dimensions')
        for index, crop in enumerate(crops):
            x = gap + index * (width + gap)
            pasted = montage.crop((x, header, x + width, header + height)).convert('RGB')
            require(np.array_equal(np.asarray(pasted), np.asarray(crop)),
                    'Review image resampled or altered a source crop')
    return data, {'source_figure': 1, 'crop_box_native_pixels': list(box),
                  'source_rasterization_dpi': DPI, 'source_pixel_to_review_pixel': '1:1',
                  'source_crops_pixel_identical_in_montage': True,
                  'review_canvas_pixels': [canvas_width, canvas_height]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', default='3e8d987', help='Previous figure export Git revision')
    parser.add_argument('--expected-label-px', type=float, default=40)
    parser.add_argument('--expected-tick-px', type=float, default=32)
    args = parser.parse_args()
    save_report({'status': 'RUNNING'})
    try:
        reference = subprocess.check_output(['git', 'rev-parse', args.reference], cwd=REPOSITORY, text=True).strip()
        figures, images, current_digests = [], [], {}
        for number in (1, 2, 3):
            filename = f'figures/figure_{number:02d}_gsea_synthetic.pdf'
            previous, before_image = measure_pdf(git_bytes(reference, filename))
            current_data = (REPOSITORY / filename).read_bytes()
            current, after_image = measure_pdf(current_data)
            current_digests[filename] = digest(current_data)
            for role, title in current['axis_titles'].items():
                require(abs(title['full_size_native_px'] - args.expected_label_px) < .01,
                        f'Figure {number} {role}: expected {args.expected_label_px} px title')
            for role, ticks in current['numeric_ticks'].items():
                require(all(abs(tick['full_size_native_px'] - args.expected_tick_px) < .01 for tick in ticks),
                        f'Figure {number} {role}: expected {args.expected_tick_px} px ticks')
            figures.append({'figure': number, 'previous': previous, 'current': current})
            if number == 1:
                images = [before_image, after_image]
        unchanged = verify_unchanged_data(reference)
        png, crop = comparison_image(*images, figures[0]['previous'], figures[0]['current'])
        for filename, checksum in current_digests.items():
            require(digest((REPOSITORY / filename).read_bytes()) == checksum,
                    f'Export changed during measurement: {filename}')
        atomic_write(OUTPUT / 'axis_type_before_after.png', png)
        report = {'status': 'PASS', 'verified_at_utc': datetime.now(timezone.utc).isoformat(),
                  'reference_commit': reference,
                  'current_version': 'Current working-tree PDF exports identified by SHA-256, measured before publication commit.',
                  'measurement_method': 'Actual searchable PDF text spans from PyMuPDF; no renderer-reported font sizes used.',
                  'native_raster_dpi': DPI, 'pdf_points_to_native_pixels': DPI / 72,
                  'axis_title_rule': 'Largest actual span in each full label; smaller math subscripts reported separately.',
                  'numeric_tick_rule': 'Upper-left volcano rows, upper-right colorbar, bottom numeric row in each GSEA half; q-value rows excluded.',
                  'all_measured_text_fonts': 'Nimbus Sans',
                  'expected_axis_label_native_px': args.expected_label_px,
                  'expected_numeric_tick_native_px': args.expected_tick_px,
                  'figures': figures, 'analysis_data_unchanged': True,
                  'unchanged_analysis_file_count': len(unchanged), 'unchanged_analysis_sha256': unchanged,
                  'review_image': {**crop, 'file': 'axis_type_before_after.png', 'sha256': digest(png)},
                  'verifier_sha256': digest(Path(__file__).read_bytes())}
        save_report(report)
        print(f'PASS: actual PDF titles {args.expected_label_px:g}px and numeric ticks {args.expected_tick_px:g}px in all three figures; '
              f'{len(unchanged)} analysis files byte-identical to {args.reference}.')
        print(f'Native-scale comparison: {OUTPUT / "axis_type_before_after.png"}')
    except BaseException as error:
        save_report({'status': 'FAILED', 'error': f'{type(error).__name__}: {error}'})
        raise


if __name__ == '__main__':
    main()
