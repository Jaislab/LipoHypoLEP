#!/usr/bin/env python3
"""Render and independently inspect supplied-gene figures, with 32px Nimbus Sans."""
from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from analysis import analyze_male_2026_09_30 as analysis
from figures.synthetic_demo import render_figure as style
from figures.synthetic_demo import verify_export_fonts as fonts_audit
import fitz
from fontTools.ttLib import TTFont
import matplotlib
from matplotlib import font_manager as fm
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

OUTPUT = ROOT / 'figures/male_2026_09_30'
SUFFIXES = ('.png', '_300dpi.png', '.pdf', '.svg', '_editable.svg')


def state(data_dir):
    manifest = analysis.validate_results(data_dir)
    paths = [Path(__file__), Path(analysis.__file__), Path(style.__file__), Path(fonts_audit.__file__),
             data_dir / 'analysis_manifest.json']
    paths += [data_dir / name for name in manifest['artifact_sha256']]
    paths += [fonts_audit.SOURCE / f'NimbusSans-{face}.otf' for face in fonts_audit.STYLES]
    return {analysis.saved_path(p): analysis.digest(p) for p in paths}


def audit_typography(pdf_path):
    with fitz.open(pdf_path) as pdf:
        analysis.require(len(pdf) == 1, 'Figure must have one PDF page')
        analysis.require(abs(pdf[0].rect.width / .72 - 1536) < .01
                         and abs(pdf[0].rect.height / .72 - 1024) < .01, 'Wrong PDF canvas')
        spans = [s for b in pdf[0].get_text('dict')['blocks'] for line in b.get('lines', [])
                 for s in line['spans'] if s['text'].strip()]
        analysis.require(spans and all(s['font'].startswith('NimbusSans-')
                         and abs(s['size'] / .72 - 32) < .01 for s in spans), 'PDF font fallback or mixed text sizes')
        return {'every_span_32_native_px': True, 'nonempty_spans': len(spans),
                'fonts': sorted({s['font'] for s in spans})}


def export(fig, number, font_paths, output_dir):
    stem = f'figure_{number:02d}_male_2026_09_30'
    for _ in range(2):
        for text in fig.findobj(Text):
            bold = text.get_fontweight() in ('bold', 700) or text.get_fontproperties().get_file() == str(font_paths['Bold'])
            text.set_fontproperties(fm.FontProperties(fname=str(font_paths['Bold' if bold else 'Regular']), size=style.font_points()))
            text.set_parse_math(False)
        fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            analysis.require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= 1536 and box.y1 <= 1024,
                             f'Text outside canvas: {text.get_text()}')
    metadata = {'Title': f'Figure {number} | Male DE exports | 2026-09-30',
                'Description': 'Confirmed filtered male DE lists. See separate legend for export scope, cutoffs and missing enrichment inputs.'}
    fig.savefig(output_dir / f'{stem}.png', dpi=100, metadata=metadata)
    fig.savefig(output_dir / f'{stem}_300dpi.png', dpi=300, metadata=metadata)
    fig.savefig(output_dir / f'{stem}.pdf', metadata={'Title': metadata['Title'], 'Subject': metadata['Description']})
    fig.savefig(output_dir / f'{stem}.svg', metadata=metadata)
    with matplotlib.rc_context({'svg.fonttype': 'none'}):
        fig.savefig(output_dir / f'{stem}_editable.svg', metadata=metadata)
    css = []
    for face, weight, slant in [('Regular', 400, 'normal'), ('Bold', 700, 'normal'), ('Italic', 400, 'italic')]:
        encoded = base64.b64encode(font_paths[face].read_bytes()).decode()
        css.append(f'@font-face {{ font-family: "Nimbus Sans"; font-weight: {weight}; font-style: {slant}; src: url("data:font/ttf;base64,{encoded}") format("truetype"); }}')
    path = output_dir / f'{stem}_editable.svg'
    path.write_text(path.read_text().replace('<defs>', '<defs>\n<style type="text/css"><![CDATA[\n' + '\n'.join(css) + '\n]]></style>', 1))
    for suffix in ('.svg', '_editable.svg'):
        path = output_dir / f'{stem}{suffix}'
        path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
    style.plt.close(fig)
    return stem


def draw(number, first, second, table, output_dir, font_paths):
    fig = style.plt.figure(figsize=(15.36, 10.24), dpi=100)
    colors = style.COMPARISON_COLORS[number]
    for letter, left in [('a', 40), ('b', 788)]:
        fig.text(left / 1536, 1 - 27 / 1024, letter, ha='left', va='top', weight='bold')
    fig.text(400 / 1536, 1 - 64 / 1024, f'{first} vs {second}', ha='center', va='baseline')
    fig.text(1180 / 1536, 1 - 64 / 1024, 'Genes', ha='center', va='baseline')
    a = fig.add_axes([130 / 1536, 184 / 1024, 540 / 1536, 700 / 1024])
    hit = analysis.joint_hits(table)
    values = -np.log10(table.padj)
    a.scatter(table.loc[~hit, 'log2FoldChange'], values[~hit], s=24, c=style.NS, alpha=.75, edgecolors='none')
    for direction, color in [(-1, colors['lower']), (1, colors['higher'])]:
        mask = hit & (np.sign(table.log2FoldChange) == direction)
        a.scatter(table.loc[mask, 'log2FoldChange'], values[mask], s=28, c=color, alpha=.9, edgecolors='none')
    a.axhline(-np.log10(.05), color=style.AXIS, lw=.8, ls=(0, (4, 4)))
    limit = max(1., float(np.ceil(table.log2FoldChange.abs().max())))
    a.set(xlim=(-limit, limit), ylim=(0, float(np.ceil(values.max() / 2) * 2)),
          xlabel='log2 fold change', ylabel='−log10(adjusted P)')
    xticks = MaxNLocator(nbins=5, symmetric=True).tick_values(-limit, limit)
    a.set_xticks(xticks[(xticks >= -limit) & (xticks <= limit)])
    upper = a.get_ylim()[1]
    yticks = MaxNLocator(nbins=5, integer=True).tick_values(0, upper)
    a.set_yticks(yticks[(yticks >= 0) & (yticks <= upper)])
    selected = analysis.selected_genes(table)
    b = fig.add_axes([1110 / 1536, 184 / 1024, 330 / 1536, 700 / 1024])
    names = selected.symbol.tolist()
    if len(set(names)) != len(names):
        names = [f'{r.symbol} [..{r.gene_id[-6:]}]' for r in selected.itertuples()]
    b.scatter(selected.log2FoldChange, np.arange(len(selected)), s=110,
              c=[colors['higher'] if x > 0 else colors['lower'] for x in selected.log2FoldChange], edgecolors='none')
    b.set_yticks(np.arange(len(selected)), labels=names)
    b.set_ylim(max(len(selected) - .5, .5), -.5)
    effect_limit = max(.4, float(np.ceil(selected.log2FoldChange.abs().max() * 1.12 / .2) * .2))
    b.set(xlim=(-effect_limit, effect_limit), xlabel='log2 fold change')
    ticks = MaxNLocator(nbins=4, symmetric=True).tick_values(-effect_limit, effect_limit)
    b.set_xticks(ticks[(ticks >= -effect_limit) & (ticks <= effect_limit)])
    b.spines['left'].set_visible(False)
    b.tick_params(axis='y', length=0, pad=14)
    for ax in (a, b):
        ax.axvline(0, color=style.GRID, lw=.65, ls=(0, (3, 4)), zorder=0)
        for cutoff in (-.2, .2):
            ax.axvline(cutoff, color=style.AXIS, lw=.8, ls=(0, (4, 4)), zorder=0)
        ax.tick_params(direction='out', pad=6)
        ax.grid(False)
    fig.canvas.draw()
    for label in b.get_yticklabels():
        analysis.require(label.get_window_extent(fig.canvas.get_renderer()).x0 >= 788, 'Gene label crosses panels')
    stem = export(fig, number, font_paths, output_dir)
    return {'figure': number, 'file_stem': stem, 'comparison': f'{first} vs {second}',
            'supplied_rows_only': True, 'enrichment_performed': False, 'sample_expression_inferred': False,
            'volcano_absolute_x_limit': limit, 'gene_effect_absolute_x_limit': effect_limit,
            'selected_gene_ids': selected.gene_id.tolist(), 'selected_effects': selected.log2FoldChange.tolist()}


def check_report(data_dir=analysis.OUTPUT, output_dir=OUTPUT):
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    report = json.loads((output_dir / 'figure_validation.json').read_text())
    analysis.require(report.get('complete') is True and report.get('status') == 'PASS', 'Figure verification is incomplete')
    analysis.require(state(data_dir) == report['input_sha256'], 'Figure inputs, code or fonts changed')
    analysis.require(len(report.get('figures', [])) == 3 and len(report.get('artifact_sha256', {})) == 15, 'Incomplete export inventory')
    for name, checksum in report['artifact_sha256'].items():
        analysis.require(analysis.digest(output_dir / name) == checksum, f'Figure changed: {name}')
    print('PASS: current male figure exports, fonts, code and analysis match saved verification.')
    return report


def run(data_dir=analysis.OUTPUT, output_dir=OUTPUT):
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / 'figure_validation.json'
    analysis.write_json(report_path, {'complete': False})
    before = state(data_dir)
    style.configure_theme()
    font_paths = style.prepare_nimbus()
    expected = {}
    for face in fonts_audit.STYLES:
        font = TTFont(fonts_audit.SOURCE / f'NimbusSans-{face}.otf')
        expected[face] = (font, fonts_audit.converted_glyphs(font))
    figures = []
    for number, first, second, _, _ in analysis.SPECS:
        table = pd.read_csv(data_dir / f'figure_{number:02d}_supplied.tsv', sep='\t', float_precision='round_trip')
        record = draw(number, first, second, table, output_dir, font_paths)
        stem = record['file_stem']
        pdf = output_dir / f'{stem}.pdf'
        record.update({'typography': audit_typography(pdf), 'pdf': fonts_audit.audit_pdf(pdf, expected),
                       'outlined_svg': fonts_audit.audit_outlined_svg(output_dir / f'{stem}.svg'),
                       'editable_svg': fonts_audit.audit_editable_svg(output_dir / f'{stem}_editable.svg', expected),
                       'png': fonts_audit.audit_png(output_dir / f'{stem}.png', pdf),
                       'high_resolution_png': fonts_audit.audit_png(output_dir / f'{stem}_300dpi.png', pdf, (4608, 3072))})
        figures.append(record)
    for font, _ in expected.values():
        font.close()
    analysis.require(before == state(data_dir), 'Figure input changed during export verification')
    report = {'complete': True, 'status': 'PASS', 'font_pixels': 32, 'canvas_pixels': [1536, 1024],
              'verified_at_utc': datetime.now(timezone.utc).isoformat(), 'input_sha256': before,
              'figures': figures,
              'artifact_sha256': {f'{r["file_stem"]}{s}': analysis.digest(output_dir / f'{r["file_stem"]}{s}')
                                 for r in figures for s in SUFFIXES}}
    analysis.write_json(report_path, report)
    check_report(data_dir, output_dir)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=analysis.OUTPUT)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    parser.add_argument('--check-report', action='store_true')
    args = parser.parse_args()
    (check_report if args.check_report else run)(args.data_dir, args.output_dir)
