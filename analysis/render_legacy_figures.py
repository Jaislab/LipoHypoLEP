#!/usr/bin/env python3
"""Render the legacy workbook's supplied estimates and exploratory enrichment.

No sample expression or confidence intervals are inferred. Methods, filtering
limitations and the proxy GSEA ranking belong in the separate legacy legends.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import json
from pathlib import Path
import re
import sys
import tempfile

REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))

from figures.synthetic_demo import render_figure as style

import fitz
import matplotlib
from matplotlib import font_manager as fm
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from PIL import Image

W, H, DPI = style.W, style.H, style.DPI
PRIMARY_RANK_METRIC = 'proxy_signed_log10_p'
SPECS = tuple((number, question, first, second, f"figure_{number:02d}_deg.tsv.gz")
              for number, question, first, second, _ in style.SPECS)


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


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', encoding='utf-8', dir=path.parent,
                                     suffix='.json', delete=False) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write('\n')
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(path)


def read_contrast(path):
    table = pd.read_csv(path, sep='\t', dtype={'gene_id': str, 'supplied_symbol': str},
                        float_precision='round_trip')
    required = {'gene_id', 'supplied_symbol', 'log2FoldChange', 'pvalue', 'padj'}
    require(required.issubset(table), f"Missing legacy contrast columns in {path}: {sorted(required-set(table))}")
    for column in ('gene_id', 'supplied_symbol'):
        require(table[column].notna().all() and table[column].map(lambda x: bool(x.strip())).all(),
                f"Missing {column} in {path}")
    require(len(table) and table.gene_id.is_unique, f"Empty or duplicate legacy Ensembl rows in {path}")
    for column in ('log2FoldChange', 'pvalue', 'padj'):
        require(pd.api.types.is_numeric_dtype(table[column]) and np.isfinite(table[column]).all(),
                f"Invalid legacy {column} in {path}")
    require(table.pvalue.between(0, 1).all() and table.padj.between(0, 1).all(),
            f"Probabilities outside [0,1] in {path}")
    return table


def select_gene_effects(table, per_direction=3):
    """Select existing point estimates; never construct intervals from legacy SEs."""
    significant = table[(table.padj < style.ADJUSTED_P_CUTOFF)
                        & (table.log2FoldChange.abs() >= style.LOG2FC_CUTOFF)]
    ordered = significant.sort_values(['padj', 'gene_id'], kind='stable')
    selected = pd.concat([ordered[ordered.log2FoldChange > 0].head(per_direction),
                          ordered[ordered.log2FoldChange < 0].head(per_direction)])
    return selected.sort_values(['log2FoldChange', 'gene_id'], ascending=[False, True], kind='stable').copy()


def gene_effect_labels(selected):
    """Keep symbols readable while identifying distinct Ensembl rows sharing one."""
    duplicated = selected.supplied_symbol.duplicated(keep=False)
    labels = []
    for row, ambiguous in zip(selected.itertuples(index=False), duplicated):
        labels.append(f"{row.supplied_symbol} [..{row.gene_id[-6:]}]" if ambiguous else row.supplied_symbol)
    require(len(set(labels)) == len(labels), 'Selected gene-label suffixes are ambiguous; review full Ensembl IDs')
    return labels


def effect_limit(selected):
    largest = selected.log2FoldChange.abs().max() if len(selected) else 0
    return max(.4, float(np.ceil(largest * 1.12 / .2) * .2))


def gsea_display(table):
    """Honor the analysis pipeline's audited selection, including scarce-direction fallback."""
    style.displayed_gsea(table)  # Reuse finite-score, probability and term-label checks.
    require(len(table) <= 4, 'Legacy GSEA display files must contain at most four audited terms')
    return table.sort_values('NES', ascending=False, kind='stable')


def export_figure(fig, number, question, fonts, output_dir, details):
    stem = f'figure_{number:02d}_legacy'
    for _ in range(2):
        for text in fig.findobj(Text):
            bold = text.get_fontweight() in ['bold', 700] or text.get_fontproperties().get_file() == str(fonts['Bold'])
            text.set_fontproperties(fm.FontProperties(fname=str(fonts['Bold' if bold else 'Regular']),
                                                     size=style.font_points()))
            text.set_parse_math(False)
        fig.canvas.draw()
    labels = [{'text': text.get_text(), 'size_px': round(text.get_fontsize() * DPI / 72, 4),
               'font': 'Nimbus Sans Bold' if text.get_fontproperties().get_file() == str(fonts['Bold']) else 'Nimbus Sans Regular'}
              for text in fig.findobj(Text) if text.get_visible() and text.get_text()]
    require(all(abs(record['size_px'] - style.FONT_PIXELS) < 1e-8 for record in labels), 'Mixed legacy font sizes')
    style.validate_gsea_layout(fig)
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        require(box.x0 >= 0 and box.y0 >= 0 and box.x1 <= W and box.y1 <= H,
                f"Legacy label outside canvas: {text.get_text()!r}, {box.bounds}")
    for label in getattr(fig, '_legacy_gene_labels', []):
        require(label.get_window_extent(renderer).x0 >= 788,
                f"Gene label crosses into panel a: {label.get_text()!r}")
    metadata = {'Title': f'FIGURE {number} | {question} | Legacy exploratory analysis',
                'Description': 'Supplied legacy differential-expression estimates; exploratory enrichment. See separate legend.'}
    fig.savefig(output_dir / f'{stem}.png', dpi=100, metadata=metadata)
    fig.savefig(output_dir / f'{stem}_300dpi.png', dpi=300, metadata=metadata)
    fig.savefig(output_dir / f'{stem}.pdf', metadata={'Title': metadata['Title'], 'Subject': metadata['Description']})
    fig.savefig(output_dir / f'{stem}.svg', metadata=metadata)
    outlined = output_dir / f'{stem}.svg'
    outlined.write_text('\n'.join(line.rstrip() for line in outlined.read_text().splitlines()) + '\n')
    require('<text' not in outlined.read_text(), 'Outlined legacy SVG contains live text')
    with matplotlib.rc_context({'svg.fonttype': 'none'}):
        fig.savefig(output_dir / f'{stem}_editable.svg', metadata=metadata)
    editable = output_dir / f'{stem}_editable.svg'
    css = []
    for name, weight, fontstyle in [('Regular', 400, 'normal'), ('Bold', 700, 'normal'), ('Italic', 400, 'italic')]:
        encoded = base64.b64encode(fonts[name].read_bytes()).decode()
        css.append(f'@font-face {{ font-family: "Nimbus Sans"; font-weight: {weight}; font-style: {fontstyle}; src: url("data:font/ttf;base64,{encoded}") format("truetype"); }}')
    code = editable.read_text().replace('<defs>', '<defs>\n<style type="text/css"><![CDATA[\n' + '\n'.join(css) + '\n]]></style>', 1)
    editable.write_text('\n'.join(line.rstrip() for line in code.splitlines()) + '\n')
    style.plt.close(fig)
    with Image.open(output_dir / f'{stem}.png') as image:
        require(image.size == (W, H), 'Wrong legacy canvas dimensions')
    with fitz.open(output_dir / f'{stem}.pdf') as document:
        spans = [span for page in document for block in page.get_text('dict')['blocks']
                 for line in block.get('lines', []) for span in line['spans'] if span['text'].strip()]
        require(spans and all(span['font'].startswith('NimbusSans-') for span in spans), 'Legacy PDF font fallback')
        require(all(abs(span['size'] * DPI / 72 - style.FONT_PIXELS) < .01 for span in spans),
                'Actual legacy PDF text has mixed sizes')
    exports = [output_dir / f'{stem}{suffix}' for suffix in ('.png', '_300dpi.png', '.pdf', '.svg', '_editable.svg')]
    return {'figure': number, 'file_stem': stem, 'question': question, 'canvas_pixels': [W, H],
            'font_pixels': style.FONT_PIXELS, 'all_visible_typography': 'Nimbus Sans', 'labels': labels,
            'pdf_nonempty_text_spans_checked': len(spans),
            'export_sha256': {path.name: digest(path) for path in exports}, **details}


def render_question(spec, table, gsea_tables, fonts, output_dir, xmax, ymax, nes_limit):
    number, question, first, second, _ = spec
    colors = style.COMPARISON_COLORS[number]
    lower, higher = colors['lower'], colors['higher']
    fig = style.plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)

    def axes(left, top, width, height):
        return fig.add_axes([left / W, 1 - (top + height) / H, width / W, height / H])

    def panel(letter, left, top):
        fig.text(left / W, 1 - top / H, letter, ha='left', va='top',
                 fontproperties=fm.FontProperties(fname=str(fonts['Bold']), size=style.font_points()))

    def theme(ax):
        ax.tick_params(direction='out', pad=5)
        ax.grid(False)
        ax.set_axisbelow(True)

    panel('a', 40, 27)
    ax = axes(130, 88, 540, 310)
    comparison = f'{first} vs {second}'
    ax.set_title(comparison, pad=20)
    x = table.log2FoldChange.to_numpy()
    y = style.volcano_ordinates(table.padj)
    selected = (table.padj.to_numpy() < style.ADJUSTED_P_CUTOFF) & (np.abs(x) >= style.LOG2FC_CUTOFF)
    ax.scatter(x[~selected], y[~selected], s=8, c=style.NS, alpha=.55, edgecolors='none', rasterized=True)
    for mask, color in [(selected & (x < 0), lower), (selected & (x > 0), higher)]:
        ax.scatter(x[mask], y[mask], s=11, c=color, alpha=.82, edgecolors='none', rasterized=True)
    ax.axvline(0, color=style.GRID, lw=.65, ls=(0, (3, 4)), zorder=0)
    for cutoff in (-style.LOG2FC_CUTOFF, style.LOG2FC_CUTOFF):
        ax.axvline(cutoff, color=style.AXIS, lw=.8, ls=(0, (4, 4)), zorder=0)
    ax.axhline(-np.log10(style.ADJUSTED_P_CUTOFF), color=style.GRID, lw=.65, ls=(0, (3, 4)), zorder=0)
    ax.set(xlim=(-xmax, xmax), ylim=(-.08, ymax), xlabel='log2 fold change', ylabel='−log10(adjusted P)')
    yticks = MaxNLocator(nbins=5, integer=True, steps=[1, 2, 2.5, 5, 10]).tick_values(0, ymax)
    ax.set_yticks(yticks[(yticks >= 0) & (yticks <= ymax)])
    xticks = MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10], symmetric=True).tick_values(-xmax, xmax)
    ax.set_xticks(xticks[(xticks >= -xmax) & (xticks <= xmax)])
    theme(ax)

    panel('b', 788, 27)
    fig.text(1150 / W, 1 - 58 / H, 'Genes', ha='center', va='baseline', fontsize=style.font_points())
    effects = select_gene_effects(table)
    names = gene_effect_labels(effects)
    ax = axes(1120, 110, 290, 288)
    rows = np.arange(len(effects))
    if len(effects):
        ax.scatter(effects.log2FoldChange, rows, s=90,
                   c=[higher if value > 0 else lower for value in effects.log2FoldChange], edgecolors='none', zorder=3)
        ax.set_yticks(rows, labels=names)
        ax.set_ylim(len(effects) - .5, -.5)
    else:
        ax.set_yticks([])
        ax.text(.5, .5, 'No genes meet cutoffs', transform=ax.transAxes, ha='center', va='center')
        ax.set_ylim(.5, -.5)
    fig._legacy_gene_labels = ax.get_yticklabels()
    ax.axvline(0, color=style.GRID, lw=.65, ls=(0, (3, 4)), zorder=0)
    for cutoff in (-style.LOG2FC_CUTOFF, style.LOG2FC_CUTOFF):
        ax.axvline(cutoff, color=style.AXIS, lw=.65, ls=(0, (4, 4)), zorder=0)
    limit = effect_limit(effects)
    ax.set_xlim(-limit, limit)
    ticks = MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10], symmetric=True).tick_values(-limit, limit)
    ax.set_xticks(ticks[(ticks >= -limit) & (ticks <= limit)])
    ax.set_xlabel('log2 fold change', labelpad=8)
    ax.spines['left'].set_visible(False)
    theme(ax)
    ax.tick_params(axis='y', length=0, pad=12)

    gsea_details = {}
    fig._gsea_layout_regions = []
    for library, letter, term_left, left, q_left, title in [
            ('kegg', 'c', 40, 516, 662, 'KEGG'), ('go_bp', 'd', 808, 1284, 1430, 'GO BP')]:
        shown = gsea_display(gsea_tables[library])
        panel(letter, 40 if library == 'kegg' else 788, 494)
        fig.text((390 if library == 'kegg' else 1160) / W, 1 - 524 / H, title,
                 ha='center', va='baseline', fontsize=style.font_points())
        top, height, width = 552, 376, 130
        ax = axes(left, top, width, height)
        labels = [style.wrap_term(re.sub(r'^(?:GO:\d+|mmu\d+)\s+', '', term).strip(),
                                  left - term_left - 16, fonts['Regular']) for term in shown.term]
        positions = style.pathway_row_positions(labels, height)
        ax.axvline(0, color=style.GRID, lw=.75, ls=(0, (3, 4)), zorder=0)
        term_labels, q_labels = [], []
        for row, (_, result), label in zip(positions, shown.iterrows(), labels):
            color = higher if result.NES > 0 else lower
            ax.scatter(result.NES, row, s=90, facecolors=color if result.gsea_fdr_q < .05 else 'white',
                       edgecolors=color, linewidths=1.4, zorder=3)
            term_labels.append(ax.text((term_left-left) / width, row, label, transform=ax.get_yaxis_transform(),
                                       ha='left', va='center', fontsize=style.font_points(), linespacing=1.1))
            q_labels.append(ax.text((q_left-left) / width, row, style.q_label(result.gsea_fdr_q),
                                    transform=ax.get_yaxis_transform(), ha='left', va='center', fontsize=style.font_points()))
        ax.set_yticks([])
        if shown.empty:
            ax.text(.5, .5, 'No eligible sets', transform=ax.transAxes, ha='center', va='center')
        limit = max(nes_limit, float(shown.NES.abs().max()) * 1.1) if len(shown) else nes_limit
        ax.set(xlim=(-limit, limit), ylim=(height, 0))
        if limit == 2.6:
            ax.set_xticks([-2, 0, 2])
        else:
            ticks = MaxNLocator(nbins=3, steps=[1, 2, 2.5, 5, 10], symmetric=True).tick_values(-limit, limit)
            ax.set_xticks(ticks[(ticks >= -limit) & (ticks <= limit)])
        ax.set_xlabel('NES', labelpad=5)
        ax.spines['left'].set_visible(False)
        theme(ax)
        fig.text((q_left+30) / W, 1 - 540 / H, 'q', ha='center', va='bottom', fontsize=style.font_points())
        fig._gsea_layout_regions.append({'axes': ax, 'collection': title, 'q_labels': q_labels,
                                         'term_labels': term_labels, 'left': term_left,
                                         'right': 768 if library == 'kegg' else W})
        gsea_details[library] = {'terms': shown.term.tolist(), 'NES': shown.NES.tolist(),
                                 'gsea_fdr_q': shown.gsea_fdr_q.tolist(), 'absolute_nes_limit': limit,
                                 'selection': 'audited display file: up to two strongest NES per direction; scarce-direction remainder by absolute NES; no significance filter',
                                 'rank_metric': PRIMARY_RANK_METRIC,
                                 'rank_metric_description': 'sign(log2FC) * -log10(raw P); exploratory proxy, not an original test statistic'}
    details = {'comparison': comparison, 'panels': ['volcano', 'gene-effect points', 'KEGG preranked GSEA', 'GO BP preranked GSEA'],
               'legacy_exploratory': True, 'sample_expression_available': False, 'heatmap_panel': False,
               'pca_panel': False, 'volcano_genes': len(table),
               'volcano_limits': {'absolute_log2fc': xmax, 'negative_log10_adjusted_p': ymax},
               'absolute_log2fc_cutoff': style.LOG2FC_CUTOFF, 'adjusted_p_cutoff': style.ADJUSTED_P_CUTOFF,
               'joint_DEG_genes': int(selected.sum()), 'higher_at_joint_cutoffs': int((selected & (x > 0)).sum()),
               'lower_at_joint_cutoffs': int((selected & (x < 0)).sum()), 'comparison_colors': colors,
               'gene_effects': {'gene_ids': effects.gene_id.tolist(), 'display_labels': names,
                                'log2FoldChange': effects.log2FoldChange.tolist(), 'padj': effects.padj.tolist(),
                                'absolute_axis_limit': effect_limit(effects), 'confidence_intervals': None,
                                'selection': 'up to three smallest adjusted-P hits per direction passing both display cutoffs'},
               'gsea': gsea_details}
    return export_figure(fig, number, question, fonts, output_dir, details)


def run(data_dir, gsea_dir, output_dir):
    from analysis.legacy_gsea import validate_completed_run
    from analysis.legacy_workbook import validate_completed_import

    data_dir, gsea_dir, output_dir = map(lambda path: Path(path).resolve(), (data_dir, gsea_dir, output_dir))
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / 'legacy_plot_audit.json'
    write_json(report_path, {'complete': False})
    provenance_path = data_dir / 'legacy_provenance.json'
    provenance = validate_completed_import(data_dir)
    require(provenance.get('importer_sha256') == digest(REPOSITORY / 'analysis/legacy_workbook.py'),
            'Legacy importer changed; regenerate the imported results before rendering')
    manifest = validate_completed_run(gsea_dir, data_dir)
    require(manifest.get('rank_metric') == PRIMARY_RANK_METRIC,
            'Render the primary signed-log-P ranking; sensitivity results belong in their separate report')
    for spec, source in zip(SPECS, provenance['figures']):
        require((spec[2], spec[3]) == (source['numerator'], source['denominator']),
                f'Legacy group direction differs for figure {spec[0]}')
    paths = [Path(__file__), Path(style.__file__), provenance_path, gsea_dir / 'gsea_manifest.json',
             REPOSITORY / 'analysis/legacy_workbook.py', REPOSITORY / 'analysis/legacy_gsea.py']
    paths += [data_dir / spec[4] for spec in SPECS]
    paths += [gsea_dir / f'figure_{number:02d}_{library}_display.tsv' for number in (1, 2, 3) for library in ('kegg', 'go_bp')]
    input_hashes = {saved_path(path): digest(path) for path in paths}
    tables = [read_contrast(data_dir / spec[4]) for spec in SPECS]
    displays = [{library: style.read_table(gsea_dir / f'figure_{number:02d}_{library}_display.tsv')
                 for library in ('kegg', 'go_bp')} for number in (1, 2, 3)]
    xmax, ymax = style.volcano_limits(tables)
    nes_limit = max([2.6] + [float(gsea_display(table).NES.abs().max()) * 1.1
                           for group in displays for table in group.values() if len(table)])
    fonts = style.prepare_nimbus()
    style.configure_theme()
    figures = [render_question(spec, table, group, fonts, output_dir, xmax, ymax, nes_limit)
               for spec, table, group in zip(SPECS, tables, displays)]
    for path in paths:
        require(digest(path) == input_hashes[saved_path(path)], f'Legacy render input changed: {path.name}')
    validate_completed_run(gsea_dir, data_dir)
    # Preserve a compact, explicit per-panel data table for downstream review.
    for figure in figures:
        effect = figure['gene_effects']
        effect_path = output_dir / f"figure_{figure['figure']:02d}_gene_effects.tsv"
        pd.DataFrame({'gene_id': effect['gene_ids'], 'display_label': effect['display_labels'],
                      'log2FoldChange': effect['log2FoldChange'], 'padj': effect['padj']}).to_csv(
            effect_path, sep='\t', index=False, float_format='%.17g')
        figure['plot_data_file'] = saved_path(effect_path)
        figure['plot_data_sha256'] = digest(effect_path)
    report = {'complete': True, 'legacy_exploratory': True, 'font_pixels': style.FONT_PIXELS,
              'input_sha256': input_hashes, 'figures': figures,
              'gsea_manifest_sha256': digest(gsea_dir / 'gsea_manifest.json'),
              'gsea_primary_analyses': len(manifest.get('analyses', []))}
    write_json(report_path, report)
    print(json.dumps([{'figure': figure['figure'], 'comparison': figure['comparison'],
                       'joint_DEG_genes': figure['joint_DEG_genes'], 'font_pixels': figure['font_pixels']}
                      for figure in figures], indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=REPOSITORY / 'analysis/legacy_results')
    parser.add_argument('--gsea-dir', type=Path, help='Defaults to DATA_DIR/gsea')
    parser.add_argument('--output-dir', type=Path, default=REPOSITORY / 'figures/legacy')
    args = parser.parse_args()
    run(args.data_dir, args.gsea_dir or args.data_dir / 'gsea', args.output_dir)


if __name__ == '__main__':
    main()
