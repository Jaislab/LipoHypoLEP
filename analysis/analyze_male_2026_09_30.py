#!/usr/bin/env python3
"""Preserve and summarize the confirmed filtered male DE exports dated 2026-09-30.

This workflow does not refit DE, read sample expression, recompute adjusted
probabilities, or perform enrichment on a selected-gene background.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import pandas as pd

SOURCE = ROOT / 'data/male_2026_09_30/source'
OUTPUT = ROOT / 'analysis/male_2026_09_30_results'
SPECS = (
    (1, 'LD SALINE', 'NON-LD control', 'DEGmale_LDsaline_WT.csv', 1),
    (2, 'LD LEPTIN', 'LD SALINE', 'DEGmale_LDsaline_LDleptin.csv', -1),
    (3, 'LD LEPTIN', 'NON-LD control', 'DEGmale_LDleptin_WT.csv', 1),
)
FIELDS = ['ENSEMBL', 'SYMBOL', 'log2FoldChange', 'pvalue', 'padj']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def saved_path(path):
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def resolve(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def write_json(path, value):
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False, encoding='utf-8') as handle:
        pending = Path(handle.name)
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
    pending.replace(path)


def read_source(payload, filename, multiplier):
    require(multiplier in (-1, 1), 'Contrast multiplier must be +1 or -1')
    reader = csv.reader(io.StringIO(payload.decode('utf-8-sig'), newline=''), strict=True)
    require(next(reader, None) == FIELDS, f'Unexpected or duplicate CSV headers: {filename}')
    rows, seen = [], set()
    for line, row in enumerate(reader, 2):
        require(len(row) == len(FIELDS), f'Wrong row width: {filename}:{line}')
        gene, symbol, *numeric = row
        require(re.fullmatch(r'ENSMUSG[0-9]{11}(?:\.[0-9]+)?', gene), f'Invalid mouse Ensembl ID: {gene}')
        require(gene not in seen, f'Duplicate Ensembl ID: {gene}')
        require(symbol and symbol.strip() == symbol, f'Missing or padded symbol: {gene}')
        seen.add(gene)
        try:
            effect, raw_p, adjusted_p = map(float, numeric)
        except ValueError as error:
            raise ValueError(f'Invalid numeric field: {filename}:{line}') from error
        require(all(math.isfinite(v) for v in (effect, raw_p, adjusted_p)), f'Nonfinite statistics: {gene}')
        require(all(Decimal(raw) == 0 or value != 0 for raw, value in zip(numeric, (effect, raw_p, adjusted_p))),
                f'Numeric underflow would change a nonzero source value to zero: {gene}')
        require(0 < raw_p <= 1 and 0 < adjusted_p <= 1, f'Probabilities must be in (0,1]: {gene}')
        rows.append({'gene_id': gene, 'symbol': symbol, 'log2FoldChange': effect * multiplier,
                     'pvalue': raw_p, 'padj': adjusted_p, 'source_log2FoldChange': effect,
                     'source_file': filename})
    require(rows, f'Empty source: {filename}')
    return pd.DataFrame(rows)


def joint_hits(table):
    return (table.padj < .05) & (table.log2FoldChange.abs() >= .2)


def selected_genes(table, count=10):
    return (table.loc[joint_hits(table)].sort_values(['padj', 'gene_id'], kind='stable').head(count)
            .sort_values(['log2FoldChange', 'gene_id'], ascending=[False, True], kind='stable').copy())


def validate_results(output_dir=OUTPUT):
    output_dir = Path(output_dir)
    manifest_path = output_dir / 'analysis_manifest.json'
    payload = manifest_path.read_bytes()
    manifest = json.loads(payload)
    require(manifest.get('complete') is True, 'Male analysis is incomplete')
    require(manifest.get('analyzer_sha256') == digest(Path(__file__)), 'Analyzer changed; regenerate results')
    require(manifest.get('filtered_lists_user_confirmed') is True, 'Filtered-export status is missing')
    require(manifest.get('first_minus_second_user_confirmed') is True, 'Source direction confirmation is missing')
    require(manifest.get('enrichment', {}).get('performed') is False, 'Filtered exports cannot supply this enrichment workflow')
    require(len(manifest.get('figures', [])) == 3, 'Three contrasts are required')
    expected = {'baseline_posttreatment_shared_hits.tsv'}
    for item, (number, first, second, filename, multiplier) in zip(manifest['figures'], SPECS):
        require((item['figure'], item['numerator'], item['denominator'], item['source_file'], item['multiplier'])
                == (number, first, second, filename, multiplier), 'Contrast direction changed')
        expected.update({f'figure_{number:02d}_supplied.tsv', f'figure_{number:02d}_joint_hits.tsv',
                         f'figure_{number:02d}_selected_genes.tsv'})
    require(set(manifest.get('artifact_sha256', {})) == expected, 'Incomplete result artifact inventory')
    require(len(manifest.get('source_sha256', {})) == 3, 'Three source files are required')
    require({Path(p).name for p in manifest['source_sha256']} == {spec[3] for spec in SPECS},
            'Unexpected source file inventory')
    for name, checksum in manifest['source_sha256'].items():
        require(digest(resolve(name)) == checksum, f'Source changed: {name}')
    for name, checksum in manifest['artifact_sha256'].items():
        require(digest(output_dir / name) == checksum, f'Result changed: {name}')
    require(manifest_path.read_bytes() == payload, 'Manifest changed during validation')
    return manifest


def run(source_dir=SOURCE, output_dir=OUTPUT):
    source_dir, output_dir = Path(source_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / 'analysis_manifest.json', {'complete': False})
    code_hash = digest(Path(__file__))
    sources = {source_dir / spec[3]: (source_dir / spec[3]).read_bytes() for spec in SPECS}
    tables = {number: read_source(sources[source_dir / filename], filename, multiplier)
              for number, _, _, filename, multiplier in SPECS}
    summary = {'complete': False, 'data_date': '2026-09-30', 'sex': 'male (source filenames)',
               'filtered_lists_user_confirmed': True, 'first_minus_second_user_confirmed': True,
               'source_sha256': {saved_path(p): hashlib.sha256(data).hexdigest() for p, data in sources.items()},
               'analyzer_sha256': code_hash,
               'display_rule': {'adjusted_p_strictly_below': .05, 'absolute_log2fc_at_least': .2},
               'gene_selection': 'Up to ten smallest adjusted P among joint-cutoff hits; Ensembl tie-break; effect ordering for display.',
               'preservation': 'All exported rows, identifiers, raw and adjusted P values retained. Figure 2 effects negated; original effects retained. No DE refit, effect rescaling or new BH correction.',
               'sample_expression_available': False, 'signed_test_statistic_available': False,
               'sample_inclusion': 'No expression columns or mouse-to-library mapping supplied; replicate counts not assigned from the project animal register.',
               'enrichment': {'performed': False, 'status': 'pending_complete_tested_gene_tables',
                              'reason': 'User confirmed filtered lists. Full tested-gene rankings/backgrounds are unavailable. No GSEA or overrepresentation analysis performed; no legacy background substituted.'},
               'figures': []}
    with tempfile.TemporaryDirectory(prefix='.male-analysis-', dir=output_dir) as temporary:
        stage = Path(temporary)
        for number, first, second, filename, multiplier in SPECS:
            table = tables[number]
            hits = joint_hits(table)
            for suffix, data in [('supplied', table), ('joint_hits', table.loc[hits]), ('selected_genes', selected_genes(table))]:
                data.to_csv(stage / f'figure_{number:02d}_{suffix}.tsv', sep='\t', index=False, float_format='%.17g')
            summary['figures'].append({'figure': number, 'numerator': first, 'denominator': second,
                'source_file': filename, 'multiplier': multiplier, 'supplied_genes': len(table),
                'joint_cutoff_hits': int(hits.sum()), 'higher': int((hits & (table.log2FoldChange > 0)).sum()),
                'lower': int((hits & (table.log2FoldChange < 0)).sum()),
                'min_absolute_supplied_log2fc': float(table.log2FoldChange.abs().min()),
                'maximum_supplied_padj': float(table.padj.max())})
        first = tables[1].loc[joint_hits(tables[1])]
        third = tables[3].loc[joint_hits(tables[3])]
        shared = first[['gene_id', 'symbol', 'log2FoldChange', 'padj']].merge(
            third[['gene_id', 'symbol', 'log2FoldChange', 'padj']], on='gene_id', suffixes=('_saline', '_leptin'), validate='one_to_one')
        shared['same_direction'] = (shared.log2FoldChange_saline > 0) == (shared.log2FoldChange_leptin > 0)
        shared.sort_values('gene_id').to_csv(stage / 'baseline_posttreatment_shared_hits.tsv', sep='\t', index=False, float_format='%.17g')
        summary['shared_baseline_posttreatment_hits'] = {
            'both_pass': len(shared), 'same_direction': int(shared.same_direction.sum()),
            'opposite_direction': int((~shared.same_direction).sum()),
            'interpretation': 'Intersection of supplied joint-cutoff hits. An absent gene has unknown status; absence or nonsignificance is not normalization.'}
        require(digest(Path(__file__)) == code_hash and all(p.read_bytes() == data for p, data in sources.items()),
                'An analysis input changed during calculation')
        summary['artifact_sha256'] = {p.name: digest(p) for p in sorted(stage.iterdir())}
        for path in stage.iterdir():
            path.replace(output_dir / path.name)
    summary['complete'] = True
    write_json(output_dir / 'analysis_manifest.json', summary)
    validate_results(output_dir)
    print(json.dumps(summary['figures'], indent=2))
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, default=SOURCE)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT)
    args = parser.parse_args()
    run(args.source_dir, args.output_dir)
