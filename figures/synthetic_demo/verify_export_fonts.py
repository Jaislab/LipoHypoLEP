#!/usr/bin/env python3
"""Independently verify actual exported Nimbus Sans font programs and glyphs.

Run after every render. Reads the exports, not the renderer's own audit claims.
Requires PyMuPDF, FontTools, Pillow and NumPy; pdffonts is optional evidence.
"""
from __future__ import annotations

import base64
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import fitz
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.ttLib import TTFont
import numpy as np
from PIL import Image, ImageFilter

ROOT = Path(__file__).resolve().parent
OUT = ROOT.parent
SOURCE = Path(os.environ.get('NIMBUS_SANS_FONT_DIR', '/usr/share/fonts/opentype/urw-base35'))
NS = {'s': 'http://www.w3.org/2000/svg'}
XLINK = '{http://www.w3.org/1999/xlink}href'
STYLES = ('Regular', 'Bold', 'Italic')


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def names(font):
    return {str(i): sorted({r.toUnicode() for r in font['name'].names if r.nameID == i})
            for i in (1, 2, 3, 4, 6)}


def converted_glyphs(font):
    """Recompute expected outlines from independently opened source OTF."""
    glyphset = font.getGlyphSet()
    result = {}
    for name in font.getGlyphOrder():
        pen = TTGlyphPen(glyphset)
        glyphset[name].draw(Cu2QuPen(pen, max_err=0.5, reverse_direction=True))
        result[name] = pen.glyph()
    return result


def outline_signature(glyph, glyphset):
    coordinates, ends, flags = glyph.getCoordinates(glyphset)
    return (tuple(map(tuple, coordinates)), tuple(ends), tuple(flags))


def audit_program(data, expected):
    font = TTFont(io.BytesIO(data))
    n = names(font)
    assert n['1'] == ['Nimbus Sans'], f'Unexpected embedded family: {n}'
    assert len(n['2']) == 1 and n['2'][0] in STYLES, f'Unexpected font style: {n}'
    style = n['2'][0]
    reference, expected_glyphs = expected[style]
    assert font['head'].unitsPerEm == reference['head'].unitsPerEm
    cmap = font.getBestCmap()
    assert cmap, 'Embedded font has no Unicode character map'
    reference_cmap = reference.getBestCmap()
    mismatch = []
    for codepoint, glyph in cmap.items():
        source_glyph = reference_cmap.get(codepoint)
        if source_glyph is None:
            mismatch.append(f'U+{codepoint:04X}: absent from source Nimbus')
            continue
        if font['hmtx'].metrics[glyph] != reference['hmtx'].metrics[source_glyph]:
            mismatch.append(f'U+{codepoint:04X}: horizontal metrics')
        actual = outline_signature(font['glyf'][glyph], font['glyf'])
        trusted = outline_signature(expected_glyphs[source_glyph], expected_glyphs)
        if actual != trusted:
            mismatch.append(f'U+{codepoint:04X}: glyph outline')
    assert not mismatch, f'{style} source mismatch: {mismatch[:12]}'
    return {'style': style, 'names': n, 'sha256': sha256(data),
            'bytes': len(data), 'units_per_em': font['head'].unitsPerEm,
            'unicode_glyphs_compared': len(cmap),
            'advance_widths_and_sidebearings_match_source': True,
            'outlines_match_source_after_cubic_to_quadratic_conversion': True,
            'conversion_max_error_font_units': 0.5}


def audit_pdf(path, expected):
    with fitz.open(path) as doc:
        entries = {entry[0]: entry for page in doc for entry in page.get_fonts(full=True)}
        programs = []
        for xref, entry in entries.items():
            name, extension, kind, content = doc.extract_font(xref)
            assert extension == 'ttf' and content, f'Font is not embedded TrueType: {entry}'
            item = audit_program(content, expected)
            item.update({'pdf_resource_name': name, 'pdf_font_type': kind, 'xref': xref})
            programs.append(item)
        assert programs, 'PDF contains no embedded fonts'
        spans = [span for page in doc for block in page.get_text('dict')['blocks']
                 if 'lines' in block for line in block['lines'] for span in line['spans']]
        span_fonts = sorted({span['font'] for span in spans if span['text'].strip()})
        assert spans and all(font.startswith('NimbusSans-') for font in span_fonts), span_fonts
        result = {'font_programs': programs, 'text_span_fonts': span_fonts,
                  'text_spans_checked': len(spans), 'pages': len(doc),
                  'fallback_font_found': False}
    if shutil.which('pdffonts'):
        evidence = subprocess.run(['pdffonts', str(path)], text=True, capture_output=True, check=True)
        result['pdffonts_stdout'] = evidence.stdout
    return result


def audit_outlined_svg(path):
    tree = ET.parse(path)
    assert not tree.findall('.//s:text', NS), 'Outlined SVG still contains live text'
    refs = []
    for group in tree.findall('.//s:g', NS):
        if group.get('id', '').startswith('text_'):
            refs.extend(x.get(XLINK, '') for x in group.findall('.//s:use', NS))
    assert refs, 'No glyph references found in outlined SVG'
    assert all(re.match(r'#NimbusSans-(Regular|Bold|Italic)-[0-9a-fA-F]+$', x) for x in refs), refs
    glyphs = {element.get('id') for element in tree.findall('.//s:path', NS)
              if element.get('id', '').startswith('NimbusSans-')}
    assert all(ref[1:] in glyphs for ref in refs), 'Missing outlined glyph definition'
    return {'live_text_elements': 0, 'glyph_uses_checked': len(refs),
            'glyph_prefixes': sorted({ref.rsplit('-', 1)[0].lstrip('#') for ref in refs}),
            'all_text_glyph_references_resolve': True,
            'note': 'Outlined artwork cannot substitute a font; identity is traced through the Nimbus glyph definitions and companion PDF font-program checks.'}


def audit_editable_svg(path, expected):
    code = path.read_text()
    programs = re.findall(r'data:font/ttf;base64,([A-Za-z0-9+/=]+)', code)
    assert len(programs) == 3, f'Expected embedded Regular, Bold, Italic fonts, found {len(programs)}'
    audits = [audit_program(base64.b64decode(program), expected) for program in programs]
    assert sorted(item['style'] for item in audits) == sorted(STYLES)
    tree = ET.fromstring(code)
    text = tree.findall('.//s:text', NS)
    assert text, 'Editable SVG contains no text'
    declarations = []
    for element in tree.iter():
        style = element.get('style', '')
        families = re.findall(r'font-family:\s*([^;]+)', style)
        for family in families:
            assert family.strip().strip('"\'') == 'Nimbus Sans', f'Unexpected SVG family: {family}'
            declarations.append(family)
    assert declarations, 'Editable SVG has no font-family declarations'
    return {'embedded_font_programs': audits, 'live_text_elements': len(text),
            'family_declarations_checked': len(declarations), 'fallback_font_found': False,
            'viewer_limitation': 'A viewer must honor embedded @font-face fonts. The outlined SVG and embedded-font PDF avoid that dependency.'}


def audit_png(path, pdf_path):
    with Image.open(path) as image:
        image.load()
        size = image.size
        info = {k: v for k, v in image.info.items() if isinstance(v, (str, int, float, tuple))}
        actual = image.convert('RGB')
    assert size == (1536, 1024), f'Unexpected image size: {size}'
    with fitz.open(pdf_path) as pdf:
        pixmap = pdf[0].get_pixmap(matrix=fitz.Matrix(100 / 72, 100 / 72), alpha=False)
        reference = Image.frombytes('RGB', (pixmap.width, pixmap.height), pixmap.samples)
    assert reference.size == actual.size
    # Different rasterizers have different antialiasing; slight blur makes this a
    # layout/content consistency measurement, never a claim of pixel identity.
    first = np.asarray(actual.filter(ImageFilter.GaussianBlur(1)), dtype=np.float32)
    second = np.asarray(reference.filter(ImageFilter.GaussianBlur(1)), dtype=np.float32)
    mae = float(np.abs(first - second).mean())
    assert mae < 5.0, f'PNG and current PDF appear inconsistent: blurred RGB MAE {mae:.3f}/255'
    return {'pixels': list(size), 'metadata': info,
            'pdf_raster_comparison': {'gaussian_blur_radius_px': 1,
                                      'mean_absolute_rgb_difference_0_to_255': round(mae, 5)},
            'font_identity_limit': 'Pixels do not encode a font name. Identity is supported by verified companion exports, explicit Nimbus source selection in the renderer, and matching rendered content.'}


def main():
    expected = {}
    sources = []
    for style in STYLES:
        path = SOURCE / f'NimbusSans-{style}.otf'
        assert path.is_file(), f'Trusted source font missing: {path}'
        font = TTFont(path)
        assert 'Nimbus Sans' in names(font)['1']
        expected[style] = (font, converted_glyphs(font))
        sources.append({'style': style, 'path': str(path), 'sha256': sha256(path.read_bytes()),
                        'names': names(font)})
    figures = []
    for number in range(1, 4):
        stem = f'figure_{number:02d}_gsea_synthetic'
        paths = {ext: OUT / f'{stem}{ext}' for ext in ('.pdf', '.svg', '_editable.svg', '.png', '_300dpi.png')}
        assert all(path.is_file() for path in paths.values()), f'Missing export in {stem}'
        record = {'figure': number, 'file_stem': stem,
                  'export_sha256': {path.name: sha256(path.read_bytes()) for path in paths.values()},
                  'pdf': audit_pdf(paths['.pdf'], expected),
                  'outlined_svg': audit_outlined_svg(paths['.svg']),
                  'editable_svg': audit_editable_svg(paths['_editable.svg'], expected),
                  'png': audit_png(paths['.png'], paths['.pdf'])}
        with Image.open(paths['_300dpi.png']) as image:
            assert image.size == (4608, 3072)
            record['high_resolution_png_pixels'] = list(image.size)
        figures.append(record)
    renderer = ROOT / 'render_figure.py'
    report = {'verified_at_utc': datetime.now(timezone.utc).isoformat(), 'status': 'PASS',
              'scope': 'Current exported files identified by SHA-256; rerun after every export.',
              'verification_method': 'Actual PDF and SVG font programs independently compared to installed URW Nimbus Sans source metrics and outlines.',
              'source_fonts': sources, 'renderer_sha256': sha256(renderer.read_bytes()),
              'verifier_sha256': sha256(Path(__file__).read_bytes()), 'figures': figures}
    output = ROOT / 'export_font_verification.json'
    output.write_text(json.dumps(report, indent=2) + '\n')
    for record in figures:
        print(f"Figure {record['figure']}: PASS; PDF={','.join(record['pdf']['text_span_fonts'])}; "
              f"editable SVG={len(record['editable_svg']['embedded_font_programs'])} verified font programs; "
              f"PNG/PDF blurred MAE={record['png']['pdf_raster_comparison']['mean_absolute_rgb_difference_0_to_255']}/255")
    print(f'Report: {output}')


if __name__ == '__main__':
    main()
