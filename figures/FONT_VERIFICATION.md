# Nimbus Sans verification

The three `figure_0*_gsea_synthetic` export sets use **Nimbus Sans**. The check examines the exported font programs, rather than trusting a font setting or filename. The latest file hashes, source-font hashes and detailed results are recorded in [the verification report](synthetic_demo/export_font_verification.json).

- **PDF:** every embedded TrueType font is extracted and its family/style names, Unicode character map, advance widths, sidebearings and glyph outlines are checked against the installed URW Nimbus Sans source fonts. All text-span font names are checked for substitutions. `pdffonts` output supplies a second embedding check. Regular and Bold are used in the current artwork; no fallback font was found.
- **Outlined SVG:** all text is converted into glyph paths. Every text-glyph reference resolves to a Nimbus Sans glyph definition, so a viewer cannot substitute another typeface.
- **Editable SVG:** the embedded Regular, Bold and Italic font programs undergo the same source-metric and glyph-outline comparison; all explicit family declarations specify Nimbus Sans. A viewer must support the embedded `@font-face` fonts. Use the PDF or outlined SVG when exchanging final artwork.
- **PNG:** a bitmap does not contain a recoverable font identity. The renderer explicitly selects Nimbus Sans font files and forbids missing-source substitution; a raster comparison checks that each PNG matches the verified companion PDF, allowing differences in antialiasing. Both the 1536 × 1024 PNG and 4608 × 3072 export are hashed in the report.

The source fonts are URW OpenType/CFF Nimbus Sans. For reliable PDF embedding, the renderer converts their cubic curves into TrueType quadratic curves with a maximum conversion error of **0.5 font units out of 1000 units per em**. Character spacing is retained. The audit independently recomputes that conversion from the source fonts and compares the actual exported outlines; a renamed substitute would fail.

Recheck after every figure export:

```bash
python figures/synthetic_demo/verify_export_fonts.py
```

The command exits on a mismatch and writes `synthetic_demo/export_font_verification.json` only after all three export sets pass. A report applies to its recorded SHA-256 hashes; rerendering requires a fresh check. Dependencies are PyMuPDF, FontTools, Pillow and NumPy; `pdffonts` is used when available. If fonts are installed elsewhere, set `NIMBUS_SANS_FONT_DIR` to the directory containing the genuine `NimbusSans-{Regular,Bold,Italic}.otf` files.
