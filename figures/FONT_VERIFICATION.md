# Nimbus Sans verification

The three `figure_0*_gsea_synthetic` export sets use **Nimbus Sans**. The check examines the exported font programs, rather than trusting a font setting or filename. The latest file hashes, source-font hashes and detailed results are recorded in [the verification report](synthetic_demo/export_font_verification.json).

The current renderer applies a single **1.20 typography multiplier** to all text roles, preserving the size hierarchy. Main labels range from 24 px (q values) through 28.8–30 px (axis/comparison and collection labels), with 25.2 px ticks/pathway terms and 38.4 px panel letters on the 1536 px canvas. At 180 mm placement width this corresponds to approximately 8–10 pt main labels and 12.8 pt panel letters. The plot audit records `typography_scale` and label sizes; artwork dimensions stay 1536 × 1024. Mathematical subscripts retain their normal relative sizing.

- **PDF:** every embedded TrueType font is extracted and its family/style names, Unicode character map, advance widths, sidebearings and glyph outlines are checked against the installed URW Nimbus Sans source fonts. All text-span font names are checked for substitutions. `pdffonts` output supplies a second embedding check. Regular and Bold are used in the current artwork; no fallback font was found.
- **Outlined SVG:** all text is converted into glyph paths. Every text-glyph reference resolves to a Nimbus Sans glyph definition, so a viewer cannot substitute another typeface.
- **Editable SVG:** the embedded Regular, Bold and Italic font programs undergo the same source-metric and glyph-outline comparison; all explicit family declarations specify Nimbus Sans. A viewer must support the embedded `@font-face` fonts. Use the PDF or outlined SVG when exchanging final artwork.
- **PNG:** a bitmap does not contain a recoverable font identity. The renderer explicitly selects Nimbus Sans font files and forbids missing-source substitution; a raster comparison checks that each PNG matches the verified companion PDF, allowing differences in antialiasing. Both the 1536 × 1024 PNG and 4608 × 3072 export are hashed in the report.

The source fonts are URW OpenType/CFF Nimbus Sans. For reliable PDF embedding, the renderer converts their cubic curves into TrueType quadratic curves with a maximum conversion error of **0.5 font units out of 1000 units per em**. Character spacing is retained. The audit independently recomputes that conversion from the source fonts and compares the actual exported outlines; a renamed substitute would fail.

Recheck after every figure export:

```bash
python figures/synthetic_demo/verify_export_fonts.py
```

The command exits on a mismatch and records RUNNING, FAILED or PASS atomically, so a failed rerun cannot retain an earlier PASS. `python figures/synthetic_demo/verify_export_fonts.py --check-report` checks freshness against exports, source fonts and scripts. Validation remains active under `python -O`. A report applies to its recorded SHA-256 hashes; rerendering requires a fresh check. Dependencies are PyMuPDF, FontTools, Pillow and NumPy; `pdffonts` is used when available. If fonts are installed elsewhere, set `NIMBUS_SANS_FONT_DIR` to the directory containing the genuine `NimbusSans-{Regular,Bold,Italic}.otf` files. Converted fonts preserve source attribution and license name records; Nimbus is supplied by [URW Base 35](https://github.com/ArtifexSoftware/urw-base35-fonts) under its AGPL font exception.
