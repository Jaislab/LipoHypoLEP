# Nimbus Sans verification

The three `figure_0*_gsea_synthetic` export sets use **Nimbus Sans**. The check examines the exported font programs, rather than trusting a font setting or filename. The latest file hashes, source-font hashes and detailed results are recorded in [the verification report](synthetic_demo/export_font_verification.json).

**Every text element uses the same 32 px Nimbus Sans size.** This includes axes, ticks, comparison/group names, KEGG/GO pathway labels, q values, collection headings, colorbar labels and panel letters. Panel letters use bold weight at the same size. The renderer has one font-size constant; there are no role-specific size overrides or automatic font shrinking. At 180 mm placement width, 32 px corresponds to approximately 10.6 pt. Artwork dimensions stay 1536 × 1024.

The distinctly named `figure_0*_uniform_type_synthetic.png` and `.pdf` files used by the README are byte-identical copies of the verified corresponding canonical exports. Cross-figure validation checks these copies and their hashes. Superseded previews remain available in Git history.

An [independent check of every actual PDF text span](typography/uniform_type_measurements.json) verifies the single size across all three figures, including GSEA. It fails on any size mismatch or fallback family. The [same-scale before/after comparison](typography/uniform_type_before_after.png) shows the GSEA change. Plain `log2` and `−log10` axis notation keeps all characters at the same size without unsupported Unicode subscripts or smaller mathematical font runs.

- **PDF:** every embedded TrueType font is extracted and its family/style names, Unicode character map, advance widths, sidebearings and glyph outlines are checked against the installed URW Nimbus Sans source fonts. All text-span font names are checked for substitutions. `pdffonts` output supplies a second embedding check. Regular and Bold are used in the current artwork; no fallback font was found.
- **Outlined SVG:** all text is converted into glyph paths. Every text-glyph reference resolves to a Nimbus Sans glyph definition, so a viewer cannot substitute another typeface.
- **Editable SVG:** the embedded Regular, Bold and Italic font programs undergo the same source-metric and glyph-outline comparison; all explicit family declarations specify Nimbus Sans. A viewer must support the embedded `@font-face` fonts. Use the PDF or outlined SVG when exchanging final artwork.
- **PNG:** a bitmap does not contain a recoverable font identity. The renderer explicitly selects Nimbus Sans font files and forbids missing-source substitution; a raster comparison checks that each PNG matches the verified companion PDF, allowing differences in antialiasing. Both the 1536 × 1024 PNG and 4608 × 3072 export are hashed in the report.

The source fonts are URW OpenType/CFF Nimbus Sans. For reliable PDF embedding, the renderer converts their cubic curves into TrueType quadratic curves with a maximum conversion error of **0.5 font units out of 1000 units per em**. Character spacing is retained. The audit independently recomputes that conversion from the source fonts and compares the actual exported outlines; a renamed substitute would fail.

Recheck after every figure export:

```bash
python figures/synthetic_demo/verify_export_fonts.py
python figures/synthetic_demo/verify_uniform_typography.py
```

The command exits on a mismatch and records RUNNING, FAILED or PASS atomically, so a failed rerun cannot retain an earlier PASS. `python figures/synthetic_demo/verify_export_fonts.py --check-report` checks freshness against exports, source fonts and scripts. Validation remains active under `python -O`. A report applies to its recorded SHA-256 hashes; rerendering requires a fresh check. Dependencies are PyMuPDF, FontTools, Pillow and NumPy; `pdffonts` is used when available. If fonts are installed elsewhere, set `NIMBUS_SANS_FONT_DIR` to the directory containing the genuine `NimbusSans-{Regular,Bold,Italic}.otf` files. Converted fonts preserve source attribution and license name records; Nimbus is supplied by [URW Base 35](https://github.com/ArtifexSoftware/urw-base35-fonts) under its AGPL font exception.

The uniform-typography command uses Git revision `46faeea` to build its historical comparison; run it from a checkout with that history. Its `--check-report` mode checks the supplied evidence without requiring Git history.
