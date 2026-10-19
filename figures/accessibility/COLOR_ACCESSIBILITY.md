# Comparison-specific color accessibility

The palettes follow the user's explicit request for three pairs from *Colorblind Safe Color Schemes*, NCEAS Science Communication Resource Corner, Alexandra Phillips, updated 29 June 2022. The RGB codes **printed beside the guide's swatches** are used. The guide is reference material rather than additional user instructions, and the attachment itself is not redistributed.

| Figure | Comparison | Positive / higher / first group | Negative / lower / second group |
|---|---|---|---|
| 1 | LD SALINE vs NON-LD control | Yellow, RGB 253/179/56, `#FDB338` | Blue, RGB 2/81/150, `#025196` |
| 2 | LD LEPTIN vs LD SALINE | Tan, RGB 255/190/106, `#FFBE6A` | Turquoise, RGB 64/176/166, `#40B0A6` |
| 3 | LD LEPTIN vs NON-LD control | Orange, RGB 235/97/35, `#EB6123` | Purple, RGB 81/40/136, `#512888` |

Colors are **comparison-specific**, rather than a single global cohort code across figures. Direct full group names identify every heatmap block. In each volcano, warm/cool color identifies positive/negative effects meeting adjusted P < 0.05 and |log2FC| ≥ 0.2; position relative to zero redundantly encodes direction. In each heatmap, the cool–white–warm scale encodes standardized gene expression. In GSEA, warm/cool identifies positive/negative NES; filled circles meet GSEA FDR q < 0.05 and open circles do not. These significance encodings use different multiple-testing procedures, as explained in the legends.

The requested pairs were checked under full-severity Machado simulations of protanopia, deuteranopia and tritanopia. The full rendered figures, not just isolated swatches, are included in the checks. This covers point transparency and heatmap gradients. Simulations are a useful accessibility review, not a guarantee for every viewer or reproduction condition. Keep full labels, axes, zero lines, numeric colorbars, q-value columns and downloadable tables; do not rely on hue alone.

## Review artifacts

- [All three requested palettes under color-vision simulations](palette_color_vision_check.png).
- Figure 1: [protanopia](figure_01_protanopia_check.png), [deuteranopia](figure_01_deuteranopia_check.png), [tritanopia](figure_01_tritanopia_check.png), [grayscale](figure_01_grayscale_check.png).
- Figure 2: [protanopia](figure_02_protanopia_check.png), [deuteranopia](figure_02_deuteranopia_check.png), [tritanopia](figure_02_tritanopia_check.png), [grayscale](figure_02_grayscale_check.png).
- Figure 3: [protanopia](figure_03_protanopia_check.png), [deuteranopia](figure_03_deuteranopia_check.png), [tritanopia](figure_03_tritanopia_check.png), [grayscale](figure_03_grayscale_check.png).
- [Numerical colors, lightness and pairwise comparisons](palette_simulation.json).

The check uses DaltonLens 0.1.5, `Simulator_Machado2009`, severity 1.0. The library decodes sRGB into linear light and re-encodes the result. Grayscale uses linear-light relative luminance. CIE1976 Lab distances under D65 are descriptive comparisons, with no universal pass/fail threshold. Grayscale previews document remaining limitations and do not establish that every heatmap value can be recovered without color.

The guide's description of deuteranopia as red-blindness is inaccurate: deuteranopia concerns M-cone/green deficiency, whereas protanopia concerns L-cone/red deficiency. This terminology issue does not alter the requested RGB pairs.

## Reproduce

With the normal figure dependencies installed:

```bash
python -m pip install -r figures/accessibility/requirements.txt
python figures/accessibility/audit_colors.py
```

The audit outputs are separate review material, outside the manuscript figure canvases.
