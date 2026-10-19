# Muted color accessibility check

The original palette was **not uniformly colorblind safe**. The expression colors and categorical group colors serve different purposes and were checked separately. This review used the user's supplied *Colorblind Safe Color Schemes*, NCEAS Science Communication Resource Corner, Alexandra Phillips, updated 29 June 2022, as reference material, not as additional user instructions.

The guide recommends tested palettes, color-vision simulation, and information encoded through labels or position as well as color. Those principles informed this check. The attachment itself is not redistributed. Its description of deuteranopia as red-blindness is inaccurate: deuteranopia concerns M-cone/green deficiency; protanopia concerns L-cone/red deficiency.

## Findings and adjustment

The original group colors had nearly identical CIE lightness, L* approximately 64. They became almost the same gray. In full-severity Machado simulations, the original LD SALINE/LD LEPTIN pair became close under protanopia, and LD LEPTIN/NON-LD control became close under tritanopia. This does not mean that every reader would confuse them, but it prevents an unqualified accessibility claim.

The muted aesthetic is preserved with a modest lightness adjustment:

| Group | Previous | Revised | Revised L* |
|---|---|---|---|
| LD SALINE | `#C28E85` | **`#C28E85`** | 63.76 |
| LD LEPTIN | `#80A596` | **`#527D6B`** | 48.91 |
| NON-LD control | `#8B9FB3` | **`#A8BED2`** | 75.95 |

Full group names and a physical gap between sample blocks remain essential. They identify groups independently of hue. The revised group swatches are also more distinguishable in grayscale.

The muted **blue `#779AAC` / coral `#C28A80` effect colors** are retained. Their hue separation remains visible under the three modeled deficiencies. In volcano and effect-estimate plots, direction is also encoded by horizontal position relative to zero. These effect colors have similar gray values, so grayscale printing is not equivalent to the color version, particularly for the diverging expression heatmap. A numeric colorbar and the supplied expression/effect tables remain available; the heatmap should be distributed in color.

## Evidence and limitations

- [Before/after palette simulations](palette_color_vision_check.png).
- Full revised Figure 1 under [protanopia](figure_01_protanopia_check.png), [deuteranopia](figure_01_deuteranopia_check.png), [tritanopia](figure_01_tritanopia_check.png), and [grayscale](figure_01_grayscale_check.png).
- [Numerical color results and method metadata](palette_simulation.json).

The check uses DaltonLens 0.1.5, `Simulator_Machado2009`, severity 1.0. The library decodes sRGB into linear light for the simulation and encodes it afterward. Grayscale uses linear-light relative luminance. CIE1976 Lab distances under D65 are descriptive comparisons, not a universal accessibility cutoff. The full rendered Figure 1, including point transparency and the entire heatmap ramp, was simulated rather than evaluating endpoints alone.

These are modeled accessibility checks, not certification for every type or degree of color-vision difference, visual impairment, display, print process, or viewing condition.

## Reproduce

With the normal figure dependencies installed:

```bash
python -m pip install -r figures/accessibility/requirements.txt
python figures/accessibility/audit_colors.py
```

The audit outputs are review material, separate from the manuscript figures.
