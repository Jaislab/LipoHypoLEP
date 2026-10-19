# Visual reference review

Reviewed on 19 October 2026, before rebuilding the analytical figure. The earlier slide-like diagrams and placeholder plots have been replaced. There is no single universal "Nature / Cell style" or publisher-specific font. The design decisions below combine direct observations of author-archived research figures with the user's explicit instructions: large Nimbus Sans type, exact full group names, and interpretive text in the external legend.

## Actual research figures inspected

### Nature Biotechnology: Milo

Dann et al., *Differential abundance testing on single-cell data using k-nearest neighbor graphs*, Nature Biotechnology, [DOI 10.1038/s41587-021-01033-z](https://doi.org/10.1038/s41587-021-01033-z).

- [Author manuscript figure assembly](https://github.com/MarioniLab/milo_analysis_2020/blob/main/ms/figures/milo_figures.pdf): pages 1 and 5 were rendered and visually inspected.
- [Author repository](https://github.com/MarioniLab/milo_analysis_2020) identifies the Milo manuscript analyses; [miloR](https://github.com/MarioniLab/miloR) supplies the publication attribution.
- [Figure 5 analysis source](https://github.com/MarioniLab/milo_analysis_2020/blob/main/notebooks/Fig5_liver_cirrhosis.Rmd) was inspected for plotting conventions.

Observed: a data-rich panel arrangement; short analytical axis labels; gray background points with color reserved for measured effects; restrained axes and compact legends; direct cell-type labels where useful; a heatmap annotation strip and small continuous colorbar. The manuscript's interpretive title and extensive caption are below the graphical panels, not a banner in the data field. The source uses classic plotting themes and explicitly removes panel gridlines. Some manuscript microtext is smaller than desirable here, so its density was not copied.

### Nature Methods: Cellpose and Cellpose3

- [Cellpose author Figure 1 preview](https://github.com/MouseLand/cellpose/blob/main/docs/_static/fig1.PNG), with [documentation establishing the figure's paper provenance](https://github.com/MouseLand/cellpose/blob/main/docs/index.rst).
- [Cellpose3 author figure utilities](https://github.com/MouseLand/cellpose/blob/main/paper/3.0/fig_utils.py).

Observed: white field, simple panel identifiers, short local technical labels, compact arrangement, and an external caption. The plotting utilities explicitly choose Arial, suppress top/right spines, use normal-weight axes titles, and define panel-letter styling. We retain **Nimbus Sans**, as required by the user; Arial in that one author's code is not treated as a journal mandate. The Cellpose author preview uses lowercase letters, whereas the Milo manuscript uses uppercase letters; neither case is universal.

### Cell: Dynamo

Qiu et al., *Mapping transcriptomic vector fields of single cells*, Cell (2022), [publisher article](https://www.sciencedirect.com/science/article/pii/S0092867421015774).

- Author-archived final-paper components: [Figure 6B](https://raw.githubusercontent.com/aristoteleo/dynamo-release/master/docs/introduction/dynamo_paper_figures/fig6_b.png), [Figure 6C](https://raw.githubusercontent.com/aristoteleo/dynamo-release/master/docs/introduction/dynamo_paper_figures/fig6_c.png), [Figure 7A](https://raw.githubusercontent.com/aristoteleo/dynamo-release/master/docs/introduction/dynamo_paper_figures/fig7_a.png).
- [Author documentation attributing those panels to the paper](https://github.com/aristoteleo/dynamo-release/blob/master/docs/introduction/index_cellfate.md).
- Complete [Figure 2](https://raw.githubusercontent.com/aristoteleo/dynamo_submission/master/dynamo-main-si-figures/Fig2.pdf) and [Figure 6](https://raw.githubusercontent.com/aristoteleo/dynamo_submission/master/dynamo-main-si-figures/Fig6.pdf) from the authors' **2020 submission archive** were also rendered and inspected. Their numbering and typesetting are not claimed to be the final 2022 journal versions.

Observed: coordinated repeated axes, large data area, thin axes or matrix rules, compact keys for actual encodings, and category or variable labels rather than editorial design claims. Some schematic panels do have scientific headings; the evidence does not imply that Cell prohibits all headings. Here, the user's stronger instruction determines the sparse analytical labeling.

## Publisher-guidance access limitation

Direct access to the [Nature final-submission guidance](https://www.nature.com/nature/for-authors/final-submission), [Nature formatting guide](https://www.nature.com/nature/for-authors/formatting-guide), and [Cell figure guidelines](https://www.cell.com/figureguidelines) was blocked by the environment proxy with `Tunnel connection failed: 403 Forbidden`. No current numeric specifications or mandatory font sizes are attributed to unread guidance.

A [secondary translated archive of Cell author instructions](https://github.com/ixxmu/mp_duty/issues/5964), dated 19 November 2024, was accessible. It puts figure legends in a separate manuscript section, with a short figure title and panel descriptions. This is corroborating secondary evidence, not a verified current publisher policy. Its 12–16 pt Arial passage concerns **graphical abstracts**, so it was not applied to research plots.

No copyrighted paper figures or downloaded article PDFs are redistributed in this repository. The links document the references actually examined.

## Decisions implemented

- Three separate question-led figures, each with four populated analytical panels: a volcano plot, a two-group expression heatmap, native-mouse KEGG preranked GSEA, and native-mouse GO Biological Process preranked GSEA. Gene-effect interval tables remain available outside the canvas. All expression statistics are synthetic; enrichment uses real memberships with randomly assigned demonstration gene symbols.
- Figure 1 compares LD SALINE with NON-LD control; Figure 2 compares LD LEPTIN with LD SALINE; Figure 3 compares LD LEPTIN with NON-LD control. PCA is omitted.
- All three figures retain the original synthetic dataset without effect-size rescaling. It is not empirically calibrated to brain tissue and does not predict experimental effects or power. The approximate log-scale OLS model is limited to this demonstration; a production reader/count-model pipeline for real data has not been implemented or tested.
- Volcano highlights and heatmap gene selection require both adjusted *P* < 0.05 and absolute log2 fold change ≥ 0.2. Vertical lines at ±0.2 mark effect cutoffs; zero remains a direction reference. GSEA uses the complete 7,647-gene signed-statistic ranking without those filters.
- No page title, subtitle, promotional label, explanatory footer, interpretation box, cohort badge, or generic RNA-seq icon inside the artwork. Questions, figure numbers, interpretation, and full methods appear in external legends and repository documentation.
- Panel letters, full contrast names, axes, heatmap group annotations, pathway names, q values, and the color scale provide the required analytical labeling. Thin gray left/bottom axes; no decorative chart containers or background grids.
- Comparison-specific warm/cool pairs follow the explicit RGB values in the supplied color-scheme guide: Figure 1 yellow `#FDB338` / blue `#025196`; Figure 2 tan `#FFBE6A` / turquoise `#40B0A6`; Figure 3 orange `#EB6123` / purple `#512888`. Heatmaps run cool–white–warm, and group strips use the first group's warm and second group's cool color. Full names identify cohorts because colors vary across comparisons. See the [color-accessibility review](accessibility/COLOR_ACCESSIBILITY.md).
- GSEA panels show two strongest NES values in each direction regardless of q. Filled versus hollow markers distinguish empirical GSEA FDR q < 0.05 from q ≥ 0.05; full result tables remain available. No gene sets, memberships, or random seeds are chosen to support a biological narrative.
- Nimbus Sans throughout, retaining large labels and panel identifiers with smaller mathematical subscripts. These are the project's design choices, not claimed journal requirements. [Export verification](FONT_VERIFICATION.md) compares embedded font metrics and glyph outlines to source Nimbus Sans fonts, rather than inferring identity from family settings.
- Matching volcano scales and repeated panel types support comparison across figures. Each heatmap uses its own contrast's DEG selection and only the two relevant groups. The GSEA calculations retain all tested genes and use documented native-mouse annotation versions.
- Synthetic-data provenance, sample counts, assumed distributions, statistical methods, cutoffs, annotation limitations, uncertainty, and interpretation appear in the [shared methods](FIGURE_LEGEND.md), [Figure 1 legend](FIGURE_1_LEGEND.md), [Figure 2 legend](FIGURE_2_LEGEND.md), [Figure 3 legend](FIGURE_3_LEGEND.md), and repository documentation.
- The post-treatment figure describes differences present after treatment, without assigning all of them to uncorrected LD effects or treating nonsignificance as normalization.
