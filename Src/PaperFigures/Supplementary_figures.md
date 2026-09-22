# Supplementary figures

All panels produced by the four figure notebooks that are **not** part of main Figure 4,
organised into the two requested groups. Numbering is a single sequence (S01-S17) so it
can be cited directly, with the group recorded in both the number block and the filename:

- **Group 1 — Data quality control:** S01-S03, filenames `FigS0N_QC_*`
- **Group 2 — Complementary to the presented analysis:** S04-S17, filenames `FigSNN_Comp_*`

Filenames are assigned centrally by `PAPER_PANELS` in `_figutils.py`; the notebooks keep
their own internal panel names, so re-running any notebook writes the published names
automatically. Shared methods (cell counts, normalisation, state definitions, pseudotime
construction) are in `Figure4_legends.md` and are not repeated here.

Six panels the notebooks used to produce are **not** published: the plain pseudotime
UMAP and the CDKN1A UMAP, the pseudotime density plot, both fixed-width-window versions of
the Doxo1-vs-pseudotime curves, and the mixed-model forest plot. Their plotting code has
been removed from the notebooks, so every figure a notebook now draws is written as a PDF.
Their `PAPER_PANELS` entries are kept at `None` so that restoring any of that code cannot
silently emit an unnumbered file. The numbering below is therefore contiguous and stable.

Supplementary figures S1-S15 support main Figure 4; S16-S17 support Figure 5.
The NEUROG1+SIM1 senescence-gene panel that was S16 has been promoted to a main
panel of Figure 5 and is written as `Fig5G_senescence_interaction_NEUROG1_SIM1_day10.pdf` (Figure 5G)
(legend in `Figure5_legends.md`).

Three of the published panels are the day04 counterparts of day10 main panels and are
marked **[day04 of 4X]**.

---

## Group 1 — Data quality control

### Figure S1. Sequencing depth and gene complexity per cell
`FigS01_QC_sequencing_depth.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

Distribution of total UMIs per cell (left) and genes detected per cell (right) for all
291,301 cells. Each histogram has a log-scaled y axis and carries a horizontal boxplot
above it: box = interquartile range, central line = median, diamond = mean, whiskers =
1.5x IQR; solid and dashed vertical lines mark the median and mean. Median 5,388 UMIs per
cell (mean 6,247; IQR 3,701-7,779) and 2,141 genes per cell (mean 2,266; IQR 1,611-2,782).

### Figure S2. Library complexity, cell recovery and cell-cycle composition
`FigS02_QC_library_composition_panels.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

Six quality-control panels. **(a)** UMIs versus genes detected per cell (hexbin, log cell
density). **(b)** Cells recovered per time point (day04 146,083; day10 145,218).
**(c)** Cells carrying one versus two guide assignments per time point (day04 91,472 versus
54,611; day10 85,242 versus 59,976), i.e. single- versus double-knockout recovery.
**(d)** Cells per transcriptional state (NEPC-N 181,079; NEPC-A&N-1 81,545; NEPC-A&N-2
7,427; NEPC-A&N-3 2,341; NEPC-A-1 15,321; NEPC-A-2 3,588). **(e)** UMIs per cell by state:
median depth spans only 5,028-5,787 across the six states, so sequencing depth does not
track state and cannot explain the state assignment. **(f)** Cell-cycle phase composition
per state, split by time point (left bar day04, right bar day10).

### Figure S3. Cells recovered per perturbation
`FigS03_QC_cells_per_perturbation.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

Cells recovered for each of the 92 perturbation groups, one panel per time point,
log-scaled. Labels use the collapsed convention in which a single knockout carrying a
non-targeting partner guide (`X+NTC`) counts as `X`. Bars are coloured by class (control,
single knockout, double knockout) and ordered by total abundance across both time points,
so the panels share an x ordering. NTC is the largest group (7,550 cells at day04; 7,001 at
day10); the smallest double knockouts retain roughly 230-250 cells per time point.

---

## Group 2 — Complementary to the presented analysis

### Figure S4. Cell-cycle phase across the transcriptional landscape
`FigS04_Comp_umap_cellcycle_phase.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

UMAP coloured by inferred cell-cycle phase (scanpy `score_genes_cell_cycle`), showing that
the transcriptional states of Figure 4E are not organised by proliferative state. Across
all cells, 43.6% are assigned G2M, 31.0% G1 and 25.4% S — a heavily cycling population, as
expected for this line, so phase is a covariate to rule out rather than a quality metric.

### Figure S5. Trajectory signature scores
`FigS05_Comp_umap_state_signature_scores.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

UMAPs of the seven signature scores carried in the object: one per trajectory state
(`Neuroendocrine_score`, `Intermediate-1/-2/-3_score`, `Differentiated-1/-2_score`) plus
`Doxo-1-differentiated_score`. Colour scales are clipped at the 1st and 99th percentiles.
Each signature localises to the region of the landscape occupied by its corresponding
state. Note these are the upstream signature scores; the `Doxo1program_score` used in
Figures 4B, 4F and 4J is scored separately in `00_Prepare_data.ipynb`.

### Figure S6. Marker gene expression per cell state
`FigS06_Comp_state_marker_dotplot.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

Dot plot of neuroendocrine markers (ASCL1, INSM1, CHGA, CHGB, SYP, NEUROD1), adeno-like and
epithelial markers (KRT8, KRT18, EPCAM, VIM, TWIST1) and CDKN1A across the six states in
trajectory order. Dot size is the fraction of expressing cells, colour the mean expression
scaled per gene. Establishes state identity from marker expression independently of the
Doxo1 signature. **NEUROG1 is requested by the code but is absent from the expression
matrix** (see author note 8), so 12 of the 13 intended markers are drawn.

### Figure S7. Pseudotime distribution per cell state
`FigS07_Comp_pseudotime_by_state.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

Violin plots of diffusion pseudotime within each state, drawn in the manuscript state
order. The states occupy overlapping but distinguishable pseudotime ranges. Median
pseudotime is **not** monotone in the manuscript ordering: NEPC-A&N-3 0.068, NEPC-N 0.077,
NEPC-A&N-1 0.137, NEPC-A&N-2 0.151, NEPC-A-2 0.174, NEPC-A-1 0.442. Two departures are
visible — NEPC-A&N-3 sits marginally *before* the neuroendocrine root, and NEPC-A-1 is far
later than NEPC-A-2 despite its lower index. See author note 7.

### Figure S8. CDKN1A expression along pseudotime
`FigS08_Comp_CDKN1A_vs_pseudotime.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

Log-normalised CDKN1A (p21) expression as a function of diffusion pseudotime, separately
for day04 and day10, as an independent readout of cell-cycle arrest and differentiation.
Line = local mean over a sliding window of consecutive pseudotime-ordered cells (window =
the larger of 10% of the group's cells or 50 cells); band = +/- 1 SD of single-cell values in
the window, i.e. cell-to-cell spread, not a confidence interval. Background shading marks
approximate state territories (see author note 7).

### Figure S9. Cell-state composition of control versus perturbed cells
`FigS09_Comp_composition_control_vs_perturbed.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

For each time point, the percentage of control cells (no target-gene guide) and of perturbed
cells (at least one target knockout) residing in each state; within each group and time
point the bars sum to 100%. Asterisks mark a two-sided Fisher exact test of perturbed versus
control occupancy of that state, Benjamini-Hochberg corrected across all 12 state x time
point tests (* FDR < 0.05, ** < 0.01, *** < 0.001, ns otherwise). One state per time point
reaches FDR < 0.05 (minimum FDR 0.024 at day04, 0.005 at day10), so the perturbed pool as a
whole is only mildly redistributed relative to control — the strong shifts are
perturbation-specific and appear in Figures 4G and 4H.

### Figure S10. Time point and state composition
`FigS10_Comp_timepoint_state_composition.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

Left: cell-state composition within each time point (each bar sums to 100%). Right: the
day04 and day10 share of each state (each bar sums to 100%). Together these show how the
landscape is populated over time and that no state is exclusive to one time point.

### Figure S11. Doxo1 score along pseudotime, day04
`FigS11_Comp_doxo1_vs_pseudotime_day04.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb — **[day04 of 4J]**

As Figure 4J but for day04. Perturbation set selected by the same rule (enriched in
NEPC-A-1/-2 at that time point, one-sided Fisher, odds ratio > 1, FDR < 0.1) plus
NEUROG1+SIM1, included on prior experimental validation rather than by the enrichment
criterion. Eleven perturbations qualify at day04 versus 19 at day10; per-curve n =
413-1,031 cells against NTC n = 7,550. Smoother and state shading as in Figure 4J; no
uncertainty band is drawn (see author note 4).

### Figure S12. Decomposition of double-knockout effects on the Doxo1 score, day04
`FigS12_Comp_interaction_heatmap_day04.pdf` · Fig4K_S12-S14__Interaction_Doxo1.ipynb — **[day04 of 4K]**

As Figure 4K but for day04: for each of the 11 day04 double knockouts, the main effect of
each constituent gene, their interaction, and the total double-knockout effect. Main and
interaction terms come from `Doxo1 ~ g1 + g2 + g1:g2` fit on NTC plus cells whose only
target guides are g1 and/or g2, restricted to the non-neuroendocrine compartment; the total
comes from a separate `Doxo1 ~ combo` model. Double-knockout cells per combination range
from 152 to 422. Two-sided Wald p-values, Benjamini-Hochberg corrected across combinations
within each column; * FDR < 0.1, ** < 0.01, *** < 0.001. At day04, 10 of 11 interaction
terms and 11 of 11 total effects reach FDR < 0.1, while no single-gene main effect does.

### Figure S13, S14. Perturbation to cell state to Doxo1: two routes to a higher score
`FigS13_Comp_state_vs_doxo1_day04.pdf`, `FigS14_Comp_state_vs_doxo1_day10.pdf`
Fig4K_S12-S14__Interaction_Doxo1.ipynb

For each double knockout and each state, two regressions fit fresh per displayed cell.
Left: **composition shift**, logistic `I(cell in state) ~ combo`, colour = log-odds of a
combination cell occupying that state versus NTC. Right: **within-state effect**, ordinary
least squares `Doxo1 ~ combo` restricted to that state, colour = change in Doxo1 score
versus NTC. One-sided tests (combination > NTC), Benjamini-Hochberg corrected across the
displayed grid; * marks FDR < 0.1, a centred dot marks a cell that cannot be fit (fewer
than 3 combination cells, or logistic separation) and is shaded grey. NEUROG1+SIM1 is boxed.

The two routes dissociate. Composition shifts are common — 13 of 66 cells at day04 and 21
of 114 at day10 reach FDR < 0.1, almost all of them entries into NEPC-A-1 (log-odds +0.39
to +0.75) — whereas within-state increases are almost absent: **none** at day04 and exactly
**one** at day10, NEUROG1+SIM1 in NEPC-A-2 (beta = +0.073, FDR = 0.035, n = 9 combination
cells versus 67 NTC). NEUROG1+SIM1 itself shifts composition into NEPC-A-1 at day04
(log-odds +0.68, FDR = 0.002) but shows no significant composition shift at day10, making it
the one perturbation that raises the score within a state rather than by relocating cells.

### Figure S15. Gene-level effects of the double knockouts, day04
`FigS15_Comp_perturbation_gene_heatmap_day04.pdf` · Fig4L_Fig5_S15__Perturbation_Gene_Effects.ipynb
**[day04 of 4L]**

As Figure 4L but for day04: log2 fold-change versus NTC for the 11 day04 double knockouts
across 84 genes, from per-perturbation differential expression with FDR recomputed by
Benjamini-Hochberg against a fixed test universe of 8,000 genes. Entries with FDR >= 0.1 are
set to exactly zero before display and clustering. Twenty-four senescence/SASP genes are
present (red column strip and red bold labels). All surviving changes at day04 are
increases (log2 fold-change +0.16 to +3.92). Rows and columns are clustered by Ward linkage
on euclidean distances of the thresholded matrix (see author note 9).

### Related to Figure 5

### Figure S16. Doxo1 score along pseudotime for Wilcoxon-significant perturbations
`FigS16_Comp_doxo1_vs_pseudotime_wilcoxon_hits.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb — **[related to Figure 5]**

Doxo1 score against diffusion pseudotime for the perturbations significant by a single-cell
one-sided Wilcoxon test (perturbation > NTC) in at least one time point x cell state
stratum, Benjamini-Hochberg corrected within each stratum. Both time points are drawn on
one axis (colour = perturbation, line style = time point). This panel supports **Figure 5**
rather than Figure 4.

### Figure S17. NEUROG1 x SIM1 interaction on the Doxo1 score, per cell state and timepoint
`FigS17_Comp_NEUROG1_SIM1_interaction.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb — **[related to Figure 5]**

Coefficients of the saturated model `Doxo1 ~ NEUROG1 + SIM1 + NEUROG1:SIM1`, fit on NTC
cells plus cells whose only target guides are NEUROG1 and/or SIM1 (NTC = reference). The
model is fit **separately within each cell state**, and in the final row on the pooled
adeno-like compartment (NEPC-A-1 + NEPC-A-2). One heatmap per timepoint; columns are the
individual effect of each gene, their interaction, and the total double-knockout effect
(the sum of the three, tested as a linear combination). Row labels give the number of
double-knockout cells contributing to that fit. Rows with fewer than 5 double-knockout
cells, or fewer than 5 NTC or single-knockout reference cells, are not identifiable and are
marked `n.d.` on grey.

Stars are Benjamini-Hochberg FDR across all fitted cells in the figure (* FDR < 0.1,
** < 0.01, *** < 0.001). The interaction is significant in exactly one stratum —
**NEPC-A-2 at day10** (beta = +0.072, p = 0.007, FDR = 0.078, n = 9 double-knockout cells)
— matching the within-state estimate in S14; every other state is null, and the pooled
adeno-like fit is not significant (day10 beta = +0.023, FDR = 0.39; day04 beta = -0.005).
This is why the effect is diluted in panels that pool states, e.g. Figure 4K
(beta = +0.005 over five states). Per-stratum statistics, including cell counts, are
exported to `figures/FigS17_NEUROG1_SIM1_interaction_by_state.csv`.


---

## Author notes — NOT for publication

1. **Grouping.** S4 (cell-cycle phase UMAP) sits in group 2: it supports the cell-state
   definitions rather than reporting data quality. Group 1 is therefore the three panels
   that report library, depth and recovery statistics only. S7 (pseudotime per state) is
   likewise a method validation and stays in group 2. Moving a figure between groups is a
   one-line change in `PAPER_PANELS`.

2. **S17 now stars by FDR.** It previously showed raw two-sided p-values while every
   other statistical panel used Benjamini-Hochberg FDR. Expanding it to a per-state grid
   multiplied the number of tests, so it now BH-corrects across the 22 fitted strata within
   the panel. The NEPC-A-2 day10 interaction survives that correction (FDR = 0.078);
   nothing else does. Correction is within the panel, not across panels.

3. **S9 could not previously be regenerated.** The cell that produces it called
   `multipletests` without importing it, so it raised a `NameError` on any clean kernel
   run; the existing PDF came from a session where the name happened to be defined. I
   added the import to that cell, so it now regenerates.

4. **Nothing now shows where 4J's curves are supported.** The two fixed-width-window
   panels were dropped, so no published figure shows that 34-56% of the pseudotime axis
   holds fewer than 3 cells per perturbation while 4J draws a continuous line across it.
   The plotting code is no longer in the notebook; regenerating it means restoring the
   `fixed_window` smoother over `fu.pseudotime_trend_fixedwindow`, which is still in
   `_figutils.py`, and giving both `fixedpt` entries in `PAPER_PANELS` a name.

5. **The mixed-model result is now unpublished.** `Doxo1 ~ C(perturbation) + (1 | state)`
   found 0 of 91 perturbations significant at either time point (largest coefficients
   +0.0029 at day04, +0.0035 at day10; smallest FDR 0.98 and 0.999). That is the
   statistical support for the "redistribution, not within-state elevation" framing used
   in the Figure 4 title and the 4J legend. The model is still fit in the Pseudotime
   notebook and still writes `figures/Fig2_doxo1_mixedlm_state.csv`; only its forest plot
   was removed, so the check stays reproducible without producing an unpublished panel.

6. **Run order.** The Pseudotime notebook
   (`Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb`) must run before
   `Fig4K_S12-S14__Interaction_Doxo1.ipynb` and
   `Fig4L_Fig5_S15__Perturbation_Gene_Effects.ipynb`, which read
   `figures/Fig2_state_enrichment_{day}.csv` to define their perturbation sets. This
   affects S12-S15 as well as 4K/4L and the Figure 5 senescence panel.

7. **The shaded "state territories" are not a valid partition.** `shade_state_regions`
   places band edges at the midpoints between the median pseudotimes of *consecutive
   states in the manuscript ordering*, but those medians are not monotone in that ordering
   (S7). The resulting edges are [0.000, 0.107, 0.144, 0.110, 0.255, 0.308, 1.000] — the
   third is greater than the fourth, so the NEPC-A&N-3 band is drawn backwards and the
   bands are not a partition of the axis. This affects the background shading of Figures
   4F, 4J, S8 and S11. The fix is to order the bands by observed median pseudotime
   (NEPC-A&N-3 < NEPC-N < NEPC-A&N-1 < NEPC-A&N-2 < NEPC-A-2 < NEPC-A-1) rather than by
   `STATE_ORDER`; the palette comment in `_figutils.py` already acknowledges this true
   ordering for the intermediate states. Not yet changed, since it alters main panels.

8. **NEUROG1, SIM1 and VSX1 are absent from the expression matrix.** The processed object
   carries 23,423 genes and these three screened targets are not among them (filtered
   upstream, before `ComboScreen.h5ad`). Consequences: S6 silently drops NEUROG1 from the
   dot plot; knockdown of these three targets cannot be confirmed from expression in this
   object; and they can never appear as columns in 4L or S15. The perturbation analyses are
   unaffected, since they identify perturbations from guide one-hot columns rather than
   from expression — but the headline combination NEUROG1+SIM1 is one whose target genes
   are both unmeasurable here, which is worth stating if a reviewer asks for knockdown
   validation.

9. **4L and S15 cluster thresholded values.** Non-significant log2 fold-changes are set to
   exactly zero before Ward clustering, so the dendrograms group perturbations by which
   genes passed FDR as much as by effect size. Describe the clustering as being on
   thresholded values.

## ATAC (chromatin) supplementary figures — Supplementary Figure S5

Generated by the ATAC notebooks in `Src/PaperFigures/` (data under `Src/ComboScreen/ATAC/`).
These currently use an independent **S5A–H** panel scheme, which collides with the RNA **S05**;
renumber into the single S-sequence when finalising.

### Figure S5A. scATAC cell calling
`FigS5A_QC_CellCalling.pdf` · FigS5_QC_ATAC.ipynb
TSS enrichment vs. log10 unique fragments over all 8.17M barcodes (2-D density); dashed lines mark
the retained-set lower bounds. High-quality nuclei (upper right) separate from empty/ambient
barcodes (lower left); 280,177 cells pass.

### Figure S5B. Per-cell ATAC QC distributions
`FigS5B_QC_Distributions.pdf` · FigS5_QC_ATAC.ipynb
TSS enrichment, log10 unique fragments, and FRIP (fraction of reads in peaks) over the retained
cells — medians 11.8, 6,734, and 0.29, respectively.

### Figure S5C. Aggregate TSS-enrichment profile
`FigS5C_QC_TSSprofile.pdf` · FigS5_QC_ATAC.ipynb
Insertion signal ± 2 kb around transcription start sites (fold-enrichment over flanking
background). The sharp central peak and +1-nucleosome shoulder indicate high-quality,
promoter-focused signal. (Fragment-size/nucleosome QC is omitted: the width-1 insertion-site input
makes those ArchR metrics degenerate.)

### Figure S5D. Pseudotime segments across modalities
`FigS5D_PseudotimeSegments_UMAP.pdf` · FigS5D_PseudotimeSegments_UMAP.ipynb
The 16 adaptive even-trajectory pseudotime segments (as used for peak calling) shown on the
ATAC-profile UMAP (left) and the RNA UMAP (right); segments are defined once from pseudotime and
colored identically, so matching layouts indicate the trajectory bins occupy concordant
territories in chromatin and transcriptome.

### Figure S5E. TF motif accessibility across pseudotime — full heatmap
`FigS5E_PseudotimeTFActivity_heatmap.pdf` · FigS5E_PseudotimeTFActivity_heatmap.ipynb
All TFs significant (FDR ≤ 0.1) in ≥ 1 segment, ordered by Spearman trend (rising top → falling
bottom); columns root NE → differentiated tip. Asterisks, FDR ≤ 0.1. Main Fig 5D shows the top
rising/falling TFs as line traces.

### Figure S5F. TF motif accessibility per genetic perturbation (pooled)
`FigS5F_PertTFActivity_heatmap.pdf` · FigS5F_PertTFActivity_heatmap.ipynb
Each perturbation vs. pooled non-targeting-control cells (no state stratification), pseudotime
peak set; per-peak log2 fold-change regressed on the motif matrix. Non-zero FDR-significant
(≤ 0.1) coefficients (TFs × perturbations), rows/columns hierarchically clustered, cell values =
coefficients. Few perturbations reach significance at this marginal, chromatin-motif readout.

### Figure S5G. TF motif accessibility per cell state × perturbation
`FigS5G_StateXPert_TFActivity.pdf` · FigS5G_StateXPert_TFActivity.ipynb
Each perturbation vs. same-state NTC; per-peak log2 fold-change regressed on the motif matrix.
Non-zero FDR-significant (≤ 0.1) coefficients (all conditions with ≥ 1 significant TF): columns
grouped by cell state (trajectory order, annotated by the colored strip/labels above; perturbation
on the x-axis), TF rows hierarchically clustered.

### Figure S5H. TF motif accessibility vs. RNA expression across pseudotime
`FigS5H_TFActivity_vs_mRNA.pdf` · FigS5H_TFActivity_vs_mRNA.ipynb
Remaining individual TFs (the four examples in main Fig 5E, and motif groups, are excluded;
SNAI2/IKZF1/LYL1 are absent from the RNA reference). TF motif accessibility coefficient (vs. NE
root) vs. mean mRNA, one point per pseudotime segment (colored NE root → differentiated tip);
Spearman ρ per TF. Concordance is TF-specific — some TFs track their expression while others
(repressors, or factors whose motif is bound by other family members) move oppositely.
