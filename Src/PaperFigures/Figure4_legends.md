# Figure 4 — panel legends

Draft legends for the ten panels currently assigned to Figure 4. Each entry gives the
file, the notebook that writes it, what is plotted, and the statistics as actually
implemented. Numbers were read off the saved result tables, not estimated.

Author notes that need reconciling with the manuscript text are in the final section
and are **not** intended for publication.

---

## Figure title

**Figure 4. Combinatorial transcription factor knockouts shift NEPC cells along a
neuroendocrine-to-adeno-like differentiation trajectory, with non-additive interactions
between gene pairs and induction of senescence-associated genes.**

Shorter alternative, if the journal caps title length:
**Figure 4. Paired transcription factor knockouts drive neuroendocrine-to-adeno-like
state transitions in NEPC through non-additive gene interactions.**

The title deliberately claims *redistribution along the trajectory* rather than
elevation of the Doxo1 score — see author note 2. "Non-additive" covers the 14 of 19
significant interaction terms in 4K; "senescence-associated" covers the 31 senescence
genes in 4L.

---

## Shared methods (for the figure preamble or Methods)

Combinatorial CRISPR knockout screen of 13 transcription factors plus non-targeting
controls (NTC) in NEPC cells, assayed by scRNA-seq at two time points. After removing
cells with 3 or more guide assignments, **291,301 cells** were analysed (day04
n = 146,083; day10 n = 145,218), spanning **92 perturbation groups** (13 single
knockouts, 78 double knockouts, NTC). Counts were normalised to 20,000 UMIs per cell
and log1p-transformed; 3,000 highly variable genes, scaled (max 9), 50 principal
components, k = 8 (euclidean) neighbours, UMAP.

Six transcriptional states were defined from the whole transcriptome, independently of
any differentiation signature, and are named along the neuroendocrine-to-adeno axis:
**NEPC-N** (n = 181,079), **NEPC-A&N-1** (81,545), **NEPC-A&N-2** (7,427),
**NEPC-A&N-3** (2,341), **NEPC-A-1** (15,321), **NEPC-A-2** (3,588).

The **Doxo1 differentiation score** is a mean-centred signature score (scanpy
`score_genes`, control-gene binned) over the Doxo-1 differentiated DE gene set.
**Diffusion pseudotime** was computed on 15 diffusion components (10 used for DPT),
rooted at the diffusion-space medoid of the cells in the lowest 1% of Doxo1 score and
oriented to increase from NEPC-N to NEPC-A; values span 0 to 1.

---

## 4A
Not assigned from the analysis notebooks — presumably a schematic or experimental
design panel. (Panel letter I is skipped, per convention, to avoid confusion with 1.)

## 4B — Doxo1 differentiation score across the transcriptional landscape
`Fig4B_umap_doxo1.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

UMAP of all 291,301 cells coloured by the Doxo1 differentiation score (magma; colour
scale clipped at the 1st and 99th percentiles). The score increases across the
landscape from the neuroendocrine compartment toward the adeno-like compartment.

## 4C — Differentiation trajectory backbone
`Fig4C_pseudotime_backbone.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

UMAP coloured by diffusion pseudotime (viridis), with the trajectory backbone overlaid.
The backbone is constructed by dividing cells into 30 equal-occupancy pseudotime
windows, taking the median UMAP coordinate of each window, and joining these in
pseudotime order (centred 3-window rolling mean for smoothing). Arrowheads indicate
increasing pseudotime; the filled circle marks the root cell and the star the terminal
end of the trajectory.

## 4D — Time point
`Fig4D_umap_timepoint.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

UMAP coloured by time point (day04, blue; day10, orange), showing that both time points
populate the full landscape.

## 4E — Transcriptional cell states
`Fig4E_umap_cellstate.pdf` · Fig4DE_S01-S06_S09__Overview_UMAP.ipynb

UMAP coloured by the six transcriptional states (colour key), ordered along the
neuroendocrine-to-adeno axis. State assignment derives from whole-transcriptome
clustering and is independent of the Doxo1 signature used in 4B, 4F and 4J.

## 4F — Doxo1 score along pseudotime, by time point
`Fig4F_doxo1_vs_pseudotime.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

Doxo1 differentiation score as a function of diffusion pseudotime, separately for day04
and day10. Each line is a local mean over a sliding window of consecutive
pseudotime-ordered cells (window = the larger of 10% of that group's cells or 50 cells).
The shaded band is **+/- 1 standard deviation of single-cell scores within the window**,
i.e. the dispersion of cells, not a confidence interval on the mean. Background shading
marks approximate state territories, bounded by the midpoints between the median
pseudotimes of consecutive states.

## 4G, 4H — Perturbation enrichment across cell states (4G day04; 4H day10)
`Fig4G_state_enrichment_heatmap_day04.pdf`, `Fig4H_state_enrichment_heatmap_day10.pdf`
Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

For every perturbation with at least 30 cells at that time point and every cell state,
a one-sided Fisher exact test (enrichment, odds ratio > 1) compares that perturbation's
cells in the state against NTC cells in the state: 546 tests per time point
(91 perturbations x 6 states), Benjamini-Hochberg corrected within the time point.
Colour is the log2 odds ratio (clipped at +/- 6); asterisks mark FDR < 0.1. Rows are the
perturbations significant in at least one state (11 at day04, 20 at day10), ordered by
their maximum log2 odds ratio across the two adeno-like states NEPC-A-1 and NEPC-A-2.

## 4J — Doxo1 score along pseudotime for differentiated-state-enriched double knockouts (day10)
`Fig4J_doxo1_vs_pseudotime_rollcells_day10.pdf` · Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb

Doxo1 score against diffusion pseudotime at day10 for NTC (black, n = 7,001) and for the
19 double knockouts enriched in the adeno-like states NEPC-A-1/-2 (one-sided Fisher,
odds ratio > 1, FDR < 0.1; the set shown in 4H), plus **NEUROG1+SIM1** (red), which is
included on the basis of prior experimental validation rather than by the enrichment
criterion. Each line is a local mean over a sliding window of consecutive
pseudotime-ordered cells (window = the larger of 15% of that group's cells or 50 cells);
per-perturbation n ranges from 245 to 1,031 cells. No uncertainty band is shown.
Background shading as in 4F.

## 4K — Decomposition of double-knockout effects on the Doxo1 score (day10)
`Fig4K_interaction_heatmap_day10.pdf` · Fig4K_S12-S14__Interaction_Doxo1.ipynb

For each of the 19 double knockouts (rows), the effect on the Doxo1 score is decomposed
into the main effect of each constituent gene, their interaction, and the total
double-knockout effect (columns). Main and interaction effects come from the ordinary
least squares model `Doxo1 ~ g1 + g2 + g1:g2`, fit on NTC cells plus cells whose only
target guides are g1 and/or g2, restricted to the non-neuroendocrine compartment
(NEPC-A&N-1/-2/-3 and NEPC-A-1/-2); the total effect comes from a separate model
`Doxo1 ~ combo` on NTC cells plus that combination's cells. Combinations required at
least 10 double-knockout cells (observed range 103-397). Cell values are regression
coefficients (change in Doxo1 score relative to NTC); p-values are two-sided Wald tests,
Benjamini-Hochberg corrected across the 19 combinations separately within each column.
Asterisks: * FDR < 0.1, ** FDR < 0.01, *** FDR < 0.001. Rows are ordered by total effect.
At day10, 14 of 19 interaction terms and 11 of 19 total effects are significant at
FDR < 0.1; no single-gene main effect is.

## 4L — Gene-level effects of the double knockouts (day10)
`Fig4L_perturbation_gene_heatmap_day10.pdf` · Fig4L_Fig5_S15__Perturbation_Gene_Effects.ipynb

Log2 fold-change relative to NTC for the 19 double knockouts (rows) across 86 genes
(columns), from per-perturbation differential expression with FDR recomputed by
Benjamini-Hochberg against a fixed test universe of 8,000 genes. **Entries with
FDR >= 0.1 are set to exactly zero** before display and clustering, so colour reflects
significant changes only (observed range -0.57 to +3.49 log2 fold-change). Genes were
selected as significant (FDR < 0.1) in at least 2 perturbations, taking the top 60 ranked
by the number of perturbations affected and then by smallest FDR, plus any
senescence/SASP-associated gene significant in at least one perturbation (31 such genes
are present), plus CDKN1A, included by design. Rows and columns are clustered by Ward
linkage on euclidean distances of the thresholded matrix. The red column strip and red
bold gene labels mark senescence/SASP-associated genes.

---

## Author notes — NOT for publication

1. **4K vs the NEUROG1+SIM1 claim.** In the model plotted in 4K (pooled across the
   non-NE compartment), NEUROG1+SIM1 has interaction beta = +0.005, FDR = 0.50 and total
   effect = +0.005, FDR = 0.52 — not significant. The +0.073 effect quoted elsewhere is a
   different model: within NEPC-A-2 only, 9 combination cells vs 67 NTC cells
   (Mann-Whitney p = 3.7e-4). If the text cites +0.073, it must say it is the
   within-NEPC-A-2 estimate and give n = 9, otherwise it will appear to contradict 4K.

2. **4J wording.** With no uncertainty band, the panel invites the reading that these
   perturbations raise the Doxo1 score. At matched pseudotime (equal-occupancy deciles,
   95% CI on single-cell means) the shifts versus NTC are -0.003 to +0.002 and 0-2 of 10
   deciles separate (chance expectation ~0.5). The defensible claim is that these
   perturbations **redistribute cells along the trajectory** — consistent with the
   composition enrichment in 4H and the interaction effects in 4K — rather than raising
   the score at matched pseudotime. The legend above is worded accordingly.

3. **4F band.** Stated as +/- 1 SD above, because that is what the code draws. If a
   confidence interval is preferred, `_figutils.plot_pseudotime_trend` produces
   mean +/- 1.96 SE with the same smoother.

4. **Circularity to disclose.** Pseudotime is rooted on the Doxo1 score itself, so the
   monotone rise seen in 4B, 4C and 4F is partly definitional. Between-group comparisons
   at matched pseudotime (4F, 4J) remain informative; the rise itself should not be
   presented as independent validation.

5. **In-figure titles carry stale cross-references.** 4K is titled "... Fig2 e1
   perturbations ..." and 4L "Fig3 double KOs x top N genes ...", both referring to the
   old notebook numbering; the panel they point to is now 4J. These strings are baked
   into the PDFs and need updating or stripping before submission.

6. **Multiplicity scope.** 4G/4H FDR is computed across the 546 tests within each time
   point, not across both time points; 4K and 4L FDRs are within-column and
   within-perturbation respectively. Worth stating if a reviewer asks.

7. **Regeneration order.** The Pseudotime notebook
   (`Fig4BCFGHJ_S07-S08_S10-S11_S16-S17__Pseudotime.ipynb`) must be run before the
   Interaction_Doxo1 and Perturbation_Gene_Effects notebooks: those two read
   `figures/Fig2_state_enrichment_{day}.csv` to define the perturbation set for 4K
   and 4L. Nothing enforces this in code.
