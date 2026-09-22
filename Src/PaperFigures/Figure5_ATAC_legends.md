# Figure 5 — figure legends

Draft legends for the chromatin (scATAC) transcription-factor motif-accessibility figure. Statistics
(ARI, ρ, counts) are the values produced by the panel notebooks; update if the analysis is rerun.

## Figure 5

**Figure 5. Transcription-factor motif accessibility across the NEPC differentiation landscape from single-cell ATAC-seq.**

**(A)** **Cell embedding computed from ATAC (chromatin-accessibility) profiles, colored by cell
states obtained from the paired RNA modality.** Cells were embedded by iterative LSI on a
genome-wide 500-bp tile matrix (features independent of the RNA annotation) and projected with
UMAP; each cell is colored by its RNA-defined state. RNA-defined states occupy coherent
territories of the chromatin embedding (adjusted Rand index = 0.68 between unsupervised ATAC
clusters and RNA states), showing that the transcriptional cell states are recapitulated at the
level of chromatin accessibility.

**(B)** The same ATAC-profile embedding colored by RNA diffusion pseudotime (`dpt_pseudotime`).
The black line is the RNA differentiation trajectory drawn in chromatin space — the
running-median ATAC-UMAP position across 80 pseudotime bins — from the neuroendocrine root (dot)
to the differentiated tip (arrow). Pseudotime forms a smooth gradient along the ATAC embedding
(|Spearman ρ| = 0.71 with UMAP position), i.e. the differentiation trajectory is also visible in
chromatin.

**(C)** TF motif accessibility per cell state. For each state, per-peak log2 fold-change (state vs. all other
cells) was regressed on the binary HOCOMOCO motif matrix; the OLS coefficient is the TF motif accessibility
(coefficient > 0, motif-bearing peaks open → more accessible). Rows are TFs significant (FDR ≤ 0.1)
in ≥ 1 state, ordered by the state where activity peaks (neuroendocrine-active at top →
differentiated-active at bottom); columns are states in trajectory order. Asterisks, FDR ≤ 0.1.

**(D)** TF motif accessibility along the differentiation trajectory (line view). Per-peak log2 fold-change
of each pseudotime segment vs. the neuroendocrine root segment was regressed on the motif matrix;
the coefficient is the TF motif accessibility of that segment relative to the root. Trajectories are shown
for the top 10 rising (solid) and top 10 falling (dashed) TFs by Spearman trend across segments
(root NE → differentiated tip). The full heatmap of all significant TFs is Fig S5F.

**(E)** TF motif accessibility vs. RNA expression (examples). Motif-accessibility coefficient
(vs. NE root) against mean mRNA, one point per pseudotime segment (colored NE root →
differentiated tip). Motif accessibility reports opening/closing of a TF's target sites, not the
TF's own transcription, so the sign of the correlation reflects the TF's regulatory role.
**GRHL2** (epithelial pioneer/activator) is concordant — sites and mRNA rise together (ρ = +0.89).
The discordant TFs are transcriptional repressors whose target sites open as the factor is lost:
**ZEB1** (EMT repressor; recruits CtBP/HDAC) ρ = −0.84; **NFIC** (NFI-family, context-dependent
repressor) ρ = −0.70; **ZBTB14** (BTB/POZ zinc-finger repressor) ρ = −0.44. Full set in Fig S5H.

## Supplementary Figure S5

**scATAC quality control (A–C).**

**(A)** Cell calling. TSS enrichment vs. log10 unique fragments over all 8.17M barcodes (2-D
density); dashed lines mark the retained-set lower bounds. High-quality nuclei (upper right)
separate from empty/ambient barcodes (lower left); 280,177 cells pass.

**(B)** Per-cell QC distributions of the retained cells: TSS enrichment, log10 unique fragments,
and FRIP (fraction of reads in peaks) — medians 11.8, 6,734, and 0.29, respectively.

**(C)** Aggregate TSS-enrichment profile (± 2 kb around transcription start sites; fold-enrichment
over flanking background). The sharp central peak and +1-nucleosome shoulder indicate
high-quality, promoter-focused signal. *(Fragment-size/nucleosome QC is omitted: the width-1
insertion-site input makes those ArchR metrics degenerate.)*

**Analysis supplements (D–H).**

**(D)** Pseudotime segments across modalities. The 16 adaptive even-trajectory pseudotime
segments (as used for peak calling) shown on the ATAC-profile UMAP (left) and the RNA UMAP
(right); segments are defined once from pseudotime and colored identically, so matching layouts
indicate the trajectory bins occupy concordant territories in chromatin and transcriptome.

**(E)** TF motif accessibility across pseudotime — full heatmap. All TFs significant (FDR ≤ 0.1) in ≥ 1
segment, ordered by Spearman trend (rising top → falling bottom); columns root NE →
differentiated tip. Asterisks, FDR ≤ 0.1. (Main Fig 5D shows the top rising/falling TFs as line
traces.)

**(F)** TF motif accessibility per genetic perturbation (pooled). Each perturbation vs. pooled
non-targeting-control cells (no state stratification), pseudotime peak set; per-peak log2
fold-change regressed on the motif matrix. Non-zero FDR-significant (≤ 0.1) coefficients (TFs ×
perturbations), rows/columns hierarchically clustered, cell values = coefficients. Few
perturbations reach significance at this marginal, chromatin-motif readout.

**(G)** TF motif accessibility per cell state × perturbation. Each perturbation vs. same-state NTC; per-peak
log2 fold-change regressed on the motif matrix. Non-zero FDR-significant (≤ 0.1) coefficients
(all conditions with ≥ 1 significant TF): columns grouped by cell state (trajectory order,
annotated by the colored strip/labels above; perturbation on the x-axis), TF rows hierarchically
clustered.

**(H)** TF motif accessibility vs. RNA expression across pseudotime — remaining individual TFs (the four
examples in main Fig 5E, and motif groups, are excluded; SNAI2/IKZF1/LYL1 are absent from the RNA
reference). TF motif accessibility coefficient (vs. NE root) vs. mean mRNA, one point per pseudotime
segment (colored NE root → differentiated tip); Spearman ρ per TF. Concordance is TF-specific —
some TFs track their expression while others (repressors, or factors whose motif is bound by other
family members) move oppositely.
