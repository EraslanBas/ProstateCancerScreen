# Figure 5 — panel legends

Panels assigned to main Figure 5. Panel letters come from `PAPER_PANELS`
(`_figutils.py`), which is the single place a published name is edited; the notebook
filenames and this document follow from it.

Shared methods (cell counts, normalisation, state definitions, pseudotime construction) are
in `Figure4_legends.md`.

---

### Figure 5F. NEUROG1+SIM1 effects on senescence genes, day10
`Fig5F_senescence_interaction_NEUROG1_SIM1_day10.pdf`
Fig4L_Fig5F_S15__Perturbation_Gene_Effects.ipynb

For each senescence/SASP gene displayed in Figure 4L, the effect of NEUROG1+SIM1 on that
gene's expression decomposed as in Figure 4K: individual effect of NEUROG1, individual
effect of SIM1, their interaction, and the total double-knockout effect, from
`expression ~ g1 + g2 + g1:g2` plus a separate total model, fit on all day10 cells.
Rows are ordered by total effect; annotations give the coefficient with stars at
FDR < 0.1, < 0.01 and < 0.001.

---

## Supplementary figures related to Figure 5

- **S16** — Doxo1 score along pseudotime for Wilcoxon-significant perturbations
  (`FigS16_Comp_doxo1_vs_pseudotime_wilcoxon_hits.pdf`)
- **S17** — NEUROG1 x SIM1 interaction on the Doxo1 score in differentiated cells
  (`FigS17_Comp_NEUROG1_SIM1_interaction.pdf`)

Legends for both are in `Supplementary_figures.md`.

---

## Author notes — NOT for publication

1. **This panel currently covers one combination.** The notebook hard-codes
   `combos = ['NEUROG1+SIM1']` in the senescence cell; the commented-out line would
   produce one panel per day10 double knockout. If Figure 5 needs more than
   NEUROG1+SIM1, uncomment it and add the new filenames to `PAPER_PANELS` (unmapped
   names are written unchanged).

2. **Panel letter assigned.** This panel is Figure 5F, written as
   `Fig5F_senescence_interaction_NEUROG1_SIM1_day10.pdf`. Panels 5A to 5E are not
   produced by these notebooks.
