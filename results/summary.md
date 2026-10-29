# Results summary

Cells: 24673 loaded, 24463 after QC (188 flagged as doublets). Subsample: None.

## Integration check (kNN mixing, 1 = perfectly mixed)
| key | mixing_before | integrated | mixing_after |
|---|---|---|---|
| patient | 0.73 | True | 0.785 |
| condition | 0.159 | True | 0.792 |

## Annotation agreement with the authors' labels
Exact label agreement 90.7%, adjusted Rand index 0.853.

## Interferon response (ISG score, stim - ctrl, mean over patients)
| cell_type | mean_delta | sd_delta | n_patients |
|---|---|---|---|
| Dendritic cells | 2.175 | 0.119 | 8 |
| CD14+ Monocytes | 2.037 | 0.129 | 8 |
| FCGR3A+ Monocytes | 1.833 | 0.137 | 8 |
| Plasmacytoid DCs | 1.659 | 0.21 | 8 |
| B cells | 1.473 | 0.064 | 8 |
| NK cells | 1.465 | 0.09 | 8 |
| CD4 T cells | 1.2 | 0.057 | 8 |
| CD8 T cells | 1.198 | 0.106 | 8 |
| Megakaryocytes | 1.134 | 0.175 | 8 |
| Erythrocytes | 0.821 | 0.398 | 8 |

## Pseudobulk DE (PyDESeq2, ~ patient + condition)
Significant = padj < 0.05 and |log2FC| >= 1.0.

| cell_type | n_patients | tested | n_cells | genes_tested | n_up | n_down | n_de |
|---|---|---|---|---|---|---|---|
| CD14+ Monocytes | 8 | True | 5542.0 | 9915.0 | 1169.0 | 911.0 | 2080.0 |
| FCGR3A+ Monocytes | 8 | True | 1042.0 | 7291.0 | 511.0 | 323.0 | 834.0 |
| Dendritic cells | 7 | True | 434.0 | 6335.0 | 469.0 | 244.0 | 713.0 |
| B cells | 8 | True | 2610.0 | 8820.0 | 390.0 | 151.0 | 541.0 |
| CD4 T cells | 8 | True | 9997.0 | 10908.0 | 379.0 | 118.0 | 497.0 |
| NK cells | 8 | True | 2276.0 | 8253.0 | 270.0 | 79.0 | 349.0 |
| CD8 T cells | 8 | True | 1844.0 | 7665.0 | 193.0 | 26.0 | 219.0 |
| Megakaryocytes | 6 | True | 374.0 | 4279.0 | 77.0 | 8.0 | 85.0 |

**Shared** genes (DE in >= 75% of tested cell types): 166. Top: ISG20, TNFSF10, PSMB9, SAT1, ISG15, GBP1, MX1, MT2A, DRAP1, IRF7, IFI16, NMI, IFI35, UBE2L6, DYNLT1, HERC5, IFI6, IFITM2, OASL, SP110

**Cell-type specific** genes (DE in exactly one cell type): 1690.

| cell_types | example specific genes |
|---|---|
| B cells | COX5A, LINC00926, GSTP1, HLA-DMA, EAF2, C16orf74, LSP1, LACTB |
| CD14+ Monocytes | H3F3B, STX11, C19orf59, HLA-A, COTL1, ANXA4, TRAM1, ZFP36 |
| CD4 T cells | PSMA3, NELL2, RGS1, UBQLN2, CD27, BLMH, KLF5, PTPN1 |
| CD8 T cells | NEAT1, GZMM, HM13, TPM4, PLEKHF1, GPR108, LMAN1, TFG |
| Dendritic cells | IL1R2, NACA, YBX1, NPM1, HSPH1, LGALS3, EEF2, EEF1B2 |
| FCGR3A+ Monocytes | PFN1, CD44, CYP1B1, FCGR3B, CD300E, HNRNPM, HMGN1, DDX5 |
| Megakaryocytes | FTH1, FTL, CFP, LGALS1 |
| NK cells | HSPA5, DNAJB9, RAMP1, HERPUD1, SMAGP, TAGLN2, THEMIS2, IFNG |
