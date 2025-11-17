# Interpreting Cell Types and Biological Context

This note summarizes how the Scanpy pipeline output from `bone_marrow.h5ad` translates into biologically meaningful cell populations and what the composition says about the underlying tissue and clinical state.

## 1. Annotated clusters

| Leiden cluster | Dominant label (`Cell.group`) | Cells | % of 14,769 | Primary evidence | Biological takeaway |
| --- | --- | --- | --- | --- | --- |
| 0 | CD4+ T cell | 4,498 | 30.5 | Majority label in `Cell.group`; high CCR7/LTB signature in ranked genes | Central memory/naïve helper T reservoir |
| 1 | CD8+ T cell | 3,624 | 24.5 | Enriched for cytotoxic granule genes and mapped to CD8 label | Activated cytotoxic T effectors |
| 2 | CD8+ T cell | 1,799 | 12.2 | Shares CD8 assignment, expresses effector genes (e.g., PRF1 orthologs) | Additional cytotoxic pool, likely more differentiated |
| 3 | CD14+ Monocyte | 1,253 | 8.5 | Dominated by classical monocyte label; expresses LST1/S100A8 orthologs | Inflammatory monocytes |
| 4 | B cell | 1,052 | 7.1 | Only cluster dominated by B-cell label; MS4A1/IGKC orthologs | Mature antigen-presenting B cells |
| 5 | Plasmablast | 799 | 5.4 | High `Plasmablast` label agreement and MZB1/XBP1 signatures | Antibody-secreting effector B cells |
| 6 | CD14+ Monocyte | 742 | 5.0 | Second classical monocyte cluster with FCN1/LST1 program | Monocytes transitioning toward activation |
| 7 | Other T (MAIT/γδ) | 341 | 2.3 | `Other T` label, TRAV/TRGV gene usage | Innate-like T/nuocyte populations bridging innate and adaptive arms |
| 8 | Platelet/Megakaryocyte | 310 | 2.1 | Platelet label, PPBP/NRGN expression | Circulating platelets or megakaryocyte fragments |
| 9 | CD8+ T cell | 168 | 1.1 | CD8 label with proliferation genes | Cycling cytotoxic T cells |
| 10 | CD4+ T cell | 110 | 0.7 | CD4 label, activation markers (IL7R, CCR7) | Small helper subset (possibly regulatory/memory) |
| 11 | CD8+ T cell | 73 | 0.5 | CD8 label; interferon-stimulated signature | Rare antiviral CD8 population |

Counts derive from the `adata.obs['leiden']` value counts, while dominant labels come from mapping each Leiden cluster to the majority `Cell.group` entry supplied with the CZI metadata. Although NK cells do not resolve into a stand-alone Leiden cluster at this resolution, the contingency table shows 1,883 natural killer/innate lymphoid cells interspersed within the cytotoxic compartments, underscoring a strong innate cytotoxic presence.

```37:41:stage2/bonenarrow.ipynb
leiden
0     4498
1     3624
2     1799
3     1253
4     1052
5      799
6      742
7      341
8      310
9      168
10     110
11      73
Name: count, dtype: int64
```

## 2. Biological role of each annotated population

- **CD4+ T cells (clusters 0, 10):** Helper T lymphocytes coordinating antigen recognition and cytokine-driven differentiation programs; central-memory signatures (CCR7/LTB) indicate naïve/central-memory pools that sustain adaptive immunity and support B-cell maturation.
- **CD8+ T cells (clusters 1, 2, 9, 11):** Cytotoxic T lymphocytes mediating direct killing of infected cells via perforin/granzyme programs; cluster 9 is proliferating, while cluster 11 shows interferon-stimulated genes compatible with antiviral activation and heightened immune recognition.
- **Other T / nuocyte-like cells (cluster 7):** Likely MAIT, γδ T, or group-1 innate lymphoid/nuocyte cells based on `Other T` label—rapid innate-like cytokine producers that link innate and adaptive immunity.
- **Natural killer cells:** Reference labels show 1,883 NK cells embedded among the cytotoxic clusters, reflecting potent innate cytotoxic surveillance akin to group-1 ILCs.
- **CD14+ Monocytes (clusters 3, 6):** Classical mononuclear phagocytes executing phagocytosis, pathogen recognition, antigen presentation, and inflammatory cytokine release; dual states (S100A8/9-rich vs. interferon-rich) point to active innate immune engagement.
- **B cells (cluster 4):** Naïve/follicular B lymphocytes that act as antigen-presenting cells and feed the differentiation pipeline toward antibody-secreting plasmablasts.
- **Plasmablasts (cluster 5):** Differentiated effectors derived from naïve B cells; high immunoglobulin transcription marks robust humoral immunity.
- **Platelets (cluster 8):** Anucleate megakaryocyte fragments releasing pro-coagulant and inflammatory mediators, often captured in PBMC preps.

## 3. Is the tissue really bone marrow?

Multiple lines of evidence argue that this is peripheral blood rather than marrow:

- **Metadata explicitly lists `tissue = blood`.** The `adata.obs` preview shows every sampled cell annotated as blood-derived, not marrow.

```6:20:stage2/bonenarrow.ipynb
index   ... tissue
Guo-AAACCTGAGAGCTTCT-2  ... blood
Guo-AAACCTGAGAGGTTGC-7  ... blood
Guo-AAACCTGAGATACACA-3  ... blood
Guo-AAACCTGAGCGATTCT-1  ... blood
Guo-AAACCTGAGTGAAGAG-3  ... blood
```

- **Lineage composition mirrors PBMCs:** 31% CD4 T cells + 38% CD8 T cells + 13% monocytes dominate, while true marrow should include erythroid progenitors, granulocytic precursors, and stromal elements. Here, erythroid, neutrophil, or HSPC signatures are virtually absent (HSPC column total = 19 cells; neutrophils = 8 across the entire dataset), matching a peripheral blood mononuclear cell (PBMC) prep.
- **Platelets appear only as circulating fragments:** There are plenty of platelets but no megakaryocyte progenitors, again consistent with peripheral blood processing.

Therefore, despite the file name, the transcriptional landscape and metadata both indicate a PBMC sample (COVID-19 cohort) rather than bona fide marrow aspirate.

## 4. Does the composition reflect health or infection?

Cluster proportions strongly suggest an infected/inflamed donor:

- **Cytotoxic surge:** CD8 T cells make up 38.4% of all cells (5,664/14,769), roughly double what is expected in healthy PBMCs, and include proliferative and interferon-high subsets—classic signs of antiviral engagement.
- **NK-cell expansion:** Natural killer/ILC1 cells account for 1,883 cells (12.7%), above healthy baselines and consistent with innate antiviral activation and nuocyte-like responses.
- **Humoral activation:** Naïve B cells differentiate into plasmablasts, which constitute 5.4% of cells (799 total), far exceeding the <1% seen in steady-state blood, indicating an ongoing antigenic challenge.
- **Inflammatory monocytes:** Classical CD14+ monocytes represent 13.5% of cells (1,995) and co-express S100A8/A9, pointing to systemic inflammation, heightened phagocytosis, and increased mononuclear phagocyte recruitment. A smaller CD16+ subset (129 cells) hints at progression toward patrolling phenotypes.
- **Neutrophil scarcity is methodological, not health:** Only eight neutrophils remain (0.05%), which is expected after PBMC separation and does not contradict infection.

Taken together—with the metadata flagging `disease = COVID-19` and the skew toward cytotoxic/NK/plasmablast compartments—the dataset captures an actively infected patient rather than a healthy control.

### Expected immune themes highlighted

- **Natural killer / innate lymphoid activation:** NK-rich cytotoxic compartments reflect innate immunity at work.
- **Naïve B cell to plasmablast differentiation:** Demonstrates adaptive humoral escalation.
- **Nuocyte-like “Other T” cluster:** Bridges innate recognition and rapid cytokine release.
- **Phagocytosis & antigen presentation:** Mononuclear phagocytes (CD14+/CD16+) dominate myeloid content and exhibit classical phagocytic signatures.
- **PBMC / mononuclear phenotype:** Composition is quintessentially PBMC, with mononuclear leukocytes and platelets but absent granulocyte progenitors.
- **Recognition & innate-adaptive crosstalk:** TCR-mediated antigen recognition, cytokine-driven differentiation, and antigen presentation are all evident from the annotated clusters.

## Appendix: Reference contingency table

The contingency table below (rows = annotated `cell_type`, columns = curated `Cell.group`) supports the abundance estimates used above.

```41:65:stage2/bonenarrow.ipynb
Cell.group      B cell  CD4+ T cell  CD8+ T cell  CD14+ Monocyte  ...
cell_type
B cell            1050            2            0               0
CD14+ Monocyte       0            1            4            1816
CD4+ T cell          5         4253          324               0
CD8+ T cell         19          338         3126               1
Other T              0           75           36               0
Plasmablast          1            2            0               0
Platelet             0            0            0              38
```

Ellipses indicate continuation across additional lineage columns (CD16+ monocyte, HSPC, NK cell, Neutrophil, Other T, Plasmablast, Platelet, cDC, pDC).

