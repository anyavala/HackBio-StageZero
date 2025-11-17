---
marp: true
paginate: true
backgroundColor: #0d0d0d
color: #f5f5f5
---

# Bone Marrow (PBMC) Snapshot  
### Scanpy UMAP & Immune Landscape

14,769 single cells • COVID-19 patient PBMC • Leiden clustering @ resolution 0.5

---

## Slide 1 · UMAP Overview

- 2D UMAP derived from 40 PCs, k=10 neighbors  
- Colors = Leiden clusters mapped to immune phenotypes  
- Clear separation between helper T, cytotoxic T, monocyte, B/plasmablast, and platelet islands  
- NK/ILC populations intermingle with cytotoxic T territory, indicating shared effector programs

---

<!-- _color: #000000 -->
<!-- _backgroundColor: #f5f5f5 -->

## Slide 2 · Major Cell Types

| Cell Type | % of cells | Role |
| --- | --- | --- |
| CD8⁺ T (clusters 1/2/9/11) | 38% | Cytotoxic clearance of infected cells |
| CD4⁺ T (clusters 0/10) | 31% | Helper coordination & differentiation cues |
| CD14⁺ Monocytes (3/6) | 14% | Phagocytosis, antigen presentation, cytokines |
| B cells (4) | 7% | Antigen recognition & presentation |
| Plasmablasts (5) | 5% | Antibody secretion |
| NK / Other T / Platelets | 5% | Innate cytotoxicity & hemostasis |

---

## Slide 3 · Biological Interpretation

- **Innate alarm:** Elevated CD14⁺ monocytes with S100A8/A9 signature → inflammatory recruitment & phagocytosis  
- **Adaptive surge:** Naïve B cells differentiating into plasmablasts; CD4⁺ helpers retain CCR7/LTB memory phenotype  
- **Cytotoxic dominance:** Expanded CD8⁺ and NK/ILC1 cells with interferon-stimulated genes → antiviral state  
- **PBMC context:** Absence of erythroid/granulocyte progenitors confirms peripheral blood mononuclear prep

---

## Slide 4 · Key Insight

> **Nuocyte / NK collaboration:** Reference labels reveal 1,883 NK/innate lymphoid cells interlaced within CD8⁺ clusters, illustrating how innate cytotoxic cells and T cells co-occupy transcriptional space during severe COVID-19.  
> This blending hints at convergent effector programs—perforin/granzyme, IFN-response, and tissue-homing cues—that could be exploited when tracking disease severity or therapeutic response.

