# What Single-Cell Data is Teaching Us About Cancer Evolution

---

## Introduction

Clonal evolution (CE) model describes tumor progression through sequential genetic mutations in individual cells, analogous to Darwinian evolution, where diversification and selection based on cellular fitness lead to the emergence of dominant, more malignant subpopulations. The model has advanced in parallel with improvements in sequencing technologies, deepening our ability to observe genetic and transcriptomic heterogeneity in tumours [1].

RNA sequencing methods have progressed rapidly from bulk RNA-seq, to laser-capture micro-dissection, single-cell RNA-seq, and more recently spatial transcriptomics. This evolution of technology has revealed that cancer does not evolve by a single linear clone. Instead, multiple subclones often emerge and coexist, sometimes due to similar fitness or spatial separation, and can even complement each other’s survival. Clonal architecture varies substantially between individual tumours of the same type. In response, computational and experimental methods continue evolving to better model clonal dynamics and inform treatment strategies [1,3].

---

## Transition from Bulk RNA-seq to scRNA-seq

Bulk RNA-seq has been foundational in cancer genomics, helping to classify molecular subtypes, measure gene-expression changes associated with mutations, and discover biomarkers. However, bulk RNA-seq averages the transcriptomes of all cells in a sample, masking the contribution of rare or distinct cell populations.

To overcome this limitation, **single-cell RNA sequencing (scRNA-seq)** was developed. scRNA-seq enables transcriptomic profiling at the resolution of individual cells, which allows detection of rare cell types, intermediate states, and cellular dynamics. RNA sequencing has progressed to single-molecular, single-cell, and spatial transcriptome approaches, improving accuracy of cell-level resolution and contributing to our understanding of cancer heterogeneity and evolution [2,3].

---

## Uncovering Tumour Heterogeneity with scRNA-seq

Cancer exhibits extensive transcriptional heterogeneity, reflecting differences in proliferation, differentiation, immune evasion, and drug response [4,5]. This variability arises from genetic mutations as well as non-genetic mechanisms, including epigenetic modifications, stochastic gene expression, and signals from the tumour microenvironment such as hypoxia, cytokines, or interactions with stromal and immune cells [6,7,11]. For instance, connective tissue cells and infiltrating immune cells can modulate gene expression in cancer cells, promoting invasion and resistance [8].

Single-cell RNA sequencing (scRNA-seq) is essential for dissecting these layers of complexity. By profiling individual cells, scRNA-seq uncovers rare cell populations, intermediate states, and developmental trajectories that are obscured in bulk analyses [9,10,13]. Tumour cells from the same patient may display divergent transcriptional programs, including proliferative, quiescent, epithelial-to-mesenchymal transition (EMT), stress-response, or drug-tolerant states [12,14].

Trajectory analyses, pseudotime reconstruction, and dimensionality-reduction techniques (e.g., UMAP, t-SNE) reveal transitions between states, including EMT and drug-tolerant phenotypes [12,14]. When combined with spatial transcriptomics, scRNA-seq additionally maps how transcriptional programs are organized within the tumour and influenced by local microenvironmental interactions [15]. This approach highlights the contributions of both genetic and non-genetic heterogeneity to tumour progression, plasticity, and therapy response.

---

## Tracking Cancer Evolution with scRNA-seq

scRNA-seq has also enabled the reconstruction of tumour evolutionary dynamics at single-cell resolution [14]. By profiling thousands of cells simultaneously, researchers can trace differentiation hierarchies, identify stem-like and drug-resistant subpopulations, and map transcriptional trajectories associated with metastasis or therapy adaptation. Unlike bulk RNA-seq, which averages signals across millions of cells, scRNA-seq reveals coexisting subclones and dynamic state transitions, providing a more accurate depiction of tumour evolution [12,14].

Genetic inference from scRNA-seq data such as detection of copy-number alterations (CNAs), expressed somatic mutations, and mitochondrial DNA variants allows subclonal reconstruction and phylogenetic mapping [16,17,18]. Computational tools like **inferCNV**, **CopyKAT**, and **SCEVAN** assign cells to likely genomic lineages, while experimental lineage-tracing approaches, including CRISPR-based barcoding, enhance clonal tracking in model systems [16].

Integration with complementary single-cell modalities, including single-cell DNA sequencing, ATAC-seq, and spatial transcriptomics, links genotype to phenotype, connecting transcriptional reprogramming to genetic diversification and microenvironmental adaptation [12,15]. Together, these methods provide a dynamic view of cancer as an evolving ecosystem, allowing identification of drug-resistant subclones before they expand and informing precision therapies designed to anticipate and disrupt tumour evolution [17].

---

## Connecting scRNA-seq to Therapy Resistance

Because genetic and epigenetic differences between tumours of the same type in different patients can lead to drug resistance and therapeutic failure [16], integrating drug response data with scRNA-seq enables analysis at single-cell resolution and the identification of therapeutic clusters. Before **Beyondcell**, there was no computational method capable of inferring drug sensitivity from single-cell gene expression combined with high-throughput drug screening. With the Beyondcell approach, tumour responses to treatment can be monitored across different cellular subpopulations using scRNA-seq [17]. Studies have shown that scRNA-seq contributes to more effective therapies in various cancers, including glioblastoma, which previously yielded less accurate results [18].

---

## Challenges and Future Perspectives

A major limitation of single-cell RNA sequencing (scRNA-seq) in solid tumours is that samples must be dissociated into viable single cells before sequencing. The mechanical and enzymatic steps used for dissociation often trigger stress responses that alter transcription and add technical noise, making it hard to tell whether gene-expression differences come from real tumour heterogeneity or from sample handling [19].

Because fresh material is usually required, collecting multiple or serial biopsies from the same patient is also difficult, which limits studies of tumour evolution over time. Some groups work with organoids or cell lines to avoid this problem, but those models do not fully reproduce interactions with the tumour microenvironment. Even with improvements such as single-nucleus RNA-seq and fixation-compatible workflows, technical noise and sampling bias still restrict the accuracy and clinical use of scRNA-seq in cancer research [20,21]. Optimising fixation and storage protocols could greatly expand the clinical application of scRNA-seq, allowing more flexible sample handling and improving the study of tumour heterogeneity. Widespread implementation of these methods may enhance our ability to track cancer evolution and guide therapeutic strategies.

---

## References

1. Wang, Y., et al. (2020). *Changing technologies of RNA sequencing and their applications in clinical oncology.* Frontiers in Oncology, 10, Article 447. https://doi.org/10.3389/fonc.2020.00447
2. Hong, M., et al. (2020). *RNA sequencing: new technologies and applications in cancer research.* Journal of Hematology & Oncology, 13, Article 166. https://doi.org/10.1186/s13045-020-01005-x
3. Jiang, L., et al. (2016). *GiniClust: detecting rare cell types from single‑cell gene expression data with Gini index.* Genome Biology, 17, Article 144. https://doi.org/10.1186/s13059-016-1010-4 
4. Patel, A. P., et al. (2014). *Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma.* Science, 344(6190), 1396–1401. https://doi.org/10.1126/science.1254257
5. Tirosh, I., et al. (2016). *Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq.* Science, 352(6282), 189–196. https://doi.org/10.1126/science.aad0501
6. Neftel, C., et al. (2019). *An integrative model of cellular states, plasticity, and genetics for glioblastoma.* Cell, 178(4), 835–849.e21. https://doi.org/10.1016/j.cell.2019.06.024
7. Zhang, L., et al. (2021). *Single-cell transcriptomics reveals the landscape of intratumoral heterogeneity and transcriptional plasticity in liver cancer.* Nature Communications, 12, 7279.  
8. Keren, L., et al. (2018). *A structured tumor-immune microenvironment in triple negative breast cancer revealed by multiplexed ion beam imaging.* Cell, 174(6), 1373–1387.e19. https://doi.org/10.1016/j.cell.2018.08.039
9. Stuart, T., & Satija, R. (2019). *Integrative single-cell analysis.* Nature Reviews Genetics, 20(5), 257–272.  
10. Suvà, M. L., & Tirosh, I. (2023). *Cancer cell state plasticity as a driver of tumor evolution and therapy resistance.* Cancer Cell, 41(5), 563–577.
11. Beyes, S., Bediaga, N. G., & Zippo, A. (2021). *An Epigenetic Perspective on Intra-Tumour Heterogeneity: Novel Insights and New Challenges from Multiple Fields.* Cancers, 13(19), 4969. https://doi.org/10.3390/cancers13194969
12. Chang, X., Zheng, Y., & Xu, K. (2024). *Single-Cell RNA Sequencing: Technological Progress and Biomedical Application in Cancer Research.* Molecular biotechnology, 66(7), 1497–1519. https://doi.org/10.1007/s12033-023-00777-0
13. Sharma, A., Merritt, E., Hu, X., Cruz, A., Jiang, C., Sarkodie, H., Zhou, Z., Malhotra, J., Riedlinger, G. M., & De, S. (2019). *Non-Genetic Intra-Tumor Heterogeneity Is a Major Predictor of Phenotypic Heterogeneity and Ongoing Evolutionary Dynamics in Lung Tumors.* Cell Reports, 29(8), 2164-2174.e5. https://doi.org/10.1016/j.celrep.2019.10.045
14. Wu, F., Fan, J., He, Y., Xiong, A., Yu, J., Li, Y., Zhang, Y., Zhao, W., Zhou, F., Li, W., Zhang, J., Zhang, X., Qiao, M., Gao, G., Chen, S., Chen, X., Li, X., Hou, L., Wu, C., & Su, C. (2021). *Single-cell profiling of tumor heterogeneity and the microenvironment in advanced non-small cell lung cancer.* Nature Communications, 12(1), 2540. https://doi.org/10.1038/s41467-021-22801-0
15. Zhang, Y., Yang, C., Chen, X., Wu, L., Yuan, Z., Zhang, F., & Qian, B. (2025). *Cancer therapy resistance from a spatial‐omics perspective.* Clinical and Translational Medicine, 15(7). https://doi.org/10.1002/ctm2.70396
16. Turajlic S, Sottoriva A, Graham T, Swanton C. *Resolving genetic heterogeneity in cancer.* Nat Rev Genet. 2019;20(7):404–16. https://doi.org/10.1038/s41576-019-0114-6.
17. Fustero-Torre, C., Jiménez-Santos, M. J., García-Martín, S., Carretero-Puche, C., García-Jimeno, L., Ivanchuk, V., Di Domenico, T., Gómez-López, G., & Al-Shahrour, F. (2021). *Beyondcell: Targeting cancer therapeutic heterogeneity in single-cell RNA-seq data.* Genome Medicine, 13(1). https://doi.org/10.1186/s13073-021-01001-x.
18. Wu, H., Guo, C., Wang, C., Xu, J., Zheng, S., Duan, J., Li, Y., Bai, H., Xu, Q., Ning, F., Wang, F., & Yang, Q. (2023). *Single‐cell rna sequencing reveals tumor heterogeneity, microenvironment, and drug‐resistance mechanisms of recurrent glioblastoma.* Cancer Science, 114(6), 2609–2621. https://doi.org/10.1111/cas.15773.
19. Tung, P. Y., Blischak, J. D., Hsiao, C. J., Knowles, D. A., Burnett, J. E., Pritchard, J. K., & Gilad, Y. (2017). *Batch effects and the effective design of single-cell gene expression studies.* Scientific reports, 7(1), 39921.
20. Grün, D., Kester, L., & van Oudenaarden, A. (2014). *Validation of noise models for single-cell transcriptomics.* Nature Methods, 11(6), 637–640. https://doi.org/10.1038/nmeth.2930.
21. González-Silva, L., Quevedo, L., & Varela, I. (2021). *Tumor functional heterogeneity unraveled by scRNA-Seq Technologies.* Trends in Cancer, 7(3), 265. https://doi.org/10.1016/j.trecan.2021.02.001.
22. Tirosh, Itay et al. *Cancer cell states: Lessons from ten years of single-cell RNA-sequencing of human tumors.* Cancer Cell, Volume 42, Issue 9, 1497 - 1506.
23. Ludwig L. S. Lareau C. A., Ulirsch J.C.,  Buenrostro J. D., Regev A.,  Sankaran V. G. *Cancer cell states: Lessons from ten years of single-cell RNA-sequencing of human tumors.* Volume 42, Issue 9P1497-1506 September 09, 2024
24. Yang D, Jones MG, Naranjo S, Rideout WM 3rd, Min KHJ, Ho R, Wu W, Replogle JM, Page JL, Quinn JJ, Horns F, Qiu X, Chen MZ, Freed-Pastor WA, McGinnis CS, Patterson DM, Gartner ZJ, Chow ED, Bivona TG, Chan MM, Yosef N, Jacks T, Weissman JS. *Lineage tracing reveals the phylodynamics, plasticity, and paths of tumor evolution.* Cell. 2022 May 26;185(11):1905-1923.e25. doi: 10.1016/j.cell.2022.04.015. Epub 2022 May 5. PMID: 35523183; PMCID: PMC9452598. 


## Contributors

- Özlem Kalkan  
- Oni David Sunday  
- Ayşenur Akcan  
- Daria Kriuchkova  
- Noor Ul Ain Amir  
- Minenhle Mayisela
