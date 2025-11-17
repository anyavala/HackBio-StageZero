# Bone Marrow scRNA-seq Analysis Pipeline

This repository contains scripts and notebooks for processing, analyzing, and visualizing single-cell RNA sequencing (scRNA-seq) data from bone marrow samples using **Scanpy** and **AnnData**.

## 📌 Overview

The analysis follows a standard scRNA-seq preprocessing and clustering workflow, including:

* Loading the dataset
* Quality control and filtering
* Normalization
* Log-transformation
* Identifying highly variable genes
* Scaling
* Dimensionality reduction (PCA)
* Neighborhood graph construction
* UMAP visualization
* Leiden clustering
* Marker gene detection
* Cell type annotation

All steps are implemented in **`bone_narrow.ipynb`** and **`bonemarrow_stage2.py`**.

---

## 📂 Files Included

| File                     | Description                                                |
| ------------------------ | ---------------------------------------------------------- |
| **bone_narrow.ipynb**    | Jupyter notebook containing the full exploratory analysis. |
| **bonemarrow.md**        | Notes and documentation describing the analysis.           |
| **bonemarrow_stage2.py** | Python script version of the full Scanpy workflow.         |

---

## 🧬 scRNA-seq Analysis Steps

Below is a summary of the workflow implemented in the code.

### **1. Load the Data**

```python
adata = sc.read_h5ad("bone_marrow.h5ad")
```

### **2. Quality Control**

* Make variable names unique
* Calculate mitochondrial gene percentage

### **3. Filtering**

Cells are filtered based on:

* number of genes detected
* mitochondrial percentage

### **4. Normalization**

Total-count normalization to 10,000 counts per cell.

### **5. Log Transformation**

```python
sc.pp.log1p(adata)
```

### **6. Highly Variable Genes**

Used to select informative genes before PCA.

### **7. Data Scaling**

```python
sc.pp.scale(adata, max_value=10)
```

### **8. PCA**

Dimensionality reduction using 40 principal components.

### **9. Neighborhood Graph**

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
```

### **10. UMAP Visualization**

Used for low-dimensional embedding.

### **11. Leiden Clustering**

Community detection to define clusters.

### **12. Marker Gene Identification**

```python
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
```

### **13. Cell Type Annotation**

Clusters are annotated using majority-vote mapping from `Cell.group`:

```python
cluster2celltype = adata.obs.groupby('leiden')["Cell.group"] \
    .agg(lambda x: x.value_counts().idxmax()).to_dict()
```

---

## 📊 Outputs

* PCA variance plots
* UMAP plots colored by QC metrics, clusters, and cell types
* Marker gene visualizations
* Contingency table of predicted vs. known cell labels

---

## 🛠 Requirements

* Python 3.10+
* Scanpy
* AnnData
* Pandas

Install dependencies:

```bash
pip install scanpy anndata pandas
```

---

## 📜 License

This project is for educational and research purposes.

---

## 👩‍💻 Author

Created by **anyavala** as part of the HackBio Stage Zero project.
