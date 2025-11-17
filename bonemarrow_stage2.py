import scanpy as sc
import anndata as ad
import pandas as pd

#1.LOAD DATA
"""
Load the single-cell RNA-seq dataset in AnnData format.

This step reads a preprocessed `.h5ad` file containing:
- adata.X: gene expression matrix (cells x genes)
- adata.obs: cell metadata
- adata.var: gene metadata

After loading, basic properties of the dataset are displayed to verify:
- Number of cells and genes
- First few genes and cell metadata
- Conversion to DataFrame (optional) for quick inspection
"""
adata = sc.read_h5ad("bone_marrow.h5ad")

adata.shape
adata.X
adata.var.head()
adata.obs.head()
adata.to_df()

#2.QUALITY CONTROL
"""
Perform initial quality control (QC) on the single-cell dataset.

Steps included:
1. Ensure unique gene and cell names to avoid conflicts during analysis.
2. Identify mitochondrial genes by checking if gene names start with 'MT-'.
3. Compute QC metrics for each cell, including:
   - Total counts
   - Number of genes detected
   - Percentage of counts from mitochondrial genes

High mitochondrial content can indicate stressed or dying cells, which may be removed in filtering steps.
"""
adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
adata.obs.head()

#3.QUALITY CONTROL FILTERING
"""
Filter out low-quality or potentially problematic cells based on QC metrics.

Filtering criteria:
1. Remove cells with fewer than 300 genes detected – likely empty droplets or dead cells.
2. Remove cells with more than 4000 genes detected – potential doublets or multiplets.
3. Remove cells with >10% of counts from mitochondrial genes – indicates stressed or dying cells.

After filtering, the dataset contains higher-quality cells suitable for downstream analysis.
"""

adata = adata[adata.obs.n_genes_by_counts < 4000, :]
adata = adata[adata.obs.n_genes_by_counts > 300, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]
adata.obs.head()

#4.NORMALIZATION
"""
Normalize gene expression counts for each cell to account for differences
in sequencing depth.

- `target_sum=1e4` scales each cell's total counts to 10,000.
- This ensures comparability across cells.
"""
sc.pp.normalize_total(adata, target_sum=1e4)

#5.LOGARITHMIC TRANSFORMATION
"""
Apply log1p transformation (log(x + 1)) to normalized counts.

- Reduces the effect of extreme values and makes the data more normally distributed.
- Facilitates downstream analyses like PCA and clustering.
"""
sc.pp.log1p(adata)

#6.IDENTIFY HIGHLY VARIABLE GENES   
"""
Select genes with high variability across cells for downstream analysis.

- Parameters:
  - min_mean=0.0125: minimum mean expression threshold
  - max_mean=3: maximum mean expression threshold
  - min_disp=0.5: minimum dispersion threshold
- Only highly variable genes are retained in the dataset for PCA and clustering,
  reducing noise from lowly expressed or uninformative genes.
"""
sc.pp
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
adata.var.head()

#7.SCALING THE DATa
"""
Scale each gene to unit variance and zero mean.

- `max_value=10` caps extreme values to reduce the influence of outliers.
- Standardization ensures that all genes contribute equally to downstream analyses
  like PCA and clustering.
"""
sc.pp.scale(adata, max_value=10)

#8.PCA
"""
Perform PCA on the scaled, highly variable genes to reduce dimensionality.

- `svd_solver="arpack"`: efficient solver for PCA
- Captures major axes of variation in the data
- The resulting principal components are used for neighborhood graph construction
- Visualize variance explained by each PC to decide how many PCs to use for downstream analyses
"""
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, log=True)

#9.COMPUTE NEIGHBORHOOD GRAPH
"""
Construct a K-nearest neighbors graph based on PCA-reduced dimensions.

- `n_neighbors=10`: number of nearest neighbors to consider for each cell
- `n_pcs=40`: number of principal components used
- This graph is the foundation for clustering and UMAP visualization
"""
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#10.UMAP EMBEDDING FOR VISUALIZATION
"""
Compute a 2D UMAP embedding of the cells for visualization.

- Preserves both local and global structure of the data
- Colors can highlight QC metrics or cluster assignments
"""
sc.tl.umap(adata)
sc.pl.umap(adata, color=["n_genes_by_counts", "pct_counts_mt"])

#11.LEIDEN CLUSTERING
"""
Cluster cells using the Leiden algorithm on the neighborhood graph.

- `resolution=0.5`: controls cluster granularity (higher → more clusters)
- Clustering identifies transcriptionally distinct groups of cells
- Visualize clusters on the UMAP embedding
"""
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=["leiden"])

#12.FIND MARKER GENES FOR EACH CLUSTER
"""
Identify genes that are differentially expressed between clusters.

- `groupby="leiden"`: tests each Leiden cluster against all others
- `method="wilcoxon"`: non-parametric statistical test for marker detection
- Plot top marker genes to help interpret cluster identities
"""
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

#13.ANNOTATE CELL TYPES BASED ON MARKER GENES
"""
Map clusters to known cell types using the most frequent label in 'Cell.group'.

Steps:
1. For each Leiden cluster, determine the most frequent Cell.group label
2. Map cluster IDs to these cell type labels
3. Store results in adata.obs['cell_type']
4. Visualize annotated cell types on the UMAP
5. Optional: cross-tabulate to compare cluster assignments with original labels
"""
clus
# For each Leiden cluster, pick the most frequent Cell.group
cluster2celltype = adata.obs.groupby('leiden')['Cell.group'] \
                            .agg(lambda x: x.value_counts().idxmax()) \
                            .to_dict()
print(cluster2celltype)
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster2celltype).astype('category')
sc.pl.umap(adata, color='cell_type', legend_loc='on data')

pd.crosstab(adata.obs['cell_type'], adata.obs['Cell.group'])
