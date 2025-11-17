import scanpy as sc
import anndata as ad
import pandas as pd

#1.LOAD DATA
adata = sc.read_h5ad("bone_marrow.h5ad")

adata.shape
adata.X
adata.var.head()
adata.obs.head()
adata.to_df()

#2.QUALITY CONTROL

adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
adata.obs.head()

#3.QUALITY CONTROL FILTERING

adata = adata[adata.obs.n_genes_by_counts < 4000, :]
adata = adata[adata.obs.n_genes_by_counts > 300, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]
adata.obs.head()

#4.NORMALIZATION
sc.pp.normalize_total(adata, target_sum=1e4)

#5.LOGARITHMIC TRANSFORMATION
sc.pp.log1p(adata)

#6.IDENTIFY HIGHLY VARIABLE GENES   
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
adata.var.head()

#7.SCALING THE DAT
sc.pp.scale(adata, max_value=10)

#8.PCA
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, log=True)

#9.COMPUTE NEIGHBORHOOD GRAPH
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

#10.UMAP EMBEDDING FOR VISUALIZATION
sc.tl.umap(adata)
sc.pl.umap(adata, color=["n_genes_by_counts", "pct_counts_mt"])

#11.LEIDEN CLUSTERING
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=["leiden"])

#12.FIND MARKER GENES FOR EACH CLUSTER
sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

#13.ANNOTATE CELL TYPES BASED ON MARKER GENES
# For each Leiden cluster, pick the most frequent Cell.group
cluster2celltype = adata.obs.groupby('leiden')['Cell.group'] \
                            .agg(lambda x: x.value_counts().idxmax()) \
                            .to_dict()
print(cluster2celltype)
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster2celltype).astype('category')
sc.pl.umap(adata, color='cell_type', legend_loc='on data')

pd.crosstab(adata.obs['cell_type'], adata.obs['Cell.group'])