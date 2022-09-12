####generate colon epithelial cells using R
###I did use the integrated seurate object and SCT transfromed for normalization
##converting the colon seurat to anndata h5ad file

anndata <- import("anndata", convert = FALSE)
adata1 <- anndata$AnnData(
  X = t(as.matrix(GetAssayData(object = SCT_Seurat_colon,assay = "RNA", slot = "data"))),
  obs = data.frame(SCT_Seurat_colon@meta.data),
  var = GetAssay(SCT_Seurat_colon,assay = "RNA")[[]],
  obsm  = list(
    "X_pca" = Embeddings(SCT_Seurat_colon[["pca"]]),
    "X_umap" = Embeddings(SCT_Seurat_colon[["umap"]])
  )
)
anndata$AnnData$write(adata1, 'colon_1.h5ad')


###I used velocyto to generate the loom files from the cell ranger I did batch script.
import os
import loompy
import igraph
import desc
import scanpy
import scvelo as scv
import anndata

###combinging all the loom files into one files
files =["BLM.loom","BLM2.loom", "Min.loom", "MSH2.loom"]
loompy.combine(files,"merged_files.loom")

###reading the colon data
colon_adata = desc.read_h5ad('colon_1.h5ad')
##reading loom file
loom_file = scv.read('merged_files.loom', cache=True)

###merging loom file with the orginal anndata
adata= scv.utils.merge(colon_adata, loom_file)

###the file should have the embided umap and pca and also should have spliced
##and unspliced RNA count
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

###processing the data

scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

###calculate the RNA velocity

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

###visualizing the adata
###to get the stream
clusters = clusters = 'integrated_snn_res.0.2'
scv.pl.velocity_embedding_stream(adata, basis='umap', color=cluster, figsize=(10,10)
,dpi=300, save=("veloStream.png"),size =50)

###doing the arrow
scv.pl.velocity_embedding(adata, arrow_length=4, arrow_size=4, dpi=300,color=clusters, figsize=(10,10)
,dpi=300, save=("veloArrow.png"),size =50)

adata.write("adata.h5ad")
adata = desc.read_h5ad('adata.h5ad')


###trying to split the integrated colon into the 4 samples
split = 'orig.ident'
scv.pl.velocity_embedding(adata, arrow_length=4, arrow_size=4, dpi=300,color=split , figsize=(10,10)
, save=("veloArrowsplit.png"),size =50)

####test some gene markers
scv.pl.velocity(adata, ['Lgr5',  'Muc2', 'Guca2a', 'Espn'],
    save="markergenes.png",ncols=2, dpi=300, figsize=(10,10))
