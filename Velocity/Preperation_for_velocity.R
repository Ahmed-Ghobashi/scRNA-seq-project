library("Seurat")
library(SeuratWrappers)
library(reticulate)
library(anndata)
library(loomR)
library(hdf5r)
library("tidyverse")
library(SeuratData)
library(SeuratDisk)

###this colon section is done by substracting from the original data,
##without reclustring the cells
SCT_Seurat_colon<-readRDS("SCT_Seurat_colon.RDS")

DefaultAssay(SCT_Seurat_colon)<-"RNA"

###sperate the integrated surate object into the 4 samples
##BLM,BLM2,Min,MSH2
colon_samples <- setNames(unique(SCT_Seurat_colon$orig.ident),
                          unique(SCT_Seurat_colon$orig.ident))
SCT_Seurat_colon1 <- map(colon_samples, function(x) {
  x <- subset(SCT_Seurat_colon, subset = orig.ident == x)
  return(x)
})

saveRDS(SCT_Seurat_colon1$BLM2 ,"Seurat_BLM2.RDS")
saveRDS(SCT_Seurat_colon1$BLM ,"Seurat_BLM.RDS")
saveRDS(SCT_Seurat_colon1$Min ,"Seurat_Min.RDS")
saveRDS(SCT_Seurat_colon1$MSH2 ,"Seurat_MSH2.RDS")
###Loop to save the 4 seurat objects into RDS
for (ii in names(SCT_Seurat_colon1)){
  filename <- paste(ii, ".RDS", sep="")
  saveRDS(SCT_Seurat_colon1[[ii]], filename)
}
###convert them into h5ad

anndata <- import("anndata", convert = FALSE)
adata_Min <- anndata$AnnData(
  X = t(as.matrix(GetAssayData(object = SCT_Seurat_colon1$Min,assay = "RNA", slot = "data"))),
  obs = data.frame(SCT_Seurat_colon1$Min@meta.data),
  var = GetAssay(SCT_Seurat_colon1$Min,assay = "RNA")[[]],
  obsm  = list(
    "X_pca" = Embeddings(SCT_Seurat_colon1$Min[["pca"]]),
    "X_umap" = Embeddings(SCT_Seurat_colon1$Min[["umap"]])
  )
)
anndata$AnnData$write(adata_Min, 'Min.h5ad')
##
adata_BLM <- anndata$AnnData(
  X = t(as.matrix(GetAssayData(object = SCT_Seurat_colon1$BLM,assay = "RNA", slot = "data"))),
  obs = data.frame(SCT_Seurat_colon1$BLM@meta.data),
  var = GetAssay(SCT_Seurat_colon1$BLM,assay = "RNA")[[]],
  obsm  = list(
    "X_pca" = Embeddings(SCT_Seurat_colon1$BLM[["pca"]]),
    "X_umap" = Embeddings(SCT_Seurat_colon1$BLM[["umap"]])
  )
)
anndata$AnnData$write(adata_BLM, 'BLM.h5ad')
##
adata_BLM2 <- anndata$AnnData(
  X = t(as.matrix(GetAssayData(object = SCT_Seurat_colon1$BLM2,assay = "RNA", slot = "data"))),
  obs = data.frame(SCT_Seurat_colon1$BLM2@meta.data),
  var = GetAssay(SCT_Seurat_colon1$BLM2,assay = "RNA")[[]],
  obsm  = list(
    "X_pca" = Embeddings(SCT_Seurat_colon1$BLM2[["pca"]]),
    "X_umap" = Embeddings(SCT_Seurat_colon1$BLM2[["umap"]])
  )
)
anndata$AnnData$write(adata_BLM2, 'BLM2.h5ad')
##

adata_MSH2 <- anndata$AnnData(
  X = t(as.matrix(GetAssayData(object = SCT_Seurat_colon1$MSH2,assay = "RNA", slot = "data"))),
  obs = data.frame(SCT_Seurat_colon1$MSH2@meta.data),
  var = GetAssay(SCT_Seurat_colon1$MSH2,assay = "RNA")[[]],
  obsm  = list(
    "X_pca" = Embeddings(SCT_Seurat_colon1$MSH2[["pca"]]),
    "X_umap" = Embeddings(SCT_Seurat_colon1$MSH2[["umap"]])
  )
)
anndata$AnnData$write(adata_MSH2, 'MSH2.h5ad')

