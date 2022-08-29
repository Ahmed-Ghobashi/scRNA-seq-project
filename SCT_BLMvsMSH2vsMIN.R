#####all of this script is performed on the supercomputer 
###/N/slate/aghobash/R-project

library("Seurat")
library(SeuratWrappers)
library(sctransform)
library(glmGamPoi)
library(htmltools)

library(tidyverse)
library(DoubletFinder)
library(ROCR)
library(SCENT)
library(clustree)
library(leiden)
library(clustree)
library(leiden)
library(htmltools)

library(SeuratData)
library(SeuratDisk)
library(sceasy)
library(reticulate)
library(velocyto.R)
library(monocle3)

BLM.Object1<-readRDS("BLM.Object1.RDS")
Min.Object1<-readRDS("Min.Object1.RDS")
BLM2.Object1<-readRDS("BLM2.Object1.RDS")
MSH2.Object1<-readRDS("MSH2.Object1.RDS")

#****filtering
#We filter cells that have unique feature counts over 10000 or less than 200
#We filter cells that have >20% mitochondrial counts

BLM.Object <- subset(BLM.Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
                       percent.mt <= 20)

Min.Object <- subset(Min.Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
                       percent.mt <= 20)

BLM2.Object <- subset(BLM2.Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
                        percent.mt <= 20)

MSH2.Object <- subset(MSH2.Object1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 &
                        percent.mt <= 20)

###use Using sctransform as  a normalizating method
BLM.Object<-readRDS("BLM.Object.RDS")
Min.Object<-readRDS("Min.Object.RDS")
BLM2.Object<-readRDS("BLM2.Object.RDS")
MSH2.Object<-readRDS("MSH2.Object.RDS")

###adding a unique prefix to avoid duplicates
Min.Object<-RenameCells(Min.Object, add.cell.id="Min")
head(Min.Object)

BLM.Object<-RenameCells(BLM.Object, add.cell.id="BLM")
head(BLM.Object)

BLM2.Object<-RenameCells(BLM2.Object, add.cell.id="BLM2")
head(BLM2.Object)

MSH2.Object<-RenameCells(MSH2.Object, add.cell.id="MSH2")
head(MSH2.Object)




##combing Min, BLM, BLM2 and MSH2 together

Min_BLM_MSH2 <- list(Min.Object,BLM.Object,BLM2.Object,MSH2.Object)
names(Min_BLM_MSH2)<-c("Min","BLM","BLM2","MSH2")

Min_BLM_MSH2 <- map(Min_BLM_MSH2, SCTransform)# do SCT normalization
saveRDS(Min_BLM_MSH2,"SCT_Min_BLM_MSH2.RDS")

##select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
SCT_Min_BLM_MSH2<-readRDS("SCT_Min_BLM_MSH2.RDS")

#### Prepare for integration. using supercomputer

features <- SelectIntegrationFeatures(object.list = SCT_Min_BLM_MSH2, nfeatures = 3000)
SCT_Min_BLM_MSH2 <- PrepSCTIntegration(object.list = SCT_Min_BLM_MSH2, 
                                       anchor.features = features)

saveRDS(SCT_Min_BLM_MSH2,"SCT_Min_BLM_MSH2.RDS")

## perform integration doing it on Cluster computer

SCT_Min_BLM_MSH2<-readRDS("SCT_Min_BLM_MSH2.RDS")

reference_dataset <- which(names(SCT_Min_BLM_MSH2) %in% c("Min"))

SCT_Min_BLM_MSH2.anchors <- FindIntegrationAnchors(object.list = SCT_Min_BLM_MSH2, 
                                          normalization.method = "SCT",
                        anchor.features = features,
                        reference = reference_dataset)

saveRDS(SCT_Min_BLM_MSH2.anchors,"SCT_Min_BLM_MSH2.anchors.RDS")


SCT_Min_BLM_MSH2.anchors<-readRDS("SCT_Min_BLM_MSH2.anchors.RDS")
SCT_Min_BLM_MSH2.combined <- IntegrateData(anchorset = SCT_Min_BLM_MSH2.anchors, 
                                      normalization.method = "SCT", dims = 1:30)


SCT_Min_BLM_MSH2.combined <- RunPCA(Min_BLM_MSH2.combined, npcs = 100)

saveRDS(SCT_Min_BLM_MSH2.combined,"SCT_Min_BLM_MSH2.combined.RDS")
###doing the cluster computer
##doing the elbowplot

SCT_Min_BLM_MSH2.combined<-readRDS(file="SCT_Min_BLM_MSH2.combined.RDS")

p <- ElbowPlot(SCT_Min_BLM_MSH2.combined, ndims = 100)#the variation is high till almost 75

pdf("ElbowPlot.pdf")
ElbowPlot(SCT_Min_BLM_MSH2.combined, ndims = 100)
dev.off()

## UMAP dimension reduction for visualization.

SCT_Min_BLM_MSH2.combined <- RunUMAP(SCT_Min_BLM_MSH2.combined, reduction = "pca", dims = 1:30)





## Clustering the data.###doing on cluster computer on the cluster computer

SCT_Min_BLM_MSH2.combined  <- FindNeighbors(SCT_Min_BLM_MSH2.combined , dims = 1:75)

#Determine the clusters for various resolutions 
SCT_Min_BLM_MSH2.combined  <- FindClusters(SCT_Min_BLM_MSH2.combined ,
                                     resolution = c(0.2,0.4, 0.6, 0.8, 1.0)),
                                     algorithm = 4, 
                                     method = "igraph"
                                     )




## Plotting a cluster tree later

## Plotting a cluster tree.

p2 <- clustree(SCT_Min_BLM_MSH2.combined, prefix = "integrated_snn_res.")

pdf(file="sct_cluster-tree_4conditions.pdf", height= 12,width=9)
p2
dev.off()

# Assign identity of clusters 0.5
Idents(object = SCT_Min_BLM_MSH2.combined) <- "integrated_snn_res.0.4"

Idents(object = Min_BLM_MSH2.combined) <- "integrated_snn_res.0.4"



##visualize the data

DimPlot(SCT_Min_BLM_MSH2.combined, split.by = "orig.ident")


pdf("SCT_cluster_4condtitons.pdf",height=15, width = 20)
DimPlot(SCT_Min_BLM_MSH2.combined, split.by = "orig.ident")
dev.off()

saveRDS(SCT_Min_BLM_MSH2.combined ,"SCT_Min_BLM_MSH2.combined.RDS")

#featureplot
#use marker to annotate each cluster manually
##using cluster computer
DefaultAssay(SCT_Min_BLM_MSH2.combined)<-"SCT"

goblet<-c("Tff3","Reg4","Muc2","Spdef","Atoh1","muc13","Galnt12","Fcgbp","Muc4",
          "Rep15","Guca2a")

stem<-c("Lgr5","Ascl2","Smoc2",
        "Axin2","Ephb2","Kcne3","Hmgn2","Id1")

EEC<-c("Chga","Dll1","Pax4","Neurog3","Insm1","Neurod1","Cpe",
       "Isl1","Pax6")

Tuft_brush<-c("Cdhr2","Espn")

mast_cells<-c("Mcpt1", "Mcpt2")

MHC<-c("H2-Aa","H2-Eb1","H2-Ab1")
monocyte<-c("Cd14","Cd52")

dentric_cell<-c("Itgax","Itgam","Cx3cr1","Cd83","Cd40",
                 "H2-Aa","H2-Eb1","H2-Ab1","Clec9a","Cd52")
NK_cells<-c("Itga2","Nkg7","Ncr1","Nkrp1c","Gzma")
T_cells<-c("Cd3e","Cd3d","Cd3g","Cd2","Satb1","Cd27","Cd28","Cd5","Cd69","Trac","Cerk","Icos",
           "Cd52","Cd160","Cd8b1","Lag3","Traf1","Sept1","Bcl2","Trbc2","Tbrc1","Ltb","Il7r",
           "Itk")
fibroblast<-c("Vim","Vtn","Col6a2","Pdgfrb")

macrophage<-c("Gpr18","Fpr2","Egr2","Mgl2","Fcgr1","Fgl2","Tgfbr1","Cd68","Itgam","C1qc",
              "Arg1","Chil3","Folr2","Adgre1","Ms4a7")
B_cells<-c("Ms4a1","Ighd","Cd22").

#doing Reg3b, Reg3a genes

P3<-FeaturePlot(SCT_Min_BLM_MSH2.combined, features = B_cells,
                pt.size = 0.3, split.by = "orig.ident"
)& theme(legend.position = c(0.0,0.2))



pdf(file=" SCT_B_cells.pdf", height= 50,width=15)
P3
dev.off()

, height= 110,width=15

##removing doublets on cluster computer

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------

###assume PK .09, Will submit batch script for it

SCT_Min_BLM_MSH2.combined<-readRDS(file="SCT_Min_BLM_MSH2.combined.RDS")

##removing doublets

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------

###assume PK .09, Will submit batch script for it

nExp <- round(ncol(SCT_Min_BLM_MSH2.combined) * 0.075)  # expect 7.5% doublets

SCT_Min_BLM_MSH2.combined_doublet <- doubletFinder_v3(SCT_Min_BLM_MSH2.combined , 
                                                      pN = 0.25, pK = 0.09, nExp = nExp,
                                                      sct = TRUE,PCs = 1:50)

# name of the DF prediction can change, so extract the correct column name.
SCT_Min_BLM_MSH2.combined_doublet.name<-colnames(SCT_Min_BLM_MSH2.combined_doublet@meta.data)[grepl("DF.classification",
                                                                                                    colnames(SCT_Min_BLM_MSH2.combined_doublet@meta.data))]

#remove the predicted doublet from my data
SCT_Seurat_integrated<- SCT_Min_BLM_MSH2.combined_doublet[, SCT_Min_BLM_MSH2.combined_doublet@meta.data[, 
                                                                                                        SCT_Min_BLM_MSH2.combined_doublet.name] == "Singlet"]

saveRDS(SCT_Seurat_integrated ,"SCT_Seurat_integrated.RDS")

###visulaization

p6<-cowplot::plot_grid(ncol = 2, DimPlot(SCT_Min_BLM_MSH2.combined_doublet, group.by = "orig.ident") + NoAxes(),
                       DimPlot(SCT_Min_BLM_MSH2.combined_doublet, 
                               group.by = SCT_Min_BLM_MSH2.combined_doublet.name) + NoAxes())

pdf(file=" SCT_doublet_.pdf", height= 50,width=15)
p6
dev.off()


p7<-VlnPlot(SCT_Min_BLM_MSH2.combined_doublet, features = "nFeature_RNA",
            group.by = SCT_Min_BLM_MSH2.combined_doublet.name, pt.size = 0.1)

pdf(file=" SCT_Vin_doublet_.pdf", height= 50,width=15)
p7
dev.off()

#remove the predicted doublet from my data
SCT_Seurat_integrated<- SCT_Min_BLM_MSH2.combined_doublet[, SCT_Min_BLM_MSH2.combined_doublet@meta.data[, 
                                                                                                        SCT_Min_BLM_MSH2.combined_doublet.name] == "Singlet"]

saveRDS(SCT_Seurat_integrated ,"SCT_Seurat_integrated.RDS")
#####end of the session

### doing the cluster tree

SCT_Seurat_integrated<-readRDS("SCT_Seurat_integrated.RDS")

## Plotting a cluster tree.

p2 <- clustree(SCT_Seurat_integrated, prefix = "integrated_snn_res.")

pdf(file="sct_cluster-tree_4conditions.pdf", height= 12,width=9)
p2
dev.off()

###genertating the marker list
# Assign identity of clusters 0.2

Idents(object = SCT_Seurat_integrated) <- "integrated_snn_res.0.2"



##visualize the data

p3<-DimPlot(SCT_Seurat_integrated, split.by = "orig.ident")

pdf(file="sct_cluster_4conditions_res0.2.pdf", height= 15,width=20)
p3
dev.off()

###generating the gene marker

###prepeare for finding marker
#seurate asked me to do this before runing findallmarker
##Given a merged object with multiple SCT models, this function uses
#inimum of the median UMI 
##alculated using the raw UMI counts) of individual objects 

SCT_Seurat_integrated<-PrepSCTFindMarkers(SCT_Seurat_integrated, 
                                          assay = "SCT", verbose = TRUE)



goblet<-c("Tff3","Reg4","Muc2","Spdef","Atoh1","muc13","Galnt12","Fcgbp","Muc4",
          "Rep15","Guca2a")

P3<-FeaturePlot(SCT_Seurat_integrated, features = goblet,
                pt.size = 0.3, split.by = "orig.ident"
)& theme(legend.position = c(0.0,0.2))



pdf(file=" SCT_test.pdf", height= 90,width=25)
P3
dev.off()

SCT_Seurat_integrated1<-readRDS("SCT_Seurat_integrated.RDS")


###genertating the marker list
# Assign identity of clusters 0.2

Idents(object = SCT_Seurat_integrated1) <- "integrated_snn_res.0.2"



###generating the gene marker

###prepeare for finding marker
#seurate asked me to do this before runing findallmarker
##Given a merged object with multiple SCT models, this function uses
#inimum of the median UMI 
##alculated using the raw UMI counts) of individual objects 

SCT_Seurat_integrated1<-PrepSCTFindMarkers(SCT_Seurat_integrated1, 
                                           assay = "SCT", verbose = TRUE)

#*******Find Markers********************



markers1 <- FindAllMarkers(SCT_Seurat_integrated1, assay = "SCT",
                           slot = "data", min.pct = 0.25,return.thresh = 0.05,
                           logfc.threshold = log(1))

saveRDS(markers1,  "markers1.RDS")

saveRDS(SCT_Seurat_integrated1,  "SCT_Seurat_integrated1.RDS")

##visualize the data
markers1<-readRDS("markers1.RDS")
tail(markers1)

str(markers1)

##convert it to data.table 
library("data.table")

setDT(markers1)#to convert it to data.table

## Load and prepare marker list.

markers <- markers1[p_val_adj < 0.05]
markers[, c("avg_log2FC", "cluster") := list(log2(exp(avg_log2FC)), str_c("cluster_", cluster))]
markers <- markers[order(cluster, -avg_log2FC)]


## Split the list based on cluster and save results.

markers<- split(markers, markers$cluster)

str(markers)


myresults_marker<-list(markers$cluster_1,markers$cluster_2,markers$cluster_3,markers$cluster_4,
           markers$cluster_5,markers$cluster_6,markers$cluster_7,markers$cluster_8,
        markers$cluster_9,markers$cluster_10,markers$cluster_11,markers$cluster_12,
            markers$cluster_13,markers$cluster_14,markers$cluster_15,markers$cluster_16,
          markers$cluster_17)
myresults2<-c("cluster_1","cluster_2","cluster_3",
              "cluster_4","cluster_5", "cluster_6",
              "cluster_7","cluster_8","cluster_9",
              "cluster_10","cluster_11","cluster_12",
              "cluster_13","cluster_14","cluster_15",
              "cluster_16", "cluster_17")

for (i in (1:length(myresults_marker))){
  n<-myresults_marker[[i]]
  write.csv(n, paste0(myresults2[i],".csv"),sep='\t', 
            quote = F, row.names = T, col.names = T) 
}

####
## Load and prepare marker list for genes having Log2FC more than 1.

markers2 <- markers1[p_val_adj < 0.05 &avg_log2FC>=1]
markers2[, c("avg_log2FC", "cluster") := list(log2(exp(avg_log2FC)), str_c("cluster_", cluster))]
markers2 <- markers2[order(cluster, -avg_log2FC)]


str(markers2)
## Split the list based on cluster and save results.

markers2<- split(markers2, markers2$cluster)

str(markers)
markers2$`1`

myresults_marker<-list(markers2$cluster_1,markers2$cluster_2,markers2$cluster_3,
                       markers2$cluster_4,markers2$cluster_5,markers2$cluster_6,
                       markers2$cluster_7,markers2$cluster_8,markers2$cluster_9,
                       markers2$cluster_10,markers2$cluster_11,markers2$cluster_12,
                       markers2$cluster_13,markers2$cluster_14,markers2$cluster_15,
                       markers2$cluster_16,markers2$cluster_17)
myresults2<-c("cluster_FC-1_1","cluster_FC-1_2","cluster_FC-1_3",
              "cluster_FC-1_4","cluster_FC-1_5", "cluster_FC-1_6",
              "cluster_FC-1_7","cluster_FC-1_8","cluster_FC-1_9",
              "cluster_FC-1_10","cluster_FC-1_11","cluster_FC-1_12",
              "cluster_FC-1_13","cluster_FC-1_14","cluster_FC-1_15",
              "cluster_FC-1_16", "cluster_FC-1_17")

for (i in (1:length(myresults_marker))){
  n<-myresults_marker[[i]]
  write.csv(n, paste0(myresults2[i],".csv"),sep='\t', 
            quote = F, row.names = T, col.names = T) 
}


####doing featureplot for some marker genes

SCT_Seurat_integrated1<-readRDS("SCT_Seurat_integrated1.RDS")

DefaultAssay(SCT_Seurat_integrated1)<-"SCT"

P3<-FeaturePlot(SCT_Seurat_integrated1, features = c("Reg3b","Reg3g"),
                pt.size = 0.3, split.by = "orig.ident"
)& theme(legend.position = c(0.0,0.2))



pdf(file=" SCT_Reg_genes.pdf", height= 30,width=20)
P3
dev.off()



##visualize the data with number on each cluster

z<-DimPlot(SCT_Seurat_integrated1, split.by = "orig.ident",label=TRUE)


pdf("SCT_cluster_4condtitons with number1.pdf",height=15, width = 20)
z
dev.off()


8888888888******************************************************************8
###reducing the resolution from 0.2 to 0.4

###genertating the marker list
# Assign identity of clusters 0.4

SCT_Seurat_integrated1<-readRDS("SCT_Seurat_integrated1.RDS")

Idents(object = SCT_Seurat_integrated1) <- "integrated_snn_res.0.4"

##visualize the data

p3<-DimPlot(SCT_Seurat_integrated1, split.by = "orig.ident",label=TRUE)

pdf(file="sct_cluster_4conditions_res0.4_number.pdf", height= 15,width=20)
p3
dev.off()

###generating marker at this resolutions
DefaultAssay(SCT_Seurat_integrated1)<-"SCT"



#*******Find Markers********************



markers_Res0.4 <- FindAllMarkers(SCT_Seurat_integrated1, assay = "SCT",
                           slot = "data", min.pct = 0.25,return.thresh = 0.05,
                           logfc.threshold = log(1),recorrect_umi=FALSE)

saveRDS(markers_Res0.4,  "markers_Res0.4.RDS")

saveRDS(SCT_Seurat_integrated1,  "SCT_Seurat_integrated1.RDS")

##visualize the data
markers_Res0.4<-readRDS("markers_Res0.4.RDS")
tail(markers_Res0.4)

str(markers_Res0.4)

##convert it to data.table 
library("data.table")
library("stringr")

setDT(markers_Res0.4)#to convert it to data.table

## Load and prepare marker list.

markers_Res0.4_1 <- markers_Res0.4[p_val_adj < 0.05]
markers_Res0.4_1[, c("avg_log2FC", "cluster") := list(log2(exp(avg_log2FC)), str_c("cluster_", cluster))]
markers_Res0.4_1 <- markers_Res0.4_1[order(cluster, -avg_log2FC)]


## Split the list based on cluster and save results.

markers_Res0.4_1<- split(markers_Res0.4_1, markers_Res0.4_1$cluster)

markers<-markers_Res0.4_1######just to make it easy

str(markers)


myresults_marker<-list(markers$cluster_0,markers$cluster_1,markers$cluster_2,markers$cluster_3,markers$cluster_4,
                       markers$cluster_5,markers$cluster_6,markers$cluster_7,markers$cluster_8,
                       markers$cluster_9,markers$cluster_10,markers$cluster_11,markers$cluster_12,
                       markers$cluster_13,markers$cluster_14,markers$cluster_15,markers$cluster_16,
                       markers$cluster_17,markers$cluster_18,markers$cluster_19,
                       markers$cluster_20,
                       markers$cluster_21)
myresults2<-c("cluster_0","cluster_1","cluster_2","cluster_3",
              "cluster_4","cluster_5", "cluster_6",
              "cluster_7","cluster_8","cluster_9",
              "cluster_10","cluster_11","cluster_12",
              "cluster_13","cluster_14","cluster_15",
              "cluster_16", "cluster_17","cluster_18","cluster_19","cluster_20",
              "cluster_21")

for (i in (1:length(myresults_marker))){
  n<-myresults_marker[[i]]
  write.csv(n, paste0(myresults2[i],".csv"),sep='\t', 
            quote = F, row.names = T, col.names = T) 
}

############I am not doing this one right now
## Load and prepare marker list for genes having Log2FC more than 0.75

markers2 <- markers_Res0.4[p_val_adj < 0.05 &avg_log2FC>=0.75]
markers2[, c("avg_log2FC", "cluster") := list(log2(exp(avg_log2FC)), str_c("cluster_", cluster))]
markers2 <- markers2[order(cluster, -avg_log2FC)]


str(markers2)
## Split the list based on cluster and save results.

markers2<- split(markers2, markers2$cluster)

str(markers2)
markers2$`1`

myresults_marker<-list(markers2$cluster_1,markers2$cluster_2,markers2$cluster_3,
                       markers2$cluster_4,markers2$cluster_5,markers2$cluster_6,
                       markers2$cluster_7,markers2$cluster_8,markers2$cluster_9,
                       markers2$cluster_10,markers2$cluster_11,markers2$cluster_12,
                       markers2$cluster_13,markers2$cluster_14,markers2$cluster_15,
                       markers2$cluster_16,markers2$cluster_17)
myresults2<-c("cluster_FC-1_1","cluster_FC-1_2","cluster_FC-1_3",
              "cluster_FC-1_4","cluster_FC-1_5", "cluster_FC-1_6",
              "cluster_FC-1_7","cluster_FC-1_8","cluster_FC-1_9",
              "cluster_FC-1_10","cluster_FC-1_11","cluster_FC-1_12",
              "cluster_FC-1_13","cluster_FC-1_14","cluster_FC-1_15",
              "cluster_FC-1_16", "cluster_FC-1_17")

for (i in (1:length(myresults_marker))){
  n<-myresults_marker[[i]]
  write.csv(n, paste0(myresults2[i],".csv"),sep='\t', 
            quote = F, row.names = T, col.names = T) 
}


#### removing all of the clusters except colon and converting to loom file
SCT_Seurat_colon<-subset(SCT_Seurat_integrated,idents = c(0,1,3,4,5,7,9))
saveRDS(SCT_Seurat_colon,  "SCT_Seurat_colon.RDS")

colon.loom <- as.loom(SCT_Seurat_colon, filename = "colon.loom", verbose = FALSE)
colon.loom$close_all()

p3<-DimPlot(SCT_Seurat_colon, split.by = "orig.ident",label=TRUE)

pdf(file="colon_res0.2_number1.pdf", height= 15,width=20)
p3
dev.off()

seurat_velocity <- RunVelocity( SCT_Seurat_colon, ambiguous = "ambiguous", ncores = 1,
  deltaT = 1, kCells = 25, fit.quantile = 0.02)


###converting seurate object to AnnData for velocity
library(reticulate)
library(anndata)
library(loomR)
SCT_Seurat_colon<-readRDS("SCT_Seurat_colon.RDS")
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




###convert seurat object to loom file
colon.loom <- as.loom(SCT_Seurat_colon, filename = "colon.loom", verbose = FALSE)
colon.loom$close_all()



###calculating differentiation potential
library(devtools)
devtools::install_github("aet21/SCENT")
