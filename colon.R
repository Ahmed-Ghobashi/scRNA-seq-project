
library(reticulate)
library("Seurat")
library(SeuratWrappers)
library(sctransform)
library(glmGamPoi)
library(htmltools)
library(tidyverse)
library(scCATCH)
library(leiden)
library(clustree)
###taking only the epithelial cells and re doing the analysis,

SCT_Seurat_integrated1<-readRDS("SCT_Seurat_integrated1.RDS")##Duplicate removed
###get the colon cells
Seurat_colon<-subset(SCT_Seurat_integrated1,idents = c(0,1,3,4,5,7,9))
Seurat_fibroblast<-subset(SCT_Seurat_integrated1,idents = c(2,14))
Seurat_immuncells<-subset(SCT_Seurat_integrated1,idents = c(2,14,0,1,3,4,5,7,9),
                          invert=TRUE)

###in order to do reclustering, I will split the samples and redo the preprocessing 
###before integration

###1##split the samples into list
colon_list <- SplitObject(Seurat_colon, split.by = "orig.ident")
fibroblast_list <- SplitObject(Seurat_fibroblast, split.by = "orig.ident")
immun_list <- SplitObject(Seurat_immuncells, split.by = "orig.ident")

##do SCT transformation
colon_list  <- lapply(colon_list, SCTransform, method = "glmGamPoi")
fibroblast_list  <- lapply(fibroblast_list, SCTransform, method = "glmGamPoi")
immun_list  <- lapply(immun_list, SCTransform, method = "glmGamPoi")

### prepare for integration
features_colon  <- SelectIntegrationFeatures(colon_list, nfeatures = 3000)
features_fibro  <- SelectIntegrationFeatures(fibroblast_list, nfeatures = 3000)
features_immun  <- SelectIntegrationFeatures(immun_list, nfeatures = 3000)

colon_list<- PrepSCTIntegration(colon_list, anchor.features = features_colon)
fibroblast_list<- PrepSCTIntegration(fibroblast_list, anchor.features = features_fibro)
immun_list<- PrepSCTIntegration(immun_list, anchor.features = features_immun)

##making the anchores

anchors_colon <- FindIntegrationAnchors(colon_list, normalization.method = "SCT", 
                                  anchor.features = features_colon)
anchors_fibro <- FindIntegrationAnchors(fibroblast_list, normalization.method = "SCT", 
                                        anchor.features = features_fibro)
anchors_immun <- FindIntegrationAnchors(immun_list, normalization.method = "SCT", 
                                        anchor.features = features_immun)

##do the integration

Seurat_colon <- IntegrateData(anchorset = anchors_colon, normalization.method = "SCT",
                                 k.weight = 20)
Seurat_fibroblast <- IntegrateData(anchorset = anchors_fibro, normalization.method = "SCT",
                                 k.weight = 20)
Seurat_immuncells <- IntegrateData(anchorset = anchors_immun, normalization.method = "SCT",
                                 k.weight = 20)
##rerun SCT trnasformation

p<- ElbowPlot(Seurat_colon, ndims = 100)#the variation is high till almost 75

pdf("ElbowPlot_fibro.pdf")
ElbowPlot(Seurat_fibroblast, ndims = 100)
dev.off()

## umap dimension reduction for visualization.
##pca
Seurat_colon <- RunPCA(Seurat_colon, verbose = FALSE,npcs = 100)
Seurat_fibroblast <- RunPCA(Seurat_fibroblast, verbose = FALSE,npcs = 100)
Seurat_immuncells <- RunPCA(Seurat_immuncells, verbose = FALSE,npcs = 100)
##umap

Seurat_colon<- RunUMAP(Seurat_colon, reduction = "pca", dims = 1:30)
Seurat_fibroblast<- RunUMAP(Seurat_fibroblast, reduction = "pca", dims = 1:30)
Seurat_immuncells<- RunUMAP(Seurat_immuncells, reduction = "pca", dims = 1:30)


##clustering
Seurat_colon<- FindNeighbors(Seurat_colon , dims = 1:50)
Seurat_fibroblast<- FindNeighbors(Seurat_fibroblast, dims = 1:50)
Seurat_immuncells<- FindNeighbors(Seurat_immuncells , dims = 1:50)
####clustering at different reslution
Seurat_colon<- FindClusters(Seurat_colon, resolution = c(0.2,0.4, 0.6, 0.8, 1.0))
Seurat_fibroblast<- FindClusters(Seurat_fibroblast, resolution = c(0.2,0.4, 0.6, 0.8, 1.0))
Seurat_immuncells<- FindClusters(Seurat_immuncells, resolution = c(0.2,0.4, 0.6, 0.8, 1.0))
                 

p2 <- clustree(Seurat_immuncells, prefix = "integrated_snn_res.")

pdf(file="immun_clustertree.pdf", height= 12,width=9)
p2
dev.off()

# Assign identity of clusters 0.08
Idents(object = Seurat_colon) <- "integrated_snn_res.0.08"
Idents(object = Seurat_fibroblast) <- "integrated_snn_res.0.2"
Idents(object = Seurat_immuncells) <- "integrated_snn_res.0.08"


##visualize the data
seurats<-list(Seurat_colon,Seurat_fibroblast,Seurat_immuncells)
names(seurats)<-c("colon","Fibroblast","Immun_cells")###to add name to the list 

####saving the serutat objects
for (ii in names(seurats)){
  filename <- paste(ii, ".RDS", sep="")
  saveRDS(seurats[[ii]], filename)
}


pdf("Fibro_cluster-split.pdf",height=15, width = 20)
DimPlot(Seurat_fibroblast, split.by = "orig.ident")
dev.off()

##feature plots
DefaultAssay(Seurat_colon)<-"SCT"

goblet<-c("Tff3","Reg4","Muc2","Spdef","Atoh1","muc13","Galnt12","Fcgbp","Muc4",
          "Rep15","Guca2a")

stem<-c("Lgr5","Ascl2","Smoc2",
        "Axin2","Ephb2","Kcne3","Hmgn2","Id1")

EEC<-c("Chga","Dll1","Pax4","Neurog3","Insm1","Neurod1","Cpe",
       "Isl1","Pax6")

Tuft_brush<-c("Cdhr2","Espn")

epithelial<-c("Epcam","Krt8","Krt18")

reg<-c("Reg3b", "Reg3a")

P3<-FeaturePlot(Seurat_colon, features =c("Lyz1"),
                pt.size = 0.3, split.by = "orig.ident"
)& theme(legend.position = c(0.0,0.2))



pdf(file=" paneth_cells.pdf", height=20,width=25)
P3
dev.off()

, height= 110,width=15


####focusing on colon epithelial cells.

Seurat_colon<-readRDS("colon.RDS")
DefaultAssay(Seurat_colon)<-"SCT"
Idents(object = Seurat_colon) <- "integrated_snn_res.0.2"

Microfold_cells<-c('Fabp5','Mtf1','Relb','Tulp4')
tuft<-c('Gng13',"Rgs2",'Klf6','Klf3','Siglecf')
endothelial<-c('Pecam1','Egfl7')
Enrocytes<-c('Lgals2','Alpi','Slc10a2')
enterocytes1<-c('Rnf186', 'Gata4','Cd63','Cldn15')
enterochromafin<-c('Tph1','Adm','Crh','Afp')
P3<-FeaturePlot(Seurat_colon, features = selected_marker_1,
                pt.size = 0.3, split.by = "orig.ident"
)& theme(legend.position = c(0.0,0.2))



pdf(file="selected_marker_1.pdf", height= 15,width=19)
P3
dev.off()
###try different resolution
Idents(object = Seurat_colon) <- "integrated_snn_res.0.2"

pdf("colon_2-split.pdf",height=10, width = 20)
DimPlot(Seurat_colon, split.by = "orig.ident",label = TRUE)
dev.off()

###the number of cluuster from the whole cell dataset
selected_marker_0<-c('Smoc2','Prdx2','Mgp')
selected_marker_1<-c('Nrg1','Pcdh7','Piezo2','Saa3','Mmp10')
selected_marker_3<-c('Nrg1','Pcdh7','Piezo2','Saa3')

###doing ScBATCH to automatically annotate the clusters
#Create scCATCH object from Seurat object 
Seurat_colon<-readRDS("colon.RDS")
DefaultAssay(Seurat_colon)<-"SCT"
Idents(object = Seurat_colon) <- "integrated_snn_res.0.2"


obj <- createscCATCH(data = Seurat_colon[['RNA']]@data,
                     cluster = as.character(Idents(Seurat_colon)))

# find highly expressed genes
obj <- findmarkergene(object = obj, species = "Mouse", 
                      marker = cellmatch, tissue = "Small intestine")

obj_1 <- findmarkergene(object = obj, species = "Mouse",
                      marker = cellmatch, tissue = "Colon")##didn't work,not enough marker

#Evidence-based score and annotation for each cluster

obj <- findcelltype(object = obj)

saveRDS(obj,'scBATCH.RDS')

###find markers in colon

##calculated using the raw UMI counts) of individual objects 

Seurat_colon_1<-PrepSCTFindMarkers(Seurat_colon, 
                                           assay = "SCT", verbose = TRUE)

colon_markers <- FindAllMarkers(Seurat_colon_1, assay = "SCT",
                           slot = "data", min.pct = 0.25,return.thresh = 0.05,
                           logfc.threshold = log(1))

saveRDS(colon_markers,  "colon_markers.RDS")


