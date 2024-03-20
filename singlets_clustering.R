library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
library(scCATCH)
library(haven)

seurat_obj_subset <- readRDS("/data/My_scRNAseq_subset_scaled.rds")

###getting doublets info
doublets_metadata <- read.csv("/data/doublet_finder_14March.csv")
colnames(doublets_metadata) <- c("cellID", "doublet_finder_res","doublet_finder")

names <- row.names(seurat_obj_subset@meta.data)
seurat_meta <- seurat_obj_subset@meta.data %>%
  mutate(cellID = names) %>%
  left_join(doublets_metadata,by="cellID")
rownames(seurat_meta) <- names
seurat_obj_subset@meta.data <- seurat_meta


### subseting singlets

seurat_obj_singlets <- subset(seurat_obj_subset, doublet_finder_res == "Singlet")


### clustering and umap vizualization

seurat_obj_singlets <- FindNeighbors(seurat_obj_singlets, dims = 1:10,graph.name = "Dim10")
seurat_obj_singlets <- FindClusters(seurat_obj_singlets , resolution = 0.3,graph.name = "Dim10")

seurat_obj_singlets <- RunUMAP(seurat_obj_singlets, dims = 1:10)
DimPlot(seurat_obj_singlets, reduction = "umap",label = TRUE)

FeaturePlot(seurat_obj_singlets,features = c("MS4A1","CD19","CD8A","CD8B","CD4","KLRB1","NKG7","NCAM1","PPBP","CD14","CD68","FCGR3A","CD1C","LILRA4"))

FeaturePlot(seurat_obj_singlets,features = c("CD1C","LILRA4"))
FeaturePlot(seurat_obj_singlets,features = c("MS4A1","CD19","CD8B","CD4"))
FeaturePlot(seurat_obj_singlets,features = c("CD14","CD68","FCGR3A"))


DotPlot(seurat_obj_singlets,features = c("MS4A1","CD19","CD3D","CD2","CD8B","CD4","KLRB1","NKG7","NCAM1","PPBP","CD14","CD68","FCGR3A","CD1C","LILRA4"))

new_cluster_ids <- c("CD4 Tcell (1)", "CD14+ Monocytes", "CD8 Tcell (1)", "CD4 Tcell (2)", "NKs", "CD4 Tcell (3)", "B cell", "CD8 Tcell (2)", "FCGR3A+ Monocytes", "DC","Megakaryocytes and platelet")
names(new_cluster_ids) <- levels(seurat_obj_singlets)
seurat_obj_singlets <- RenameIdents(seurat_obj_singlets, new_cluster_ids)

DimPlot(seurat_obj_singlets, reduction = "umap", label = TRUE) + NoLegend()

DotPlot(seurat_obj_singlets,features = c("CCR7","SELL","S100A4"))
DotPlot(seurat_obj_singlets,features = c("TESPA1","GPR183","CD28"))

### adding patietns data 
patient_meta <- read_csv("~/Projects/scRNAseq_Omniscope/IMC-F106C-101_patient_groups.csv")

names <- row.names(seurat_obj_singlets@meta.data)
seurat_meta <- seurat_obj_singlets@meta.data %>%
  mutate(Patient = names) %>%
  mutate(Patient = gsub(".*-", "", Patient)) %>%
  left_join(patient_meta %>%
              mutate(SUBJID = as.character(SUBJID)), by = c("Patient" = "SUBJID"))
rownames(seurat_meta) <- names
seurat_obj_singlets@meta.data <- seurat_meta

saveRDS(seurat_obj_singlets,file="/data/My_seurat_obj_singlets_lvl1.rds")

#### SUBCLUSTERING
seurat_obj_singlets <- readRDS("/data/My_seurat_obj_singlets_lvl1.rds")

FeaturePlot(seurat_obj_singlets,features = "nFeature_RNA")

FeaturePlot(seurat_obj_singlets,features = "nCount_RNA")

FeaturePlot(seurat_obj_singlets,features = "percent.mt")

FeaturePlot(seurat_obj_singlets,features = "percent.ribo")


### DC subclastering

seurat_obj_singlets <- FindSubCluster(seurat_obj_singlets,"DC",resolution = 0.1, graph.name = "Dim10",subcluster.name = "DCs")
DimPlot(seurat_obj_singlets, reduction = "umap",group.by="DCs",label = TRUE) + NoLegend()

seurat_obj_singlets <- SetIdent(seurat_obj_singlets,value = seurat_obj_singlets@meta.data$DCs)
DimPlot(seurat_obj_singlets, reduction = "umap",label = TRUE) + NoLegend()
DotPlot(seurat_obj_singlets,features = c("CD1C","LILRA4"))

seurat_obj_singlets <- RenameIdents(seurat_obj_singlets, "DC_0"="mDC","DC_1"="pDC")
DotPlot(seurat_obj_singlets,features = c("CD1C","LILRA4"))

### B cells clusters

seurat_obj_singlets <- FindSubCluster(seurat_obj_singlets,"B cell",resolution = 0.15, graph.name = "Dim10",subcluster.name = "Bcells")
DimPlot(seurat_obj_singlets, reduction = "umap",group.by="Bcells",label = TRUE) + NoLegend()

seurat_obj_singlets <- SetIdent(seurat_obj_singlets,value = seurat_obj_singlets@meta.data$Bcells)

DotPlot(seurat_obj_singlets,features = c("MS4A1","XBP1","CD27","CD38","TCL1A","IGHG1","IGHA1"))
seurat_obj_singlets <- RenameIdents(seurat_obj_singlets, "B cell_0"="B Naive","B cell_1"="B Memory","B cell_2"="Plasma","B cell_3"="B Naive")
DotPlot(seurat_obj_singlets,features = c("MS4A1","XBP1","CD27","CD38","TCL1A","IGHG1","IGHA1"))

DimPlot(seurat_obj_singlets, reduction = "umap",label = TRUE) + NoLegend()

### NK clusters

seurat_obj_singlets <- FindSubCluster(seurat_obj_singlets,"NKs",resolution = 0.1, graph.name = "Dim10",subcluster.name = "NK")
DimPlot(seurat_obj_singlets, reduction = "umap",group.by="NK",label = TRUE) + NoLegend()
seurat_obj_singlets <- SetIdent(seurat_obj_singlets,value = seurat_obj_singlets@meta.data$NK)
DimPlot(seurat_obj_singlets, reduction = "umap",label = TRUE) + NoLegend()

DotPlot(seurat_obj_singlets,features = c("KLRB1","NKG7","NCAM1","PPBP","FCGR3A","CD3D","CD3E","CD2"))

NK1_markers <- FindMarkers(seurat_obj_singlets,ident.1="NKs_1",ident.2="NKs_0", only.pos = TRUE)
head(NK1_markers,n=20)

NK2_markers <- FindMarkers(seurat_obj_singlets,ident.1="NKs_2",ident.2="NKs_0", only.pos = TRUE)
head(NK2_markers,n=20)

seurat_obj_singlets <- RenameIdents(seurat_obj_singlets, "NKs_0"="NKs","NKs_1"="NKT","NKs_2"="PPBP+ NKs")
DotPlot(seurat_obj_singlets,features = c("KLRB1","NKG7","NCAM1","PPBP","FCGR3A","CD3D","CD3E","CD2"))

## Monocyte clustering

DimPlot(seurat_obj_singlets, reduction = "umap",label = TRUE,split.by = "group") + NoLegend()

seurat_obj_singlets <- FindSubCluster(seurat_obj_singlets,"CD14+ Monocytes",resolution = 0.1, graph.name = "Dim10",subcluster.name = "C14Mono")
DimPlot(seurat_obj_singlets, reduction = "umap",group.by="C14Mono",label = TRUE,split.by = "group") + NoLegend()
DimPlot(seurat_obj_singlets, reduction = "umap",group.by="C14Mono",label = TRUE) + NoLegend()

seurat_obj_singlets <- SetIdent(seurat_obj_singlets,value = seurat_obj_singlets@meta.data$C14Mono)

CD14Mono_Normal_markers <- FindMarkers(seurat_obj_singlets,ident.1="CD14+ Monocytes_2",ident.2=c("CD14+ Monocytes_0","CD14+ Monocytes_1"), only.pos = TRUE)
head(CD14Mono_Normal_markers,n=20)

DotPlot(seurat_obj_singlets,features = c("S100A8","S100A9","CD14","FCGR3A","LYZ"))

dontknow_markers <- FindMarkers(seurat_obj_singlets,ident.1="Dont Know",ident.2="CD14+ Monocytes_2", only.pos = TRUE)
dontknow_markers_v2 <- FindMarkers(seurat_obj_singlets,ident.1="Dont Know", only.pos = TRUE)

write.csv(dontknow_markers_v2,file="~/Projects/scRNAseq_Omniscope/NotKnown_cluster_vs_all_s100a8_s100a9.csv")
saveRDS(seurat_obj_singlets,file="/data/seurat_obj_singlets_lvl2.rds")

seurat_obj_singlets <- readRDS("/data/seurat_obj_singlets_lvl2.rds")

### subseting T cells

seurat_obj_Tcells <- subset(seurat_obj_singlets, NK %in% c("CD4 Tcell (1)", "CD4 Tcell (2)","CD4 Tcell (3)","CD8 Tcell (1)","CD8 Tcell (2)"))
saveRDS(seurat_obj_Tcells,file="/data/seurat_obj_Tcells.rds")

##some plots
FeaturePlot(seurat_obj_singlets,features = c("TESPA1"),split.by = "group")
FeaturePlot(seurat_obj_singlets,features = c("GPR183"),split.by = "group")
FeaturePlot(seurat_obj_singlets,features = c("CD28"),split.by = "group")


DimPlot(seurat_obj_singlets, reduction = "umap",split.by="group",label = TRUE) + NoLegend()
DotPlot(seurat_obj_singlets,features=c("TESPA1","GPR183","CD28"),split.by = "group",cols=c("black","red","green"))


VlnPlot(seurat_obj_singlets,features="TESPA1",log=TRUE,pt.size=0,split.by = "group")
VlnPlot(seurat_obj_singlets,features="GPR183",log=TRUE,pt.size=0,split.by = "group")
VlnPlot(seurat_obj_singlets,features="CD28",log=TRUE,pt.size=0,split.by = "group")


