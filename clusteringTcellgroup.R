library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
library(scCATCH)
library(haven)

seurat_obj_Tcells <- readRDS(file="/data/seurat_obj_Tcells.rds")

seurat_obj_Tcells  <- FindNeighbors(seurat_obj_Tcells , dims = 1:15,graph.name = "Dim15")
seurat_obj_Tcells  <- FindClusters(seurat_obj_Tcells  , resolution = 0.8,graph.name = "Dim15")

seurat_obj_Tcells  <- RunUMAP(seurat_obj_Tcells , dims = 1:15)
#DimPlot(seurat_obj_Tcells , reduction = "umap",label=TRUE, split.by="group")
DimPlot(seurat_obj_Tcells , reduction = "umap",label=TRUE)

## plots

#FeaturePlot(seurat_obj_Tcells,features = c("CD4","CD8A","CD8B","CD34"))
#FeaturePlot(seurat_obj_Tcells,features = c("TESPA1","GPR183","CD28"))
#FeaturePlot(seurat_obj_Tcells,features = c("GZMK","GZMB"))
#FeaturePlot(seurat_obj_Tcells,features = c("CCR7","SELL","S100A4"))

VlnPlot(seurat_obj_Tcells,features=c("PCLAF","MKI67","STMN1"),log=TRUE,pt.size=0)

#DotPlot(seurat_obj_Tcells,features = c("CD4","CD8A","CD8B","CD34"))
DotPlot(seurat_obj_Tcells,features = c("PCLAF","MKI67","STMN1","ATXN1","SESN3","TESPA1"))


####resolution 0.65

all_markers <- FindAllMarkers(seurat_obj_Tcells,only.pos = TRUE)
saveRDS(all_markers,file="/data/all_markers_Tcell_16March.rds")
  
all_markers <- readRDS("/data/all_markers_Tcell_16March.rds")

all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >1) %>%
  slice_head(n=10) %>%
  ungroup() -> top10


DoHeatmap(subset(seurat_obj_Tcells,downsample=100),features=top10$gene)

new_cluster_ids <- c("CD4_act_1","CD8 CTL_1","CD4_act_2","CD4_act_2","CD8 Naive","CD4 Naive (1)","CD4 Naive (1)","CD4 Naive (2)","CD8 EM","CD4 CM","MAIT","CD4 act_3","CD8 CM","PPBP+ CD4","CD8 CTL_2","CD8(?)")
names(new_cluster_ids) <- levels(seurat_obj_Tcells)
seurat_obj_Tcells <- RenameIdents(seurat_obj_Tcells, new_cluster_ids)

DimPlot(seurat_obj_Tcells , reduction = "umap",label=TRUE) + NoLegend()

## plots

VlnPlot(seurat_obj_Tcells,features=c("CD4","CD8B","CD8A","CD3"),log=TRUE,pt.size=0)

VlnPlot(seurat_obj_Tcells,features=c("CCR7","SELL","S100A4"),log=TRUE,pt.size=0, idents=c("CD4_act_1","CD4_act_2","CD4 Naive (1)","CD4 Naive (2)","CD4 CM","CD4 act_3","PPBP+ CD4"))
VlnPlot(seurat_obj_Tcells,features=c("CCR7","SELL","S100A4"),log=TRUE,pt.size=0, idents=c("CD8 CTL_1","CD8 Naive","CD8 EM","MAIT","CD8 CM","CD8 CTL_2","CD8(?)"))
VlnPlot(seurat_obj_Tcells,features=c("KLRB1","GZMK","GZMB","NCAM1"),log=TRUE,pt.size=0)
VlnPlot(seurat_obj_Tcells,features=c("ITGB1","AHNAK","CRIP1"),log=TRUE,pt.size=0, idents=c("CD4_act_1","CD4_act_2","CD4 Naive (1)","CD4 Naive (2)","CD4 CM","CD4 act_3","PPBP+ CD4"))

VlnPlot(seurat_obj_Tcells,features=c("MKI67","PPBP"),log=TRUE,pt.size=0)

markers_qustion <- FindMarkers(seurat_obj_Tcells,ident.1="CD8(?)",only.pos = TRUE)
FeaturePlot(seurat_obj_Tcells,features = c("KIT","CD34"))

##finding HSPC cluster
seurat_obj_Tcells <- FindSubCluster(seurat_obj_Tcells,"PPBP+ CD4",resolution = 0.1, graph.name = "Dim15",subcluster.name = "PPBP_CD4")
DimPlot(seurat_obj_Tcells, reduction = "umap",group.by="PPBP_CD4",label = TRUE) + NoLegend()
seurat_obj_Tcells <- SetIdent(seurat_obj_Tcells,value = seurat_obj_Tcells@meta.data$PPBP_CD4)
seurat_obj_Tcells <- RenameIdents(seurat_obj_Tcells, "PPBP+ CD4_0"="PPBP+ CD4","PPBP+ CD4_1"="HPSC")

VlnPlot(seurat_obj_Tcells,features=c("KIT","KLRB1","GATA3","CD8A","CD4","NCAM1"),log=TRUE,pt.size=0)
VlnPlot(seurat_obj_Tcells,features=c("CD69"),log=TRUE,pt.size=0,idents=c("CD4_act_1","CD4_act_2","CD4 Naive (1)","CD4 Naive (2)","CD4 CM","CD4 act_3","PPBP+ CD4"))
seurat_obj_Tcells <- RenameIdents(seurat_obj_Tcells, "CD8(?)"="DN T cells")

saveRDS(seurat_obj_Tcells,file="/data/seurat_obj_Tcells.rds")

FeaturePlot(seurat_obj_Tcells,features = c("FOXP3"))

seurat_obj_Tcells <- RenameIdents(seurat_obj_Tcells, "CD4_act_2"="CD4_act_1")

seurat_obj_Tcells <- FindSubCluster(seurat_obj_Tcells,"CD4_act_1",resolution = 0.5, graph.name = "Dim15",subcluster.name = "CD4_act_withTregs")
DimPlot(seurat_obj_Tcells, reduction = "umap",group.by="CD4_act_withTregs",label = TRUE) + NoLegend()

VlnPlot(seurat_obj_Tcells,features=c("FOXP3"),log=TRUE,pt.size=0,group.by="CD4_act_withTregs")
DotPlot(seurat_obj_Tcells,features = c("FOXP3"),group.by="CD4_act_withTregs")

DimPlot(subset(seurat_obj_Tcells, CD4_act_withTregs %in% c("CD4_act_1_5","CD4_act_1_7","CD4_act_1_0","CD4 Naive (1)","CD8 CTL_1")), reduction = "umap",group.by="CD4_act_withTregs",label = TRUE) + NoLegend()
DimPlot(subset(seurat_obj_Tcells, CD4_act_withTregs %in% c("CD4_act_1_3","CD4_act_1_4","CD4 Naive (1)","CD8 CTL_1")), reduction = "umap",group.by="CD4_act_withTregs",label = TRUE) + NoLegend()
DimPlot(subset(seurat_obj_Tcells, CD4_act_withTregs %in% c("CD4_act_1_1","CD4_act_1_2","CD4_act_1_7","CD4 Naive (1)","CD8 CTL_1")), reduction = "umap",group.by="CD4_act_withTregs",label = TRUE) + NoLegend()
DimPlot(subset(seurat_obj_Tcells, CD4_act_withTregs %in% c("CD4_act_1_6","CD4_act_1_0","CD4_act_1_7","CD4 Naive (1)","CD8 CTL_1")), reduction = "umap",group.by="CD4_act_withTregs",label = TRUE) + NoLegend()

seurat_obj_Tcells <- SetIdent(seurat_obj_Tcells,value = seurat_obj_Tcells@meta.data$CD4_act_withTregs)

seurat_obj_Tcells <- RenameIdents(seurat_obj_Tcells, "CD4_act_1_5"="FOXP3+ Tregs","CD4_act_1_0"="CD4_act_1",
                                  "CD4_act_1_1"="CD4_act_1",
                                  "CD4_act_1_2"="CD4_act_1",
                                  "CD4_act_1_3"="CD4_act_1",
                                  "CD4_act_1_4"="CD4_act_1",
                                  "CD4_act_1_6"="CD4_act_1",
                                  "CD4_act_1_7"="CD4_act_1")

DotPlot(seurat_obj_Tcells,features = c("FOXP3","TRDC","TYMS"))

DimPlot(seurat_obj_Tcells, reduction = "umap",label = TRUE) + NoLegend()
VlnPlot(seurat_obj_Tcells,features=c("CCR7","SELL","S100A4"),log=TRUE,pt.size=0,idents = c("CD4 Naive (1)","CD4 Naive (2)","CD4 CM","CD4 act_1","CD4 act_3","FOXP3+ Tregs","PPBP+ CD4"))

FeaturePlot(subset(seurat_obj_Tcells,idents=c("FOXP3+ Tregs","CD4 Naive (1)","CD8 CTL_1")),features = c("FOXP3"))

DimPlot(seurat_obj_Tcells, reduction = "umap",label = TRUE,split.by="group") + NoLegend()
DotPlot(seurat_obj_Tcells,features=c("TESPA1","GPR183","CD28"))

VlnPlot(seurat_obj_Tcells,features="TESPA1",log=TRUE,pt.size=0,split.by = "group")
VlnPlot(seurat_obj_Tcells,features="GPR183",log=TRUE,pt.size=0,split.by = "group")
VlnPlot(seurat_obj_Tcells,features="CD28",log=TRUE,pt.size=0,split.by = "group")

FeaturePlot(seurat_obj_Tcells,features = c("TESPA1"),split.by = "group")
FeaturePlot(seurat_obj_Tcells,features = c("GPR183"),split.by = "group")
FeaturePlot(seurat_obj_Tcells,features = c("CD28"),split.by = "group")


