library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
library(scCATCH)
library(haven)

#Cluster the cells

seurat_obj_subset <- readRDS("/data/My_scRNAseq_subset_scaled.rds")
seurat_obj_subset <- FindNeighbors(seurat_obj_subset, dims = 1:15,graph.name = "Dim15")
seurat_obj_subset <- FindClusters(seurat_obj_subset , resolution = 0.3,graph.name = "Dim15")

seurat_obj_subset  <- RunUMAP(seurat_obj_subset , dims = 1:15)
DimPlot(seurat_obj_subset , reduction = "umap",label=TRUE)

#Add metadata


patient_meta <- read_csv("~/Projects/scRNAseq_Omniscope/IMC-F106C-101_patient_groups.csv")

names <- row.names(seurat_obj_subset@meta.data)
seurat_meta <- seurat_obj_subset@meta.data %>%
  mutate(Patient = names) %>%
  mutate(Patient = gsub(".*-", "", Patient)) %>%
  left_join(patient_meta %>%
              mutate(SUBJID = as.character(SUBJID)), by = c("Patient" = "SUBJID"))
rownames(seurat_meta) <- names
seurat_obj_subset@meta.data <- seurat_meta

#Try DoubletFinder

seurat_obj_subset.split <- SplitObject(seurat_obj_subset, split.by = "Patient") 

# loop through samples to find doublets
for (i in 1:length(seurat_obj_subset.split)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  pbmc.sample <- NormalizeData(seurat_obj_subset.split[[i]])
  pbmc.sample <- FindVariableFeatures(pbmc.sample)
  pbmc.sample <- ScaleData(pbmc.sample)
  pbmc.sample <- RunPCA(pbmc.sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- pbmc.sample[["pca"]]@stdev
  sum.stdv <- sum(pbmc.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  pbmc.sample <- RunUMAP(pbmc.sample, dims = 1:min.pc)
  pbmc.sample <- FindNeighbors(object = pbmc.sample, dims = 1:min.pc)              
  pbmc.sample <- FindClusters(object = pbmc.sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(pbmc.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- pbmc.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(pbmc.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  pbmc.sample <- doubletFinder(seu = pbmc.sample, 
                               PCs = 1:min.pc, 
                               pK = optimal.pk,
                               nExp = nExp.poi.adj)
  metadata <- pbmc.sample@meta.data
  colnames(metadata)[15] <- "doublet_finder"
  pbmc.sample@meta.data <- metadata 
  
  # subset and save
  #  pbmc.singlets <- subset(pbmc.sample, doublet_finder == "Singlet")
  seurat_obj_subset.split[[i]] <- pbmc.sample
  #  remove(pbmc.singlets)
}


for (i in 1:length(seurat_obj_subset.split)) {
  metadata <- seurat_obj_subset.split[[i]]@meta.data
  colnames(metadata)[16] <- "doublet_finder"
  seurat_obj_subset.split[[i]]@meta.data <- metadata
}


# converge seurat_obj_subset.split
pbmc.singlets <- merge(x = seurat_obj_subset.split[[1]],
                       y = c(seurat_obj_subset.split[[2]], 
                             seurat_obj_subset.split[[3]],
                             seurat_obj_subset.split[[4]],
                             seurat_obj_subset.split[[5]],
                             seurat_obj_subset.split[[6]],
                             seurat_obj_subset.split[[7]],
                             seurat_obj_subset.split[[8]],
                             seurat_obj_subset.split[[9]], 
                             seurat_obj_subset.split[[10]],
                             seurat_obj_subset.split[[11]],
                             seurat_obj_subset.split[[12]],
                             seurat_obj_subset.split[[13]],
                             seurat_obj_subset.split[[14]],
                             seurat_obj_subset.split[[15]],
                             seurat_obj_subset.split[[16]], 
                             seurat_obj_subset.split[[17]],
                             seurat_obj_subset.split[[18]],
                             seurat_obj_subset.split[[19]],
                             seurat_obj_subset.split[[20]],
                             seurat_obj_subset.split[[21]],
                             seurat_obj_subset.split[[22]],
                             seurat_obj_subset.split[[23]],
                             seurat_obj_subset.split[[24]],
                             seurat_obj_subset.split[[25]],
                             seurat_obj_subset.split[[26]], 
                             seurat_obj_subset.split[[27]],
                             seurat_obj_subset.split[[28]],
                             seurat_obj_subset.split[[29]]),
                       project = "Omniscope scRNAseq")
pbmc.singlets

doublets_metadata <- pbmc.singlets@meta.data[,c("doublet_finder","doublet_finder_res")]

write.csv(doublets_metadata,file="/data/doublet_finder_14March.csv")

#seurat_obj_subset <- readRDS("/data/scRNAseq_subset_scaled.rds")
#seurat_obj_subset <- FindNeighbors(seurat_obj_subset, dims = 1:10,graph.name = "Dim10")
#seurat_obj_subset <- FindClusters(seurat_obj_subset , resolution = 0.3,graph.name = "Dim10")


#seurat_obj_subset <- RunUMAP(seurat_obj_subset, dims = 1:10)
#DimPlot(seurat_obj_subset, reduction = "umap")

doublets_metadata <- read.csv("/data/doublet_finder_14March.csv")
colnames(doublets_metadata) <- c("cellID", "doublet_finder","doublet_finder_res")

names <- row.names(seurat_obj_subset@meta.data)
seurat_meta <- seurat_obj_subset@meta.data %>%
  mutate(cellID = names) %>%
  left_join(doublets_metadata,by="cellID")
rownames(seurat_meta) <- names
seurat_obj_subset@meta.data <- seurat_meta
