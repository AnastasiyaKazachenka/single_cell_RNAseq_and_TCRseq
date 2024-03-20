library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
library(scCATCH)
library(haven)


###load 10x data

data_dir <- "/data/raw_data/results/"
samples <- list.files("/data/raw_data/results/")
sample_files <- paste0(samples,"/GEX/filtered_feature_bc_matrix")

data_10x <- sapply(sample_files, function(i) {
  print(paste0("Reading in ",i))
  d10x <- Read10X(file.path(data_dir, i)) # matrix/features/barcodes inside folder named by sample ID
  
  sample <- gsub("_sample.*", "", i)
  # Get colnames as cellID-sampleID (Edit code as necessary)
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),
                                          split = "-"), "[[", 1L), sample, sep = "-")
  d10x
})


all_data <- do.call("cbind", data_10x) # Need to bind all of the matrices from each sample together
seurat_obj <- CreateSeuratObject(all_data, project = "t_cell_fitness", min.cells = 3, min.features = 3) # Generate a seurat object
rm(all_data, data_10x)

##pre-processing

# Calculating the percentage of mitochondrial genes within each cell
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Calculating the percentage of ribosomal genes per cell
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^rps|^rpl")

# Visualize QC metrics

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0)

FeatureScatter(seurat_obj,feature1 = "nCount_RNA", feature2="nFeature_RNA") & geom_hline(yintercept=600)

VlnPlot(seurat_obj, features = c("nFeature_RNA"), pt.size =0) & geom_hline(yintercept=600)

seurat_obj[["QC"]] <- ifelse(seurat_obj@meta.data$nFeature_RNA>700 & seurat_obj@meta.data$nFeature_RNA<4000,"Pass","Low/High Feature")
seurat_obj[["QC"]] <- ifelse(seurat_obj@meta.data$percent.mt>10 & seurat_obj@meta.data$QC == "Pass","High MT",seurat_obj@meta.data$QC )

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0,group.by = "QC")

seurat_obj_subset <- subset(seurat_obj, subset = nFeature_RNA > 700 & nFeature_RNA < 4000 & percent.mt < 10)


### Normalize data
seurat_obj_subset <- NormalizeData(seurat_obj_subset, normalization.method = "LogNormalize", scale.factor = 10000)

### Selection of highly variable data

seurat_obj_subset <- FindVariableFeatures(seurat_obj_subset, selection.method = "vst", nfeatures = 2000)

head(VariableFeatures(seurat_obj_subset), 10)
VariableFeaturePlot(seurat_obj_subset)
LabelPoints(VariableFeaturePlot(seurat_obj_subset), points = head(VariableFeatures(seurat_obj_subset), 10), repel = T)

### Scaling data

seurat_obj_subset <- ScaleData(seurat_obj_subset, features = rownames(x = seurat_obj_subset), vars.to.regress = c("percent.mt", "nCount_RNA", "percent.ribo"))



### Dimentional reduction

seurat_obj_subset <- RunPCA(seurat_obj_subset, features = VariableFeatures(object = seurat_obj_subset))

ElbowPlot(seurat_obj_subset)

saveRDS(seurat_obj_subset, file = "/data/My_scRNAseq_subset_scaled.rds")
