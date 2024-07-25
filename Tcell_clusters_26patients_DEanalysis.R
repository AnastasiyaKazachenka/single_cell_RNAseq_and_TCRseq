library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
library(scCATCH)
library(haven)
library(scRepertoire)
library(org.Hs.eg.db)
library(ReactomePA)
library(openxlsx)
library(biomaRt)
library(pheatmap)
library(ggpubr)
library(ggfortify)
library(clusterProfiler)

seurat_26patients_Tcells <- readRDS("/data/Seurat_object_26patients_Tcells.rds")

seurat_idents <- seurat_26patients_Tcells@active.ident %>% as.data.frame()
colnames(seurat_idents) <- c("Seurat_idents")
seurat_metadata <- seurat_26patients_Tcells@meta.data

seurat_metadata <- merge(seurat_metadata,seurat_idents,by=0)

seurat_metadata <- seurat_metadata %>%
  select(Row.names,orig.ident,nCount_RNA,nFeature_RNA,percent.mt,percent.ribo,QC,cellID,
         Patient,group,BOR,ORIGDIAG,HISTCLAS,barcode,CTgene,CTnt,CTaa,CTstrict,Frequency,
         cloneType,Seurat_idents)

seurat_metadata <- seurat_metadata %>%
  unite(Seurat_idents_group,Seurat_idents,group,remove=FALSE)

row.names(seurat_metadata) <- seurat_metadata$Row.names
seurat_metadata$Row.names <- NULL

seurat_26patients_Tcells@meta.data <- seurat_metadata

Idents(seurat_26patients_Tcells) <- "Seurat_idents_group"
cd4_naive_1_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD4 Naive (1)_High", ident.2= "CD4 Naive (1)_Low")
cd4_naive_2_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD4 Naive (2)_High", ident.2= "CD4 Naive (2)_Low")

cd4_EM_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD4 EM_High", ident.2= "CD4 EM_Low")
cd4_CM_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD4 CM_High", ident.2= "CD4 CM_Low")

cd8_naive_1_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD8 Naive (1)_High", ident.2= "CD8 Naive (1)_Low")
cd8_naive_2_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD8 Naive (2)_High", ident.2= "CD8 Naive (2)_Low")

cd8_CM_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD8 CM_High", ident.2= "CD8 CM_Low")
cd8_EM_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD8 EM_High", ident.2= "CD8 EM_Low")

cd8_CTL_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD8 CTL_High", ident.2= "CD8 CTL_Low")
MAIT_DEgenes <- FindMarkers(seurat_26patients_Tcells, ident.1 = "MAIT_High", ident.2= "MAIT_Low")

save.image("/data/Tcell_clusters_26patients_DEanalysis.RData")

load("/data/Tcell_clusters_26patients_DEanalysis.RData")

levels(seurat_26patients_Tcells) <- c("CD4 Naive (1)_High","CD4 Naive (1)_Low","CD4 Naive (1)_Normal",
                                      "CD4 Naive (2)_High","CD4 Naive (2)_Low","CD4 Naive (2)_Normal",
                                      "CD4 CM_High","CD4 CM_Low","CD4 CM_Normal",
                                      "CD4 EM_High","CD4 EM_Low","CD4 EM_Normal",
                                      "CD8 Naive (1)_High","CD8 Naive (1)_Low","CD8 Naive (1)_Normal",
                                      "CD8 Naive (2)_High","CD8 Naive (2)_Low","CD8 Naive (2)_Normal",
                                      "CD8 CM_High","CD8 CM_Low","CD8 CM_Normal",
                                      "CD8 EM_High","CD8 EM_Low","CD8 EM_Normal",
                                      "CD8 CTL_High","CD8 CTL_Low","CD8 CTL_Normal",
                                      "MAIT_High","MAIT_Low","MAIT_Normal",
                                      "Proliferating T cells_High","Proliferating T cells_Low","Proliferating T cells_Normal")

DotPlot(seurat_26patients_Tcells,features=c("CD28","TESPA1","GPR183","IFI27","TCF7","IL7R"),cols="RdYlBu",idents = c("CD4 Naive (1)_Low","CD4 Naive (1)_High","CD4 Naive (2)_Low","CD4 Naive (2)_High","CD4 EM_Low","CD4 EM_High","CD4 CM_Low","CD4 CM_High"))
DotPlot(seurat_26patients_Tcells,features=c("CD28","TESPA1","GPR183","IFI27","TCF7","IL7R"),cols="RdYlBu",idents = c("CD8 Naive (1)_Low","CD8 Naive (1)_High","CD8 Naive (2)_Low","CD8 Naive (2)_High","CD8 EM_Low","CD8 EM_High","CD8 CM_Low","CD8 CM_High","CD8 CTL_Low","CD8 CTL_High"))

#### Pathway enrichment analysis
#get ensembl ids
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_ensembl_id <-  getBM(attributes = c('ensembl_gene_id_version', 'external_gene_name', 'hgnc_symbol'),
                          mart = ensembl)
gene_ensembl_id <- gene_ensembl_id %>% mutate(ensembl_gene_id_version = gsub("\\.[0-9]*$", "", ensembl_gene_id_version))



#Enrichment analysis function
enrichment <- function(genes, analysis_type){
  
  gene_list <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") # first the gene names must be converted to entrez IDs
  
  if(analysis_type == "GO"){
    GO <- enrichGO(gene = gene_list$ENTREZID, # This determines signficiantly enriched pathways 
                   OrgDb = "org.Hs.eg.db", 
                   pvalueCutoff=0.05, 
                   qvalueCutoff=0.05, ont="BP")
    p_go <- dotplot(GO, showCategory = 5) + ggtitle("GO") # View enriched pathways as a dotplot
    print(p_go)
    
  } else if (analysis_type == "KEGG"){
    
    KEGG <- enrichKEGG(gene = gene_list$ENTREZID, 
                       organism = "hsa", 
                       pvalueCutoff=0.05,
                       qvalueCutoff=0.05)
    p_kegg <- dotplot(KEGG, showCategory = 5) + ggtitle("KEGG")
    print(p_kegg)
    
  } else {
    
    reactome <- enrichPathway(gene_list$ENTREZID, 
                              pvalueCutoff = 0.05, 
                              pAdjustMethod = 'BH')
    p_reactome <- dotplot(reactome, showCategory = 5) + ggtitle("Reactome")
    print(p_reactome)
  }
  
}

##heatmap

genes_toplot <- cd4_naive_2_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol)

DoHeatmap(subset(seurat_26patients_Tcells,idents=c("CD4 Naive (2)_High","CD4 Naive (2)_Low","CD4 Naive (1)_High","CD4 Naive (1)_Low"),downsample=1000),features=genes_toplot$hgnc_symbol) + scale_fill_gradientn(colors = c("blue","white","red"))

#upregulated genes DE_Cd4 Naive 2
up_genes_DE_cd4naive2 <- cd4_naive_2_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(up_genes_DE_cd4naive2$hgnc_symbol, "REACTOME")
enrichment(up_genes_DE_cd4naive2$hgnc_symbol, "GO")
enrichment(up_genes_DE_cd4naive2$hgnc_symbol, "KEGG")

DotPlot(seurat_26patients_Tcells,features=up_genes_DE_cd4naive2$hgnc_symbol,idents = c("CD4 Naive (2)_Low","CD4 Naive (2)_High"))

down_genes_DE_cd4naive2 <- cd4_naive_2_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC<0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(down_genes_DE_cd4naive2$hgnc_symbol, "REACTOME")
enrichment(down_genes_DE_cd4naive2$hgnc_symbol, "GO")
enrichment(down_genes_DE_cd4naive2$hgnc_symbol, "KEGG")


all_genes_DE_cd4naive2 <- cd4_naive_2_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(all_genes_DE_cd4naive2$hgnc_symbol, "REACTOME")
enrichment(all_genes_DE_cd4naive2$hgnc_symbol, "GO")
enrichment(all_genes_DE_cd4naive2$hgnc_symbol, "KEGG")

#upregulated genes DE_Cd4 Naive 1
up_genes_DE_cd4naive1 <- cd4_naive_1_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(up_genes_DE_cd4naive1$hgnc_symbol, "REACTOME")
enrichment(up_genes_DE_cd4naive1$hgnc_symbol, "GO")
enrichment(up_genes_DE_cd4naive1$hgnc_symbol, "KEGG")

down_genes_DE_cd4naive1 <- cd4_naive_1_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC<0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(down_genes_DE_cd4naive1$hgnc_symbol, "REACTOME")
enrichment(down_genes_DE_cd4naive1$hgnc_symbol, "GO")
enrichment(down_genes_DE_cd4naive1$hgnc_symbol, "KEGG")


all_genes_DE_cd4naive1 <- cd4_naive_1_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(all_genes_DE_cd4naive1$hgnc_symbol, "REACTOME")
enrichment(all_genes_DE_cd4naive1$hgnc_symbol, "GO")
enrichment(all_genes_DE_cd4naive1$hgnc_symbol, "KEGG")

#upregulated genes cd4_CM_DEgenes
up_genes_DE_cd4_CM <- cd4_CM_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(up_genes_DE_cd4_CM$hgnc_symbol, "REACTOME")
enrichment(up_genes_DE_cd4_CM$hgnc_symbol, "GO")
enrichment(up_genes_DE_cd4_CM$hgnc_symbol, "KEGG")

down_genes_DE_cd4_CM <- cd4_CM_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC<0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(down_genes_DE_cd4_CM$hgnc_symbol, "REACTOME")
enrichment(down_genes_DE_cd4_CM$hgnc_symbol, "GO")
enrichment(down_genes_DE_cd4_CM$hgnc_symbol, "KEGG")


#upregulated genes cd4_EM_DEgenes
up_genes_DE_cd4_EM <- cd4_EM_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(up_genes_DE_cd4_EM$hgnc_symbol, "REACTOME")
enrichment(up_genes_DE_cd4_EM$hgnc_symbol, "GO")
enrichment(up_genes_DE_cd4_EM$hgnc_symbol, "KEGG")

down_genes_DE_cd4_EM <- cd4_EM_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC<0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(down_genes_DE_cd4_EM$hgnc_symbol, "REACTOME")
enrichment(down_genes_DE_cd4_EM$hgnc_symbol, "GO")
enrichment(down_genes_DE_cd4_EM$hgnc_symbol, "KEGG")


#upregulated genes cd8_naive2_DEgenes
up_genes_DE_cd8_naive2 <- cd8_naive_2_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(up_genes_DE_cd8_naive2$hgnc_symbol, "REACTOME")
enrichment(up_genes_DE_cd8_naive2$hgnc_symbol, "GO")
enrichment(up_genes_DE_cd8_naive2$hgnc_symbol, "KEGG")

down_genes_DE_cd8_naive2 <- cd8_naive_2_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC<0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(down_genes_DE_cd8_naive2$hgnc_symbol, "REACTOME")
enrichment(down_genes_DE_cd8_naive2$hgnc_symbol, "GO")
enrichment(down_genes_DE_cd8_naive2$hgnc_symbol, "KEGG")

#upregulated genes cd8_CTL_DEgenes
up_genes_DE_cd8_CTL <- cd8_CTL_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(up_genes_DE_cd8_CTL$hgnc_symbol, "REACTOME")
enrichment(up_genes_DE_cd8_CTL$hgnc_symbol, "GO")
enrichment(up_genes_DE_cd8_CTL$hgnc_symbol, "KEGG")

down_genes_DE_cd8_CTL <- cd8_CTL_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC<0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(down_genes_DE_cd8_CTL$hgnc_symbol, "REACTOME")
enrichment(down_genes_DE_cd8_CTL$hgnc_symbol, "GO")
enrichment(down_genes_DE_cd8_CTL$hgnc_symbol, "KEGG")



#upregulated genes cd8_EM_DEgenes
up_genes_DE_cd8_EM <- cd8_EM_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC>0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(up_genes_DE_cd8_EM$hgnc_symbol, "REACTOME")
enrichment(up_genes_DE_cd8_EM$hgnc_symbol, "GO")
enrichment(up_genes_DE_cd8_EM$hgnc_symbol, "KEGG")

down_genes_DE_cd8_EM <- cd8_EM_DEgenes %>%
  rownames_to_column(var="hgnc_symbol") %>%
  filter(p_val_adj < 0.05,
         avg_log2FC<0) %>%
  arrange(p_val_adj) %>%
  dplyr::select(hgnc_symbol, avg_log2FC, p_val_adj)

enrichment(down_genes_DE_cd8_EM$hgnc_symbol, "REACTOME")
enrichment(down_genes_DE_cd8_EM$hgnc_symbol, "GO")
enrichment(down_genes_DE_cd8_EM$hgnc_symbol, "KEGG")


#### getting signif

cd4_naive_1_DEgenes <- cd4_naive_1_DEgenes %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")
cd4_naive_2_DEgenes <- cd4_naive_2_DEgenes %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")


cd4_EM_DEgenes <- cd4_EM_DEgenes %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")
cd4_CM_DEgenes <- cd4_CM_DEgenes %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")

cd8_naive_1_DEgenes <- cd8_naive_1_DEgenes  %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")
cd8_naive_2_DEgenes <- cd8_naive_2_DEgenes  %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")

cd8_CM_DEgenes <- cd8_CM_DEgenes  %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")
cd8_EM_DEgenes <- cd8_EM_DEgenes  %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")

cd8_CTL_DEgenes <- cd8_CTL_DEgenes  %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")
MAIT_DEgenes <- MAIT_DEgenes  %>%
  filter(p_val_adj < 0.05) %>%
  rownames_to_column(var="hgnc_symbol")


#### getting TCF signature

DE_res <- list(cd4_naive_1_DEgenes,cd4_naive_2_DEgenes,cd4_CM_DEgenes,
               cd4_EM_DEgenes,cd8_naive_1_DEgenes,cd8_naive_2_DEgenes,
               cd8_CM_DEgenes,cd8_EM_DEgenes,cd8_CTL_DEgenes,MAIT_DEgenes)

TCF_sig_names <- c("CD4_Naive_1","CD4_Naive_2","CD4_CM","CD4_EM",
                   "CD8_Naive_1","CD8_Naive_2","CD8_CM","CD8_EM",
                   "CD8_CTL","MAIT")

TCF_sig <- list()
for (x in 1:10) {
  assign(TCF_sig_names[x],DE_res[[x]][c("TESPA1","CD28","GPR183","TCF7","IL7R","IFI27","CD2"),])
  TCF_sig <- append(TCF_sig,list(get(TCF_sig_names[x])))
}

names(TCF_sig) <- TCF_sig_names


##Vulcano plot

deGenes <- cd4_EM_DEgenes %>%
    filter(p_val_adj < 0.05) %>%
    rownames_to_column(var="Gene.stable.ID") %>%
    dplyr::select(Gene.stable.ID)

geneUniverse <- cd4_EM_DEgenes %>%
  rownames_to_column(var="Gene.stable.ID") %>%
  dplyr::select(Gene.stable.ID)

EnhancedVolcano::EnhancedVolcano(cd4_EM_DEgenes,lab = rownames(cd4_EM_DEgenes),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0,
#                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)

EnhancedVolcano::EnhancedVolcano(cd4_CM_DEgenes,lab = rownames(cd4_CM_DEgenes),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0,
#                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)

EnhancedVolcano::EnhancedVolcano(cd4_naive_1_DEgenes,lab = rownames(cd4_naive_1_DEgenes),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 2,
#                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)

EnhancedVolcano::EnhancedVolcano(cd4_naive_2_DEgenes,lab = rownames(cd4_naive_2_DEgenes),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0,
                                 #                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)

EnhancedVolcano::EnhancedVolcano(cd8_naive_2_DEgenes,lab = rownames(cd8_naive_2_DEgenes),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0,
                                 #                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)


EnhancedVolcano::EnhancedVolcano(cd8_CTL_DEgenes,lab = rownames(cd8_CTL_DEgenes),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0,
                                 #                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)


EnhancedVolcano::EnhancedVolcano(cd8_EM_DEgenes,lab = rownames(cd8_EM_DEgenes),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0,
                                 #                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)

##### GO genes

gene_list <- bitr(up_genes_DE_cd4naive2$hgnc_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") # first the gene names must be converted to entrez IDs

GO <- enrichGO(gene = gene_list$ENTREZID, # This determines signficiantly enriched pathways 
                 OrgDb = "org.Hs.eg.db",
                 pvalueCutoff=0.05, 
                 qvalueCutoff=0.05, ont="BP",
               readable=TRUE)

GO_PW <- head(GO,n=5) %>% as.data.frame() %>% select(Description,geneID)

all_de_list <- list(up_genes_DE_cd4naive2$hgnc_symbol,up_genes_DE_cd4naive1$hgnc_symbol,up_genes_DE_cd4_CM$hgnc_symbol,up_genes_DE_cd4_EM$hgnc_symbol,
             up_genes_DE_cd8_naive2$hgnc_symbol,up_genes_DE_cd8_CTL$hgnc_symbol,up_genes_DE_cd4_EM$hgnc_symbol,down_genes_DE_cd4naive2$hgnc_symbol,down_genes_DE_cd4naive1$hgnc_symbol,down_genes_DE_cd4_CM$hgnc_symbol,down_genes_DE_cd4_EM$hgnc_symbol,
             down_genes_DE_cd8_naive2$hgnc_symbol,down_genes_DE_cd8_CTL$hgnc_symbol,down_genes_DE_cd4_EM$hgnc_symbol)


for (x in (2:length(all_de_list))) {
  gene_list <- bitr(all_de_list[[x]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
  GO <- enrichGO(gene = gene_list$ENTREZID, 
                 OrgDb = "org.Hs.eg.db",
                 pvalueCutoff=0.05, 
                 qvalueCutoff=0.05, ont="BP",
                 readable=TRUE)
  GO_PW_tmp <- head(GO,n=5) %>% as.data.frame() %>% select(Description,geneID)
  GO_PW <- rbind(GO_PW,GO_PW_tmp)
  
}

GO_PW <- GO_PW %>% distinct()


## same for KEGG

gene_list <- bitr(up_genes_DE_cd4naive2$hgnc_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") # first the gene names must be converted to entrez IDs

KEGG <- enrichKEGG(gene = gene_list$ENTREZID, 
                   organism = "hsa", 
                   pvalueCutoff=0.05,
                   qvalueCutoff=0.05)
KEGG <-setReadable(KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
KEGG_PW <- head(KEGG,n=5) %>% as.data.frame() %>% select(Description,geneID)

for (x in (2:length(all_de_list))) {
  gene_list <- bitr(all_de_list[[x]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") 
  KEGG <- enrichGO(gene = gene_list$ENTREZID, 
                 OrgDb = "org.Hs.eg.db",
                 pvalueCutoff=0.05, 
                 qvalueCutoff=0.05, ont="BP",
                 readable=TRUE)
  KEGG <-setReadable(KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  KEGG_PW_tmp <- head(KEGG,n=5) %>% as.data.frame() %>% select(Description,geneID)
  KEGG_PW <- rbind(KEGG_PW,KEGG_PW_tmp)
  
}

KEGG_PW <- KEGG_PW %>% distinct()

write.csv(GO_PW,file="~/Projects/scRNAseq_Omniscope/GO_PW_genes_annot.csv")
write.csv(KEGG_PW,file="~/Projects/scRNAseq_Omniscope/KEGG_PW_genes_annot.csv")



##### Surat plots
VlnPlot(seurat_26patients_Tcells,features = c("CD2"),pt.size = 0,idents = c("CD8 Naive (1)_Low","CD8 Naive (1)_High","CD8 Naive (2)_Low","CD8 Naive (2)_High","CD8 EM_Low","CD8 EM_High","CD8 CM_Low","CD8 CM_High","CD8 CTL_Low","CD8 CTL_High")) + NoLegend()
DotPlot(seurat_26patients_Tcells,features = c("CD2"),idents = c("CD4 Naive (1)_Low","CD4 Naive (1)_High","CD4 Naive (2)_Low","CD4 Naive (2)_High","CD4 EM_Low","CD4 EM_High","CD4 CM_Low","CD4 CM_High","CD8 Naive (1)_Low","CD8 Naive (1)_High","CD8 Naive (2)_Low","CD8 Naive (2)_High","CD8 EM_Low","CD8 EM_High","CD8 CM_Low","CD8 CM_High","CD8 CTL_Low","CD8 CTL_High")) + NoLegend()

FeaturePlot(seurat_26patients_Tcells,features = c("CD2"),pt.size = 0)

## DE genes plots CD4
up_genes <- DEgenes_in_allCD4clsuters %>%
  filter(up_down=="up") %>%
  select(gene)
DotPlot(seurat_26patients_Tcells,features = up_genes$gene,
        idents = c("CD4 Naive (1)_Low","CD4 Naive (1)_High",
                   "CD4 Naive (2)_Low","CD4 Naive (2)_High",
                   "CD4 EM_Low","CD4 EM_High","CD4 CM_Low","CD4 CM_High")) +
  RotatedAxis()

down_genes <- DEgenes_in_allCD4clsuters %>%
  filter(up_down=="down") %>%
  select(gene)
DotPlot(seurat_26patients_Tcells,features = down_genes$gene,
        idents = c("CD4 Naive (1)_Low","CD4 Naive (1)_High",
                   "CD4 Naive (2)_Low","CD4 Naive (2)_High",
                   "CD4 EM_Low","CD4 EM_High","CD4 CM_Low","CD4 CM_High")) +
  RotatedAxis()


DotPlot(seurat_26patients_Tcells,features = c("FOSB","FOS","JUN","JUNB"),
        idents = c("CD4 Naive (1)_Low","CD4 Naive (1)_High",
                   "CD4 Naive (2)_Low","CD4 Naive (2)_High",
                   "CD4 EM_Low","CD4 EM_High","CD4 CM_Low","CD4 CM_High")) +
  RotatedAxis()

VlnPlot(seurat_26patients_Tcells,features = c("FOSB","JUN"),pt.size = 0,
        idents = c("CD4 Naive (1)_Low","CD4 Naive (1)_High",
                   "CD4 Naive (2)_Low","CD4 Naive (2)_High",
                   "CD4 EM_Low","CD4 EM_High","CD4 CM_Low","CD4 CM_High")) +
  RotatedAxis()


VlnPlot(seurat_26patients_Tcells,features = c("SOCS1","SOCS2","CISH"),pt.size = 0,
        idents = c("CD4 Naive (1)_Low","CD4 Naive (1)_High",
                   "CD4 Naive (2)_Low","CD4 Naive (2)_High",
                   "CD4 EM_Low","CD4 EM_High","CD4 CM_Low","CD4 CM_High")) +
  RotatedAxis()



## DE genes all CD8 clusters

CD8_DE_results <- c("cd8_naive_1_DEgenes","cd8_naive_2_DEgenes","cd8_EM_DEgenes","cd8_CTL_DEgenes","cd8_CM_DEgenes")

cd8_genes_up <- c()

for (x in 1:5) {
  genes_up <- get(CD8_DE_results[x]) %>%
    filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    rownames_to_column(var="genes") %>%
    select(genes)
  cd8_genes_up <- c(cd8_genes_up,genes_up$genes)
}

cd8_genes_up <- table(cd8_genes_up) %>% as.data.frame()
cd8_genes_up$up_down <- rep("up",nrow(cd8_genes_up))
colnames(cd8_genes_up) <- c("genes","Freq","up_down")



cd8_genes_down <- c()

for (x in 1:5) {
  genes_down <- get(CD8_DE_results[x]) %>%
    filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
    rownames_to_column(var="genes") %>%
    select(genes)
  cd8_genes_down <- c(cd8_genes_down,genes_down$genes)
}

cd8_genes_down <- table(cd8_genes_down) %>% as.data.frame()
cd8_genes_down$up_down <- rep("down",nrow(cd8_genes_down))
colnames(cd8_genes_down) <- c("genes","Freq","up_down")



DE_genes_inall_CD8clusters <- rbind(cd8_genes_up,cd8_genes_down) %>%
  filter(Freq >2)



##### FInd markers all CD4 High vs Low

seurat_26patients_Tcells <- RenameIdents(object=seurat_26patients_Tcells,
                                         `CD4 Naive (1)_Low` = "CD4_Low",
                                         `CD4 Naive (2)_Low` = "CD4_Low",
                                         `CD4 CM_Low` = "CD4_Low",
                                         `CD4 EM_Low` = "CD4_Low",
                                         `CD4 Naive (1)_High` = "CD4_High",
                                         `CD4 Naive (2)_High` = "CD4_High",
                                         `CD4 CM_High` = "CD4_High",
                                         `CD4 EM_High` = "CD4_High")


cd4_combined_HvsL <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD4_High", ident.2= "CD4_Low")

EnhancedVolcano::EnhancedVolcano(cd4_combined_HvsL,lab = rownames(cd4_combined_HvsL),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0.585,
                                 #                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)

cd4_combined_HvsL %>%
  filter(p_val_adj<0.05) %>%
  filter(avg_log2FC>0.585 | avg_log2FC < -0.585) %>%
  View()

seurat_26patients_Tcells <- RenameIdents(object=seurat_26patients_Tcells,
                                         `CD8 Naive (1)_Low` = "CD8_Low",
                                         `CD8 Naive (2)_Low` = "CD8_Low",
                                         `CD8 CM_Low` = "CD8_Low",
                                         `CD8 EM_Low` = "CD8_Low",
                                         `CD8 CTL_Low` = "CD8_Low",
                                         `CD8 Naive (1)_High` = "CD8_High",
                                         `CD8 Naive (2)_High` = "CD8_High",
                                         `CD8 CM_High` = "CD8_High",
                                         `CD8 EM_High` = "CD8_High",
                                         `CD8 CTL_High` = "CD8_High")

cd8_combined_HvsL <- FindMarkers(seurat_26patients_Tcells, ident.1 = "CD8_High", ident.2= "CD8_Low")

EnhancedVolcano::EnhancedVolcano(cd8_combined_HvsL,lab = rownames(cd8_combined_HvsL),
                                 x="avg_log2FC",y="p_val_adj",
                                 pCutoff = 0.05,FCcutoff = 0.585,
                                 #                                 selectLab = c("TESPA1"),
                                 shape=20,pointSize = 1)

save(cd4_combined_HvsL,cd8_combined_HvsL,file="/data/Tcell_clusters_26patients_DEanalysis_combined_cd4_cd8.RData")
