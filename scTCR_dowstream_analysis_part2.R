seurat_obj_26patients <- readRDS("/data/Seurat_object_26patients_complete.rds")
combined_TCR <- readRDS("~/Projects/scRNAseq_Omniscope/combined_scTCRseq_26patients_withoutcontig3.rds")
patient_meta <- read_csv("~/Projects/scRNAseq_Omniscope/IMC-F106C-101_patient_groups.csv")
samples_name <- c("10040003","10040008","10050104","10070075","10090003","10090008",
                  "10090014","10090016","10090018","10090029","10090040","10120062",
                  "20010031","20010043","20010066","20010096","DFI_030","KOI_017",
                  "LZI_011","PKI_025","SXI_029","WVY_018","XMI_017","XPI_012",
                  "YCI_008","YZI_003")
patient_meta <- patient_meta %>%
  filter(SUBJID %in% samples_name) %>%
  arrange(match(SUBJID,samples_name))

seurat_metadata <- seurat_obj_26patients@meta.data

seurat_metadata <- seurat_metadata %>%
  select(orig.ident,cellID,Patient,group,BOR,ORIGDIAG,HISTCLAS,Cell_type) %>%
  mutate(Cell_type = recode(Cell_type,"CD4 Naive (1)"="CD4 Naive",
                            "CD4 Naive (2)"="CD4 SCM",
                            "CD8 Naive (1)"="CD8 Naive",
                            "CD8 Naive (2)"="CD8 SCM")) %>%
  unite(Cell_type_Patient,Cell_type,Patient,remove=FALSE) %>%
  mutate(Cell_type2 = case_when(Cell_type %in% c("CD14+ Mono","FCGR3A+ Mono") ~ "Monocytes",
                                Cell_type %in% c("CD4 Naive","CD4 CM","CD4 EM","CD4 SCM") ~ "CD4",
                                Cell_type %in% c("CD8 Naive","CD8 CM","CD8 EM","CD8 SCM", "CD8 CTL") ~ "CD8",
                                Cell_type %in% c("NK") ~ "NK",
                                Cell_type %in% c("B Memory","B Naive") ~ "B cells",
                                Cell_type %in% c("PPBP+ CD14+ Mono","PPBP+ T cells","Platelet") ~ "PPBP+",
                                Cell_type %in% c("HPSC") ~ "HPSC",
                                Cell_type %in% c("MAIT") ~ "MAIT",
                                Cell_type %in% c("Proliferating T cells") ~ "Proliferating T cells",
                                Cell_type %in% c("Mix B and T markers") ~ "Mix B and T markers",
                                Cell_type %in% c("Plasma") ~ "Plasma",
                                Cell_type %in% c("pDC","mDC") ~ "DC")) %>%
seurat_metadata <- seurat_metadata %>%  unite(Cell_type_group,Cell_type,group,remove=FALSE) %>%
  unite(Cell_type2_group,Cell_type2,group,remove=FALSE)

seurat_obj_26patients@meta.data <- seurat_metadata

new_seurat <- combineExpression(
  combined_TCR,
  seurat_obj_26patients,
  cloneCall = "strict",
  chain = "both",
  group.by = "Patient",
  proportion = TRUE,
  filterNA = FALSE,
  addLabel = FALSE
)

newList <- expression2List(new_seurat, split.by = "Cell_type_Patient")

clonalHomeostasis(newList,cloneCall="strict",group.by="Cell_type2")
quantContig(newList,cloneCall="strict",group.by="Cell_type2")

Tcells_List <- subsetContig(newList,name="Cell_type2",variables = c("CD4","CD8","MAIT"))
clonalHomeostasis(Tcells_List,cloneCall="strict",group.by="Cell_type2_group")

saveRDS(Tcells_List,file="~/Projects/scRNAseq_Omniscope/combined_scTCRseq_26patients_CD4_CD8_MAIT.rds")
