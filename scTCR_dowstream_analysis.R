library(tidyverse)
library(Seurat)
library(scRepertoire)
library(SingleR)
library(celldex)
library(scCATCH)
library(haven)
library(circlize)
library(scales)

### create combinedTCR object
data_dir <- "/data/raw_data/results/"
samples <- list.files("/data/raw_data/results/")
sample_files <- paste0(samples,"/TCR/filtered_contig_annotations.csv")

samples_name <- gsub("_sample","",samples)
#samples_name <- paste("s_",samples_name,sep="")

for (x in 1:length(samples_name)) {
  my_file <- read_csv(file.path(data_dir,sample_files[x]))
  contig3 <- my_file %>%
    filter(grepl("_contig_3",contig_id)==T)
  my_file <- my_file %>%
    filter(!(barcode %in% contig3$barcode))
  assign(samples_name[x], my_file)
}

contig_list <- list(`10040003`,`10040008`,`10050104`,`10070075`,`10090003`,`10090008`,
                    `10090014`,`10090016`,`10090018`,`10090029`,`10090040`,`10120062`,
                    `20010031`,`20010043`,`20010066`,`20010096`,DFI_030,KOI_017,
                    LZI_011,PKI_025,SXI_029,WVY_018,XMI_017,XPI_012,YCI_008,YZI_003)

samples_name <- c("10040003","10040008","10050104","10070075","10090003","10090008",
                  "10090014","10090016","10090018","10090029","10090040","10120062",
                  "20010031","20010043","20010066","20010096","DFI_030","KOI_017",
                  "LZI_011","PKI_025","SXI_029","WVY_018","XMI_017","XPI_012",
                  "YCI_008","YZI_003")

patient_meta <- read_csv("~/Projects/scRNAseq_Omniscope/IMC-F106C-101_patient_groups.csv")

patient_meta <- patient_meta %>%
  filter(SUBJID %in% samples_name) %>%
  arrange(match(SUBJID,samples_name))

combined_TCR <- combineTCR(contig_list,
                       samples=samples_name,
                       removeNA=TRUE)

combined_TCR <- addVariable(combined_TCR,
                            variable.name = "group",
                            variables =patient_meta$group)

combined_TCR <- addVariable(combined_TCR,
                            variable.name = "BOR",
                            variables =patient_meta$BOR)

combined_TCR <- addVariable(combined_TCR,
                            variable.name = "ORIGDIAG",
                            variables =patient_meta$ORIGDIAG)

combined_TCR <- addVariable(combined_TCR,
                            variable.name = "HISTCLAS",
                            variables =patient_meta$HISTCLAS)

combined_TCR <- addVariable(combined_TCR,
                            variable.name = "SUBJID",
                            variables =patient_meta$SUBJID)


saveRDS(combined_TCR,file="~/Projects/scRNAseq_Omniscope/combined_scTCRseq_26patients_withoutcontig3.rds")

### diversity indexes

diversity_indexes <- clonalDiversity(combined_TCR,cloneCall="strict",chain="both",exportTable = TRUE)
diversity_indexes <- merge(patient_meta,diversity_indexes,by.x="SUBJID",by.y="Group")

write.csv(diversity_indexes,file="~/Projects/scRNAseq_Omniscope/Diversity_indexes_scTCRseq_26patients_woithoutcontig3.csv")
