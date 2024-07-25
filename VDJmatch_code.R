setwd("~/Projects/scRNAseq_Omniscope/VDJmatch")

library(pheatmap)
library(readr)

patient_meta <- read_csv("~/Projects/scRNAseq_Omniscope/IMC-F106C-101_patient_groups.csv")
patient_meta <- as.data.frame(patient_meta)
rownames(patient_meta) <- patient_meta$SUBJID
patient_meta <- patient_meta[order(patient_meta$group),]

samples_name <- c("10040003","10040008","10050104","10070075","10090003","10090008",
                  "10090014","10090016","10090018","10090029","10090040","10120062",
                  "20010031","20010043","20010066","20010096","DFI_030","KOI_017",
                  "LZI_011","PKI_025","SXI_029","WVY_018","XMI_017","XPI_012",
                  "YCI_008","YZI_003")

combined_TCR <- readRDS("~/Projects/scRNAseq_Omniscope/combined_scTCRseq_26patients.rds")

sample_TRA <- combined_TCR$`YZI_003` %>% 
  as.data.frame() %>%
  select(TCR1,cdr3_aa1,cdr3_nt1) %>%
  filter(grepl(";",TCR1)==FALSE) %>%
  separate(TCR1,c("V","J","C"),sep="\\.") %>%
  rename("cdr3_aa1"="CDR3aa",
         "cdr3_nt1"="CDR3nt") %>%
  group_by(CDR3aa,CDR3nt,V,J) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  mutate(D=rep(".",nrow(.))) %>%
  mutate(freq=count/sum(count)) %>%
  select(count,freq,CDR3nt,CDR3aa,V,D,J)

sample_TRB <- combined_TCR$`YZI_003` %>% 
  as.data.frame() %>%
  select(TCR2,cdr3_aa2,cdr3_nt2) %>%
  filter(grepl(";",TCR2)==FALSE) %>%
  separate(TCR2,c("V","D","J","C"),sep="\\.") %>%
  rename("cdr3_aa2"="CDR3aa",
         "cdr3_nt2"="CDR3nt") %>%
  group_by(CDR3aa,CDR3nt,V,J) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  mutate(D=rep(".",nrow(.))) %>%
  mutate(freq=count/sum(count)) %>%
  select(count,freq,CDR3nt,CDR3aa,V,D,J)

write_delim(sample_TRA,"sample_YZI_003_TRA.txt",delim="\t")
write_delim(sample_TRB,"sample_YZI_003_TRB.txt",delim="\t")



### running vdjmatch

for (x in 1:26) {
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRA --min-epi-size 1 --search-scope 0,0,0 --v-match --j-match sample_",samples_name[x],"_TRA.txt vdjmatch.000",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRB --min-epi-size 1 --search-scope 0,0,0 --v-match --j-match sample_",samples_name[x],"_TRB.txt vdjmatch.000",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRA --min-epi-size 1 --search-scope 2,1,2 --v-match --j-match sample_",samples_name[x],"_TRA.txt vdjmatch.212",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRB --min-epi-size 1 --search-scope 2,1,2 --v-match --j-match sample_",samples_name[x],"_TRB.txt vdjmatch.212",sep=""))
}

### analysing results

filter_10x <- "https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#"

TRA_export <- read_delim(paste("vdjmatch.000.sample_",samples_name[1],"_TRA.txt",sep=""), show_col_types = FALSE)
TRA_export <- TRA_export %>% filter(reference.id != filter_10x) %>%
  unite(antigene.full, antigen.gene, antigen.species, sep="_",remove=FALSE) 
# %>%
#  filter(vdjdb.score >0)
antigenes <- unique(TRA_export$antigene.full)
TRA_antigenes<- data.frame(antigenes , rep(1,length(antigenes)))
names(TRA_antigenes) <- c("antigenes",samples_name[1])


for (x in 2:26) {
  TRA_export <- read_delim(paste("vdjmatch.000.sample_",samples_name[x],"_TRA.txt",sep=""), show_col_types = FALSE)
  TRA_export <- TRA_export %>% filter(reference.id != filter_10x) %>%
    unite(antigene.full, antigen.gene, antigen.species, sep="_",remove=FALSE)
#  %>%
#    filter(vdjdb.score >0)
  antigenes <- unique(TRA_export$antigene.full)
  TRA_antigenes_part <- data.frame(antigenes , rep(1,length(antigenes)))
  names(TRA_antigenes_part) <- c("antigenes",samples_name[x])
  TRA_antigenes <- full_join(TRA_antigenes,TRA_antigenes_part)
}

TRA_antigenes[is.na(TRA_antigenes)] <- 0
row.names(TRA_antigenes) <- TRA_antigenes$antigenes


TRB_export <- read_delim(paste("vdjmatch.000.sample_",samples_name[1],"_TRB.txt",sep=""), show_col_types = FALSE)
TRB_export <- TRB_export %>% filter(reference.id != filter_10x) %>%
  unite(antigene.full, antigen.gene, antigen.species, sep="_",remove=FALSE) 
# %>%
#  filter(vdjdb.score >0)
antigenes <- unique(TRB_export$antigene.full)
TRB_antigenes<- data.frame(antigenes , rep(1,length(antigenes)))
names(TRB_antigenes) <- c("antigenes",samples_name[1])


for (x in 2:26) {
  TRB_export <- read_delim(paste("vdjmatch.000.sample_",samples_name[x],"_TRB.txt",sep=""), show_col_types = FALSE)
  TRB_export <- TRB_export %>% filter(reference.id != filter_10x) %>%
    unite(antigene.full, antigen.gene, antigen.species, sep="_",remove=FALSE) 
#  %>%
#    filter(vdjdb.score >0)
  antigenes <- unique(TRB_export$antigene.full)
  TRB_antigenes_part <- data.frame(antigenes , rep(1,length(antigenes)))
  names(TRB_antigenes_part) <- c("antigenes",samples_name[x])
  TRB_antigenes <- full_join(TRB_antigenes,TRB_antigenes_part)
}

TRB_antigenes[is.na(TRB_antigenes)] <- 0
row.names(TRB_antigenes) <- TRB_antigenes$antigenes

pheatmap(TRB_antigenes[,c(2:27)],scale="none",annotation_col = patient_meta[,c("group","ORIGDIAG")])
pheatmap(TRA_antigenes[,c(2:27)],scale="none",annotation_col = patient_meta[,c("group","ORIGDIAG")])

TRB_antigenes %>%
  select(`10050104`,`10090014`,`10090016`,`10090018`,
  `10090029`,`20010031`,`10040003`,`10040008`,`10070075`,
  `10090003`,`10090008`,`10090040`,`20010096`) %>%
#  filter_all(any_vars(. != 0)) %>%
#  filter_all(any_vars(. != 1)) %>%
  pheatmap(scale="none",
         annotation_col = patient_meta[c("10050104","10090014","10090016","10090018",
                                         "10090029","20010031","10040003",
                                         "10040008","10070075","10090003","10090008","10090040","20010096"),
                                       c("group","ORIGDIAG","HISTCLAS")],cluster_cols = F)


