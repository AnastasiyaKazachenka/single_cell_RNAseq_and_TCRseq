setwd("~/Projects/scRNAseq_Omniscope/VDJmatch")

library(readr)

patient_meta <- read_csv("~/Projects/scRNAseq_Omniscope/IMC-F106C-101_patient_groups.csv")
patient_meta <- as.data.frame(patient_meta)
rownames(patient_meta) <- patient_meta$SUBJID
patient_meta <- patient_meta[order(patient_meta$group),]

combined_TCR <- readRDS("~/Projects/scRNAseq_Omniscope/combined_scTCRseq_26patients_CD4_CD8_MAIT.rds")

#cd4_YCI_008 <- c("CD4 CM_YCI_008","CD4 EM_YCI_008","CD4 Naive_YCI_008","CD4 SCM_YCI_008")
#cd8_YCI_008 <- c("CD8 CM_YCI_008","CD8 EM_YCI_008","CD8 Naive_YCI_008","CD8 SCM_YCI_008","CD8 CTL_YCI_008")

#CD4

CD4_YCI_008 <- combined_TCR["CD4 CM_YCI_008"][[1]][c("CTgene","CTnt","CTaa")]
CD4_YCI_008 <- rbind(CD4_YCI_008,combined_TCR["CD4 EM_YCI_008"][[1]][c("CTgene","CTnt","CTaa")])
CD4_YCI_008 <- rbind(CD4_YCI_008,combined_TCR["CD4 SCM_YCI_008"][[1]][c("CTgene","CTnt","CTaa")])
CD4_YCI_008 <- rbind(CD4_YCI_008,combined_TCR["CD4 Naive_YCI_008"][[1]][c("CTgene","CTnt","CTaa")])


CD4_TRA <- CD4_YCI_008 %>%
  separate(CTgene, c("TRA_gene","TRB_gene"),sep="_") %>%
  separate(CTnt, c("TRA_nt","TRB_nt"),sep="_") %>%
  separate(CTaa, c("TRA_aa","TRB_aa"),sep="_") %>% 
  select(TRA_gene,TRA_nt,TRA_aa) %>%
#  filter(grepl(";",TRA_gene)==FALSE) %>%
  separate(TRA_gene,c("V","J","C"),sep="\\.") %>%
  rename("CDR3aa" = "TRA_aa",
         "CDR3nt" = "TRA_nt") %>%
  group_by(CDR3aa,CDR3nt,V,J) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  mutate(D=rep(".",nrow(.))) %>%
  mutate(freq=count/sum(count)) %>%
  select(count,freq,CDR3nt,CDR3aa,V,D,J)

CD4_TRB <- CD4_YCI_008 %>%
  separate(CTgene, c("TRA_gene","TRB_gene"),sep="_") %>%
  separate(CTnt, c("TRA_nt","TRB_nt"),sep="_") %>%
  separate(CTaa, c("TRA_aa","TRB_aa"),sep="_") %>% 
  select(TRB_gene,TRB_nt,TRB_aa) %>%
#  filter(grepl(";",TRB_gene)==FALSE) %>%
  separate(TRB_gene,c("V","D","J","C"),sep="\\.") %>%
  rename("CDR3aa" = "TRB_aa",
         "CDR3nt" = "TRB_nt") %>%
  group_by(CDR3aa,CDR3nt,V,J) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  mutate(D=rep(".",nrow(.))) %>%
  mutate(freq=count/sum(count)) %>%
  select(count,freq,CDR3nt,CDR3aa,V,D,J)

write_delim(CD4_TRA,"./CD4_TCR/data/sample_YCI_008_CD4_TRA.txt",delim="\t")
write_delim(CD4_TRB,"./CD4_TCR/data/sample_YCI_008_CD4_TRB.txt",delim="\t")


#CD8

CD8_YCI_008 <- combined_TCR["CD8 CM_YCI_008"][[1]][c("CTgene","CTnt","CTaa")]
CD8_YCI_008 <- rbind(CD8_YCI_008,combined_TCR["CD8 EM_YCI_008"][[1]][c("CTgene","CTnt","CTaa")])
CD8_YCI_008 <- rbind(CD8_YCI_008,combined_TCR["CD8 SCM_YCI_008"][[1]][c("CTgene","CTnt","CTaa")])
CD8_YCI_008 <- rbind(CD8_YCI_008,combined_TCR["CD8 Naive_YCI_008"][[1]][c("CTgene","CTnt","CTaa")])
CD8_YCI_008 <- rbind(CD8_YCI_008,combined_TCR["CD8 CTL_YCI_008"][[1]][c("CTgene","CTnt","CTaa")])


CD8_TRA <- CD8_YCI_008 %>%
  separate(CTgene, c("TRA_gene","TRB_gene"),sep="_") %>%
  separate(CTnt, c("TRA_nt","TRB_nt"),sep="_") %>%
  separate(CTaa, c("TRA_aa","TRB_aa"),sep="_") %>% 
  select(TRA_gene,TRA_nt,TRA_aa) %>%
  #  filter(grepl(";",TRA_gene)==FALSE) %>%
  separate(TRA_gene,c("V","J","C"),sep="\\.") %>%
  rename("CDR3aa" = "TRA_aa",
         "CDR3nt" = "TRA_nt") %>%
  group_by(CDR3aa,CDR3nt,V,J) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  mutate(D=rep(".",nrow(.))) %>%
  mutate(freq=count/sum(count)) %>%
  select(count,freq,CDR3nt,CDR3aa,V,D,J)

CD8_TRB <- CD8_YCI_008 %>%
  separate(CTgene, c("TRA_gene","TRB_gene"),sep="_") %>%
  separate(CTnt, c("TRA_nt","TRB_nt"),sep="_") %>%
  separate(CTaa, c("TRA_aa","TRB_aa"),sep="_") %>% 
  select(TRB_gene,TRB_nt,TRB_aa) %>%
  #  filter(grepl(";",TRB_gene)==FALSE) %>%
  separate(TRB_gene,c("V","D","J","C"),sep="\\.") %>%
  rename("CDR3aa" = "TRB_aa",
         "CDR3nt" = "TRB_nt") %>%
  group_by(CDR3aa,CDR3nt,V,J) %>%
  summarise(count=n()) %>%
  as.data.frame() %>%
  mutate(D=rep(".",nrow(.))) %>%
  mutate(freq=count/sum(count)) %>%
  select(count,freq,CDR3nt,CDR3aa,V,D,J)

write_delim(CD8_TRA,"./CD8_TCR/data/sample_YCI_008_CD8_TRA.txt",delim="\t")
write_delim(CD8_TRB,"./CD8_TCR/data/sample_YCI_008_CD8_TRB.txt",delim="\t")


### Runing VDJmatch

samples_name <- c("10040003","10040008","10050104","10070075","10090003","10090008",
                  "10090014","10090016","10090018","10090029","10090040","10120062",
                  "20010031","20010043","20010066","20010096","DFI_030","KOI_017",
                  "LZI_011","PKI_025","SXI_029","WVY_018","XMI_017","XPI_012",
                  "YCI_008","YZI_003")

for (x in 1:26) {
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRA --min-epi-size 1 --search-scope 0,0,0 --v-match --j-match ./CD4_TCR/data/sample_",samples_name[x],"_CD4_TRA.txt ./CD4_TCR/VDJmatch_results/vdjmatch.000",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRB --min-epi-size 1 --search-scope 0,0,0 --v-match --j-match ./CD4_TCR/data/sample_",samples_name[x],"_CD4_TRB.txt ./CD4_TCR/VDJmatch_results/vdjmatch.000",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRA --min-epi-size 1 --search-scope 2,1,2 --v-match --j-match ./CD4_TCR/data/sample_",samples_name[x],"_CD4_TRA.txt ./CD4_TCR/VDJmatch_results/vdjmatch.212",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRB --min-epi-size 1 --search-scope 2,1,2 --v-match --j-match ./CD4_TCR/data/sample_",samples_name[x],"_CD4_TRB.txt ./CD4_TCR/VDJmatch_results/vdjmatch.212",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRA --min-epi-size 1 --search-scope 0,0,0 --v-match --j-match ./CD8_TCR/data/sample_",samples_name[x],"_CD8_TRA.txt ./CD8_TCR/VDJmatch_results/vdjmatch.000",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRB --min-epi-size 1 --search-scope 0,0,0 --v-match --j-match ./CD8_TCR/data/sample_",samples_name[x],"_CD8_TRB.txt ./CD8_TCR/VDJmatch_results/vdjmatch.000",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRA --min-epi-size 1 --search-scope 2,1,2 --v-match --j-match ./CD8_TCR/data/sample_",samples_name[x],"_CD8_TRA.txt ./CD8_TCR/VDJmatch_results/vdjmatch.212",sep=""))
  system(paste("java -jar /nethome/anastasiya.kazachenk/vdjmatch/vdjmatch-1.3.1/vdjmatch-1.3.1.jar match -S human -R TRB --min-epi-size 1 --search-scope 2,1,2 --v-match --j-match ./CD8_TCR/data/sample_",samples_name[x],"_CD8_TRB.txt ./CD8_TCR/VDJmatch_results/vdjmatch.212",sep=""))
}