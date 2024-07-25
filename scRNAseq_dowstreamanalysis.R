library(tidyverse)
library(Seurat)
library(SingleR)
library(celldex)
library(scCATCH)
library(haven)
library(pheatmap)
library(ggpubr)


### ALL Patients

all_clusters_annot <- read.csv("~/Projects/scRNAseq_Omniscope/ALLclusters_annot.csv",row.names="Row.names")
all_clusters_annot$X <- NULL

Tcell_annot <- read.csv("~/Projects/scRNAseq_Omniscope/Tcell_clusters_annot.csv",row.names="Row.names")
Tcell_annot$X <- NULL

full_annot <- merge(all_clusters_annot,Tcell_annot,by=0,all=TRUE)
full_annot <- full_annot %>%
  select(-Patient.y,-BOR.y,-group.y) %>%
  mutate(Cell_type=coalesce(Cell_type,seurat_obj_singlets.active.ident)) %>%
  select(-seurat_obj_singlets.active.ident)
rownames(full_annot) <- full_annot$Row.names
full_annot$Row.names <- NULL
colnames(full_annot) <- c("Patient","BOR","group","Cell_type")

cell_types_numbers <- full_annot %>% group_by(Patient,Cell_type,BOR,group) %>%
  summarize(count=n())

cell_types_numbers %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=count/sum(count)) %>% 
#  filter(group=="Low") %>%
  ggplot(aes(x=Patient,y=freq,fill=Cell_type)) +
  geom_bar(stat="identity") +
#  scale_fill_manual(values=c("#000000","#003B7F","#0064D2","#DC6DFF","#DE2064","#00B097","#8F76FF","#ED83A9","#751135",
#                             "#B3D7FF","#006153","#FFA78F","#00FAD6","#65AEFF","#1584FF","#D6D6CC","#8906BE","#FF6133")) +
  facet_grid(~group, scales="free_x") +
#  ylab("% of T cells",) +
#  xlab("") +
  theme_bw()

cell_types_numbers %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=count/sum(count)) %>% 
  #  filter(group=="Low") %>%
  ggplot(aes(x=Patient,y=count,fill=Cell_type)) +
  geom_bar(stat="identity") +
  #  scale_fill_manual(values=c("#000000","#003B7F","#0064D2","#DC6DFF","#DE2064","#00B097","#8F76FF","#ED83A9","#751135",
  #                             "#B3D7FF","#006153","#FFA78F","#00FAD6","#65AEFF","#1584FF","#D6D6CC","#8906BE","#FF6133")) +
  facet_grid(~group, scales="free_x") +
  #  ylab("% of T cells",) +
  #  xlab("") +
  theme_bw()

                           
cell_type_freq <- cell_types_numbers %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>%
  as.data.frame() %>%
  select(Patient,Cell_type,freq) %>%
  pivot_wider(names_from = Patient,values_from = freq,values_fill = 0) %>%
  as.data.frame()

Patient_annot <- full_annot[!(duplicated(full_annot[,c("Patient","BOR","group")])),c("Patient","BOR","group")]
row.names(Patient_annot) <- Patient_annot$Patient
Patient_annot$Patient <- NULL

row.names(cell_type_freq) <- cell_type_freq$Cell_type
cell_type_freq$Cell_type <- NULL

pheatmap(cell_type_freq,annotation_col = Patient_annot)

mycomparisons <- list(c("High","Low"),c("High","Normal"),c("Low","Normal"))

cell_types_numbers %>%
  filter(group != "Normal") %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>% 
  filter(grepl("CD4",Cell_type)==TRUE | Cell_type=="FOXP3+ Tregs") %>%
  ggplot(aes(x=factor(Cell_type,levels=c("CD4 Naive (1)","CD4 Naive (2)","CD4 Naive (3)","CD4 CM","CD4_act_1","CD4 act_2","FOXP3+ Tregs","PPBP+ CD4")),y=freq,fill=group,color=group)) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  geom_jitter(position=position_dodge(0.7)) +
#  scale_fill_manual(values = alpha(c("#003B7F","#1857F4","#FDBD13"),0.5)) +
#  scale_color_manual(values = c("#003B7F","#1857F4","#FDBD13")) +
  scale_fill_manual(values = alpha(c("#003B7F","#FDBD13"),0.5)) +
  scale_color_manual(values = c("#003B7F","#FDBD13")) +
  stat_compare_means(aes(group = group),label="p.signif") +
#  facet_wrap(~Cell_type,scales="free_x") +
  ylab("% of all cells") +
  xlab("") +
  theme_bw()

cell_types_numbers %>%
  filter(group != "Normal") %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>% 
  filter(grepl("CD8",Cell_type)==TRUE | Cell_type=="MAIT") %>%
  ggplot(aes(x=factor(Cell_type,levels=c("CD8 Naive (0)","CD8 Naive (1)","CD8 Naive (2)","CD8 CM","CD8 EM","MAIT","CD8 CTL_1","CD8 CTL_2")),y=freq,fill=group,color=group)) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  geom_jitter(position=position_dodge(0.7)) +
  scale_fill_manual(values = alpha(c("#003B7F","#FDBD13"),0.5)) +
  scale_color_manual(values = c("#003B7F","#FDBD13")) +
  stat_compare_means(aes(group = group),label="p.signif") +
  #  facet_wrap(~Cell_type,scales="free_x") +
  ylab("% of all cells") +
  xlab("") +
  #  facet_wrap(~Cell_type,scales="free_x") +
  theme_bw()

cell_types_numbers %>%
  filter(group != "Normal") %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>% 
  filter(grepl("CD4",Cell_type)==FALSE & grepl("CD8",Cell_type)==FALSE & Cell_type != "FOXP3+ Tregs" & Cell_type!="MAIT") %>%
  ggplot(aes(x=Cell_type,y=freq,fill=group,color=group)) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  geom_jitter(position=position_dodge(0.7)) +
  scale_fill_manual(values = alpha(c("#003B7F","#FDBD13"),0.5)) +
  scale_color_manual(values = c("#003B7F","#FDBD13")) +
  stat_compare_means(aes(group = group),label="p.signif") +
  #  facet_wrap(~Cell_type,scales="free_x") +
  ylab("% of all cells") +
  xlab("") +
  #  facet_wrap(~Cell_type,scales="free_x") +
  theme_bw()


##### 26 patients


all_clusters_annot <- read.csv("~/Projects/scRNAseq_Omniscope/ALLclusters_annot_26P.csv",row.names="Row.names")
all_clusters_annot$X <- NULL

Tcell_annot <- read.csv("~/Projects/scRNAseq_Omniscope/Tcellclusters_annot_26P.csv",row.names="Row.names")
Tcell_annot$X <- NULL

full_annot <- merge(all_clusters_annot,Tcell_annot,by=0,all=TRUE)
full_annot <- full_annot %>%
  dplyr::select(-Patient.y,-BOR.y,-group.y) %>%
  mutate(Cell_type=coalesce(seurat_26patients_Tcells.active.ident,seurat_obj_26patients.active.ident)) %>%
  dplyr::select(-seurat_obj_26patients.active.ident,-seurat_26patients_Tcells.active.ident)
rownames(full_annot) <- full_annot$Row.names
full_annot$Row.names <- NULL
colnames(full_annot) <- c("Patient","BOR","group","Cell_type")

cell_types_numbers <- full_annot %>% group_by(Patient,Cell_type,BOR,group) %>%
  summarize(count=n())

cell_types_numbers %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=count/sum(count)) %>% 
  #  filter(group=="Low") %>%
  ggplot(aes(x=Patient,y=freq,fill=Cell_type
             )) +
  geom_bar(stat="identity") +
  #  scale_fill_manual(values=c("#000000","#003B7F","#0064D2","#DC6DFF","#DE2064","#00B097","#8F76FF","#ED83A9","#751135",
  #                             "#B3D7FF","#006153","#FFA78F","#00FAD6","#65AEFF","#1584FF","#D6D6CC","#8906BE","#FF6133")) +
  facet_grid(~group, scales="free_x") +
  #  ylab("% of T cells",) +
  #  xlab("") +
  theme_bw()

cell_types_numbers %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=count/sum(count)) %>% 
  #  filter(group=="Low") %>%
  ggplot(aes(x=Patient,y=count,fill=Cell_type)) +
  geom_bar(stat="identity") +
  #  scale_fill_manual(values=c("#000000","#003B7F","#0064D2","#DC6DFF","#DE2064","#00B097","#8F76FF","#ED83A9","#751135",
  #                             "#B3D7FF","#006153","#FFA78F","#00FAD6","#65AEFF","#1584FF","#D6D6CC","#8906BE","#FF6133")) +
  facet_grid(~group, scales="free_x") +
  #  ylab("% of T cells",) +
  #  xlab("") +
  theme_bw()


cell_type_freq <- cell_types_numbers %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>%
  as.data.frame() %>%
  dplyr::select(Patient,Cell_type,freq) %>%
  pivot_wider(names_from = Patient,values_from = freq,values_fill = 0) %>%
  as.data.frame()

Patient_annot <- full_annot[!(duplicated(full_annot[,c("Patient","BOR","group")])),c("Patient","BOR","group")]
row.names(Patient_annot) <- Patient_annot$Patient
Patient_annot$Patient <- NULL

row.names(cell_type_freq) <- cell_type_freq$Cell_type
cell_type_freq$Cell_type <- NULL

pheatmap(cell_type_freq,annotation_col = Patient_annot)

mycomparisons <- list(c("High","Low"),c("High","Normal"),c("Low","Normal"))

cell_types_numbers %>%
  filter(group != "Normal") %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>% 
  filter(grepl("CD4",Cell_type)==TRUE | Cell_type=="Proliferating T cells") %>%
  ggplot(aes(x=factor(Cell_type,levels=c("CD4 Naive (1)","CD4 Naive (2)","CD4 CM","CD4 EM", "Proliferating T cells")),y=freq,fill=group,color=group)) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  geom_jitter(position=position_dodge(0.7)) +
  #  scale_fill_manual(values = alpha(c("#003B7F","#1857F4","#FDBD13"),0.5)) +
  #  scale_color_manual(values = c("#003B7F","#1857F4","#FDBD13")) +
  scale_fill_manual(values = alpha(c("#003B7F","#FDBD13"),0.5)) +
  scale_color_manual(values = c("#003B7F","#FDBD13")) +
  stat_compare_means(aes(group = group),label="p.signif") +
  #  facet_wrap(~Cell_type,scales="free_x") +
  ylab("% of all cells") +
  xlab("") +
  theme_bw()

cell_types_numbers %>%
  filter(group != "Normal") %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>% 
  filter(grepl("CD8",Cell_type)==TRUE | Cell_type=="MAIT") %>%
  ggplot(aes(x=factor(Cell_type,levels=c("CD8 Naive (1)","CD8 Naive (2)","CD8 CM","CD8 EM","CD8 CTL","MAIT")),y=freq,fill=group,color=group)) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  geom_jitter(position=position_dodge(0.7)) +
  scale_fill_manual(values = alpha(c("#003B7F","#FDBD13"),0.5)) +
  scale_color_manual(values = c("#003B7F","#FDBD13")) +
  stat_compare_means(aes(group = group),label="p.signif") +
  #  facet_wrap(~Cell_type,scales="free_x") +
  ylab("% of all cells") +
  xlab("") +
  #  facet_wrap(~Cell_type,scales="free_x") +
  theme_bw()

cell_types_numbers %>%
#  filter(group != "Normal") %>%
  group_by(Patient,BOR,group) %>%
  mutate(freq=100*(count/sum(count))) %>% 
  filter(grepl("CD4",Cell_type)==FALSE & grepl("CD8",Cell_type)==FALSE & Cell_type != "Proliferating T cells" & Cell_type!="MAIT") %>%
  ggplot(aes(x=Cell_type,y=freq,fill=group,color=group)) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  geom_jitter(position=position_dodge(0.7)) +
  scale_fill_manual(values = alpha(c("#003B7F","#0064D2","#FDBD13"),0.5)) +
  scale_color_manual(values = c("#003B7F","#0064D2","#FDBD13")) +
  #  facet_wrap(~Cell_type,scales="free_x") +
  ylab("% of all cells") +
  xlab("") +
  facet_wrap(~Cell_type,scales="free_x") +
  stat_compare_means(aes(group = group),label="p.signif") +
  theme_bw()


cell_types <- unique(cell_types_numbers$Cell_type)


cell_type_plot <- function(cell_type_name) {
  cell_types_numbers %>%
    #  filter(group != "Normal") %>%
    group_by(Patient,BOR,group) %>%
    mutate(freq=100*(count/sum(count))) %>% 
    filter(Cell_type==cell_type_name) %>%
    ggplot(aes(x=group,y=freq,fill=group,color=group)) +
    geom_boxplot(position=position_dodge(0.7),width=0.5) +
    geom_jitter(position=position_dodge(0.7)) +
    scale_fill_manual(values = alpha(c("#003B7F","#0064D2","#FDBD13"),0.5)) +
    scale_color_manual(values = c("#003B7F","#0064D2","#FDBD13")) +
    stat_compare_means(comparisons=list(c("High","Low"),c("Low","Normal"),c("High","Normal")),label="p.format") +
    expand_limits(y=0) +
    ylab("% of all cells") +
    xlab("") +
    ggtitle(cell_type_name) +
    theme(legend.position="none") +
    theme_bw()
  
#  ggsave(my_plot,file=paste0(cell_type_name,".svg"))
}

pdf("Cell_type_percent_plots_26patients.pdf")
for (i in 1:24) {
  print(cell_type_plot(cell_types[i]))
}
dev.off()


cell_type_plot(cell_types[1])


### remove not Melanoma samples

### adding patietns data 
patient_meta <- read_csv("~/Projects/scRNAseq_Omniscope/IMC-F106C-101_patient_groups.csv")


cell_type_plot_v2 <- function(cell_type_name) {
  cell_types_numbers %>%
  rename("SUBJID"="Patient") %>%
  left_join(patient_meta) %>%
  group_by(SUBJID,BOR,group,ORIGDIAG,HISTCLAS) %>%
  mutate(freq=100*(count/sum(count))) %>%
  as.data.frame() %>%
  dplyr::select(SUBJID,ORIGDIAG,HISTCLAS,group,Cell_type,freq) %>%
  filter(ORIGDIAG=="Melanoma" | ORIGDIAG=="Normal") %>%
  filter(Cell_type==cell_type_name) %>%
  ggplot(aes(x=group,y=freq,fill=group,color=group)) +
  geom_boxplot(position=position_dodge(0.7),width=0.5) +
  geom_jitter(position=position_dodge(0.7)) +
  scale_fill_manual(values = alpha(c("#003B7F","#0064D2","#FDBD13"),0.5)) +
  scale_color_manual(values = c("#003B7F","#0064D2","#FDBD13")) +
  stat_compare_means(comparisons=list(c("High","Low"),c("Low","Normal"),c("High","Normal")),label="p.format") +
  expand_limits(y=0) +
  ylab("% of all cells") +
  xlab("") +
  ggtitle(cell_type_name) +
  theme(legend.position="none") +
  theme_bw()
}
  

pdf("Cell_type_percent_plots_26patients_onlyMelanomaandNormal.pdf")
for (i in 1:24) {
  print(cell_type_plot_v2(cell_types[i]))
}
dev.off()



## TCF score high Low

tcf_scores <- read.csv("~/Projects/scRNAseq_Omniscope/TCS_scores.csv")

cell_types_numbers_v2 <- merge(cell_types_numbers,tcf_scores, by.x="Patient",by.y="SUBJID")

cell_types_numbers_v2 %>%
  dplyr::select(-group.y,-Log.reduction) %>%
  group_by(Patient,BOR,log2.TPM) %>%
  mutate(freq=100*(count/sum(count))) %>%
  filter(Cell_type %in% c("B Naive","B Memory","Plasma","Mix B and T markers")) %>%
  ggplot(aes(x=log2.TPM, y=freq,group=Cell_type,color=Cell_type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#003B7F","#FDBD13","#00C3A6","#8F76FF")) +
  ylab("% of all cells") +
  xlab("TCF sccore") +
  theme_bw()
  
cell_types_numbers_v2 %>%
  dplyr::select(-group.y,-Log.reduction) %>%
  group_by(Patient,BOR,log2.TPM) %>%
  mutate(freq=100*(count/sum(count))) %>%
  filter(Cell_type %in% c("CD14+ Mono","FCGR3A+ Mono","mDC","pDC")) %>%
  ggplot(aes(x=log2.TPM, y=freq,group=Cell_type,color=Cell_type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#003B7F","#FDBD13","#00C3A6","#8F76FF")) +
  ylab("% of all cells") +
  xlab("TCF sccore") +
  theme_bw()

cell_types_numbers_v2 %>%
  dplyr::select(-group.y,-Log.reduction) %>%
  group_by(Patient,BOR,log2.TPM) %>%
  mutate(freq=100*(count/sum(count))) %>%
  filter(Cell_type %in% c("PPBP+ CD14+ Mono","PPBP+ T cells","Platelet")) %>%
  ggplot(aes(x=log2.TPM, y=freq,group=Cell_type,color=Cell_type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#003B7F","#FDBD13","#00C3A6","#8F76FF")) +
  ylab("% of all cells") +
  xlab("TCF sccore") +
  theme_bw()

cell_types_numbers_v2 %>%
  dplyr::select(-group.y,-Log.reduction) %>%
  group_by(Patient,BOR,log2.TPM) %>%
  mutate(freq=100*(count/sum(count))) %>%
  filter(Cell_type %in% c("NK","HPSC","Proliferating T cells","MAIT")) %>%
  ggplot(aes(x=log2.TPM, y=freq,group=Cell_type,color=Cell_type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#003B7F","#FDBD13","#00C3A6","#8F76FF")) +
  ylab("% of all cells") +
  xlab("TCF sccore") +
  theme_bw()


cell_types_numbers_v2 %>%
  dplyr::select(-group.y,-Log.reduction) %>%
  group_by(Patient,BOR,log2.TPM) %>%
  mutate(freq=100*(count/sum(count))) %>%
  filter(Cell_type %in% c("CD4 Naive (1)","CD4 Naive (2)","CD4 EM","CD4 CM")) %>%
  ggplot(aes(x=log2.TPM, y=freq,group=Cell_type,color=Cell_type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#003B7F","#FDBD13","#00C3A6","#8F76FF")) +
  ylab("% of all cells") +
  xlab("TCF sccore") +
  theme_bw()

cell_types_numbers_v2 %>%
  dplyr::select(-group.y,-Log.reduction) %>%
  group_by(Patient,BOR,log2.TPM) %>%
  mutate(freq=100*(count/sum(count))) %>%
  filter(Cell_type %in% c("CD8 Naive (1)","CD8 Naive (2)","CD8 EM","CD8 CM","CD8 CTL")) %>%
  ggplot(aes(x=log2.TPM, y=freq,group=Cell_type,color=Cell_type)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("#003B7F","#FDBD13","#00C3A6","#8F76FF","#CE155C")) +
  ylab("% of all cells") +
  xlab("TCF sccore") +
  theme_bw()