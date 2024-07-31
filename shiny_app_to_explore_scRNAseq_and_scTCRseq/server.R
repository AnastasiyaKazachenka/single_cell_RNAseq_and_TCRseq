expression_all <- readRDS("mean_Expression_for_Allgenes_for_allclusters.rds")
diversity_indexes <-read.csv("Diversity_indexes_scTCRseq_26patients_woithoutcontig3.csv",row.names = 1)
combined_TCR <- readRDS("combined_scTCRseq_26patients_withoutcontig3.rds")
Tcell_TCR <-readRDS("combined_scTCRseq_26patients_CD4_CD8_MAIT.rds")

function(input,output,session) {
  
  updateSelectizeInput(session,"gene",choices=colnames(expression_all[,c(10:27808)]),selected=c("CD28","TESPA1","GPR183"),server=TRUE)
  updateSelectizeInput(session,"cluster",choices=unique(expression_all$Cell_type),selected=c("CD4 Naive (2)"),server=TRUE)
  
  ### intro. scRNAseq results
  
  mytable_genes <- reactive({
    mytable <- expression_all %>%
      filter(BOR %in% input$bor) %>%
      filter(group %in% input$group) %>%
      filter(ORIGDIAG %in% input$diagnose) %>%
      filter(HISTCLAS %in% input$histology) %>%
      filter(Cell_type %in% input$cluster) 
    mytable <- mytable[,c("Patient","BOR","group","ORIGDIAG","HISTCLAS","Cell_type",input$gene)]
    mytable$mean_val <- rowMeans(mytable[,input$gene])
    mytable
  })
  
  mytable_freq <- reactive({
    mytable <- expression_all %>%
      filter(BOR %in% input$bor) %>%
      filter(group %in% input$group) %>%
      filter(ORIGDIAG %in% input$diagnose) %>%
      filter(HISTCLAS %in% input$histology) %>%
      filter(Cell_type %in% input$cluster) 
    mytable <- mytable[,c("Patient","BOR","group","ORIGDIAG","HISTCLAS","Cell_type","count","freq")]
  })
  
  
  ### TAB 1. Mean gene expression Boxplots
  output$Table1 <- DT::renderDataTable(server=FALSE,{
    DT::datatable(mytable_genes(), rownames = FALSE, extensions = "Buttons", 
                  options = list(dom = 'Bfrtip', scrollX = TRUE,
                                 buttons = list(
                                   list(extend = "copy", text = "Copy table", 
                                        # filename = "Table.csv",
                                        exportOptions = list(
                                          modifier = list(page = "all")
                                        )
                                   )
                                 )))
  })
  
  output$Boxplot1 <- renderPlot({
    mylevels <- c("B Naive","B Memory","Plasma","Mix B and T markers",
                  "CD14+ Mono","FCGR3A+ Mono","pDC","mDC",
                  "Platelet","PPBP+ CD14+ Mono","PPBP+ T cells",
                  "NK", "HPSC","CD4 Naive (1)","CD4 Naive (2)",
                  "CD4 CM","CD4 EM","CD8 Naive (1)","CD8 Naive (2)",
                  "CD8 CM","CD8 EM","CD8 CTL","MAIT","Proliferating T cells")
    myplot <- mytable_genes() %>%
      ggplot(aes(x=factor(Cell_type,level=mylevels),y=mean_val,fill=group,color=group)) +
      geom_boxplot(position=position_dodge(0.7),width=0.5) +
      geom_jitter(position=position_dodge(0.7)) +
      scale_fill_manual(values = alpha(c("#003B7F","#0064D2","#FDBD13"),0.5)) +
      scale_color_manual(values = c("#003B7F","#0064D2","#FDBD13")) +
      expand_limits(y=0) +
      ylab("Mean expression") +
      xlab("") +
      theme_bw()
    myplot
  })
  
  output$Stat1 <- renderDataTable({
      mytable_genes() %>%
      ungroup() %>%
      select(Cell_type,group,mean_val) %>%
      mutate(group=as.factor(group)) %>%
      group_by(Cell_type) %>%
      pairwise_wilcox_test(mean_val~group,p.adjust.method="none")
  })
  
  ### TAB 2. Mean gene expression Heatmap and Dotplot
  
  output$Heatmap <- renderPlot({
    annot <- mytable_genes() %>%
      ungroup() %>%
      select(Patient,BOR,group,ORIGDIAG,Cell_type) %>%
      arrange(Cell_type,group) %>%
      mutate(row_names = paste(Patient,Cell_type)) %>%
      as.data.frame() %>%
      set_rownames(.$row_names) %>%
      select(-Patient,-row_names)
    
    mean_expression <- mytable_genes()[,c("Patient","Cell_type",input$gene)] %>%
      ungroup() %>%
      mutate(row_names = paste(Patient,Cell_type)) %>%
      as.data.frame() %>%
      set_rownames(.$row_names) %>%
      select(-Patient,-row_names,-Cell_type) %>%
      t() %>%
      as.data.frame()
    
    my_plot <- pheatmap(mean_expression[,row.names(annot)],scale=input$scaled,annotation_col = annot,cluster_cols = as.logical(input$cluster_cols))
    my_plot
  })
    
  
  output$DotPlot <- renderPlot({
    mylevels <- c("B Naive","B Memory","Plasma","Mix B and T markers",
                  "CD14+ Mono","FCGR3A+ Mono","pDC","mDC",
                  "Platelet","PPBP+ CD14+ Mono","PPBP+ T cells",
                  "NK", "HPSC","CD4 Naive (1)","CD4 Naive (2)",
                  "CD4 CM","CD4 EM","CD8 Naive (1)","CD8 Naive (2)",
                  "CD8 CM","CD8 EM","CD8 CTL","MAIT","Proliferating T cells")
    
    dotplot_df1 <- mytable_genes() %>%
      ungroup() %>%
      select(Patient,Cell_type,mean_val)
    
    my_dotplot <- mytable_freq() %>%
      ungroup() %>%
      select(Patient,Cell_type,group,freq) %>%
      mutate(Patient_group=paste(group,Patient)) %>%
      left_join(dotplot_df1) %>%
      ggplot(aes(x=factor(Cell_type,level=mylevels),y=Patient_group,color=mean_val,size=freq)) +
      geom_point() +
      scale_color_gradient(low="blue",high="red") +
      xlab("") +
      theme_bw()
    
    my_dotplot
  })
  
  ### TAB 3. Cell type frequencues
  
    output$Table2 <- DT::renderDataTable(server=FALSE,{
      DT::datatable(mytable_freq(), rownames = FALSE, extensions = "Buttons", 
                    options = list(dom = 'Bfrtip', scrollX = TRUE,
                                   buttons = list(
                                     list(extend = "copy", text = "Copy table", 
                                          # filename = "Table.csv",
                                          exportOptions = list(
                                            modifier = list(page = "all")
                                          )
                                     )
                                   )))
    })
    
    output$Boxplot2 <- renderPlot({
      mylevels <- c("B Naive","B Memory","Plasma","Mix B and T markers",
                    "CD14+ Mono","FCGR3A+ Mono","pDC","mDC",
                    "Platelet","PPBP+ CD14+ Mono","PPBP+ T cells",
                    "NK", "HPSC","CD4 Naive (1)","CD4 Naive (2)",
                    "CD4 CM","CD4 EM","CD8 Naive (1)","CD8 Naive (2)",
                    "CD8 CM","CD8 EM","CD8 CTL","MAIT","Proliferating T cells")
      myplot <- mytable_freq() %>%
        ggplot(aes(x=factor(Cell_type,level=mylevels),y=freq,fill=group,color=group)) +
        geom_boxplot(position=position_dodge(0.7),width=0.5) +
        geom_jitter(position=position_dodge(0.7)) +
        scale_fill_manual(values = alpha(c("#003B7F","#0064D2","#FDBD13"),0.5)) +
        scale_color_manual(values = c("#003B7F","#0064D2","#FDBD13")) +
        expand_limits(y=0) +
        ylab("% of all cells") +
        xlab("") +
        theme(legend.position="none") +
        theme_bw()
      myplot
    })
    
    output$Stat2 <- renderDataTable({
      mytable_freq() %>%
        ungroup() %>%
        select(Cell_type,group,freq) %>%
        mutate(group=as.factor(group)) %>%
        group_by(Cell_type) %>%
        pairwise_wilcox_test(freq~group,p.adjust.method="none")
    })
    
    
    ### intro. scTCRseq
    
    mydiversityindex_table <- reactive({
      
      new_diversity_index <- diversity_indexes[,c("SUBJID",input$group_ind,input$div_ind)]
      colnames(new_diversity_index) <- c("SUBJID","ind_group","div_ind")
      new_diversity_index <- new_diversity_index[order(new_diversity_index$div_ind,decreasing = TRUE),]
      new_diversity_index <- new_diversity_index[order(new_diversity_index$ind_group),]
      new_diversity_index$SUBJID <- factor(new_diversity_index$SUBJID,levels=new_diversity_index$SUBJID[order(new_diversity_index$ind_group)])
      new_diversity_index
      
    })
    
    subset_combinedTCR <- reactive({
      subset <- subsetClones(combined_TCR,name="SUBJID",variable= mytable_genes()$Patient)
      subset
    })
    
    subset_combinedTCR_Tcells <- reactive({
#      tcell_pop <- c("CD4","CD8","MAIT")
#      tcell_pop <- tcell_pop[tcell_pop %in% input$Tcell_pop]
      subset <- subsetClones(Tcell_TCR,name="Patient",variable= mytable_genes()$Patient)
      subset <- subsetClones(subset,name="Cell_type2",variable=input$Tcell_pop)
      subset <- subsetClones(subset,name="Cell_type",variable=input$Tcell_subpop)
      subset
    })
    
    
    ### TAB 4. scTCRseq Summary
    
    output$Diversity_plot1 <- renderPlot({
      mydiversityindex_table() %>%
        filter(SUBJID %in% mytable_genes()$Patient) %>%
        ggplot(aes(x=SUBJID,y=div_ind,fill=ind_group,color=ind_group)) +
        geom_bar(stat="identity") +
        scale_fill_manual(values = alpha(c("#003B7F","#0064D2","#7AD0FF","#C0EAFE","#FDBD13"),0.5)) +
        scale_color_manual(values = c("#003B7F","#0064D2","#7AD0FF","#C0EAFE","#FDBD13")) +
        ylab(input$div_ind) +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = "bottom")
    })
    
    output$Diversity_plot2 <- renderPlot({
      mydiversityindex_table() %>%
        filter(SUBJID %in% mytable_genes()$Patient) %>%
        ggplot(aes(x=ind_group,y=div_ind,fill=ind_group,color=ind_group)) +
        geom_boxplot(outlier.shape=NA) +
        geom_jitter(position=position_jitter(0.1))+
        stat_compare_means() +
        scale_fill_manual(values = alpha(c("#003B7F","#0064D2","#7AD0FF","#C0EAFE","#FDBD13"),0.5)) +
        scale_color_manual(values = c("#003B7F","#0064D2","#7AD0FF","#C0EAFE","#FDBD13")) +
        ylab(input$div_ind) +
        xlab("") +
        theme_bw() +
        theme(axis.text.x = element_text(angle=45,hjust=1),legend.position = "bottom")
    })
    
    output$Plot_clonalOverlap <- renderPlot({
      plot <- clonalOverlap(combined_TCR,cloneCall=input$cloneCall,chain="both",method = input$method)
      plot
    })
    
    output$Plot_Homeostasis <- renderPlot({
      plot <- clonalHomeostasis(subset_combinedTCR(),cloneCall=input$cloneCall,group.by=input$group_ind)
      plot
    })
    
    output$Plot_Proportion <- renderPlot({
      plot <- clonalProportion(subset_combinedTCR(),cloneCall=input$cloneCall,group.by=input$group_ind)
      plot
    })
    
    output$Plot_UniqueClones <- renderPlot({
      plot <- clonalQuant(subset_combinedTCR(),cloneCall=input$cloneCall,group.by=input$group_ind,chain = "both")
      plot
    })
    
    output$Plot_UniqueClones_scaled <- renderPlot({
      plot <- clonalQuant(subset_combinedTCR(),cloneCall=input$cloneCall,group.by=input$group_ind,chain = "both",scale=TRUE)
      plot
    })
    
    output$Plot_clonalLength <- renderPlot({
      plot <- clonalLength(subset_combinedTCR(),cloneCall=input$cloneCall, chain = "both",group.by=input$group_ind)
      plot
    })
    
    ### TAB 5. scTCRseq T cells
    
    output$Plot_Homeostasis_Tcells <- renderPlot({
      plot <- clonalHomeostasis(subset_combinedTCR_Tcells(),cloneCall="strict",group.by=input$group_tcell)
      plot
    })
    
    output$Table_tcell <- DT::renderDataTable(server=FALSE,{
      table_out <- clonalHomeostasis(subset_combinedTCR_Tcells(),cloneCall="strict",group.by=input$group_tcell,exportTable = TRUE)
      DT::datatable(table_out, rownames = TRUE, extensions = "Buttons", 
                    options = list(dom = 'Bfrtip', scrollX = TRUE,
                                   buttons = list(
                                     list(extend = "copy", text = "Copy table", 
                                          # filename = "Table.csv",
                                          exportOptions = list(
                                            modifier = list(page = "all")
                                          )
                                     )
                                   )))
    })
    
    output$Plot_UniqueClones_scaled_Tcell <- renderPlot({
      plot <- clonalQuant(subset_combinedTCR_Tcells(),cloneCall="strict",group.by=input$group_tcell,chain = "both",scale=TRUE)
      plot
    })
    
    ### TAB 6. scTCRseq VDJ genes
    
    output$CDR3_AAsomposition <- renderPlot({
      plot <- percentAA(subset_combinedTCR(),chain = input$chain,aa.length = 20,group.by=input$group_vdj)
      plot
    })
    
    output$VJ_usage1 <- renderPlot({
      plot <- percentGenes(subset_combinedTCR(),chain = input$chain,gene = input$vdj_gene1,group.by=input$group_vdj)
      plot
    })
    
    output$VDJgene_plot1 <- renderPlot({
      plot <- vizGenes(subset_combinedTCR(),x.axis =input$vdj_gene2,plot="heatmap",group.by = input$group_vdj)
      plot
    })
    
    output$VDJgene_plot2 <- renderPlot({
      plot <- vizGenes(subset_combinedTCR(),x.axis =input$vdj_gene2,plot="barplot",group.by = input$group_vdj)
      plot
    })
    
    output$VDJ_table <- DT::renderDataTable(server=FALSE,{
      table_vdj <-  vizGenes(subset_combinedTCR(),x.axis =input$vdj_gene2,plot="heatmap",group.by = input$group_vdj,exportTable=TRUE)
      DT::datatable(table_vdj, rownames = FALSE, extensions = "Buttons", 
                    options = list(dom = 'Bfrtip', scrollX = TRUE,
                                   buttons = list(
                                     list(extend = "copy", text = "Copy table", 
                                          # filename = "Table.csv",
                                          exportOptions = list(
                                            modifier = list(page = "all")
                                          )
                                     )
                                   )))
    })
}