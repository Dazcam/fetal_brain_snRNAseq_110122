# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA GO plots 
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figure 7
#  Issues: Region IDs encoded differently in MAGMA and LSDR analysis

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(reshape2) # For melt

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/magma_GO/"
FIG_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/figures/'
SIG_CELLS <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'GE_InN_2', 'Hipp_ExN_3', 'Hipp_ExN_5')

GO_DF <- read.table(paste0(DATA_DIR, 'magma_GO.gsa.out'), header = FALSE) %>%
  janitor::row_to_names(row_number = 1)


for (CELL_TYPE in SIG_CELLS) {
  
  if (CELL_TYPE == "GE_InN_2") {
  
  GO_DF_FILT <- filter(GO_DF, grepl(CELL_TYPE, VARIABLE)) %>%
    mutate(Term = gsub(paste0("^[^", CELL_TYPE, "_]*", CELL_TYPE,"_"), "", VARIABLE)) %>%
    select(Term, P)

  PLOT <- ggplot(data = GO_DF_FILT, aes(x = -log10(as.numeric(P)), y = Term)) +
    geom_bar(stat = "identity", color = 'black', fill = '#3CBB75FF') +
    geom_vline(xintercept=-log10(0.05/121), linetype = "dashed", color = "black") +
  #  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    ggtitle(CELL_TYPE) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 12),
          legend.position = "none") +
    xlab(expression(-log[10](P))) +
    ylab('Term') +
    xlim(0, 11.5) 
  
  assign(paste0(CELL_TYPE, '_go_plot'), PLOT)
  
  } else {
    
    GO_DF_FILT <- filter(GO_DF, grepl(CELL_TYPE, VARIABLE)) %>%
      mutate(Term = gsub(paste0("^[^", CELL_TYPE, "_]*", CELL_TYPE,"_"), "", VARIABLE)) %>%
      select(Term, P)
    
    PLOT <- ggplot(data = GO_DF_FILT, aes(x = -log10(as.numeric(P)), y = Term)) +
      geom_bar(stat = "identity", color = 'black', fill = '#CEE5FD') +
      geom_vline(xintercept=-log10(0.05/121), linetype = "dashed", color = "black") +
      #  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(CELL_TYPE) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 12),
            legend.position = "none") +
      xlab(expression(-log[10](P))) +
      ylab('Term') +
      xlim(0, 11.5) 
    
    assign(paste0(CELL_TYPE, '_go_plot'), PLOT)
    
  }
  
  
}

group_plot <- plot_grid(FC_ExN_2_go_plot, FC_ExN_3_go_plot, 
                        FC_ExN_4_go_plot, GE_InN_2_go_plot, 
                        Hipp_ExN_3_go_plot, Hipp_ExN_5_go_plot, 
                        ncol = 2, labels = 'AUTO')


# Save
tiff(paste0(FIG_DIR, "Fig_7.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
group_plot
dev.off()


# Jpegs
jpeg(paste0(FIG_DIR, "Fig_7.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
group_plot
dev.off()


# All cell types in same plot draft

# group_plot <- plot_grid(FC_ExN_2_go_plot, FC_ExN_3_go_plot, FC_ExN_4_go_plot,
#                         GE_InN_2_go_plot, Hipp_ExN_3_go_plot, Hipp_ExN_5_go_plot)
# 
# 
# PLOT_DF <- GO_DF %>%
#   separate(VARIABLE, c("Cell_type", "Term"), sep = "_GO") %>%
#   mutate(Term = paste0('GO', Term)) %>%
#   select(Term, Cell_type, P)
  

# ggplot(data = PLOT_DF, aes(x = -log10(as.numeric(P)), y = Term, fill = Cell_type)) +
#   geom_bar(stat = "identity", color = 'black', position = "dodge") +
#   geom_vline(xintercept=-log10(0.05/121), linetype = "dashed", color = "black") +
#   #  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
#   theme_bw() +
# #  ggtitle(CELL_TYPE) +
#   theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(colour = "black", size = 1),
#         plot.title = element_text(hjust = 0.5),
#         axis.title.x = element_text(colour = "#000000", size = 14),
#         axis.title.y = element_text(colour = "#000000", size = 14),
#         axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
#         axis.text.y  = element_text(colour = "#000000", size = 12)) +
#   xlab(expression(-log[10](P))) +
#   ylab('Term') +
#   xlim(0, 11.5) 

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------