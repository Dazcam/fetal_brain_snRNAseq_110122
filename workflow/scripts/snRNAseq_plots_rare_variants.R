# -------------------------------------------------------------------------------------
#
#    snRNAseq rare variant test plots 
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 8, 9
#  Issues: Region IDs encoded differently in MAGMA and LSDR analysis
#          Much tweaking needed between ggplot and cowplot for FDR < 5% line

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(ggsignif)
library(cowplot)
library(reshape2) # For melt

##  Initialise variables  -------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/rare_variants/'
FIG_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/figures/'
STUDIES <- c('schema_genes', 'asd_genes')
REGIONS <- c('Cer', 'FC', 'GE', 'Hipp', 'Thal')

# Prepare df
cat('\nPreparing rare variant data ... \n')
for (STUDY in STUDIES) {

  for (REGION in REGIONS) {
    
  wilcoxon_df <- read_tsv(paste0(DATA_DIR, STUDY, '_', REGION, '_wilcoxon.txt'))
  assign(paste0(STUDY, '_', REGION, '_wilcoxon_df'), wilcoxon_df)
  
  }
  
}



# Plot
cat('\nCreate plots ... \n')
for (REGION in REGIONS) {
  
  for (STUDY in STUDIES) {
    
    if (REGION == 'FC') {
      REGION_TITLE = 'Frontal Cortex'
    } else if (REGION == 'GE') {
      REGION_TITLE = 'Ganglionic Eminence'
    } else if (REGION == 'Hipp') {
      REGION_TITLE = 'Hippocampus'
    } else if (REGION == 'Cer') {
      REGION_TITLE = 'Cerebellum'
    } else {
      REGION_TITLE = 'Thalamus'
    }
    
    PLOT_DF <- get(paste0(STUDY, '_', REGION, '_wilcoxon_df'))
    
    WILCOXON_PLOT <- ggplot(data = PLOT_DF, aes(x = -log10(P), y = factor(cell_type, rev(levels(factor(cell_type)))))) +
      geom_bar(stat = "identity", color = 'black', fill = 'darkseagreen') +
      geom_vline(xintercept=-log10(0.00054), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(REGION_TITLE) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 13),
            legend.position = "none") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 11.5) 
    
    assign(paste0(REGION, '_', STUDY, '_wilcoxon_plot'), WILCOXON_PLOT, envir = .GlobalEnv) 
    
  }
  
}

## Add FDR significance lines
#  Schema
FC_schema_genes_wilcoxon_plot <- FC_schema_genes_wilcoxon_plot +
  geom_segment(aes(x = 11, y = 13.6, xend = 11, yend = 14.4)) + # FC-ExN-2
  annotate("text", x = 11.5, y = 13.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 8.6, xend = 11, yend = 9.4)) + # FC-InN-2
  annotate("text", x = 11.5, y = 8.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 2.6, xend = 11, yend = 3.4)) + # FC-OPC
  annotate("text", x = 11.5, y = 2.77, label = "*", size = 7) 

Hipp_schema_genes_wilcoxon_plot <- Hipp_schema_genes_wilcoxon_plot +
  geom_segment(aes(x = 11, y = 10.6, xend = 11, yend = 11.4)) + # Hipp-ExN-5
  annotate("text", x = 11.5, y = 10.75, label = "*", size = 7)

# ASD
Cer_asd_genes_wilcoxon_plot <- Cer_asd_genes_wilcoxon_plot +
  geom_segment(aes(x = 11, y = 19.6, xend = 11, yend = 20.4)) + # Cer-Endo
  annotate("text", x = 11.5, y = 19.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 18.6, xend = 11, yend = 19.4)) + # Cer-ExN-1
  annotate("text", x = 11.5, y = 18.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 17.6, xend = 11, yend = 18.4)) + # Cer-ExN-2
  annotate("text", x = 11.5, y = 17.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 16.6, xend = 11, yend = 17.4)) + # Cer-ExN-3
  annotate("text", x = 11.5, y = 16.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 10.6, xend = 11, yend = 11.4)) + # Cer-InN-4
  annotate("text", x = 11.5, y = 10.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 9.6, xend = 11, yend = 10.4)) + # Cer-InN-5
  annotate("text", x = 11.5, y = 9.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 8.6, xend = 11, yend = 9.4)) + # Cer-InN-6
  annotate("text", x = 11.5, y = 8.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 7.6, xend = 11, yend = 8.4)) + # Cer-InN-7
  annotate("text", x = 11.5, y = 7.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 6.6, xend = 11, yend = 7.4)) + # Cer-InN-8
  annotate("text", x = 11.5, y = 6.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 5.6, xend = 11, yend = 6.4)) + # Cer-IP
  annotate("text", x = 11.5, y = 5.72, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 4.6, xend = 11, yend = 5.4)) + # Cer-MG
  annotate("text", x = 11.5, y = 4.72, label = "*", size = 7) 

FC_asd_genes_wilcoxon_plot <- FC_asd_genes_wilcoxon_plot +
  geom_segment(aes(x = 11, y = 16.6, xend = 11, yend = 17.4)) + # FC-CycPro
  annotate("text", x = 11.5, y = 16.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 14.6, xend = 11, yend = 15.4)) + # FC-ExN-1
  annotate("text", x = 11.5, y = 14.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 13.6, xend = 11, yend = 14.4)) + # FC-ExN-2
  annotate("text", x = 11.5, y = 13.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 12.6, xend = 11, yend = 13.4)) + # FC-ExN-3
  annotate("text", x = 11.5, y = 12.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 11.6, xend = 11, yend = 12.4)) + # FC-ExN-4
  annotate("text", x = 11.5, y = 11.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 10.6, xend = 11, yend = 11.4)) + # FC-ExN-5
  annotate("text", x = 11.5, y = 10.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 9.6, xend = 11, yend = 10.4)) + # FC-InN-1
  annotate("text", x = 11.5, y = 9.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 8.6, xend = 11, yend = 9.4)) + # FC-InN-2
  annotate("text", x = 11.5, y = 8.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 7.6, xend = 11, yend = 8.4)) + # FC-InN-3
  annotate("text", x = 11.5, y = 7.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 3.6, xend = 11, yend = 4.4)) + # FC-N-undef
  annotate("text", x = 11.5, y = 3.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 2.6, xend = 11, yend = 3.4)) + # FC-OPC
  annotate("text", x = 11.5, y = 2.77, label = "*", size = 7) 

GE_asd_genes_wilcoxon_plot <- GE_asd_genes_wilcoxon_plot +
  geom_segment(aes(x = 11, y = 9.6, xend = 11, yend = 10.4)) + # FC-InN-1
  annotate("text", x = 11.5, y = 9.85, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 8.6, xend = 11, yend = 9.4)) + # FC-InN-2
  annotate("text", x = 11.5, y = 8.85, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 7.6, xend = 11, yend = 8.4)) + # FC-InN-3
  annotate("text", x = 11.5, y = 7.85, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 6.6, xend = 11, yend = 7.4)) + # FC-InN-4
  annotate("text", x = 11.5, y = 6.85, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 5.6, xend = 11, yend = 6.4)) + # FC-InN-5
  annotate("text", x = 11.5, y = 5.85, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 4.6, xend = 11, yend = 5.4)) + # FC-InN-6
  annotate("text", x = 11.5, y = 4.85, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 3.6, xend = 11, yend = 4.4)) + # FC-InN-7
  annotate("text", x = 11.5, y = 3.85, label = "*", size = 7)

Hipp_asd_genes_wilcoxon_plot <- Hipp_asd_genes_wilcoxon_plot +
  geom_segment(aes(x = 11, y = 18.6, xend = 11, yend = 19.4)) + # Hipp-CR-1
  annotate("text", x = 11.5, y = 18.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 15.6, xend = 11, yend = 16.4)) + # Hipp-Endo
  annotate("text", x = 11.5, y = 15.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 12.6, xend = 11, yend = 13.4)) + # Hipp-3
  annotate("text", x = 11.5, y = 12.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 10.6, xend = 11, yend = 11.4)) + # Hipp-5
  annotate("text", x = 11.5, y = 10.77, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 4.6, xend = 11, yend = 5.4)) + # Hipp-6
  annotate("text", x = 11.5, y = 4.77, label = "*", size = 7) 

Thal_asd_genes_wilcoxon_plot <- Thal_asd_genes_wilcoxon_plot +
  geom_segment(aes(x = 11, y = 21.6, xend = 11, yend = 22.4)) + # Thal-Endo
  annotate("text", x = 11.5, y = 21.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 20.6, xend = 11, yend = 21.4)) + # Thal-ExN-1
  annotate("text", x = 11.5, y = 20.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 18.6, xend = 11, yend = 19.4)) + # Thal-ExN-3
  annotate("text", x = 11.5, y = 18.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 14.6, xend = 11, yend = 15.4)) + # Thal-InN-4
  annotate("text", x = 11.5, y = 14.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 13.6, xend = 11, yend = 14.4)) + # Thal-InN-5
  annotate("text", x = 11.5, y = 13.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 12.6, xend = 11, yend = 13.4)) + # Thal-InN-6
  annotate("text", x = 11.5, y = 12.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 11.6, xend = 11, yend = 12.4)) + # Thal-InN-7
  annotate("text", x = 11.5, y = 11.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 9.6, xend = 11, yend = 10.4)) + # Thal-IP
  annotate("text", x = 11.5, y = 9.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 8.6, xend = 11, yend = 9.4)) + # Thal-MG
  annotate("text", x = 11.5, y = 8.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 3.6, xend = 11, yend = 4.4)) + # Thal-RG-4
  annotate("text", x = 11.5, y = 3.7, label = "*", size = 7) +
  geom_segment(aes(x = 11, y = 1.6, xend = 11, yend = 2.4)) + # Thal-RG-6
  annotate("text", x = 11.5, y = 1.7, label = "*", size = 7)

cat('\nCreating group plots ... \n')
for (STUDY in STUDIES) {
  
  wilcoxon_group_plot <- plot_grid(get(paste0('FC_', STUDY, '_wilcoxon_plot')),
                               get(paste0('GE_', STUDY, '_wilcoxon_plot')),
                               get(paste0('Hipp_', STUDY, '_wilcoxon_plot')), 
                               get(paste0('Thal_', STUDY, '_wilcoxon_plot')),
                               get(paste0('Cer_', STUDY, '_wilcoxon_plot')))
  
  assign(paste0('all_regions_', STUDY, '_wilcoxon_plot'), wilcoxon_group_plot, envir = .GlobalEnv)
  
}

# Save plots
# Fig 8 - SCZ
tiff(paste0(FIG_DIR, "Fig_8.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_schema_genes_wilcoxon_plot
dev.off()

# Fig 9 - ASD
tiff(paste0(FIG_DIR, "Fig_9.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_asd_genes_wilcoxon_plot
dev.off()

# Jpegs
jpeg(paste0(FIG_DIR, "Fig_8.jpg"), width = 960, height = 960, 
     units = "px", pointsize = 12, quality = 150)
all_regions_schema_genes_wilcoxon_plot
dev.off()

jpeg(paste0(FIG_DIR, "Fig_9.jpg"), width = 960, height = 960, 
     units = "px", pointsize = 12, quality = 150)
all_regions_asd_genes_wilcoxon_plot
dev.off()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


