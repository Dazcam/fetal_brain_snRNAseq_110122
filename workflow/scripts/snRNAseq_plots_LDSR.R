# -------------------------------------------------------------------------------------
#
#    snRNAseq sLSDC plots  - final
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures 3, 5 

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)

##  Initialise Categorys  -------------------------------------------------------------
REGIONS <- c("Cer", "FC", "GE", "Hipp", "Thal")
GWAS <- c('SCZ', 'HEIGHT')
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/LDSR_part_herit/"
FIG_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/figures/'

## Load Data  -------------------------------------------------------------------------
for (DISORDER in GWAS) {
  
  for (REGION in REGIONS) {
    
    # Load data
    all_top10_df <- read_tsv(paste0(DATA_DIR, 'snRNAseq_LDSC_', DISORDER, '_baseline.v1.2_top10pc.tsv')) %>%
      mutate(Coeff_P = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0))
    subset_top10_df <- as.data.frame(filter(all_top10_df, grepl(REGION, Category)))
    
    # Remove regional info from cell types on y axis
    #subset_top10_df$Category <- gsub("^.*?-", "", subset_top10_df$Category)
    
    ##  Plot  ---------------------------------------------------------------------------
    # Update region names in plot titles to new ones for FC and GE
    if (REGION == "FC") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = Coeff_P, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c("#DCBEFF", '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                                             '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                                             '#3CBB75FF', '#D078FF', '#F58231', '#CCCCCC', '#FDE725FF', 
                                             '#FF5959', '#FF5959'), color = 'black') +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Frontal Cortex') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression(-log[10](P))) +
        ylab('Cell type') +
        xlim(0, 8)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
    } else if (REGION == "Cer") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = Coeff_P, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                                             '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                                             '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                                             "#D078FF", '#F58231', '#FDE725FF', '#FF5959', '#FF5959', 
                                             '#FF5959'), color = 'black') +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Cerebellum') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression(-log[10](P))) +
        ylab('Cell type') +
        xlim(0, 8)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
    } else if (REGION == "Hipp") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = Coeff_P, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c('#B200ED', '#B200ED', "#DCBEFF", '#9A6324', '#CEE5FD',   
                                             '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                                             '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#F58231',   
                                             '#FDE725FF', '#FF5959', '#FF5959', '#FF5959'), color = 'black') +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Hippocampus') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression(-log[10](P))) +
        ylab('Cell type') +
        xlim(0, 8)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
    } else if (REGION == "GE") {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = Coeff_P, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c('#DCBEFF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                             '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#FF5959', '#FF5959', 
                                             '#FF5959'), color = 'black') +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Ganglionic Eminence') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression(-log[10](P))) +
        ylab('Cell type') +
        xlim(0, 8)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_data'), subset_top10_df, envir = .GlobalEnv)
      
      
    } else {
      
      top10Plot <- ggplot(data = subset_top10_df, aes(x = Coeff_P, y = factor(Category, rev(levels(factor(Category)))))) +
        geom_bar(stat = "identity", fill = c("#DCBEFF", '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD',  
                                             '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                                             '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', "#D078FF", '#F58231', 
                                             '#FDE725FF', '#FF5959', '#FF5959', '#FF5959', '#FF5959', 
                                             '#FF5959', '#FF5959', '#FF5959'), color = 'black') +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        geom_vline(xintercept=-log10(0.05/91), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle('Thalamus') +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 12)) +
        xlab(expression(-log[10](P))) +
        ylab('Cell type') +
        xlim(0, 8)
      
      assign(paste0(REGION, '_', DISORDER, '_ldsc_top10_plot'), top10Plot, envir = .GlobalEnv)
      assign(x = paste0(REGION, '_', DISORDER, '_ldsc_top10_df'), value = subset_top10_df, envir = .GlobalEnv)
      
    }
    
  }
  
}

cat('\nCreating group plots ... \n')
for (DISORDER in GWAS) {
  
  ldsc_top10_plot <- plot_grid(get(paste0('FC_', DISORDER, '_ldsc_top10_plot')),
                               get(paste0('GE_', DISORDER, '_ldsc_top10_plot')),
                               get(paste0('Hipp_', DISORDER, '_ldsc_top10_plot')), 
                               get(paste0('Cer_', DISORDER, '_ldsc_top10_plot')),
                               get(paste0('Thal_', DISORDER, '_ldsc_top10_plot')))
  
  assign(paste0('all_regions_', DISORDER, '_ldsc_top10_plot'), ldsc_top10_plot, envir = .GlobalEnv)
  
}

# Save plots
# Fig 3 - SCZ
tiff(paste0(FIG_DIR, "Fig_3.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_SCZ_ldsc_top10_plot
dev.off()

# Fig 5 - HEIGHT
tiff(paste0(FIG_DIR, "Fig_5.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_HEIGHT_ldsc_top10_plot
dev.off()

# Jpegs
jpeg(paste0(FIG_DIR, "Fig_3.jpg"), width = 960, height = 960, 
     units = "px", pointsize = 12, quality = 150)
all_regions_SCZ_ldsc_top10_plot
dev.off()

jpeg(paste0(FIG_DIR, "Fig_5.jpg"), width = 960, height = 960, 
     units = "px", pointsize = 12, quality = 150)
all_regions_HEIGHT_ldsc_top10_plot
dev.off()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------