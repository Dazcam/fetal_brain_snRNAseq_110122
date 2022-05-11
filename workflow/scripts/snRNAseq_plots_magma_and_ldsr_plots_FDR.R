# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA Celltyping and LDSR plots 
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures X and X
#  Issues: Region IDs encoded differently in MAGMA and LSDR analysis
#          FDR processing slightly different from BF correction 

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(reshape2) # For melt

##  Initialise variables  -------------------------------------------------------------
MAGMA_REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
LDSR_REGIONS <- c("Cer", "FC", "GE", "Hipp", "Thal")
GWAS <- c('ASD', 'SCZ', 'HEIGHT')
MAGMA_DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/magma_celltyping/"
LDSC_DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/LDSR_part_herit/"
FIG_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/figures/'
NAME_BODY <- '_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN'
SUFFIX_TOP10 <- '_top10.gsa.out'

# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (MAGMA_REGION in MAGMA_REGIONS) {
  
  for (DISORDER in GWAS) {
    
    magma_top10_df <- read.table(paste0(MAGMA_DATA_DIR, DISORDER, NAME_BODY, 
                                   '.level1.', MAGMA_REGION, SUFFIX_TOP10), header = FALSE) %>%
      janitor::row_to_names(row_number = 1) %>% 
      mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
      dplyr::rename(Category = VARIABLE) # Match LDSRs cell-type column
    
    if (MAGMA_REGION == 'cer') {
      
      assign(paste0('magma_Cer_', DISORDER, '_df'), magma_top10_df, envir = .GlobalEnv) 
      
      } else if (MAGMA_REGION == 'hip') {
      
      assign(paste0('magma_Hipp_', DISORDER, '_df'), magma_top10_df, envir = .GlobalEnv) 
      
      } else if (MAGMA_REGION == 'pfc') {
      
      assign(paste0('magma_FC_', DISORDER, '_df'), magma_top10_df, envir = .GlobalEnv) 
      
      } else if (MAGMA_REGION == 'tha') {
      
      assign(paste0('magma_Thal_', DISORDER, '_df'), magma_top10_df, envir = .GlobalEnv) 
      
      } else {
      
      assign(paste0('magma_GE_', DISORDER, '_df'), magma_top10_df, envir = .GlobalEnv) 
      
      }
 
  }
  
}


# Calculate magma FDR - note we need to first cat the dfs together to calculate FDR
for (DISORDER in GWAS) {

  magma_all_regions <- rbind(get(paste0('magma_Cer_', DISORDER, '_df')),
                             get(paste0('magma_FC_', DISORDER, '_df')),
                             get(paste0('magma_GE_', DISORDER, '_df')),
                             get(paste0('magma_Hipp_', DISORDER, '_df')),
                             get(paste0('magma_Thal_', DISORDER, '_df'))) %>%
  mutate(magma_FDR = p.adjust(P, 'BH', length(P))) %>%
  mutate(magma_FDR_log10 = -log10(as.numeric(magma_FDR))) %>%
  arrange(Category)
  
  assign('MAGMA_all_regions_with_FDR_df', magma_all_regions, envir = .GlobalEnv) 
  write_tsv(MAGMA_all_regions_with_FDR_df, paste0(FIG_DIR, DISORDER, '_Magma_FDR_table.txt'))
  
  for (LDSR_REGION in LDSR_REGIONS) {
    
    MAGMA_single_region <- magma_all_regions %>%
      filter(grepl(LDSR_REGION, Category))
    
    assign(paste0('magma_', LDSR_REGION, '_', DISORDER, '_FDR_df'), MAGMA_single_region, envir = .GlobalEnv)
    
  }


}

# LDSR - prepare df
cat('\nPreparing LDSR data ... \n')
for (LDSR_REGION in LDSR_REGIONS) {
  
  for (DISORDER in GWAS) {

    ldsr_top10_df <- read_tsv(paste0(LDSC_DATA_DIR, 'snRNAseq_LDSC_', DISORDER, '_baseline.v1.2_top10pc.tsv')) %>%
      mutate(ZscoreP = pnorm(`Coefficient_z-score`, lower.tail = FALSE))
    ldsr_top10_df <- as.data.frame(filter(ldsr_top10_df, grepl(LDSR_REGION, Category))) 
    
    assign(paste0('ldsr_', LDSR_REGION, '_', DISORDER, '_df'), ldsr_top10_df, envir = .GlobalEnv) 
    
  }
  
}

# Calculate LSDR FDR 
# Need to rbind all regions together to get FDR scores for each individual cell type
# May want to remove (or move) FDR_log10 column from LDSR_all_regions_with_FDR_df
for (DISORDER in GWAS) {

  LDSR_all_regions <- rbind(get(paste0('ldsr_Cer_', DISORDER, '_df')),
                            get(paste0('ldsr_FC_', DISORDER, '_df')),
                            get(paste0('ldsr_GE_', DISORDER, '_df')),
                            get(paste0('ldsr_Hipp_', DISORDER, '_df')),
                            get(paste0('ldsr_Thal_', DISORDER, '_df'))) %>%
  mutate(ldsr_FDR = p.adjust(ZscoreP, 'BH', length(ZscoreP))) %>%
  mutate(ldsr_FDR_log10 = if_else(`Coefficient_z-score` > 0, -log10(ldsr_FDR), 0)) 
  
  assign('LDSR_all_regions_with_FDR_df', LDSR_all_regions, envir = .GlobalEnv) 
  write_tsv(LDSR_all_regions_with_FDR_df, paste0(FIG_DIR, DISORDER, '_LDSR_FDR_table.txt'))
  
  for (LDSR_REGION in LDSR_REGIONS) {
    
    LDSR_single_region <- LDSR_all_regions %>%
      filter(grepl(LDSR_REGION, Category))
    
    assign(paste0('ldsr_', LDSR_REGION, '_', DISORDER, '_FDR_df'), LDSR_single_region, envir = .GlobalEnv)
    
    
  }

}

# Plot
cat('\nCreate plots ... \n')
for (LDSR_REGION in LDSR_REGIONS) {
  
  for (DISORDER in GWAS) {
    
    if (LDSR_REGION == 'FC') {
      REGION = 'Frontal Cortex'
    } else if (LDSR_REGION == 'GE') {
      REGION = 'Ganglionic Eminence'
    } else if (LDSR_REGION == 'Hipp') {
      REGION = 'Hippocampus'
    } else if (LDSR_REGION == 'Cer') {
      REGION = 'Cerebellum'
    } else {
      REGION = 'Thalamus'
    }
      
    
    PLOT_DF <- left_join(get(paste0('magma_', LDSR_REGION, '_', DISORDER, '_FDR_df')), 
                        get(paste0('ldsr_', LDSR_REGION, '_', DISORDER, '_FDR_df')), 
                        by = 'Category') %>%
      select(Category, magma_FDR_log10, ldsr_FDR_log10) %>% 
      rename(Magma = magma_FDR_log10,
             LDSR = ldsr_FDR_log10) %>%
      melt()
      
        
    MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))), 
                                                  fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(REGION) +
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
      
    assign(paste0(LDSR_REGION, '_', DISORDER, '_magma_ldsr_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
    
  }
  
}

legend <- get_legend(
  # create some space to the left of the legend
  ggplot(data = PLOT_DF, aes(x = value, y = Category, fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    theme(legend.key.size = unit(1.5, 'cm'),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) 
)

cat('\nCreating group plots ... \n')
for (DISORDER in GWAS) {
  
  magma_ldsr_plot <- plot_grid(get(paste0('FC_', DISORDER, '_magma_ldsr_plot')),
                               get(paste0('GE_', DISORDER, '_magma_ldsr_plot')),
                               get(paste0('Hipp_', DISORDER, '_magma_ldsr_plot')), 
                               get(paste0('Thal_', DISORDER, '_magma_ldsr_plot')),
                               get(paste0('Cer_', DISORDER, '_magma_ldsr_plot')),
                               legend, label_size = 16)
  
  assign(paste0('all_regions_', DISORDER, '_magma_ldsr_FDR_plot'), magma_ldsr_plot, envir = .GlobalEnv)
  
}

# Save plots
# Fig 2 - SCZ
tiff(paste0(FIG_DIR, "Fig_SCZ_FDR.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_SCZ_magma_ldsr_FDR_plot
dev.off()

# Fig 3 - HEIGHT
tiff(paste0(FIG_DIR, "Fig_HEIGHT_FDR.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_HEIGHT_magma_ldsr_FDR_plot
dev.off()

# Fig X - ASD
tiff(paste0(FIG_DIR, "Fig_ASD_FDR.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
all_regions_ASD_magma_ldsr_FDR_plot
dev.off()

# Jpegs
jpeg(paste0(FIG_DIR, "Fig_SCZ_FDR.jpg"), width = 960, height = 960, 
     units = "px", pointsize = 12, quality = 150)
all_regions_SCZ_magma_ldsr_FDR_plot
dev.off()

jpeg(paste0(FIG_DIR, "Fig_HEIGHT_FDR.jpg"), width = 960, height = 960, 
     units = "px", pointsize = 12, quality = 150)
all_regions_HEIGHT_magma_ldsr_FDR_plot
dev.off()

jpeg(paste0(FIG_DIR, "Fig_ASD_FDR.jpg"), width = 960, height = 960, 
     units = "px", pointsize = 12, quality = 150)
all_regions_ASD_magma_ldsr_FDR_plot
dev.off()

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


