#--------------------------------------------------------------------------------------
#
#    snRNAseq create differentially expressed gene tables - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Supplementary tables 2-6

##  Load Packages  --------------------------------------------------------------------
library(cowplot)
library(tidyverse)
library(Seurat)

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/'
FIG_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/figures/'
REGIONS <- c('pfc', 'wge', 'hip',  'tha', 'cer')

# Load Seurat objects
for (REGION in REGIONS) {
  
  OBJ <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  
  assign(paste0('seurat.', REGION), OBJ)
  
}

for (REGION in c('cer', 'hip', 'pfc', 'tha', 'wge')) {
  
  cat(paste0('\n\n\nGet the ', REGION, 
             ' top DE genes per cluster vs. all other clusters ...\n\n\n'))
  
  OBJ <- get(paste0('seurat.', REGION))
  
  OBJ.markers <- FindAllMarkers(OBJ, only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25) 
  assign(paste0(REGION, '_DE_markers'), OBJ.markers, .GlobalEnv)
  
}

# Save tables - extract only genes with BF corrected
for (INDEX in 1:length(REGIONS)) {
  
  INDEX_2 <- INDEX + 1
  
  OBJ <- get(paste0(REGIONS[INDEX], '_DE_markers')) %>%
    relocate(cluster, gene) %>%
    filter(p_val_adj < 0.05)
  print(head(OBJ))
  write.table(OBJ, paste0(FIG_DIR, "Supplementary_table_", INDEX_2, ".txt"), sep = '\t', 
              row.names = FALSE, quote = FALSE, col.names = TRUE)
  
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------