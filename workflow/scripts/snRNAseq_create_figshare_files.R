#--------------------------------------------------------------------------------------
#
#    snRNAseq extract raw GEX counts and metadata for figshare
#
#------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(Seurat)
library(data.table)
library(tidyverse)

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/'
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/figshare_files/'
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')

##  Extract Data  ---------------------------------------------------------------------
for (REGION in REGIONS) {
  
  cat('\n\n\n\nCreating raw count mat and metadata for:', REGION, '...\n')
  
  if (REGION == 'cer') {
    
    REGION_NEW <-  'Cer'
    
  } else if (REGION == 'hip') {
    
    REGION_NEW <- 'Hipp'
    
  } else if (REGION == 'pfc') {
    
    REGION_NEW <- 'FC'
    
  } else if (REGION == 'tha') {
    
    REGION_NEW <- 'Thal'
    
  } else {
    
    REGION_NEW <- 'GE'
    
  }
  
  cat('REGION_NEW =', REGION_NEW, '...\n')
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION,'.final.rds'))
  
  # Extract raw count matrix
  cat('Generating matrix ...\n')
  counts_mat <- as.matrix(GetAssayData(object = seurat.obj, slot = "counts"))  
  write.table(counts_mat, file = paste0(OUT_DIR, 'cameron_2021_snRNAseq_', REGION_NEW, 
                                        '_raw_count_gEX_matrix.txt'), quote = FALSE)
  
  # Extract metadata
  cat('Generating metadata ...\n')
  seurat.obj@meta.data %>%
    select(cellIDs, Sample) %>%
    rownames_to_column('cells') %>%
    rename(sample = Sample) %>%
    write_tsv(file = paste0(OUT_DIR, 'cameron_2021_snRNAseq_', REGION_NEW, '_metadata.txt'))
  
  cat('Generating UMAP ...\n')
  cluster_plot <- DimPlot(seurat.obj, reduction = "umap", group.by = 'cellIDs',
                          label = TRUE, label.size = 5,
                          pt.size = 0.1, repel = TRUE) + ggtitle(NULL) 
  
  assign(paste0(REGION, '_plot'), cluster_plot, .GlobalEnv)
  
  cat('Done.\n\n\n\n\n\n')
  
}

data.table::fread(paste0(OUT_DIR, 'cameron_2021_snRNAseq_', 'Thal', 
                         '_raw_count_gEX_matrix.txt'))[1:6,1:6]




#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------





