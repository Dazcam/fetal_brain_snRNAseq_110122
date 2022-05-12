#--------------------------------------------------------------------------------------
#
#    snRNAseq - Fetal cell types Schema enrichment
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Extract the cell specificity scores for all cell types from ctd objects  
#  2. Take top 32 schema genes and Q values (i.e those FDR < 0.05)
#  3. Run a wilcoxon test on specificity scores of each of 91 cell types against Qvalues  
#     of schema genes

# Note: H1-4 in schema data is encoded HIST1H1E in our data

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(biomaRt)

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/'
FIG_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/figures/'
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')

## Load Data --------------------------------------------------------------------------
# Load fetal specificty scores for all 91 cell types 
for (REGION in REGIONS) {
 
  load(paste0(DATA_DIR, 'ctd_objects/CellTypeData_', REGION, '.rda'))
  ctd_specificity <- as.data.frame(ctd[[1]]$specificity)
  assign(paste0(REGION, '_specificity'), ctd_specificity)
   
}

# Load SCHEMA genes - note that although SCHEMA genes on website show gene symbols, the csv 
# genes have ensemble encoding so need conversion file
schema <- read_csv(paste0(DATA_DIR, 'public_datasets/schema/meta_results_2022_05_12_09_17_34.csv')) %>%
  head(32) %>%
  dplyr::select(Gene, `Q meta`) %>%
  rename(Q = `Q meta`)
schema_genes <- read_tsv(paste0(DATA_DIR, 'public_datasets/schema/schema_gene_conversion.txt')) %>%
  left_join(schema) %>%
  arrange(Q)


# Run Wilcoxon tests on all 91 cell types
for (REGION in REGIONS) {
  
  # Load regional specificity scores
  specificity_DF <- get(paste0(REGION, '_specificity')) %>%
    rownames_to_column(var = 'gene')
  
  for (CELL_TYPE in colnames(specificity_DF)) {
    
    
    if (CELL_TYPE == 'gene') next
    cat('\n\nRunning wilcoxon test for:', CELL_TYPE)
    
    specificity_cell <- data.frame(gene = specificity_DF$gene,
               cell_scores = specificity_DF[[CELL_TYPE]]) %>%
      mutate(schema = ifelse(gene %in% schema_genes$symbol, 'schema', 'no_schema'))
    
    schema_scores <- specificity_cell %>%
      filter(grepl('^schema', schema)) %>%
      pull(cell_scores)
    
    no_schema_scores <- specificity_cell %>%
      filter(grepl('no_schema', schema)) %>%
      pull(cell_scores)
    
    wilcox_result <- wilcox.test(schema_scores, no_schema_scores, paired = FALSE, data = specificity_cell)
    
    if (exists('wilcoxon_df')) {
      
      wilcox_new_result <- as.data.frame(t(c(CELL_TYPE, wilcox_result$p.value)))
      colnames(wilcox_new_result) <- c('cell_type', 'p-value')
      wilcoxon_df <- rbind(wilcoxon_df, wilcox_new_result)
      
    } else {
      
      wilcoxon_df <- data.frame()
      wilcoxon_df <- rbind(wilcoxon_df, as.data.frame(t(c(CELL_TYPE, wilcox_result$p.value))))
      colnames(wilcoxon_df) <- c('cell_type', 'p-value')
      
    }
    
    cat('\nNumber of schema genes in sample:', sum(specificity_cell$schema == 'schema'))
    
  }
    
  
}

write_tsv(wilcoxon_df %>% arrange(`p-value`), paste0(FIG_DIR, 'SCHEMA_wilcoxon.txt'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
