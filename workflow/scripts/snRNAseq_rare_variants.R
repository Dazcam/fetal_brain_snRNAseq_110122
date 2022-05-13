#--------------------------------------------------------------------------------------
#
#    snRNAseq - Fetal cell types Schema enrichment
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Extract the cell specificity scores for all cell types from ctd objects  
#  2. Take top 32 schema genes (i.e those FDR < 0.05)
#  3. Run a wilcoxon test on specificity scores of schema genes in each of 91 cell types 
#     against the specificty scores for the rest of the genes in each cell type

# Note: H1-4 in schema data is encoded HIST1H1E in our data

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)

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
schema <- read_csv(paste0(DATA_DIR, 'public_datasets/rare_variants/meta_results_2022_05_12_09_17_34.csv')) %>%
  head(32) %>%
  dplyr::select(gene, `Q meta`) %>%
  rename(Q = `Q meta`)
schema_genes <- read_tsv(paste0(DATA_DIR, 'public_datasets/rare_variants/schema_gene_conversion.txt')) %>%
  left_join(schema) %>%
  arrange(Q) %>%
  rename(ensemble_gene = gene,
         gene = symbol)

asd_genes <- readxl::read_excel(paste0(DATA_DIR, 'public_datasets//rare_variants/78 FDR 5% ASD genes from Satterstrom et al.xlsx')) %>%
  dplyr::select(gene)

for (STUDY in c('schema_genes', 'asd_genes')) {
  
  cat('\n\nRunning wilcoxon test for:', STUDY)
  STUDY_DF <- get(STUDY)
  STUDY_GENES <- STUDY_DF$gene
  
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
        mutate(study_status = ifelse(gene %in% STUDY_GENES, 'in_study', 'not_in_study'))
      
      in_study_scores <- specificity_cell %>%
        filter(grepl('^in_study', study_status)) %>%
        pull(cell_scores)
      
      not_in_study_scores <- specificity_cell %>%
        filter(grepl('^not_in_study', study_status)) %>%
        pull(cell_scores)
      
      wilcox_result <- wilcox.test(in_study_scores, not_in_study_scores, paired = FALSE, data = specificity_cell)
      
      if (exists('wilcoxon_df')) {
        
        wilcox_new_result <- as.data.frame(t(c(CELL_TYPE, wilcox_result$p.value)))
        colnames(wilcox_new_result) <- c('cell_type', 'P')
        wilcoxon_df <- rbind(wilcoxon_df, wilcox_new_result)
        
      } else {
        
        wilcoxon_df <- data.frame()
        wilcoxon_df <- rbind(wilcoxon_df, as.data.frame(t(c(CELL_TYPE, wilcox_result$p.value))))
        colnames(wilcoxon_df) <- c('cell_type', 'P')
        
      }
      
      cat(paste0('\nNumber of ', STUDY, ' in sample:'), sum(specificity_cell$study_status == 'in_study'))
      
    }
      
    
  }
  
  wilcoxon_df <- wilcoxon_df %>% 
    arrange(P) %>%
    mutate(BF = p.adjust(P, 'bonferroni', length(P))) %>%
    mutate(FDR = p.adjust(P, 'BH', length(P)))
  
  assign(paste0(STUDY, '_wilcoxon_df'), wilcoxon_df)
  write_tsv(wilcoxon_df, paste0(FIG_DIR, STUDY, '_wilcoxon.txt'))
  rm(wilcoxon_df)
  
}





#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

asd_genes_wilcoxon_df %>% 
  arrange(P) %>%
  mutate(BF = p.adjust(P, 'bonferroni', length(P))) %>%
  mutate(FDR = p.adjust(P, 'BH', length(P)))
