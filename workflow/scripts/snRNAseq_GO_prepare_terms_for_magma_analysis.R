# -------------------------------------------------------------------------------------
#
#    snRNAseq - prepare gene lists for GO Term MAGMA analysis
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Intersect final list of GO Terms in GO figure with Terms in GO output for each sig. cell type
#  2. Remove rows with NA - i.e. Terms not found in GO output of that cell type
#  3. Munge data to obtain a list of genes 1 line per cell type per term and export to single file
#  4. These gene lists were then used for GO Term magma analysis
#  See: snRNAseq_magma_GO.smk


##  Load packages  --------------------------------------------------------------------
library(biomaRt)
library(tidyverse)
library(readxl)
library(janitor)

##  Initialise variables  -------------------------------------------------------------
SIG_CELLS <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'GE-InN-2', 'Hipp-ExN-3', 'Hipp-ExN-5')
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/GO/"
RESOURCES_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/resources/sheets/"

##  Load GO results data  -------------------------------------------------------------
for (CELL_TYPE in SIG_CELLS) {
  
  GENES <- read_delim(paste0(DATA_DIR, CELL_TYPE, ' GO.txt'), 
                      delim = '\t', escape_double = FALSE, col_names = TRUE,
                      trim_ws = TRUE)
  assign(paste0(CELL_TYPE, '_go'), GENES)
  
}

##  Extract GO terms from GO results summary file  ------------------------------------
GO_TERMS <- read_excel(paste0(RESOURCES_DIR, 'FinalGOterms_to_plot_for_figure6.xlsx')) %>%
  dplyr::select(...1) %>%
  janitor::row_to_names(row_number = 1) 

for (CELL_TYPE in SIG_CELLS) {
  
  CELL_TYPE_EDIT <- gsub("-", "_", CELL_TYPE)
  
  # Extract only final 26 GO terms in summary file (or GO plt) from each cell type
  INTERSECT_GENES <- left_join(GO_TERMS, get(paste0(CELL_TYPE, '_go')), 
                               by = 'Term')
  
  # Drop NAs and create DF - 1 column per term 
  INTERSECT_GENES_FILT <- INTERSECT_GENES %>% drop_na() %>%
    select(Term, Genes) %>%
    group_by(Term) %>%
    dplyr::mutate(i1 = row_number()) %>%
    spread(Term, Genes) %>%
    select(-i1)
  
  cat('\n\nChecks for: ', CELL_TYPE)
  cat('\nTerm count: ', nrow(INTERSECT_GENES))
  cat('\nGene count after NA filter: ', nrow(INTERSECT_GENES_FILT))
  
  # Extract list of genes 1 line per cell type per term and export to single file
  for (i in 1:ncol(INTERSECT_GENES_FILT)) {
    
    TERM <- colnames(INTERSECT_GENES_FILT)[i]
    GENES <- strsplit(as.character(INTERSECT_GENES_FILT[1, i]), ", ")[[1]]
    
    if (file.exists(paste0(DATA_DIR, 'GO_term_genes_for_magma.txt'))) {
      
      cat('\n\nAPPEND\n\n')
      
      GENE_LIST <- paste0(CELL_TYPE_EDIT, '_', TERM, ' ', paste(GENES, collapse = ' '))
      cat('\n\n', GENE_LIST, '\n')
      cat(GENE_LIST, '\n', file = paste0(DATA_DIR, 'GO_term_genes_for_magma.txt'), 
          append = TRUE)
      
    } else {
      
      cat('\n\nCREATE FILE\n\n')
      
      GENE_LIST <- paste0(CELL_TYPE_EDIT, '_', TERM, ' ', paste(GENES, collapse = ' '))
      cat('\n', GENE_LIST, '\n')
      cat(GENE_LIST, '\n', file = paste0(DATA_DIR, 'GO_term_genes_for_magma.txt'))
      
    }
    
  }
  
  
}
  
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
  



