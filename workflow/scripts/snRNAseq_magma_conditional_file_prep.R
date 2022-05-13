#--------------------------------------------------------------------------------------
#
#    Magma conditional analysis - file prep
#    
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(ggVennDiagram)
library(readxl)
library(biomaRt)
library(cowplot)
library(org.Hs.eg.db)

##  Info  -----------------------------------------------------------------------------
# 1. Convert fetal top decile genes (from q10_gene_lists_for_LDSC) from HGNC IDs to entrez IDs 
# 2. Obtain top decile genes for Skene cell types
# 3. Convert Skene top decile genes from MGI IDs to HGNC IDs 
# 4. Convert Skene top decile genes from HGNC IDs to entrez IDs 
# 5. These were then transposed to get a list of genes in a single line per cell type
# 6. These were then bound by row in bash (using >>) for input to MAGMA
# 7. These gene lists were then used for magma conditional analysis
#    See: snRNAseq_magma_conditional.smk

# TODO Added Bryois cell types manually need to add them bioinformatically!!!!

##  Set variables  --------------------------------------------------------------------
cat('\nDefining variables for section 1 ... \n')
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/gene_lists/q10_gene_lists/'
ENTREZ_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/results/gene_lists/ENTREZ/'
dir.create(ENTREZ_DIR)
SKENE_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/public_datasets/skene/'
CELL_TYPES <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'FC-ExN-5', 'FC-InN-1', 
                'GE-InN-1', 'GE-InN-2', 'Hipp-ExN-3', 'Hipp-ExN-5', 'Thal-ExN-1',
                'Thal-ExN-3')
ALL_CELL_TYPES <- c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 'FC-ExN-5', 'FC-InN-1', 
                    'GE-InN-1', 'GE-InN-2', 'Hipp-ExN-3', 'Hipp-ExN-5', 'Thal-ExN-1',
                    'Thal-ExN-3', 'skene_InN', 'skene_MSN', 'skene_CA1', 'skene_SS')

## Load and prep data  ----------------------------------------------------------------
# Gene lists from Q10s for our significant cell types in snRNAseq data
cat('Load the fetal brain cell type Q10 gene lists ... \n')
for (CELL_TYPE in CELL_TYPES) {
  
  GENE_LIST <- read_csv(paste0(DATA_DIR, CELL_TYPE, "_Q10_genes.tsv"), col_names = FALSE)
  CELL_TYPE <- gsub("-", "_", CELL_TYPE)
  assign(CELL_TYPE, GENE_LIST)
  
}

# Convert HGNC IDs to entrez IDS
for (CELL_TYPE in CELL_TYPES) {
  
  CELL_TYPE_EDIT <- gsub("-", "_", CELL_TYPE)
  GENE_LIST <- get(CELL_TYPE_EDIT)
  
  cat(paste0('\n\nRunning conversion for: ', CELL_TYPE_EDIT, '\n'))
  cat(paste0('Total genes before converting IDs: ', nrow(GENE_LIST), '\n'))
  
  # Convert gene IDs 
  cat('Converting human gene IDs using BiomaRt ... \n')
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ENTREZ_GENES = getBM(attributes = c("hgnc_symbol", "entrezgene_id"), 
                       filters = "hgnc_symbol", 
                       values = GENE_LIST, 
                       bmHeader = T, 
                       mart = mart)
  ENTREZ_GENES_noNA <- ENTREZ_GENES %>% drop_na()
  ENTREZ_GENES_UNIQUE <- as.data.frame(unique(ENTREZ_GENES_noNA[, 2]))
  colnames(ENTREZ_GENES_UNIQUE) <- CELL_TYPE_EDIT
  
  cat(paste0('Total genes after converting IDs: ', nrow(ENTREZ_GENES), '\n'))
  cat(paste0('Total genes after removing NAs: ', nrow(ENTREZ_GENES_noNA), '\n'))
  cat(paste0('Total genes after checking uniqueness (so final count): ', nrow(ENTREZ_GENES_UNIQUE), '\n'))
  
  # Add annotation to gene list
  ENTREZ_GENES_ANNOT <- ENTREZ_GENES_UNIQUE
  ENTREZ_GENES_ANNOT$gene_list <- rep(paste0(CELL_TYPE_EDIT), nrow(ENTREZ_GENES_ANNOT))
  
  # Create transposed df
  ENTREZ_GENES_TRANS <- as.data.frame(t(ENTREZ_GENES_UNIQUE))
  
  assign(paste0(CELL_TYPE_EDIT, '_entrezID'), ENTREZ_GENES_UNIQUE)
  assign(paste0(CELL_TYPE_EDIT, '_entrezID_annot'), ENTREZ_GENES_ANNOT)
  assign(paste0(CELL_TYPE_EDIT, '_entrezID_trans'), ENTREZ_GENES_TRANS)
  
  write_tsv(ENTREZ_GENES_UNIQUE, paste0(ENTREZ_DIR, CELL_TYPE_EDIT, '_entrez_gene_list.tsv'))
  write_tsv(ENTREZ_GENES_ANNOT, paste0(ENTREZ_DIR, CELL_TYPE_EDIT, '_entrez_gene_list_annot.tsv'), col_names = FALSE)
  write.table(ENTREZ_GENES_TRANS, paste0(ENTREZ_DIR, CELL_TYPE_EDIT, '_entrez_gene_list_trans.tsv'), col.names = FALSE, 
              quote = FALSE, sep = '\t', row.names = TRUE)
  
}

# Load Skene gene specificity data
cat('Load the Skene cell specificity scores ... \n')
skene_genes <- readxl::read_xlsx(paste0(SKENE_DIR, 'skene_gene_specificity_scores.xlsx'))
colnames(skene_genes) <- c('gene', 'skene_InN', 'skene_MSN', 'skene_CA1', 'skene_SS')


# Extract top 10% most sig genes in each Skene cell-type and convert gene IDs to human
for (CELL_TYPE in c('skene_InN', 'skene_MSN', 'skene_CA1', 'skene_SS')) {
  
  # Pull out Q10 vector
  cat(paste0('Extracting genes that have highest expression specificity in: ', CELL_TYPE, '\n'))
  GENES <- skene_genes %>%
    dplyr::select(gene, .data[[CELL_TYPE]]) %>%
    arrange(desc(.data[[CELL_TYPE]])) %>%
    top_frac(.10) %>%
    pull(gene) 
  
  cat(paste0('Total genes before converting IDs: ', length(GENES), '\n'))
  
  # Convert mouse gene IDs to human
  cat('Convert gene IDs from mouse to human using BiomaRt ... \n')
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), 
                   filters = "mgi_symbol", 
                   values = GENES, 
                   mart = mouse, 
                   attributesL = c("hgnc_symbol"), 
                   martL = human, 
                   uniqueRows=T)
  humanx <- as.data.frame(unique(genesV2[, 2]))
  colnames(humanx) <- "X1"
  
  cat(paste0('Total genes after converting IDs: ', dim(humanx)[1], '\n'))
  
  assign(CELL_TYPE, humanx)
  
}

# Convert Skene HGNC IDs to entrez IDS
for (CELL_TYPE in c('skene_InN', 'skene_MSN', 'skene_CA1', 'skene_SS')) {
  
  GENE_LIST <- get(CELL_TYPE)
  
  cat(paste0('\n\nRunning conversion for: ', CELL_TYPE, '\n'))
  cat(paste0('Total genes before converting IDs: ', nrow(GENE_LIST), '\n'))
  
  # Convert gene IDs 
  cat('Converting human gene IDs using BiomaRt ... \n')
  mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  ENTREZ_GENES = getBM(attributes = c("hgnc_symbol", "entrezgene_id"), 
                       filters = "hgnc_symbol", 
                       values = GENE_LIST, 
                       bmHeader = T, 
                       mart = mart)
  ENTREZ_GENES_noNA <- ENTREZ_GENES %>% drop_na()
  ENTREZ_GENES_UNIQUE <- as.data.frame(unique(ENTREZ_GENES_noNA[, 2]))
  colnames(ENTREZ_GENES_UNIQUE) <- CELL_TYPE
  
  cat(paste0('Total genes after converting IDs: ', nrow(ENTREZ_GENES), '\n'))
  cat(paste0('Total genes after removing NAs: ', nrow(ENTREZ_GENES_noNA), '\n'))
  cat(paste0('Total genes after checking uniqueness (so final count): ', nrow(ENTREZ_GENES_UNIQUE), '\n'))
  
  # Add annotation to gene list
  ENTREZ_GENES_ANNOT <- ENTREZ_GENES_UNIQUE
  ENTREZ_GENES_ANNOT$gene_list <- rep(paste0(CELL_TYPE), nrow(ENTREZ_GENES_ANNOT))
  
  # Create transposed df
  ENTREZ_GENES_TRANS <- as.data.frame(t(ENTREZ_GENES_UNIQUE))
  
  assign(paste0(CELL_TYPE, '_entrezID'), ENTREZ_GENES_UNIQUE)
  assign(paste0(CELL_TYPE, '_entrezID_annot'), ENTREZ_GENES_ANNOT)
  assign(paste0(CELL_TYPE, '_entrezID_trans'), ENTREZ_GENES_TRANS)
  
  write_tsv(ENTREZ_GENES_UNIQUE, paste0(ENTREZ_DIR, CELL_TYPE, '_entrez_gene_list.tsv'))
  write_tsv(ENTREZ_GENES_ANNOT, paste0(ENTREZ_DIR, CELL_TYPE, '_entrez_gene_list_annot.tsv'), col_names = FALSE)
  write.table(ENTREZ_GENES_TRANS, paste0(ENTREZ_DIR, CELL_TYPE, '_entrez_gene_list_trans.tsv'), col.names = FALSE, 
              quote = FALSE, sep = '\t', row.names = TRUE)
  
}

# Cat single line gene lists to one file for conditional analysis
for (CELL_TYPE in ALL_CELL_TYPES) {
  
  if (CELL_TYPE == 'FC-ExN-2') {
    
    #file.create(paste0(DATA_DIR, 'ALL_SIG_AND_SKENE_entrez_gene_list.tsv'))
    CELL_TYPE_EDIT <- gsub("-", "_", CELL_TYPE)
    cat(readLines(paste0(ENTREZ_DIR, CELL_TYPE_EDIT, '_entrez_gene_list.tsv')), '\n', 
        file = paste0(DATA_DIR, 'ALL_SIG_AND_SKENE_entrez_gene_list.tsv'))
    
    
  } else {
    
    # Append rest of cell types to file
    CELL_TYPE_EDIT <- gsub("-", "_", CELL_TYPE)
    cat(readLines(paste0(ENTREZ_DIR, CELL_TYPE_EDIT, '_entrez_gene_list.tsv')), '\n', 
        file = paste0(DATA_DIR, 'ALL_SIG_AND_SKENE_entrez_gene_list.tsv'), 
        append = TRUE)
    
  }
  
  # Append Bryois adult human cell types to gene list
  cat(readLines(paste0(SKENE_DIR, 'bryois_human_adult_top10pc_space_delim.txt')), '\n', 
      file = paste0(DATA_DIR, 'ALL_SIG_AND_SKENE_entrez_gene_list.tsv'), 
      append = TRUE)
 
   
}

 
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------