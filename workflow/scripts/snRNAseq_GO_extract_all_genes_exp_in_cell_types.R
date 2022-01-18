#--------------------------------------------------------------------------------------
#
#    snRNAseq - Extract all genes expressed in significant fetal cell types
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  1. Generate average expression values for each cluster
#  2. Remove genes that are not expressed from each of the 11 cell types
#  3. Convert cell type gene IDs to Entrez IDs got Go analysis
#  4. Nick ran GO analyses here:

# Input for GO:

# 

# FC-InN-4 added here as it was sig. in LDSC analysis

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(biomaRt)

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/"
OUT_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/gene_lists/Ave_expression/"
dir.create(OUT_DIR, recursive = TRUE)
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
CELL_TYPES <- c('FC_ExN_2', 'FC_ExN_3', 'FC_ExN_4', 'FC_ExN_5', 'FC_InN_1', 
                'FC_InN_4', 'GE_InN_1', 'GE_InN_2', 'Hipp_ExN_3', 'Hipp_ExN_5', 
                'Thal_ExN_1', 'Thal_ExN_3')


## Load Data --------------------------------------------------------------------------
# Seurat objects and calulate average expression for each region
for (REGION in REGIONS) { 
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  seurat.ave_exp <- AverageExpression(object = seurat.obj)
  
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  assign(paste0(REGION, '_avExp'), seurat.ave_exp, .GlobalEnv)
  
}


## Extract gene list for each cell type   ---------------------------------------------
FC_ExN_2_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-2`) %>%
  filter(`FC-ExN-2` > 0) %>%
  arrange(desc(`FC-ExN-2`))

FC_ExN_3_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-3`) %>%
  filter(`FC-ExN-3` > 0) %>%
  arrange

FC_ExN_4_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-4`) %>%
  filter(`FC-ExN-4` > 0) %>%
  arrange(desc(`FC-ExN-4`))

FC_ExN_5_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-ExN-5`) %>%
  filter(`FC-ExN-5` > 0) %>%
  arrange(desc(`FC-ExN-5`))

FC_InN_1_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-InN-1`) %>%
  filter(`FC-InN-1` > 0) %>%
  arrange(desc(`FC-InN-1`))

FC_InN_4_gene_list <- as.data.frame(pfc_avExp$RNA) %>%
  dplyr::select(`FC-InN-4`) %>%
  filter(`FC-InN-4` > 0) %>%
  arrange(desc(`FC-InN-4`))

GE_InN_1_gene_list <- as.data.frame(wge_avExp$RNA) %>%
  dplyr::select(`GE-InN-1`) %>%
  filter(`GE-InN-1` > 0) %>%
  arrange(desc(`GE-InN-1`))

GE_InN_2_gene_list <- as.data.frame(wge_avExp$RNA) %>%
  dplyr::select(`GE-InN-2`) %>%
  filter(`GE-InN-2` > 0) %>%
  arrange(desc(`GE-InN-2`))

Hipp_ExN_3_gene_list <- as.data.frame(hip_avExp$RNA) %>%
  dplyr::select(`Hipp-ExN-3`) %>%
  filter(`Hipp-ExN-3` > 0) %>%
  arrange(desc(`Hipp-ExN-3`))

Hipp_ExN_5_gene_list <- as.data.frame(hip_avExp$RNA) %>%
  dplyr::select(`Hipp-ExN-5`) %>%
  filter(`Hipp-ExN-5` > 0) %>%
  arrange(desc(`Hipp-ExN-5`))

Thal_ExN_1_gene_list <- as.data.frame(tha_avExp$RNA) %>%
  dplyr::select(`Thal-ExN-1`) %>%
  filter(`Thal-ExN-1` > 0) %>%
  arrange(desc(`Thal-ExN-1`))

Thal_ExN_3_gene_list <- as.data.frame(tha_avExp$RNA) %>%
  dplyr::select(`Thal-ExN-3`) %>%
  filter(`Thal-ExN-3` > 0) %>%
  arrange(desc(`Thal-ExN-3`))

## Convert cell type gene IDs to Entrez   ---------------------------------------------
for (CELL_TYPE in CELL_TYPES) {
  
  GENE_FILE <- get(paste0(CELL_TYPE, '_gene_list'))
  GENE_LIST <- rownames(GENE_FILE)
  
  
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
  
  assign(paste0(CELL_TYPE, '_entrezID'), ENTREZ_GENES_UNIQUE)
  
  write.table(ENTREZ_GENES_UNIQUE, paste0(OUT_DIR, CELL_TYPE, '_all_genes_expresssed.tsv'),
              quote = FALSE, row.names = FALSE,  sep = '\t', col.names = FALSE)
  
}


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------