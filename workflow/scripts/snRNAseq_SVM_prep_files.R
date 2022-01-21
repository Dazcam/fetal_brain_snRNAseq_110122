# -------------------------------------------------------------------------------------
#
#   snRNAseq - prepare files for SVM analysis
#   
# -------------------------------------------------------------------------------------

# Load packages  ----------------------------------------------------------------------
library(Seurat)
library(reshape2)
library(tidyverse)

# Set variables
DATA_DIR <- "../resources/R_objects/"
OUT_DIR <- "../resources/sheets/"
REGIONS <- c("cer", "hip", "pfc", "tha", "wge")

# Pre SVM - Prep files for run cluster (using Python)  --------------------------------
for (REGION in REGIONS) {
  
  # Get seurat object per brain region 
  cat(paste0('\nLoading ', REGION, ' seurat object ... \n'))
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  
  # Create raw count gene matrix .csv - needs to be cells x genes
  cat('Creating raw count matrix ... ')
  raw_counts <- as.data.frame(t(as.matrix(seurat.obj@assays$RNA@counts)))
  write_csv(raw_counts, 
            paste0(OUT_DIR, "snRNAseq_raw_counts_", REGION, "_final.csv"), 
            escape = "double")
  
  # Create cluster labels .csv
  cat('\nCreating cluster labels ... ')
  cluster_labels <- as.data.frame(as.vector(seurat.obj$cellIDs))
  colnames(cluster_labels) <- "Class"
  write_csv(cluster_labels, 
            paste0(OUT_DIR, "snRNAseq_cluster_labels_", REGION, "_for_SVM.csv"), 
            escape = "double")
  
  # Randomise cluster label rows for negative control and generate .csv
  cat('\nCreating randomised cluster labels ... \n')
  set.seed(42)
  rows <- sample(nrow(cluster_labels))
  cluster_labels_rand <- as.data.frame(cluster_labels[rows,])
  colnames(cluster_labels_rand) <- "Class"
  write_csv(cluster_labels_rand, 
            paste0(OUT_DIR, "snRNAseq_cluster_labels_rand_", REGION, "_for_SVM.csv"), 
            escape = "double")

  
}

cat('\nDone.')

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
