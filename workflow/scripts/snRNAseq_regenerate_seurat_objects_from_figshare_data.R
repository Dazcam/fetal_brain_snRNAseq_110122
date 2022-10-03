#--------------------------------------------------------------------------------------
#
#   Cameron et al. 2022 - Code to regenerate snRNAseq Seurat objects
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("Seurat", "tidyverse", "cowplot", "data.table"))

## Load variables  --------------------------------------------------------------------
DATA_DIR <- '~/Downloads/' # Add dir for figshare files 
REGIONS <- c('Cer', 'FC', 'GE', 'Hipp', 'Thal') 
FIG_SHARE_PATH <- 'https://figshare.com/ndownloader/files/'
FIG_SHARE_IDs <- c('Cer_metadata' = 36712263, 'Cer_raw_count_gEX_matrix' = 36712266, 
                   'FC_metadata' = 36712269, 'FC_raw_count_gEX_matrix' = 36712272, 
                   'GE_metadata' = 36712275, 'GE_raw_count_gEX_matrix' = 36712278, 
                   'Hipp_metadata' = 36712281, 'Hipp_raw_count_gEX_matrix' = 36712284, 
                   'Thal_metadata' = 36712287, 'Thal_raw_count_gEX_matrix' = 36712290)
VAR_FEAT  <- 2000       # No. variable features to use for FindVariableFeatures()
PCA_UPPER_DIM <- 17     # No. PCA dimensions to be considered for dim reduction
RESOLUTION <- 0.5       # Variable for cluster granularity

# Download and unzip files from figshare  ---------------------------------------------
for (FILE_ID in 1:length(FIG_SHARE_IDs)) {
  
  cat('\n\nDownloading:', names(FIG_SHARE_IDs[FILE_ID]), '...')
  system(paste0("wget ", FIG_SHARE_PATH, FIG_SHARE_IDs[[FILE_ID]], " -P ",  DATA_DIR), 
         ignore.stdout = T, ignore.stderr = T)
  
  cat('\nRenaming:', names(FIG_SHARE_IDs[FILE_ID]), '...')
  file.rename(paste0(DATA_DIR, FIG_SHARE_IDs[[FILE_ID]]), 
              paste0(DATA_DIR, 'cameron_2022_snRNAseq_', names(FIG_SHARE_IDs[FILE_ID]), '.txt.gz'))
  
  cat('\nUnzipping:', names(FIG_SHARE_IDs[FILE_ID]), '...')
  system(paste0("gunzip ", paste0(DATA_DIR, 'cameron_2022_snRNAseq_', names(FIG_SHARE_IDs[FILE_ID]), '.txt.gz')))
  
  cat('\nDone.')
  
}

# Run files through Seurat  -----------------------------------------------------------
for (REGION in REGIONS) {
  
   # Load Data
  cat('\n\nLoading data for:', REGION, '...\n\n')
  mat <- fread(paste0(DATA_DIR, 'cameron_2022_snRNAseq_', REGION, '_raw_count_gEX_matrix.txt'))
  meta <- as.data.frame(read_tsv(paste0(DATA_DIR, 'cameron_2022_snRNAseq_', REGION, '_metadata.txt')))
  
  # Prep data for Seurat
  cat('\nPrepping data for:', REGION, '...\n')
  rownames(meta) <- meta[,1]
  meta <- meta[, -1]
  mat_genes <- mat$V1
  mat_no_genes <- as.matrix(mat[, -1])
  rownames(mat_no_genes) <- mat_genes
  
  # Run basic Seurat pipeline
  cat('\nRunning basic Seurat pipeline for:', REGION, '...\n\n')
  seurat.obj <- CreateSeuratObject(counts = mat_no_genes, meta.data = meta)
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = VAR_FEAT) 
  seurat.obj <- ScaleData(seurat.obj) 
  seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj)) 
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:PCA_UPPER_DIM) 
  seurat.obj <- FindClusters(seurat.obj, resolution = RESOLUTION) 
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:PCA_UPPER_DIM)
  
  # Create UMAPs
  cat('\nGenerating UMAP:', REGION, '...')
  seurat.umap <- DimPlot(seurat.obj, pt.size = 0.2, reduction = "umap", group.by = 'cellIDs')
  
  assign(paste0('seurat.', REGION), seurat.obj)
  assign(paste0('seurat.', REGION, '_umap'), seurat.umap)
  
  cat('\nDone.\n\n')
  
}


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

