# -------------------------------------------------------------------------------------
#
#    Cell IDs - Adds final cell IDs to the regional clusters
#
# -------------------------------------------------------------------------------------

## Explanation  -----------------------------------------------------------------------

#  Finalising IDs and colour scheme for clusters 

## Load packages  ---------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(cowplot)


## Set variables  ---------------------------------------------------------------------
DATA_DIR <- '~/Desktop/single_cell/scRNAseq/batch2_CR5_200121/r_objects/'
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
REGIONS_PASS1 <- c('cer', 'hip', 'pfc', 'tha')

## Load Data --------------------------------------------------------------------------
# Seurat objects  ----

for (REGION in REGIONS_PASS1) { 
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'pass1/seurat.', REGION, '.pass1.rds'))
  
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  
}

# WGE in pass 2
seurat.wge <- readRDS(paste0(DATA_DIR, 'pass2/seurat.wge.pass2.rds'))

# Cell IDs  ----

for (REGION in REGIONS) { 
  
  cellIDs <- read.table(paste0(DATA_DIR, 'final/', REGION, '_final_clusterIDs_110122.txt'))
  cellIDs <- as.vector(cellIDs[, 2])
  
  assign(paste0(REGION, '_cellIDs'), cellIDs, .GlobalEnv)
  
}


## Ensure cluster IDs are up to date --------------------------------------------------
seurat.cer$pfc_idents <- Idents(seurat.cer)
seurat.hip$hip_idents <- Idents(seurat.hip)
seurat.pfc$pfc_idents <- Idents(seurat.pfc)
seurat.tha$tha_idents <- Idents(seurat.tha)
seurat.wge$wge_idents <- Idents(seurat.wge)

## Check plots
for (REGION in REGIONS) { 
  
  seurat.obj <- get(paste0('seurat.', REGION))  
  cluster_plot <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, 
                          pt.size = 0.2,repel = TRUE) + ggtitle(NULL) + NoLegend()
  assign(paste0(REGION, '_cluster_plot'), cluster_plot, .GlobalEnv)
  
}


## Save cell IDs as Idents  -----------------------------------------------------------
names(cer_cellIDs) <- levels(seurat.cer)
seurat.cer <- RenameIdents(seurat.cer, cer_cellIDs)
seurat.cer$cellIDs <- Idents(seurat.cer)

names(hip_cellIDs) <- levels(seurat.hip)
seurat.hip <- RenameIdents(seurat.hip, hip_cellIDs)
seurat.hip$cellIDs <- Idents(seurat.hip)

names(pfc_cellIDs) <- levels(seurat.pfc)
seurat.pfc <- RenameIdents(seurat.pfc, pfc_cellIDs)
seurat.pfc$cellIDs <- Idents(seurat.pfc)

names(tha_cellIDs) <- levels(seurat.tha)
seurat.tha <- RenameIdents(seurat.tha, tha_cellIDs)
seurat.tha$cellIDs <- Idents(seurat.tha)

names(wge_cellIDs) <- levels(seurat.wge)
seurat.wge <- RenameIdents(seurat.wge, wge_cellIDs)
seurat.wge$cellIDs <- Idents(seurat.wge)

## Check plots with final labels  -----------------------------------------------------
for (REGION in REGIONS) { 
  
  seurat.obj <- get(paste0('seurat.', REGION))  
  cluster_plot <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, 
                          pt.size = 0.2,repel = TRUE) + ggtitle(NULL)
  
  assign(paste0(REGION, '_cluster_final_cellIDs_plot'), cluster_plot, .GlobalEnv)
  
}

## Save RDS files  --------------------------------------------------------------------
for (REGION in REGIONS) { 
  
  seurat.obj <- get(paste0('seurat.', REGION)) 
  saveRDS(seurat.obj, paste0(DATA_DIR, 'final/seurat.', REGION, '.final.rds'))
  
}

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
