#--------------------------------------------------------------------------------------
#
#    snRNAseq stacked violin plots - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Figures for stacked violin plots (may be better done in python scanpy?)
#  Extended data figures 1-5

##  Load Packages  --------------------------------------------------------------------
library(cowplot)
library(tidyverse)
library(Seurat)

## Set variables  ---------------------------------------------------------------------
DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/"
FIG_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/figures/"
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')

## Load Data --------------------------------------------------------------------------
# Seurat objects  ----
for (REGION in REGIONS) { 
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  
}

# Set features
fc_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                 "MKI67", "C3", "ITM2A", "SST", "CALB2", 
                 "SCGN", "TLE3", "FEZF2", "CRYM", "LHX2")
ge_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                 "MKI67", "C3", "ITM2A", "LHX6", "SIX3", 
                 "PROX1", "TSHZ1", "DLX1", "SCGN")
hip_features <- c("NEUROD1", "GRIK4", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "ADARB2",
                  "GAD2", "TNC", "PROX1", "RELN", "LHX6")
tha_features <- c("EOMES", "GLI3", "OLIG1", "MKI67", "C3", 
                  "ITM2A", "SLC1A6", "LHX9", "TNC", "GAD2", 
                  "PAX6", "SLC17A6")
cer_features <- c("GAD1", "EOMES", "GLI3", "OLIG1", "MKI67", 
                  "C3", "ITM2A", "CA8", "ITPR1", "RBFOX3", 
                  "RELN")

# Set colours
fc_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                '#3CBB75FF', '#D078FF', '#F58231', '#CCCCCC', '#FDE725FF', 
                '#FF5959', '#FF5959')
ge_colours <- c('#DCBEFF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#FF5959', '#FF5959',
                '#FF5959')
hip_colours <- c('#B200ED', '#B200ED', '#DCBEFF', '#9A6324', '#CEE5FD', 
                 '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#F58231', 
                 '#FDE725FF', '#FF5959', '#FF5959', '#FF5959')
tha_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                 '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', "#D078FF", '#F58231', 
                 '#FDE725FF', '#FF5959', '#FF5959', '#FF5959', '#FF5959', 
                 '#FF5959', '#FF5959', '#FF5959')
cer_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                 '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF',
                 "#D078FF", '#F58231', '#FDE725FF', '#FF5959', '#FF5959',
                 '#FF5959')



# Order Idents
Idents(seurat.pfc) <- factor(x = Idents(seurat.pfc), levels = sort(levels(seurat.pfc)))
Idents(seurat.wge) <- factor(x = Idents(seurat.wge), levels = sort(levels(seurat.wge)))
Idents(seurat.cer) <- factor(x = Idents(seurat.cer), levels = sort(levels(seurat.cer)))
Idents(seurat.tha) <- factor(x = Idents(seurat.tha), levels = sort(levels(seurat.tha)))
Idents(seurat.hip) <- factor(x = Idents(seurat.hip), levels = sort(levels(seurat.hip)))


# Plot
FC_plot <- VlnPlot(seurat.pfc, fc_features, stack = TRUE, flip = TRUE, 
                   cols = fc_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        axis.text.x  = element_text(colour = "#000000", size = 16),
        axis.text.y  = element_text(colour = "#000000", size = 16)) +
  xlab('Cell type') +
  ggtitle("Frontal Cortex")

GE_plot <- VlnPlot(seurat.wge, ge_features, stack = TRUE, flip = TRUE, 
                   cols = ge_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        axis.text.x  = element_text(colour = "#000000", size = 16),
        axis.text.y  = element_text(colour = "#000000", size = 16)) +
  xlab('Cell type') +
  ggtitle("Ganglionic Eminence")

Hip_plot <- VlnPlot(seurat.hip, hip_features, stack = TRUE, flip = TRUE, 
                    cols = hip_colours, same.y.lims = TRUE, fill.by = 'ident')  +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        axis.text.x  = element_text(colour = "#000000", size = 16),
        axis.text.y  = element_text(colour = "#000000", size = 16)) +
  xlab('Cell type') +
  ggtitle("Hippocampus")

Tha_plot <- VlnPlot(seurat.tha, tha_features, stack = TRUE, flip = TRUE, 
                    cols = tha_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        axis.text.x  = element_text(colour = "#000000", size = 16),
        axis.text.y  = element_text(colour = "#000000", size = 16)) +
  xlab('Cell type') +
  ggtitle("Thalamus")

Cer_plot <- VlnPlot(seurat.cer, cer_features, stack = TRUE, flip = TRUE, 
                    cols = cer_colours, same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        axis.text.x  = element_text(colour = "#000000", size = 16),
        axis.text.y  = element_text(colour = "#000000", size = 16)) +
  xlab('Cell type') +
  ggtitle("Cerebellum")


# Save
tiff(paste0(FIG_DIR, "Supplementary_figure_1.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
FC_plot 
dev.off()

tiff(paste0(FIG_DIR, "Supplementary_figure_2.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
GE_plot 
dev.off()

tiff(paste0(FIG_DIR, "Supplementary_figure_3.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
Hip_plot 
dev.off()

tiff(paste0(FIG_DIR, "Supplementary_figure_4.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
Tha_plot 
dev.off()

tiff(paste0(FIG_DIR, "Supplementary_figure_5.tiff"), height = 30, width = 40, units='cm', 
     compression = "lzw", res = 300)
Cer_plot 
dev.off()


# Jpegs
jpeg(paste0(FIG_DIR, "Supplementary_figure_1.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
FC_plot 
dev.off()

jpeg(paste0(FIG_DIR, "Supplementary_figure_2.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
GE_plot 
dev.off()

jpeg(paste0(FIG_DIR, "Supplementary_figure_3.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
Hip_plot 
dev.off()

jpeg(paste0(FIG_DIR, "Supplementary_figure_4.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
Tha_plot 
dev.off()

jpeg(paste0(FIG_DIR, "Supplementary_figure_5.jpg"), width = 1440, height = 960, 
     units = "px", pointsize = 12, quality = 150)
Cer_plot 
dev.off()



#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------