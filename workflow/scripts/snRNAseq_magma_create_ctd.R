# -------------------------------------------------------------------------------------
#
#    MAGMA Celltyping - create CTD object
#
# -------------------------------------------------------------------------------------

##  Resources  -------------------------------------------------------------------------
   
    #  Paper - https://www.nature.com/articles/s41588-018-0129-5
    #  Github - https://github.com/NathanSkene/MAGMA_Celltyping

    # ECWE - for creating the ctd object
      # Github - https://github.com/neurogenomics/EWCE/
      # Vingette - https://nathanskene.github.io/EWCE/articles/EWCE.html

## Requirements  -------------------------------------------------------------------------

    # Required on Hawk before opening R
      # module load libgit2/1.1.0
      # module load R/4.0.3

    # GWAS:
      # SNP, CHR, BP as first three columns.
      # GWAS has one of: Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT
      # GWAS has all of: SNP, CHR, BP, P, A1, A2

    # Net access:
      # Uses BiomaRt for annotations need to run SLURM jobs locally

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

## Load packages  ---------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(Seurat) 
library(devtools)
library(EWCE) # using srun can't access Biomart
library(MAGMA.Celltyping) # Note the "." instead of "_" - using srun can't see magma executable
library(argparser)
library(reshape2)
library(edgeR)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read brain region for magma ... \n")
p <- add_argument(p, "region", help = "No brain region provided")
p <- add_argument(p, "seurat_obj", help = "No Seurat obj provided")
p <- add_argument(p, "MHC_gene_list", help = "No gene list (intersecting MHC region and Seurat object) provided")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
REGION <- args$region
SEURAT_OBJ <- args$seurat_obj
MHC_GENES <- args$MHC_gene_list
CTD_DIR <- '/scratch/c.c1477909/fetal_brain_snRNAseq_110122/resources/ctd_objects'

##  Load seurat objects  ------------------------------------------------------------
cat(paste0("\nLoading seurat object and MHC gene list for ", toupper(REGION), " ... \n"))
seurat.obj <- readRDS(SEURAT_OBJ)
intersect.genes <- readRDS(MHC_GENES)

# Create raw count gene matrix .csv - needs to be cells x genes
cat("Creating count by gene matrix ...\n")
raw_counts <- as.matrix(seurat.obj@assays$RNA@counts)

# Scale counts
#raw_counts_scaled = Matrix::t(Matrix::t(raw_counts)*(1/Matrix::colSums(raw_counts)))

# Normalise counts
cat("Normalising counts to 1x10-6 molecules per cell  ...\n")
raw_counts_norm <- edgeR::cpm(raw_counts)

cat('Total molecules in cells pre-CPM normalisation:\n\n ')
head(colSums(raw_counts))
cat('Total molecules in cells post-CPM normalisation:\n\n ')
head(colSums(raw_counts_norm))

# Remove genes overlapping MHC region
cat("\nRemoving genes overlapping MHC region  ...\n")
raw_counts_noMHC <- raw_counts_norm[!(rownames(raw_counts_norm) %in% intersect.genes), ]

cat(paste0('\n Genes in ', toupper(REGION), ' seurat obj pre-MHC filtering: ', dim(raw_counts)[1]))
cat(paste0('\n Genes in ', toupper(REGION), ' seurat obj post-MHC filtering: ', dim(raw_counts_noMHC)[1]))
cat('\n Genes are: \n\n', sort(intersect.genes), '\n')

# Create annotations - required “cell_id”, “level1class” and “level2class”
cat("\nCreating annotations ...\n")
annotations <- annotations <- as.data.frame(cbind(colnames(seurat.obj), 
                                                  seurat.obj[[c('cellIDs')]], 
                                                  seurat.obj[[c('cellIDs')]]))
colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
rownames(annotations) <- NULL

##  Generate CTD object  --------------------------------------------------------------
  
  # Uses EWCE package - this creates an object called ctd 
  # Note: option here to scTransform data to correct for cell size - not done this
  # drop.uninformative.genes drops genes that do not show sig variance
  # between level 2 celltypes (based on ANOVA) - not necessary as we are only using
  # 1 level
  # generate.celltype.data calculates the cell specific averages for each gene

cat("Running ANOVA to drop uninformative genes ... \n")
exp_DROPPED <- drop.uninformative.genes(exp = raw_counts_noMHC, 
                                            level2annot = annotations$level2class)
cat('Uninformative genes dropped: ', dim(raw_counts_noMHC)[1] - dim(exp_DROPPED)[1], '\n\n')

annotLevels <- list(level1class = annotations$level2class, 
                    level2class = annotations$level2class)

cat("Creating celltype data ... \n\n")
ctd <- generate.celltype.data(exp = exp_DROPPED, 
                              annotLevels = annotLevels, 
                              groupName = REGION,
                              savePath = CTD_DIR)
print(ctd)
cat(paste0("\n", toupper(REGION), " ctd ... created.\n"))

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
