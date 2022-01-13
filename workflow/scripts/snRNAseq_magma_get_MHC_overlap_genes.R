# -------------------------------------------------------------------------------------
#
#    Create genes lists of genes in seurat objects that overlap MHC region 
#
# -------------------------------------------------------------------------------------

## Explanation  -----------------------------------------------------------------------

#  These gene lists will be used in the magma_create_ctd.R to filter these genes from
#  The raw count matrices before creation of the ctd objects

## Resources  -------------------------------------------------------------------------

#  MHC locations - https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13 
#  Stack exchange - https://bioinformatics.stackexchange.com/questions/14649/
#  Biostars - Biomart - https://www.biostars.org/p/44426/

## Requirements  -------------------------------------------------------------------------

# Required on Hawk before opening R
# module load libgit2/1.1.0
# module load R/4.0.3

# GWAS:
# SNP, CHR, BP as first three columns.
# GWAS has one of: Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT
# GWAS has all of: SNP, CHR, BP, P, A1, A2

# Net access:
# Uses BiomaRt for MHC gene IDing so need to run SLURM jobs locally

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

## Load packages  ---------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(Seurat)
library(biomaRt)
library(argparser)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read brain region for magma ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
p <- add_argument(p, "seurat_obj", help = "No Seurat obj specified")
p <- add_argument(p, "out_dir", help = "No output directory specified")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
REGION <- args$region
SEURAT_OBJ <- args$seurat_obj
OUT_DIR <- args$out_dir

## Load Seurat obj --------------------------------------------------------------------
seurat.obj <- readRDS(SEURAT_OBJ)

## Get genes in each region -----------------------------------------------------------
seurat.genes <- rownames(seurat.obj)
cat(paste0('\n Genes in ', toupper(REGION), ' seurat obj pre-MHC filtering: ', dim(seurat.obj)[1]))

## Get genes in MHC -------------------------------------------------------------------
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position")
filters <- c("chromosome_name","start","end")
values <- list(chromosome="6",start="28510120",end="33480577")
all.genes <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
all.genes.unique <- unique(all.genes$hgnc_symbol)

#listAttributes(mart = mart)
#searchDatasets(mart = mart, pattern = "hsapiens") # This confirms that were using GRCh38.p13

## Get intersection of seurat obj genes and MHC genes ---------------------------------
intersect.genes <- intersect(rownames(seurat.obj), all.genes.unique)
cat(paste0('\n Genes in ', toupper(REGION), ' seurat obj that intersect MHC region: ', length(intersect.genes)))
cat('\n Genes are: \n\n', sort(intersect.genes), '\n')

## Save MHC genes in seurat object ---------------------------------------------------
cat('\n Saving RDS object ... \n')
saveRDS(intersect.genes, paste0(OUT_DIR, REGION, '_MHC_overlapping_genes.rds'))

cat('Done.')

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

