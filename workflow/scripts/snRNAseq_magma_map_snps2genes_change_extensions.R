# -------------------------------------------------------------------------------------
#
#    MAGMA Celltyping - map SNPs to genes
#
# -------------------------------------------------------------------------------------

##  Resources  -------------------------------------------------------------------------

    #  Paper - https://www.nature.com/articles/s41588-018-0129-5
    #  Github - https://github.com/NathanSkene/MAGMA_Celltyping

    # ECWE - for creating the ctd object
      # Github - https://github.com/neurogenomics/EWCE/

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
library(MAGMA.Celltyping) # Note the "." instead of "_" - using srun can't see magma executable
library(argparser)
library(reshape2)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read GWAS sumstats path for magma to map SNPs to genes ... \n")
p <- add_argument(p, "gwas", help = "No GWAS sumstats file path provided")
p <- add_argument(p, "upstream_ext", help = "No upstream gene extension provided")
p <- add_argument(p, "downstream_ext", help = "No downstream gene extension provided")
args <- parse_args(p)
print(args)

##  Define variables  -----------------------------------------------------------------
MAGMA_DIR <- '/scratch/c.c1477909/fetal_brain_snRNAseq_110122/resources/magma_celltyping/'

GENOME_DIR <- paste0(MAGMA_DIR, 'g1000_eur/g1000_eur')
GWAS_PATH <- args$gwas 
UPSTREAM_EXT <- args$upstream_ext
DOWNSTREAM_EXT <- args$downstream_ext

##  Map SNPs to genes  ----------------------------------------------------------------

    # Produces a .genes.out file 
    # This requires net access so needs to be run locally in home
    # ~10 mins on Hawk home

genesOutPath <- map.snps.to.genes(path_formatted = GWAS_PATH, 
                                  genome_ref_path = GENOME_DIR,
                                  upstream_kb = UPSTREAM_EXT,
                                  downstream_kb = DOWNSTREAM_EXT)
cat('\n\n\nGenes out path: ', genesOutPath, '\n\n')


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
