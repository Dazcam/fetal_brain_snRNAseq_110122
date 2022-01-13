# -------------------------------------------------------------------------------------
#
#    MAGMA Celltyping plots and data-analysis
#
# -------------------------------------------------------------------------------------

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(cowplot)
library(rmarkdown)
library(argparser)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("Read magma input directory ... \n")
p <- add_argument(p, "magma_dir", help = "No markdown magma dir specified")
p <- add_argument(p, "markdown_file", help = "No markdown file path specified")
p <- add_argument(p, "report_dir", help = "No report output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
args <- parse_args(p)
print(args)


##  Initialise variables  -------------------------------------------------------------
REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
GWAS <- c('ADHD', 'BPD', 'MDD', 'SCZ', 'HEIGHT')
MAGMA_DIR <- args$magma_dir
MARKDOWN_FILE <- args$markdown_file
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file
NAME_BODY <- '_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN'
SUFFIX_LINEAR <- '_linear.gsa.out'
SUFFIX_TOP10 <- '_top10.gsa.out'


# Linear results
cat('\nCreating linear plots ... \n')
for (REGION in REGIONS) {
  
  for (DISORDER in GWAS) {
    
    ##  Load Data  ----------------------------------------------------------------------
    
    #   Need to skip the first 4 columns in datafile 
    
    linearData <- read.table(paste0(MAGMA_DIR, DISORDER, NAME_BODY, '/', DISORDER, NAME_BODY, 
                                     '.level1.', REGION, SUFFIX_LINEAR), header = FALSE)
    names(linearData) <- as.matrix(linearData[1, ])
    linearData <- linearData[-1, ]
    as.numeric(linearData$P)
    -log10(as.numeric(linearData$P))
    print(head(linearData))
    
    ##  Plot  ---------------------------------------------------------------------------
    linearPlot <- ggplot(data = linearData, aes(x = -log10(as.numeric(P)), y = VARIABLE)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_bw() +
      ggtitle(toupper(REGION)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste0(REGION, '_', DISORDER, '_magma_linear_plot'), linearPlot, envir = .GlobalEnv)
    assign(paste0(REGION, '_', DISORDER, '_magma_linear_data'), linearData, envir = .GlobalEnv)
    
  }
  
}


# Top 10% results
cat('\nCreating Top 10% plots ... \n')
for (REGION in REGIONS) {
  
  for (DISORDER in GWAS) {
    
    ##  Load Data  ----------------------------------------------------------------------
    
    #   Need to skip the first 4 columns in datafile 
    
    top10Data <- read.table(paste0(MAGMA_DIR, DISORDER, NAME_BODY, '/', DISORDER, NAME_BODY, 
                                    '.level1.', REGION, SUFFIX_TOP10), header = FALSE)
    names(top10Data) <- as.matrix(top10Data[1, ])
    top10Data <- top10Data[-1, ]
    as.numeric(top10Data$P)
    -log10(as.numeric(top10Data$P))
    print(head(top10Data))
    
    ##  Plot  ---------------------------------------------------------------------------
    top10Plot <- ggplot(data = top10Data, aes(x = -log10(as.numeric(P)), y = VARIABLE)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      theme_bw() +
      ggtitle(toupper(REGION)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    assign(paste0(REGION, '_', DISORDER, '_magma_top10_plot'), top10Plot, envir = .GlobalEnv)
    assign(paste0(REGION, '_', DISORDER, '_magma_top10_data'), top10Data, envir = .GlobalEnv)
    
  }
  
}

# Create group plots
cat('\nCreating group plots ... \n')
for (DISORDER in GWAS) {

  magma_linear_plot <- plot_grid(get(paste0('cer_', DISORDER, '_magma_linear_plot')),
                                 get(paste0('hip_', DISORDER, '_magma_linear_plot')), 
                                 get(paste0('pfc_', DISORDER, '_magma_linear_plot')),
                                 get(paste0('tha_', DISORDER, '_magma_linear_plot')),
                                 get(paste0('wge_', DISORDER, '_magma_linear_plot')))
  
  assign(paste0('all_regions_', DISORDER, '_magma_linear_plot'), magma_linear_plot, envir = .GlobalEnv)
  
  magma_top10_plot <- plot_grid(get(paste0('cer_', DISORDER, '_magma_top10_plot')),
                                 get(paste0('hip_', DISORDER, '_magma_top10_plot')), 
                                 get(paste0('pfc_', DISORDER, '_magma_top10_plot')),
                                 get(paste0('tha_', DISORDER, '_magma_top10_plot')),
                                 get(paste0('wge_', DISORDER, '_magma_top10_plot')))
  
  assign(paste0('all_regions_', DISORDER, '_magma_top10_plot'), magma_top10_plot, envir = .GlobalEnv)
                                        
}

# Render markdown report
cat('\nCreating markdown report ... \n')
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

cat('\nDONE.\n')
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
