Repo for single nuceli fetal brain study using cell label version 110121

***

## Downstream analyses

### MAGMA process - Note this was run using Magma Celltyping (MCT) package as we initally ran analysis using 
the linear regression model described in Skene et al. 2019. We changed this subsequently to top 10% (which is an option in MCT) but decided 
to stick with using MCT for this part of the process rather than running magma manually.

SCRIPT: snRNAseq_magma_cellTyping.smk
INPUT:  R_objects/seurat.cer.final.rds
OUTPUT: results/magma_celltyping, results/figures/

1: get_genes_in_MHC - Identify genes in each seurat object that overlap MHC region
2: create_ctd - create ctd object (using magma celltyping package) splits seurat GeX counts into quantiles
3: map_SNPs_to_genes - (using magma cell typing which uses magma) to map SNPS to genes for magma analysis
4: magma_analysis - runs the magma analysis on the cell types of each region 
5. magma_generate_plots - creates plots for magma anlysis
6: magma_create_LDSC_geneLists - extracts q10 gene lists from CTD object to run LDSC analysis
7: Plotting  
