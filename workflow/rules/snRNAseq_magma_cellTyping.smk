# -------------------------------------------------------------------------------------
#
#    Script for MAGMA cell typing analysis
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"
localrules: get_genes_in_MHC, create_ctd, map_genes_to_snps


# ----------  SET VARIABLES  ----------

MAGMA_DIR = "../resources/magma_celltyping/"
LOG_DIR = "results/logs/"
GWAS_DIR =  MAGMA_DIR + "GWAS_for_magma/"
MAGMA_OUTDIR = GWAS_DIR + "MAGMA_Files/"
RESULTS_DIR = "results/"
GENE_LIST_OUTDIR = RESULTS_DIR + "results/q10_gene_lists_for_LDSC/"
MARKDOWN_DIR = "reports/"
RESULTS_DIR = "results/"

# -------------  RULES  ---------------

rule get_genes_in_MHC:
    input:  "../resources/R_objects/seurat.{REGION}.final.rds"
    output: "../resources/R_objects/{REGION}_MHC_overlapping_genes.rds"
    params: outdir = "../resources/R_objects/"
    log:    "../results/logs/magma/get_genes_in_MHC_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/whoami/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_magma_get_MHC_overlap_genes.R {wildcards.REGION} {input} {params.outdir} 2> {log}          

            """

rule create_ctd:
    input:  seurat_obj = "../resources/R_objects/seurat.{REGION}.final.rds",
            MHC_gene_list = "../resources/R_objects/{REGION}_MHC_overlapping_genes.rds"
    output: "../resources/ctd_objects/CellTypeData_{REGION}.rda"
    log:    "../results/logs/magma/create_ctd_{REGION}.log" 
    shell: 
            """
            
            export R_LIBS_USER=/scratch/whoami/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_magma_create_ctd.R {wildcards.REGION} {input.seurat_obj} {input.MHC_gene_list} 2> {log}

            """

rule map_genes_to_snps:
    # Requires net access to run
    input:  GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv"
    output: MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw"
    log:   "../results/logs/magma/map_snps2genes_hg19_{GWAS}.log"
    shell:
            """
       
            export R_LIBS_USER=/scratch/whoami/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_magma_map_snps2genes.R {input} 2> {log}
              
            """

rule magma_analysis:
    input:  ctd_obj = "../resources/ctd_objects/CellTypeData_{REGION}.rda",
            gene_file = MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw", 
            gwas_file = GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv" 
    output: MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.level2.{REGION}_top10.gsa.out"
    log:    "../results/logs/magma/magma_analysis_hg19_{REGION}_{GWAS}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_magma_analysis.R {wildcards.REGION} {input.ctd_obj} {input.gwas_file} 2> {log}

            """

rule cp_magma_results:
    # Need this as magma celltyping put the results in a wierd place
    input:  files = expand(MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.level2.{REGION}_top10.gsa.out", REGION = config["RNA_REGIONS"], GWAS = config["SUMSTATS"])
    output: directory("../results/magma_celltyping/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/")
    params: dir = MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/"
    log:    "../results/logs/magma/cp_magma_results_{GWAS}.log" 
    shell: 
            """

            cp -r {params.dir} {output}

            """

# Maybe don't need this - currently doing manually on Desktop
#rule magma_generate_plots:
#    input:  magma = expand(MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.level2.{REGION}_top10.gsa.out", REGION = config["RNA_REGIONS"], GWAS = config["SUMSTATS"]), 
#            markdown = MARKDOWN_DIR + "snRNAseq_magma_generate_plots.Rmd"     
#    output: RESULTS_DIR + "snRNAseq_magma_generate_plots.html"
#    params: magma_dir = MAGMA_OUTDIR,
#            report_dir = RESULTS_DIR,
#            report_file = "snRNAseq_magma_generate_plots.html",
#    log:    LOG_DIR + "magma/magma_generate_plots.log"
#    shell:
#            """

#            export R_LIBS_USER=/scratch/c.c1477909/R/library
#            module load libgit2/1.1.0
#            module load pandoc/2.7.3
#            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
#            scripts/snRNAseq_magma_generate_plots.R {params.magma_dir} {input.markdown} {params.report_dir} {params.report_file}  2> {log}

#            """

rule magma_create_q10_geneLists:
    input:  ctd_obj = "../resources/ctd_objects/CellTypeData_{REGION}.rda", 
    output: "../results/gene_lists/q10_gene_lists/all_quantiles/{REGION}_complete.file"
    log:    "../results/logs/magma/magma_ldsc_gene_lists_{REGION}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_magma_create_LDSC_geneLists.R {wildcards.REGION} {input.ctd_obj} 2> {log}

            """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------


