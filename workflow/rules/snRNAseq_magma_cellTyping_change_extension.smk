# -------------------------------------------------------------------------------------
#
#    Script for MAGMA cell typing analysis: testing different gene extensions
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"
localrules:  map_genes_to_snps_35kbUP_10kbDOWN


# ----------  SET VARIABLES  ----------

MAGMA_DIR = "../resources/magma_celltyping/"
GWAS_DIR =  MAGMA_DIR + "GWAS_for_magma/"
MAGMA_OUTDIR = GWAS_DIR + "MAGMA_Files/"

rule map_genes_to_snps_change_extensions:
    # Requires net access to run
    input:  GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv"
    output: MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN.genes.raw"
    log:   "../results/logs/magma/map_snps2genes_hg19_{MAGMA_EXT_UP}UP_{MAGMA_EXT_DOWN}DOWN_{GWAS}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/whoami/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_magma_map_snps2genes_change_extensions.R {input} {wildcards.MAGMA_EXT_UP} {wildcards.MAGMA_EXT_DOWN} 2> {log}

            """

rule magma_analysis_change_extensions:
    input:  ctd_obj = "../resources/ctd_objects/CellTypeData_{REGION}.rda",
            gene_file = MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN.genes.raw", 
            gwas_file = GWAS_DIR + "{GWAS}_hg19_magma_ready.sumstats.tsv" 
    output: MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN.level2.{REGION}_top10.gsa.out"
    log:    "../results/logs/magma/magma_analysis_hg19_{MAGMA_EXT_UP}UP_{MAGMA_EXT_DOWN}DOWN_{REGION}_{GWAS}.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_magma_analysis_change_extensions.R {wildcards.REGION} {input.ctd_obj} {input.gwas_file} {wildcards.MAGMA_EXT_UP} {wildcards.MAGMA_EXT_DOWN} 2> {log}

            """

rule cp_magma_results_change_extensions:
    # Need this as magma celltyping put the results in a wierd place
    input:  files = expand(MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN/{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN.level2.{REGION}_top10.gsa.out", REGION = config["RNA_REGIONS"], GWAS = config["SUMSTATS"], MAGMA_EXT_UP = config["MAGMA_EXT_UP"], MAGMA_EXT_DOWN = config["MAGMA_EXT_DOWN"])
    output: directory("../results/magma_celltyping/{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN/")
    params: dir = MAGMA_OUTDIR + "{GWAS}_hg19_magma_ready.sumstats.tsv.{MAGMA_EXT_UP}UP.{MAGMA_EXT_DOWN}DOWN/"
    log:    "../results/logs/magma/cp_magma_results_{MAGMA_EXT_UP}UP_{MAGMA_EXT_DOWN}DOWN_{GWAS}.log" 
    shell: 
            """

            cp -r {params.dir} {output}

            """
