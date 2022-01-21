# -------------------------------------------------------------------------------------
#
#    SVM analysis
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARMAS  ----------

configfile: "../config/config.yaml"

# -------------  RULES  ---------------

rule snRNAseq_SVM_prep:  
    # Produces raw RNAseq count and cluster label tables  for SVM analysis
    input:  "../resources/R_objects/seurat.wge.final.rds"
    output: "../resources/sheets/snRNAseq_cluster_labels_rand_wge_for_SVM.csv"
    log:    "../results/logs/SVM/snRNAseq_SVM_file_prep.log"
    shell:
            """

            export R_LIBS_USER=/scratch/c.c1477909/R/library
            module load libgit2/1.1.0
            module load pandoc/2.7.3
            /apps/languages/R/4.0.3/el7/AVX512/gnu-8.1/bin/Rscript --vanilla \
            scripts/snRNAseq_SVM_prep_files.R 2> {log}

            """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------


