# -------------------------------------------------------------------------------------
#
#
#    MAGMA conditional analyses
#
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# ----------  SET VARIABLES  ----------
MAGMA_DIR = "../results/magma_celltyping/SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/"

# -------------  RULES  ---------------

rule magma_conditional:
    input:   gene_list = "../results/q10_gene_lists_for_LDSC/ALL_SIG_AND_SKENE_entrez_gene_list.tsv",
             scz_magma = MAGMA_DIR + "SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw"
    output:  "../results/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}.gsa.out"
    params:  "../results/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}"
    message: "Running magma on all significant cell types conditioning on {wildcards.CONDITION}"
    log:     "../results/logs/magma_conditional/snRNAseq.magma.conditional.{CONDITION}.log"
    shell:
             """

             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} --model condition={wildcards.CONDITION} --out {params}

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
