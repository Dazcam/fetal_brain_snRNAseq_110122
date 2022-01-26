# -------------------------------------------------------------------------------------
#
#
#    MAGMA analyses on sig. cell type GO terms 
#
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# ----------  SET VARIABLES  ----------
MAGMA_DIR = "../results/magma_celltyping/SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN/"

# -------------  RULES  ---------------

rule magma_GO:
    input:   gene_list = "../results/GO/GO_term_genes_for_magma.txt",
             scz_magma = MAGMA_DIR + "SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw"
    output:  "../results/magma_GO/magma_GO.gsa.out"
    params:  "../results/magma_GO/magma_GO"
    message: "Running magma final GO terms for all significant cell types"
    log:     "../results/logs/magma_GO/magma_GO.log"
    shell:
             """

             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} --out {params}

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
