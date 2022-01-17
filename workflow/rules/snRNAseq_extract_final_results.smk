rule extract_q10_gene_lists_for_sig_cellTpes:
    input:  "../results/q10_gene_lists_for_LDSC/all_quantiles/wge_complete.file"
    output: "../results/q10_gene_lists_for_LDSC/{SIG_CELL_TYPE}_Q10_genes.tsv"
    params: "../results/q10_gene_lists_for_LDSC/all_quantiles/"
    log:    "../results/logs/extract_final_results/extract_q10_gene_list_{SIG_CELL_TYPE}.log"
    shell:
            """
             
            cp {params}{wildcards.SIG_CELL_TYPE}_Q10_genes.tsv {output}

            """
