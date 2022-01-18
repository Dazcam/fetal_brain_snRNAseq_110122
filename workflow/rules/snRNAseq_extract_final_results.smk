rule extract_q10_gene_lists_for_sig_cellTpes:
    input:  "../results/gene_lists/q10_gene_lists/all_quantiles/wge_complete.file"
    output: "../results/gene_lists/q10_gene_lists/{SIG_CELL_TYPE}_Q10_genes.tsv"
    params: "../results/gene_lists/q10_gene_lists/all_quantiles/"
    log:    "../results/logs/extract_final_results/extract_q10_gene_list_{SIG_CELL_TYPE}.log"
    shell:
            """
             
            cp {params}{wildcards.SIG_CELL_TYPE}_Q10_genes.tsv {output}

            """
