# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snRNA-seq data - on SCZ GWAS only at this point
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------

rule make_annot:
    input:   gene_set = "../results/q10_gene_lists_for_LDSC/{CELL_TYPE}_Q{QUANTILE}_genes.tsv",
             gene_coord = "../resources/snRNAseq_LDSC_gene_coords.tsv",
             bim_file = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim" 
    output:  "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}.annot.gz"
    conda:   "../envs/ldsc.yml"    
    message: "Creating annotation files for snRNAseq: {wildcards.CELL_TYPE} Quantile {wildcards.QUANTILE}, Chr {wildcards.CHR} "
    log:     "../results/logs/ldsc/make_annot.snRNAseq.{CELL_TYPE}.Q{QUANTILE}.Chr{CHR}.log"
    shell:
        """
        
        python ../resources/ldsc/make_annot.py \
        --gene-set-file {input.gene_set} \
        --gene-coord-file {input.gene_coord} \
        --windowsize 100000 \
        --bimfile {input.bim_file} \
        --annot-file {output} 2> {log} \
        
        """

rule ldsr:
    input:   annot = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsc/reference_files/hapmap3_snps"
    output:  "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsc.yml"
    params:  bfile = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}",
             snps = "../resources/ldsc/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.CELL_TYPE} Quantile {wildcards.QUANTILE} CHR {wildcards.CHR}" 
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.Chr{CHR}_ldsc.log"
    shell:
        "python ../resources/ldsc/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule partitioned_heritability_baseline_v12:
    input:   GWAS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{CHR}.l2.ldscore.gz", CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = '10', CHR = range(1,23))
    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_Q{QUANTILE}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.",
             out_file = "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_Q{QUANTILE}_{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} Q{wildcards.QUANTILE} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_partHerit_summary:
    # This is still optimised for multiple quantiles so creating > 100 single line files
    input:   expand("../results/LDSR_part_herit/snRNAseq_LDSC_{CELL_TYPE}_Q{QUANTILE}_{GWAS}_baseline.v1.2.results", CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = '10', GWAS = config["SUMSTATS"])
    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating summary file for {wildcards.CELL_TYPE} and {wildcards.GWAS} GWAS"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/"
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.dir}snRNAseq_LDSC_Cer-RG-1_Q1_SCZ_baseline.v1.2.results > {output}
             grep L2_1 {params.dir}snRNAseq_LDSC_{wildcards.CELL_TYPE}_Q*_{wildcards.GWAS}_baseline.v1.2.results >> {output}

             """

rule create_top_decile_tables:
    input:   expand("../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_{GWAS}_baseline.v1.2_summary.tsv", CELL_TYPE = config["RNA_CELL_TYPES"], GWAS = config["SUMSTATS"])
    output:  "../results/LDSR_part_herit/snRNAseq_LDSC_{GWAS}_baseline.v1.2_top10pc.tsv"
    message: "Creating LDSC top decile tables for {wildcards.GWAS} GWAS"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/"
    log:     "../results/logs/LDSR/snRNAseq.{GWAS}_partHerit_baseline.v1.2_top10pc_summary.log"
    shell:
             """

             head -1 {params.dir}snRNAseq_LDSC_Cer-OPC_Q10_SCZ_baseline.v1.2.results > {output}

             for file in `ls {params.dir}*Q10_{wildcards.GWAS}*`; do

             CELL_TYPE=$(echo ${{file}} | cut -d'_' -f6)
             tail -1 ${{file}} >> {output}
             sed -i "s/L2_1/${{CELL_TYPE}}/g" {output}
             sed -i '/Total time elapsed/d' {output}

             done

             """        
        
        
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
