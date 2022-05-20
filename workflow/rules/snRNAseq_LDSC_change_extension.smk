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

localrules: extend_gene_window

rule extend_gene_window:
    input: coords = "../resources/sheets/snRNAseq_LDSC_gene_coords.tsv"
    output: tsv = "../results/gene_lists/gene_coords/snRNAseq_LDSC_gene_coords_UP{EXT_UP}_DOWN{EXT_DOWN}.tsv"
    message: "Extending gene coord file {wildcards.EXT_UP}bp upstream and {wildcards.EXT_DOWN}bp downstream"
    log:     "../results/logs/ldsc/extend_gene_window_UP{EXT_UP}_DOWN{EXT_DOWN}.log" 
    run:

        genes = pd.read_csv(input.coords, delim_whitespace = True)
        genes['START'] = np.maximum(1, genes['START'] - int(wildcards.EXT_UP))
        genes['END'] = genes['END'] + int(wildcards.EXT_DOWN)
        genes.to_csv(output.tsv, sep="\t", index=False)


rule make_annot:
    input:   gene_set = "../results/gene_lists/q10_gene_lists/all_q10s/{CELL_TYPE}_Q{QUANTILE}_genes.tsv",
             gene_coord = "../results/gene_lists/gene_coords/snRNAseq_LDSC_gene_coords_UP{EXT_UP}_DOWN{EXT_DOWN}.tsv",
             bim_file = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim" 
    output:  "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.{CHR}.annot.gz"
    conda:   "../envs/ldsc.yml"    
    message: "Creating annotation files for snRNAseq: {wildcards.CELL_TYPE} Quantile {wildcards.QUANTILE}, Chr {wildcards.CHR} with extensions UP{wildcards.EXT_UP}.DOWN{wildcards.EXT_DOWN}"
    log:     "../results/logs/ldsc/make_annot.snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.Chr{CHR}.log"
    shell:
        """
        
        python ../resources/ldsc/make_annot.py \
        --gene-set-file {input.gene_set} \
        --gene-coord-file {input.gene_coord} \
        --windowsize 0 \
        --bimfile {input.bim_file} \
        --annot-file {output} 2> {log} \
        
        """

rule ldsr:
    input:   annot = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsc/reference_files/hapmap3_snps"
    output:  "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsc.yml"
    params:  bfile = "../resources/ldsc/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.{CHR}",
             snps = "../resources/ldsc/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.CELL_TYPE} Quantile {wildcards.QUANTILE} CHR {wildcards.CHR} with extensions UP{wildcards.EXT_UP}.DOWN{wildcards.EXT_DOWN}" 
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.Chr{CHR}_ldsc.log"
    shell:
        "python ../resources/ldsc/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule partitioned_heritability_baseline_v12:
    input:   GWAS = "../results/GWAS_for_ldsc/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.{CHR}.l2.ldscore.gz", CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = '10', CHR = range(1,23), EXT_UP = config["EXT_UP"], EXT_DOWN = config["EXT_DOWN"])
    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_Q{QUANTILE}_UP{EXT_UP}_DOWN{EXT_DOWN}_{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsc.yml"
    params:  weights = "../resources/ldsc/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsc/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsc/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.",
             out_file = "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_Q{QUANTILE}_UP{EXT_UP}_DOWN{EXT_DOWN}_{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} Q{wildcards.QUANTILE} with extensions UP{wildcards.EXT_UP}.DOWN{wildcards.EXT_DOWN} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.Q{QUANTILE}.UP{EXT_UP}.DOWN{EXT_DOWN}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsc/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule create_partHerit_summary:
    # This is still optimised for multiple quantiles so creating > 100 single line files
    input:   expand("../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_Q{QUANTILE}_UP{EXT_UP}_DOWN{EXT_DOWN}_{GWAS}_baseline.v1.2.results", CELL_TYPE = config["RNA_CELL_TYPES"], QUANTILE = '10', EXT_UP = config["EXT_UP"], EXT_DOWN = config["EXT_DOWN"], GWAS = config["SUMSTATS"])
    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_UP{EXT_UP}_DOWN{EXT_DOWN}_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating summary file for {wildcards.CELL_TYPE} with extensions UP{wildcards.EXT_UP}.DOWN{wildcards.EXT_DOWN} and {wildcards.GWAS} GWAS"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/"
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.UP{EXT_UP}.DOWN{EXT_DOWN}.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """

             head -1 {params.dir}snRNAseq_LDSC_Cer-RG-1_Q1_SCZ_baseline.v1.2.results > {output}
             grep L2_1 {params.dir}snRNAseq_LDSC_{wildcards.CELL_TYPE}_Q*_UP{wildcards.EXT_UP}_DOWN{wildcards.EXT_DOWN}_{wildcards.GWAS}_baseline.v1.2.results >> {output}

             """

rule create_top_decile_tables:
    input:   expand("../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_UP{EXT_UP}_DOWN{EXT_DOWN}_{GWAS}_baseline.v1.2_summary.tsv", CELL_TYPE = config["RNA_CELL_TYPES"], EXT_UP = config["EXT_UP"], EXT_DOWN = config["EXT_DOWN"], GWAS = config["SUMSTATS"])
    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_UP{EXT_UP}_DOWN{EXT_DOWN}_{GWAS}_baseline.v1.2_top10pc.tsv"
    message: "Creating LDSC top decile tables for {wildcards.GWAS} GWAS with extensions UP{wildcards.EXT_UP}.DOWN{wildcards.EXT_DOWN}"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/"
    log:     "../results/logs/LDSR/snRNAseq.{GWAS}_UP{EXT_UP}_DOWN{EXT_DOWN}_partHerit_baseline.v1.2_top10pc_summary.log"
    shell:
             """

             head -1 {params.dir}snRNAseq_LDSC_Cer-OPC_Q10_SCZ_baseline.v1.2.results > {output}

             for file in `ls {params.dir}*Q10_UP{wildcards.EXT_UP}_DOWN{wildcards.EXT_DOWN}_{wildcards.GWAS}*`; do

             CELL_TYPE=$(echo ${{file}} | cut -d'_' -f6)
             tail -1 ${{file}} >> {output}
             sed -i "s/L2_1/${{CELL_TYPE}}/g" {output}
             sed -i '/Total time elapsed/d' {output}

             done

             """        
        
        
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
