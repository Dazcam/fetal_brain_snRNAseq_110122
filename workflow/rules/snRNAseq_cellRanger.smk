# -------------------------------------------------------------------------------------
#
#
#    Script for processing snRNA-seq FASTQ files with Cell Ranger (10X)
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"


# -------------  RULES  ---------------

rule CR_cnt_RNA:
    output: "../results/snRNAseq_CR5_190121/{SAMPLE}.stamp" # Need .stamp and touch {output} in shell command as both SM and CR want to mkdir
    params: fastq_dir=config["FASTQ_DIR"],
            transcriptome=config["REFERENCE_RNA"]
    log:    "../results/logs/cell_ranger/{SAMPLE}.log"
    shell:
            """

            cellranger count --id={wildcards.SAMPLE} \
            --fastqs={params.fastq_dir} \
            --sample={wildcards.SAMPLE} \
            --transcriptome={params.transcriptome} \
            --include-introns \
            --chemistry=SC3Pv3 \
            --localcores=32 \
            --localmem=230 2> {log}
            
            touch {output}

            """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
