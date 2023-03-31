# Single-Nuclei RNA Sequencing of 5 Regions of the Human Prenatal Brain Implicates Developing Neuron Populations in Genetic Risk for Schizophrenia (2023)

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). The paper is [here](https://www.biologicalpsychiatryjournal.com/article/S0006-3223(22)01404-4/fulltext). The workflow follows the the snakemake [distribution and reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) recommendations. 

***

A snakemake pipeline to process snRNAseq data. Utilising the following packages:

+ [Snakemake 6.6.1](https://snakemake.readthedocs.io/en/stable/)
+ [CellRanger 5.0.1](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/5.0/what-is-cell-ranger)
+ [Seurat 4.0.3](https://satijalab.org/seurat/articles/get_started.html) - Ran interactively on laptop still to add code for this
+ [DropletUtils 1.10.3](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html)
+ [scDblFinder 1.4.0](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
+ [Scater 1.18.5](https://bioconductor.org/packages/release/bioc/html/scater.html)
+ [MAGMA Celltyping 1.0.0](https://github.com/neurogenomics/MAGMA_Celltyping)
+ [MAGMA 1.08](https://ctg.cncr.nl/software/magma)
+ [LD Score Regression 1.0.1](https://github.com/bulik/ldsc)
+ [clustifyr 1.2.0](https://bioconductor.org/packages/release/bioc/html/clustifyr.html)
+ [scRNAseq_Benchmark](https://github.com/tabdelaal/scRNAseq_Benchmark)
+ etc.

***

**Data**

This study:

+ FASTQ files are available through the [European Genome-Phenome Archive](https://ega-archive.org/studies/EGAS00001006537) under study accession EGAS00001006537.
+ The gene expression matrix and metadata for each assayed brain region is provided through this [figshare repository](https://figshare.com/articles/dataset/11629311).

Cluster Assignment:

+ [Nowakowski](https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_cortex_dev.rda)
+ [Polioudakis](http://solo.bmap.ucla.edu/shiny/webapp/)
+ [Aldinger](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-021-00872-y/MediaObjects/41593_2021_872_MOESM3_ESM.xlsx)

GWAS:

+ [SCZ - PGC3](https://doi.org/10.6084/m9.figshare.14672178)
+ [Height](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files)

SCHEMA:

+ [SCHEMA](https://schema.broadinstitute.org/results)

***

## Copyright and Licence Information

See the [LICENCE file](LICENCE.md).

