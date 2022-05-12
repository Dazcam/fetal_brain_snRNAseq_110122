# Cells of the prenatal brain mediating common variant genetic risk for schizophrenia (2022)

TODO: Add abstract and Authors and link to paper

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). The workflow follows the the snakemake [distribution and reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) recomendations.

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

Data

This study:

+ [FastQ]()
+ [Seurat RDS]()

Cluster Assignment:

+ [Nowakowski](https://github.com/rnabioco/clustifyrdata/raw/master/data/ref_cortex_dev.rda)
+ [Polioudakis](http://solo.bmap.ucla.edu/shiny/webapp/)
+ [Aldinger](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-021-00872-y/MediaObjects/41593_2021_872_MOESM3_ESM.xlsx)

GWAS:

+ [SCZ - PGC3](https://doi.org/10.6084/m9.figshare.14672178)
+ [Height](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files)

SCHEMA:

+ [SCHEMA](https://schema.broadinstitute.org/results)

## Installation and usage
Clone this project, or copy the files into an existing project and edit as needed. This includes:

* LICENSE
* CONTRIBUTING.md
* CODE_OF_CONDUCT.md
* and this README.md

## Contributing
Contributing issues, feature proposals or code is actively welcomed - please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for more details.

## Code of Conduct
We want to create a welcoming environment for everyone who is interested in contributing. Please see the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md) file to learn more about our commitment to an open and welcoming environment.

## Copyright and Licence Information

See the [LICENCE file](LICENCE.md).

