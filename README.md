# Fungal Pipeline (Nextflow/nf-core Demo)

## Introduction

This pipeline is strictly an exercise to refresh my nextflow knowledge and build upon it by transcribing one of my existing pipelines into nextflow.

This pipeline is designed to run on a specific SLURM cluster.

## Problem

Instead of sequencing highly conserved regions of (such as 16S rRNA), fungal amplicon sequencing targets regions of high variability which includes variability in length. In order to maximize the preformance of downstream bioinformatic tools, these sequences much be trimmed to the desired region.

Here I use `itsxpress` to target the ITS2 region for 2 different paired-end sequencing runs which are part of the same experiment.

## Pipeline Outline

The pipeline is broken into 3 steps:

1) Trim the fungal sequences to only contain the region of interest using `itsxpress`.
2) Generate a `fastqc` report for the forward and reverse reads for each sequencing run/sample id combination.
3) Using `multiqc` generate an sequencing quality report for the forward and reverse reads on each sequencing run.

## TODO

- [X] add in the `itsxpress` module to the `subworkflows/local/utils_nfcore_fungal_pipeline_pipeline/main.nf` and fill in TODO
- [X] Call `nf-core` `fastqc` and `multiqc` modules and fill in relevant info
- [ ] Update config file to run on SLURM HPC with ceres specific options.
- [ ] Write a test profile options for `itsxpress` process.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
