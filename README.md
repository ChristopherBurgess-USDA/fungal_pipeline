# Fungal Pipeline (Nextflow/nf-core Demo)

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/fungal_pipeline)
[![GitHub Actions CI Status](https://github.com/nf-core/fungal_pipeline/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/fungal_pipeline/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/fungal_pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/fungal_pipeline/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/fungal_pipeline)

## Introduction

This pipeline is strictly an exercise to refresh my nextflow knowledge and build upon it by transcribing one of my existing pipelines into nextflow.

This pipeline is designed to run on a specific SLURM cluster.

## Problem

Instead of sequencing highly conserved regions of (such as 16S rRNA), fungal amplicon sequencing targets regions of high variability which includes variability in length. In order to maximize the preformance of downstream bioinformatic tools, these sequences much be trimmed to the desired region.

Here I use `itsxpress` to target the ITS2 region for 2 different paired-end sequencing runs which are part of the same experiment.

## Input data

The input data is the results of 2 different paired end sequencing runs. Due to sequencing batch effects, the sequencing runs keep to be kept separate even though sample ids overlap across runs. There is a mapping file (`csv`) containing which sequencing run, sample id, location of forward reads, and location of reverse reads.

## Pipeline Outline

The pipeline is broken into 3 steps:

1) Trim the fungal sequences to only contain the region of interest using `itsxpress`.
2) Generate a `fastqc` report for the forward and reverse reads for each sequencing run/sample id combination.
3) Using `multiqc` generate an sequencing quality report for the forward and reverse reads on each sequencing run.


## TODO

- [X] add in the `itsxpress` module to the `subworkflows/local/utils_nfcore_fungal_pipeline_pipeline/main.nf` and fill in TODO
- [X] Call `nf-core` `fastqc` and `multiqc` modules and fill in relevant info
- [ ] Write a `config` file for the `nf-core ampliseq` pipeline for its2 amplicon sequencing
- [ ] write `R` script to parse output and format to desired ouput and set it up in the `bin/` folder

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
