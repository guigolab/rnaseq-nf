# RNA-seq pipeline

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.07.0-brightgreen.svg)](http://nextflow.io)
[![Build Status](https://travis-ci.org/guigolab/rnaseq-nf.svg?branch=master)](https://travis-ci.org/guigolab/rnaseq-nf)

An example pipeline for quantification of genomic features from short read data implemented with Nextflow.

## Requirements 

* Unix-like operating system (Linux, macOS, etc)
* Java 8
* Docker/Singularity

## Quickstart 

1. Install Nextflow (version 0.24.x or higher):
      
        curl -s https://get.nextflow.io | bash

3. Launch the pipeline execution using Docker: 

        ./nextflow run guigolab/rnaseq-nf -with-docker
        
4. When the execution completes, pipeline results can be found inside the `results` folder
	
Note: the very first time you execute it, it will take a few minutes to download the pipeline from this GitHub repository and the the associated Docker images needed to execute the pipeline.

## Pipeline results

The processing output files are stored in a folder named `results`
inside the pipeline working folder. The structure of the output directory
is as follows:

``` bash
results
├── mapping
│   ├── mouse_cns_E14_rep1_Aligned.sortedByCoord.out.bam
│   ├── mouse_cns_E14_rep1_Aligned.toTranscriptome.out.bam
│   ├── mouse_cns_E14_rep1_Log.final.out
│   ├── mouse_cns_E14_rep2_Aligned.sortedByCoord.out.bam
│   ├── mouse_cns_E14_rep2_Aligned.toTranscriptome.out.bam
│   ├── mouse_cns_E14_rep2_Log.final.out
│   ├── mouse_cns_E18_rep1_Aligned.sortedByCoord.out.bam
│   ├── mouse_cns_E18_rep1_Aligned.toTranscriptome.out.bam
│   ├── mouse_cns_E18_rep1_Log.final.out
│   ├── mouse_cns_E18_rep2_Aligned.sortedByCoord.out.bam
│   ├── mouse_cns_E18_rep2_Aligned.toTranscriptome.out.bam
│   └── mouse_cns_E18_rep2_Log.final.out
├── matrix
│   ├── mouse_cns.gene.matrix.TPM.tsv
│   └── mouse_cns.gene.matrix.expected_count.tsv
└── quantification
    ├── mouse_cns_E14_rep1.genes.results
    ├── mouse_cns_E14_rep1.isoforms.results
    ├── mouse_cns_E14_rep2.genes.results
    ├── mouse_cns_E14_rep2.isoforms.results
    ├── mouse_cns_E18_rep1.genes.results
    ├── mouse_cns_E18_rep1.isoforms.results
    ├── mouse_cns_E18_rep2.genes.results
    └── mouse_cns_E18_rep2.isoforms.results

3 directories, 22 files
```

## Cluster support

RNASeq-NF execution relies on [Nextflow](http://www.nextflow.io) framework which provides an abstraction between the pipeline functional logic and the underlying processing system.

This allows the execution of the pipeline in a single computer or in a HPC cluster without modifying it.

Currently the following resource manager platforms are supported:

  + Univa Grid Engine (UGE)
  + Platform LSF
  + SLURM
  + PBS/Torque


By default the pipeline is parallelized by spawning multiple threads in the machine where the script is launched.

To submit the execution to a UGE cluster create a file named `nextflow.config` in the directory where the pipeline is going to be executed with the following content:

    process {
      executor='uge'
      queue='<queue name>'
    }

To lean more about the avaible settings and the configuration file read the Nextflow [documentation](http://www.nextflow.io/docs/latest/config.html).


## Components 

The pipeline uses the following software: 

* [STAR](https://github.com/alexdobin/STAR/) 2.7.2d
* [RSEM](https://github.com/deweylab/RSEM/) 1.3.1

