# BTV Consensus Sequence Generation with Single Reference - Long-Read Sequence Data, Nextflow Pipeline
## Created by Tillie Dunham on 07.02.24

This pipeline is written in nextflow and generates a consensus sequence for Bluetongue virus (BTV) from long-read sequencing reads generated with Oxford Nanopore technology. 

The input to this pipeline is: 
  Concatenated fastq file of all long reads for the sample.
  BTV reference sequence for a single BTV serotype.
