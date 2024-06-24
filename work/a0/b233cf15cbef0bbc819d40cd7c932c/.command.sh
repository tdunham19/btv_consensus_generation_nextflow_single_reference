#!/bin/bash -ue
minimap2 -ax map-ont reference FABADRU_600588.fastq.gz > FABADRU_600588.sam 

cat <<-END_VERSIONS > versions.yml
    "BTV_CONSENSUS:MINIMAP2_ALIGN":
        minimap2: $(minimap2 --version 2>&1)
        samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
    END_VERSIONS
