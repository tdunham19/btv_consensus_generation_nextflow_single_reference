#!/bin/bash -ue
samtools view -F4 -S -b FABADRU_891979.sam > FABADRU_891979.unsorted.bam

 cat <<-END_VERSIONS > versions.yml
 "BTV_CONSENSUS:SAMTOOLS_VIEW":
     samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
 END_VERSIONS
