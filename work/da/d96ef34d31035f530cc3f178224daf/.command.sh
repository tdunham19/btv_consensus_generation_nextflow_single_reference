#!/bin/bash -ue
samtools index FABADRU_891979.sorted.bam

cat <<-END_VERSIONS > versions.yml
"BTV_CONSENSUS:SAMTOOLS_INDEX":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
