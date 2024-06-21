#!/bin/bash -ue
samtools sort FABADRU_600588.unsorted.bam -o FABADRU_600588.sorted.bam -O bam

   cat <<-END_VERSIONS > versions.yml
   "BTV_CONSENSUS:SAMTOOLS_SORT":
       samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
   END_VERSIONS
