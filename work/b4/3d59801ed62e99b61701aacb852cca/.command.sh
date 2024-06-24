#!/bin/bash -ue
/home/tdunham/nextflow_tutorial/btv_consensus_generation_nextflow/scripts/identify_best_segments_from_sam reference FABADRU_600588.sam > FABADRU_600588_best10_refseq.fa

cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
    END_VERSIONS
