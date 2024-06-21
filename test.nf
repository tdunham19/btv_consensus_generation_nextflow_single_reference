params.outdir = 'results'
// publishDir "$params.outdir"

// minimap data against set of BTV ref seqs 
process minimap2 {

workdir = 'fastq' 

container 'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3a70f8bc7e17b723591f6132418640cfdbc88246-0'

input:
path refseq
path input

output:
path "*.sam", emit: sam

script:
"""
minimap2 -ax map-ont $refseq $input > ${input}.sam 
"""
}

workflow {
	fastq_channel = Channel.fromFilePairs('./fastq/FABADRU_*.fastq.gz') 
	refseq_channel=Channel.fromPath("refseq/btv_refseq.fasta").collect()
	minimap2(refseq_channel, fastq_channel)
    	}