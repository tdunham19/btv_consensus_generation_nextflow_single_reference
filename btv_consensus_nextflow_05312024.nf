params.outdir = 'results'
// publishDir "$params.outdir"

// minimap data against set of BTV ref seqs 
process minimap2 {

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
	fastq_channel=Channel.fromPath("fastq/*.fastq.gz")
	refseq_channel=Channel.fromPath("refseq/btv_refseq.fasta").collect()
	minimap2(refseq_channel, fastq_channel)
    	samtoolsview(minimap2.out.sam)
    	samtoolssort(samtoolsview.out.bam)
    	samtoolsindex(samtoolssort.out.bam)
    	bestsegsfromsam(refseq_channel, minimap2.out.sam)
    	minimap2_2(bestsegsfromsam.out.fa, fastq_channel)
    	samtoolsview_2(minimap2_2.out.sam)
    	samtoolssort_2(samtoolsview_2.out.bam)
    	samtoolsindex_2(samtoolssort_2.out.bam)
    	samtoolsfaidx(bestsegsfromsam.out.fa)
    	bcftoolsmpileup(bestsegsfromsam.out.fa, samtoolssort_2.out.bam)
    	createmaskfile(bcftoolsmpileup.out.vcf)
    	bcftoolsview(bcftoolsmpileup.out.vcf)
    	bcftoolsindex(bcftoolsview.out.gz)
    	bcftoolscall(bcftoolsview.out.gz)
    	bcftoolsindex_2(bcftoolscall.out.gz)
    	bcftoolsconsensus(createmaskfile.out.mask, bestsegsfromsam.out.fa, bcftoolscall.out.gz)
}


// process minimap 2 output
process samtoolsview {

container 'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0'

input: 
path alignment

output: 
path "*_alignment.unsorted.bam", emit: bam

script: 
"""
samtools view -F4 -S -b $alignment > ${alignment}_alignment.unsorted.bam
"""
}

/*
workflow {
	alignment_channel=Channel.fromPath(alignment/*_alignment.sam)
	samtoolsview(alignment_channel)
}
*/


// process minimap 2 output
process samtoolssort {

container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'

input: 
path unsorted

output: 
path "*_alignment.sorted.bam", emit: bam

script: 
"""
samtools sort $unsorted -o ${unsorted}_alignment.sorted.bam -O bam
"""
}

/*
workflow {
	unsorted_channel=Channel.fromPath(unsorted/*_alignment.unsorted.sam)
	samtoolssort(unsorted_channel)
}
*/


// process minimap 2 output
process samtoolsindex {

container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'

input: 
path sorted

output: 
path "*_alignment.sorted.bam.bai", emit: bai

script: 
"""
samtools index ${sorted}
"""
}

/*
workflow {
	sorted_channel=Channel.fromPath(sorted/*_alignment.sorted.bam)
	samtoolsindex(sorted_channel)
}
*/


// extract new fasta file containing best aligned-to seqs for this dataset
process bestsegsfromsam {

input: 
path refseq
path alignment 

output: 
path "${alignment}_best10_refseq.fa", emit: fa

script: 
"""
identify_best_segments_from_sam ${refseq} ${alignment} > ${alignment}_best10_refseq.fa

"""

}

/* 
workflow {
	alignment_channel=Channel.fromPath(alignment/*_alignment.sam)
	refseq_channel=Channel.fromPath("btv_refseq.fasta").collect()
	data = channel.fromPath('./modules/local/identify_best_segments_from_sam')
	bestsegsfromsam(data)
	
}
*/


// re-minimap data against best 10 BTV ref seqs
process minimap2_2 {

container 'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3a70f8bc7e17b723591f6132418640cfdbc88246-0'

input:
path best10refseq
path fastq

output:
path "*.sam", emit: sam

script:
"""
minimap2 -ax map-ont $best10refseq $fastq > ${fastq}_alignment_best10.sam 

"""
}

/* 
workflow {
	fastq_channel=Channel.fromPath("fastq/*.fastq.gz")
	best10refseq_channel=Channel.fromPath("best10refseq/*_best10_refseq.fa").collect()
	minimap2(best10refseq_channel, fastq_channel)
}
*/


// process minimap2 output again
process samtoolsview_2 {

container 'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0'

input: 
path best10alignment

output: 
path "*_alignment_best10.unsorted.bam", emit: bam

script: 
"""
samtools view -F4 -S -b $best10alignment > ${best10alignment}_alignment_best10.unsorted.bam
"""
}

/*
workflow {
	best10alignment_channel=Channel.fromPath(best10alignment/*_alignment_best10.sam)
	samtoolsview(alignment_channel)
}
*/


// process minimap2 output again
process samtoolssort_2 {

container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'

input: 
path best10unsorted

output: 
path "*_alignment_best10.sorted.bam", emit: bam

script: 
"""
samtools sort $best10unsorted -o ${best10unsorted}_alignment_best10.sorted.bam -O bam
"""
}

/*
workflow {
	best10unsorted_channel=Channel.fromPath(best10unsorted/*_alignment_best10.unsorted.bam)
	samtoolssort(best10unsorted_channel)
}
*/


// process minimap2 output again
process samtoolsindex_2 {

container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'

input: 
path best10sorted

output: 
path "*_alignment_best10.sorted.bam.bai", emit: bai

script: 
"""
samtools index ${best10sorted}
"""
}

/*
workflow {
	best10sorted_channel=Channel.fromPath(best10sorted/*_alignment_best10.sorted.bam)
	samtoolsindex(best10sorted_channel)
}
*/



// commands below create "new" consensus sequences


// have to make a .fai file to make mpileup happy -> output is ${barcode_var}_alignment_best10_refseq.fa.fai
process samtoolsfaidx {

container 'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0'

input: 
path best10refseq

output: 
path "*_best10_refseq.fa.fai", emit: fai

script: 
"""
samtools faidx ${best10refseq}
"""
}

/*
workflow {
	best10refseq_channel=Channel.fromPath("best10refseq/*_best10_refseq.fa").collect()
	samtoolsfaidx(best10refseq_channel)
}
*/


// mpileup calls variants -> output is a vcf file
process bcftoolsmpileup {

container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'

input: 
path best10refseq
path best10sorted

output: 
path "*_best10.paired_sorted.bam.vcf", emit: vcf

script: 
"""
bcftools mpileup -f ${best10refseq} -Ov -o ${best10sorted}_best10.paired_sorted.bam.vcf ${best10sorted}
"""
}

/*
workflow {
	best10refseq_channel=Channel.fromPath("best10refseq/*_best10_refseq.fa").collect()
	best10sorted_channel=Channel.fromPath(best10sorted/*_alignment_best10.sorted.bam)
	bcftoolsmpileup(best10refseq_channel, best10sorted_channel)
}
*/


// this script creates a mask file
// this is necessary because otherwise bcftools consensus doesn't hanlde positions with no coverage well
process createmaskfile {

input: 
path vcf

output: 
path "${vcf}.mask", emit: mask

script: 
"""
create_mask_file ${vcf} > ${vcf}.mask
"""

}

/* 
workflow {
	vcf_channel=Channel.fromPath(vcf/*.vcf)
	data = channel.fromPath('./modules/local/create_mask_file')
	createmaskfile(data)
	
}
*/


// convert vcf -> compressed vcf to make bcftools happy
process bcftoolsview {

container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'

input: 
path vcf

output: 
path "${vcf}.gz", emit: gz

script: 
"""
bcftools view ${vcf} -O z > ${vcf}.gz 
"""
}

/*
workflow {
	vcf_channel=Channel.fromPath(vcf/*.vcf)
	bcftoolsview(vcf_channel)
}
*/


// need to make an indexed bcf file to make other bcftools commands happy
// output is *.paired_sorted.bam.vcf.gz.csi
process bcftoolsindex {

container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  
input: 
path compressedvcf

output: 
path "${compressedvcf}.csi", emit: csi

script: 
"""
bcftools index ${compressedvcf}
"""
}

/*
workflow {
	compressedvcf_channel=Channel.fromPath(compressedvcf/*.vcf.gz)
	bcftoolsindex(compressedvcf_channel)
}
*/


// bcftools call creates a new vcf file that has consensus bases 
process bcftoolscall {

container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  
input: 
path compressedvcf

output: 
path "${compressedvcf}.cons.vcf.gz", emit: gz

script: 
"""
bcftools call -c ${compressedvcf} -O b > ${compressedvcf}.cons.vcf.gz
"""
}

/*
workflow {
	compressedvcf_channel=Channel.fromPath(compressedvcf/*.vcf.gz)
	bcftoolscall(compressedvcf_channel)
}
*/


// need to make an indexed bcf file to make other bcftools commands happy
// output is *.paired_sorted.bam.cons.vcf.gz.csi
process bcftoolsindex_2 {

container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'
  
input: 
path cons

output: 
path "${cons}.csi", emit: csi

script: 
"""
bcftools index ${cons}
"""
}

/*
workflow {
	cons_channel=Channel.fromPath(cons/*.cons.vcf.gz)
	bcftoolsindex_2(cons_channel)
}
*/


/*
bcftools consensus will output a fasta file containing new draft consensus sequence based on called variants
pipe output through remove_trailing_fasta_Ns to strip N characters from beginning and ends of seqs
and then through a sed to append new_X_draft_sequence to name of fasta record
*/
process bcftoolsconsensus {

container 'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0'

publishDir "$params.outdir"
  
input: 
path mask
path best10refseq
path cons

output: 
path "${cons}_consensus_seqs.fa", emit: fa

script: 
"""
bcftools consensus -m ${mask} -f ${best10refseq} ${cons} -o ${cons}_consensus_seqs.fa
"""
}

/*
workflow {
	mask_channel=Channel.fromPath(mask/*.mask)
	best10refseq_channel=Channel.fromPath("best10refseq/*_best10_refseq.fa").collect()
	cons_channel=Channel.fromPath(cons/*.cons.vcf.gz)
	bcftoolsconsensus(mask_channel, best10refseq_channel, cons_channel)
}
*/


/*
code that still needs to be ran: 

bcftools consensus -m "${mask}" -f "${best10refseq}" "${cons}" | 
remove_trailing_fasta_Ns | 
sed "/^>/s/$/_new_${cons}_draft_sequence/" > ${cons}_new_draft_seqs.fa


Finally, re-run minimap against new draft seqs
re-minimap data against best 10 BTV ref seqs
./run_bt_align_paired $fastq_r1 $fastq_r2 ${sample_id}_new_draft_seqs_index

minimap2 -ax map-ont ${barcode_var}_alignment_best10_refseq.fa ${barcode_var} > ${barcode_var}_best10_alignment_consensus.sam

process output again
samtools view -F4 -S -b ${barcode_var}_best10_alignment_consensus.sam > ${barcode_var}_best10_alignment_consensus.unsorted.bam
samtools sort ${barcode_var}_best10_alignment_consensus.unsorted.bam -o ${barcode_var}_best10_alignment_consensus.sorted.bam -O bam
samtools index ${barcode_var}_best10_alignment_consensus.sorted.bam
*/