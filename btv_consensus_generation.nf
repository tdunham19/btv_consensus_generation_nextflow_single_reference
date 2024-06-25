include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_1   } from './modules/local/minimap2/align/main1.nf'
include { SAMTOOLS_VIEW	             		   } from './modules/local/samtools/view/main.nf'
include { SAMTOOLS_SORT               		   } from './modules/local/samtools/sort/main.nf'
include { SAMTOOLS_INDEX              		   } from './modules/local/samtools/index/main.nf'
include { IDENTIFY_BEST_SEGMENTS_FROM_SAM      } from './modules/local/identify_best_segments_from_sam/main.nf'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_2   } from './modules/local/minimap2/align/main2.nf'


workflow BTV_CONSENSUS {

  ch_versions = Channel.empty()                                               

  // fastq input files

  Channel.fromFilePairs("${params.fastq_dir}/${params.input_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
  .map{ name, reads ->
         def matcher = name =~ /FABADRU_\d+/
         def meta = [:]
         if (matcher.find()) {
             meta.id = matcher.group(0)
         } else {
             meta.id = "UNKNOWN" 
         }
         [ meta, reads ]
     }
  .set { ch_reads }
  
  // refseq input files

    Channel.fromPath("${params.reference_dir}")
    .collect()
    .map { reference ->
            def meta2 = [:]
            meta2.id = "reference"
            [meta2, reference]
        }
    .set { ch_reference }
    
  // run minimap2 on input reads
  MINIMAP2_ALIGN_1 ( ch_reads, ch_reference )
  
  // run samtools to process minimap2 alignment - view
  SAMTOOLS_VIEW ( MINIMAP2_ALIGN_1.out.sam )
  
  // run samtools to process minimap2 alignment - sort
  SAMTOOLS_SORT ( SAMTOOLS_VIEW.out.bam )
  
  // run samtools to process minimap2 alignment - index
  SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
  
  // extract new fasta file containing best aligned-to seqs for this dataset
  IDENTIFY_BEST_SEGMENTS_FROM_SAM ( MINIMAP2_ALIGN_1.out.sam, ch_reference )
  
  // re-minimap data against best 10 BTV ref seqs.
  MINIMAP2_ALIGN_2 ( IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa, ch_reads )
  
  }
  
  // specify the entry point for the workflow
workflow {
  BTV_CONSENSUS()
}
