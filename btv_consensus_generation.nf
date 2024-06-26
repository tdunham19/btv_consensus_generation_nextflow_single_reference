include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_1   } from './modules/local/minimap2/align/main.nf'
include { SAMTOOLS_VIEW	 as SAMTOOLS_VIEW_1    } from './modules/local/samtools/view/main.nf'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_1    } from './modules/local/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_1   } from './modules/local/samtools/index/main.nf'
include { IDENTIFY_BEST_SEGMENTS_FROM_SAM      } from './modules/local/identify_best_segments_from_sam/main.nf'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_2   } from './modules/local/minimap2/align/main2.nf'
include { SAMTOOLS_VIEW	 as SAMTOOLS_VIEW_2    } from './modules/local/samtools/view/main.nf'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_2    } from './modules/local/samtools/sort/main.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_2   } from './modules/local/samtools/index/main.nf'
include { SAMTOOLS_FAIDX  					   } from './modules/local/samtools/faidx/main.nf'
include { BCFTOOLS_MPILEUP 					   } from './modules/local/bcftools/mpileup/main.nf'
include { CREATE_MASK_FILE				       } from './modules/local/create_mask_file/main.nf'

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

    Channel.fromPath("${params.reference_fasta}")
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
  SAMTOOLS_VIEW_1 ( MINIMAP2_ALIGN_1.out.sam )
  
  // run samtools to process minimap2 alignment - sort
  SAMTOOLS_SORT_1 ( SAMTOOLS_VIEW_1.out.bam )
  
  // run samtools to process minimap2 alignment - index
  SAMTOOLS_INDEX_1 ( SAMTOOLS_SORT_1.out.bam )
  
  // extract new fasta file containing best aligned-to seqs for this dataset
  IDENTIFY_BEST_SEGMENTS_FROM_SAM ( MINIMAP2_ALIGN_1.out.sam, ch_reference )
  
  // re-minimap data against best 10 BTV ref seqs.
  MINIMAP2_ALIGN_2 ( IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa, ch_reads )
  
  // run samtools to process minimap2 alignment again - view
  SAMTOOLS_VIEW_2 ( MINIMAP2_ALIGN_2.out.sam )
  
  // run samtools to process minimap2 alignment again - sort
  SAMTOOLS_SORT_2 ( SAMTOOLS_VIEW_2.out.bam )
  
   // run samtools to process minimap2 alignment - index
  SAMTOOLS_INDEX_2 ( SAMTOOLS_SORT_2.out.bam )
  
  
  // commands below create "new" consensus sequences
  
  
  // have to make a .fai file to make mpileup happy
  SAMTOOLS_FAIDX ( IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa )
  
  // mpileup calls variants -> output is a vcf file
  BCFTOOLS_MPILEUP ( SAMTOOLS_SORT_2.out.bam, IDENTIFY_BEST_SEGMENTS_FROM_SAM.out.fa )
  
  // this script creates a mask file
  // this is necessary because otherwise bcftools consensus doesn't hanlde positions with no coverage well
  CREATE_MASK_FILE ( BCFTOOLS_MPILEUP.out.vcf )

  }
  
  // specify the entry point for the workflow
workflow {
  BTV_CONSENSUS()
}
