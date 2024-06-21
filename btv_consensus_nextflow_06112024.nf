include { MINIMAP2_ALIGN              } from './modules/nf-core/modules/minimap2/align/main.nf'
include { SAMTOOLS_VIEW	              } from './modules/nf-core/modules/samtools/view/main.nf'

workflow BTV_CONSENSUS {

  ch_versions = Channel.empty()                                               

  // fastq input files

  Channel.fromFilePairs("${params.fastq_dir}/${params.input_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
  .map{ name, reads ->

         // define a new empty map named meta for each sample
         // and populate it with id and single_end values
         // for compatibility with nf-core module expected parameters
         // reads are just the list of fastq
         def meta        = [:]

         // map to input_pattern
         meta.id         = "${params.input_pattern}"

         // since there is only 1 fastq file per sample, single_end is true
         meta.single_end = true

         // this last statement in the map closure is the return value
         [ meta, reads ] }

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
  MINIMAP2_ALIGN ( ch_reads, ch_reference )
  
  // run samtools to process minimap2 alignment - view
  SAMTOOLS_VIEW ( MINIMAP2_ALIGN.out.sam )
  
  }
  
  // specify the entry point for the workflow
workflow {
  BTV_CONSENSUS()
}