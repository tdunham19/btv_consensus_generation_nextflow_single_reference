process IDENTIFY_BEST_SEGMENTS_FROM_SAM {
    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
	// container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //   'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
    //    'quay.io/biocontainers/ubuntu:20.04' }"
        
	input: 
    tuple val(meta), path(input)
    tuple val(meta2), path(reference)

	output: 
	tuple val(meta), path("*.fa") , emit: fa

	script: 
	"""
	${params.script_dir}/identify_best_segments_from_sam ${reference} ${input} > ${meta.id}_best10_refseq.fa

cat <<-END_VERSIONS > versions.yml
    "!{task.process}":
    END_VERSIONS
	"""
	}