process TABIX_BGZIP {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.19.1--h81da01d_1' :
        'biocontainers/htslib:1.19.1--h81da01d_1' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz"),     emit: gz

    script:
    """
    bgzip $input > ${meta.id}.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}