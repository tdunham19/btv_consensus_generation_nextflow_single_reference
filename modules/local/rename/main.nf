process RENAME {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.gz"),     emit: gz

    script:
    """
    cp $input ${meta.id}.cons.vcf.gz
    """
}