process CONCAT_FASTQ {
    tag "$meta.id"

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}.merged.fastq.gz")

    script:
    """
    cat ${fastqs} > ${meta.id}.merged.fastq.gz
    """
}
