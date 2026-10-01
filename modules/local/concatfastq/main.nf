process CONCAT_FASTQ {
    tag "$meta.id"

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("${meta.id}.merged.fastq.gz"), emit: reads
    tuple val("${task.process}"), val('coreutils'), eval("cat --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_coreutils

    script:
    """
    cat ${fastqs} > ${meta.id}.merged.fastq.gz
    """
}
