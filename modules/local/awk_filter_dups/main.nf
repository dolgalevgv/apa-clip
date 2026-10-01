process AWK_FILTER_DUPS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.id}.nodup.bed"), emit: bed

    script:
    """
    awk \
        '\$5 == 1' \
        ${bed} \
        > ${meta.id}.nodup.bed
    """
}
