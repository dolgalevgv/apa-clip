process BED_SORT {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("${meta.id}.sorted.bed"), emit: bed

    script:
    """
    sort \
        -k1,1 -k2,2n -k3,3n -k6,6 \
        ${bed} \
        > ${meta.id}.sorted.bed
    """
}
