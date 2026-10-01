process EXTRACT_CHROM_SIZES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(star_index)

    output:
    tuple val(meta), path("chrom.sizes"), emit: chrom_sizes

    script:
    """
    cp ${star_index}/chrNameLength.txt chrom.sizes
    """
}
