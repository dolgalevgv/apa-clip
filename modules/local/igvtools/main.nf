process IGVTOOLS_TOTDF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/igvtools:2.17.3--99f2430c080b7bdb'

    input:
    tuple val(meta), path(wig)
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}.tdf"), emit: tdf

    script:
    def wig_unique = wig.filter { !it.contains('UniqueMultiple') }
    """
    igvtools toTDF $wig_unique ${meta.id}.tdf $genome
    """
}
