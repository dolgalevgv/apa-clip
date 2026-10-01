process IGVTOOLS_TOTDF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/igvtools:2.17.3--99f2430c080b7bdb'

    input:
    tuple val(meta), path(wig)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.tdf"), emit: tdf
    tuple val("${task.process}"), val('igvtools'), eval("igvtools version 2>/dev/null | sed -n 's/^IGV Version \\([^ ]*\\).*/\\1/p'"), topic: versions, emit: versions_igvtools

    script:
    """
    igvtools toTDF $wig ${meta.id}.tdf $fasta
    """
}
