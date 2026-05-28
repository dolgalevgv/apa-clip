#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from './modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_RAW } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_TRIMMED } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_ALIGN } from './modules/nf-core/multiqc/main'
include { CUTADAPT } from './modules/nf-core/cutadapt/main'
include { STAR_GENOMEGENERATE } from './modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main'
include { IGVTOOLS_TOTDF } from './modules/local/igvtools/totdf/main.nf'


def parse_samplesheet(csv_path) {
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                donor: row.donor,
                condition: row.condition,
                single_end: true
            ]
            def reads = [ file(row.fastq_1) ]
            [ meta, reads ]
        }
}


workflow {
    ch_reads = parse_samplesheet(params.samplesheet)

    FASTQC_RAW(ch_reads)
    
    ch_multiqc_files = FASTQC_RAW.out.zip
        .map { meta, zip -> zip }
        .collect()
        .map { files -> [ [id: 'ALL'], files, [], [], [], [] ] }

    MULTIQC_RAW(ch_multiqc_files)
    
    CUTADAPT(ch_reads)
    
    ch_trimmed = CUTADAPT.out.reads

    FASTQC_TRIMMED(ch_trimmed)
    
    ch_multiqc_files = FASTQC_TRIMMED.out.zip
        .map { meta, zip -> zip }
        .collect()
        .map { files -> [ [id: 'ALL'], files, [], [], [], [] ] }
    
    MULTIQC_TRIMMED(ch_multiqc_files)

    ch_fasta = Channel.value([ [ id: 'genome' ], file(params.genome_fasta) ])
    ch_gtf = Channel.value([ [ id: 'genome' ], file(params.genome_gtf) ])

    STAR_GENOMEGENERATE(ch_fasta, ch_gtf)

    STAR_ALIGN(ch_trimmed, STAR_GENOMEGENERATE.out.index, ch_gtf, false)

    ch_multiqc_files = STAR_ALIGN.out.log_final
        .map { meta, log -> log }
        .collect()
        .map { files -> [ [id: 'ALL'], files, [], [], [], [] ] }

    MULTIQC_ALIGN(ch_multiqc_files)
    
    // STAR produces two .wig files, but we only need the first one (unique alignments)
    ch_wig = STAR_ALIGN.out.wig
        .map { meta, wigs -> [ meta, wigs[0] ] }

    IGVTOOLS_TOTDF(ch_wig, ch_fasta)
}
