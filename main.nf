#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from './modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_RAW } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_TRIMMED } from './modules/nf-core/multiqc/main'
include { TRIMGALORE } from './modules/nf-core/trimgalore/main'


def parse_samplesheet(csv_path) {
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                donor: row.donor,
                condition: row.condition
            ]
            def reads = [ file(row.fastq_1), file(row.fastq_2) ]
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

    TRIMGALORE(ch_reads)
    
    ch_trimmed = TRIMGALORE.out.reads

    FASTQC_TRIMMED(ch_trimmed)
    
    ch_multiqc_files = FASTQC_TRIMMED.out.zip
        .map { meta, zip -> zip }
        .collect()
        .mix(
            TRIMGALORE.out.log
                .map { meta, log -> log }
                .collect()
        )
        .collect()
        .map { files -> [ [id: 'ALL'], files, [], [], [], [] ] }
    
    MULTIQC_TRIMMED(ch_multiqc_files)
}
