#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CONCAT_FASTQ } from './modules/local/concatfastq/main'
include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED } from './modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_RAW } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_TRIMMED } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_ALIGN } from './modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_MARKDUP } from './modules/nf-core/multiqc/main'
include { CUTADAPT } from './modules/nf-core/cutadapt/main'
include { STAR_GENOMEGENERATE } from './modules/nf-core/star/genomegenerate/main'
include { EXTRACT_CHROM_SIZES } from './modules/local/extract_chrom_sizes/main'
include { STAR_ALIGN } from './modules/nf-core/star/align/main'
include { IGVTOOLS_TOTDF } from './modules/local/igvtools/totdf/main'
//include { PICARD_MARKDUPLICATES } from './modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
include { BEDTOOLS_BAMTOBED } from './modules/nf-core/bedtools/bamtobed/main'
include { BED_SORT } from './modules/local/bed_sort/main'
include { BEDTOOLS_GROUPBY as BEDTOOLS_COLLAPSE } from './modules/nf-core/bedtools/groupby/main'
include { BEDTOOLS_BEDTOBAM } from './modules/local/bedtools_bedtobam/main'
include { MACS3_CALLPEAK } from './modules/nf-core/macs3/callpeak/main'
include { IDR } from './modules/nf-core/idr/main'


def parse_samplesheet(csv_path) {
    Channel
        .fromPath(csv_path)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id: row.sample_id,
                group: row.group,
                condition: row.condition,
                single_end: true
            ]
            def fastq = row.fastq.split(';').collect { file(it.trim()) }
            [ meta, fastq ]
        }
}


workflow {
    ch_fastq = parse_samplesheet(params.samplesheet)

    ch_fastq = CONCAT_FASTQ(ch_fastq)
    
    FASTQC_RAW(ch_fastq)
    
    ch_multiqc_files = FASTQC_RAW.out.zip
        .map { meta, zip -> zip }
        .collect()
        .map { files -> [ [id: 'ALL'], files, [], [], [], [] ] }

    MULTIQC_RAW(ch_multiqc_files)
    
    CUTADAPT(ch_fastq)
    
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
    EXTRACT_CHROM_SIZES(STAR_GENOMEGENERATE.out.index)

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
    
    // Don't need to use the reference so pass placeholder as a second argument
    //ch_dup = PICARD_MARKDUPLICATES(STAR_ALIGN.out.bam_sorted_aligned, [ '', [], [] ])
    
    //ch_multiqc_files = ch_dup.metrics
    //    .map { meta, metrics -> metrics }
    //    .collect()
    //    .map { files -> [ [id: 'ALL'], files, [], [], [], [] ] }

    //MULTIQC_MARKDUP(ch_multiqc_files)

    //ch_bam = ch_dup.bam 
    ch_bam = STAR_ALIGN.out.bam_sorted_aligned
    //    .map { meta, bam -> [ meta.group, meta.condition, bam ] }
    //    .groupTuple(by: [0, 1])
    //    .map { group, condition, bams -> [ [ id: group, condition: condition ], bams, [] ] }

    //SAMTOOLS_MERGE(ch_bam, [ '', [], [], [] ])
    
    //BEDTOOLS_BAMTOBED(SAMTOOLS_MERGE.out.bam)
    BEDTOOLS_BAMTOBED(ch_bam)
    BED_SORT(BEDTOOLS_BAMTOBED.out.bed)
    BEDTOOLS_COLLAPSE(BED_SORT.out.bed, '4')
    BEDTOOLS_BEDTOBAM(BEDTOOLS_COLLAPSE.out.bed, EXTRACT_CHROM_SIZES.out.chrom_sizes)

    //ch_bam_ip = BEDTOOLS_BEDTOBAM.out.bam
    //    .filter { meta, bed -> meta.condition == 'IP' }
    
    //ch_bam_input = BEDTOOLS_BEDTOBAM.out.bam
    //    .filter { meta, bed -> meta.condition == 'INP' }

    //ch_bam = ch_bam_ip.combine(ch_bam_input)
    //    .map { meta, ipbam, meta2, controlbam -> [ [ id: meta.id, single_end: true ], ipbam, controlbam ] }
    
    ch_bam = BEDTOOLS_BEDTOBAM.out.bam
        .map { meta, bam -> [ meta, bam, [] ] }
    MACS3_CALLPEAK(ch_bam, '1.5e+8')

    //ch_idr = MACS3_CALLPEAK.out.peak
    //    .map { meta, peaks -> peaks }
    //    .collect()
    //    .map { all_peaks -> [ [ id: 'merged' ], all_peaks, 'narrowPeak' ] }

    //IDR(ch_idr)
}
