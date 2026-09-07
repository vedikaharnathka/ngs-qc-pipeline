#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc.nf'
include { BWA_INDEX; FAIDX; ALIGN; MARK_DUPLICATES } from './modules/align.nf'
include { CALL_VARIANTS } from './modules/call_variants.nf'
include { QC_METRICS; PLOT_QC } from './modules/qc.nf'

workflow {
    ref = file(params.ref)

    read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

    FASTQC(read_pairs_ch)

    BWA_INDEX(ref)
    FAIDX(ref)
    ALIGN(read_pairs_ch, BWA_INDEX.out)
    MARK_DUPLICATES(ALIGN.out)

    CALL_VARIANTS(MARK_DUPLICATES.out[0], ref, FAIDX.out)

    qc_input_ch = MARK_DUPLICATES.out[0]
        .join(CALL_VARIANTS.out)

    QC_METRICS(qc_input_ch, MARK_DUPLICATES.out[1], ref)

    db_ch = QC_METRICS.out[0]
        .combine(read_pairs_ch.map { it[0] })
        .map { db, sample_id -> tuple(sample_id, db) }

    PLOT_QC(db_ch)
}
