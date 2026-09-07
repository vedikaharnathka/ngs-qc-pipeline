process QC_METRICS {
    tag "$sample_id"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(vcf), path(tbi)
    path dup_metrics
    path ref

    output:
    path "qc.db"
    path "${sample_id}.qc_metrics.json"

    script:
    """
    python3 ${projectDir}/scripts/parse_vcf_qc.py \
        --sample-id ${sample_id} \
        --bam ${bam} \
        --vcf ${vcf} \
        --dup-metrics ${dup_metrics} \
        --ref ${ref} \
        --db qc.db \
        --json ${sample_id}.qc_metrics.json \
        --min-mean-depth ${params.min_mean_depth} \
        --min-q30-pct ${params.min_q30_pct} \
        --titv-min ${params.titv_min} \
        --titv-max ${params.titv_max} \
        --min-variant-pass-pct ${params.min_variant_pass_pct} \
        --max-dup-rate-pct ${params.max_dup_rate_pct}
    """
}

process PLOT_QC {
    tag "$sample_id"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    tuple val(sample_id), path(db)

    output:
    path "*.png"

    script:
    """
    python3 ${projectDir}/scripts/plot_qc.py --db ${db} --sample-id ${sample_id} --outdir .
    """
}
