process BWA_INDEX {
    tag "$ref.name"

    input:
    path ref

    output:
    path "${ref}*", includeInputs: true

    script:
    """
    bwa index ${ref}
    """
}

process FAIDX {
    tag "$ref.name"

    input:
    path ref

    output:
    path "${ref}.fai"

    script:
    """
    samtools faidx ${ref}
    """
}

process ALIGN {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)
    path ref_files

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    def ref = ref_files[0]
    """
    bwa mem -t ${task.cpus} -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA" \
        ${ref} ${reads[0]} ${reads[1]} \
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    """
}

process MARK_DUPLICATES {
    tag "$sample_id"
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai")
    path "${sample_id}.dup_metrics.txt"

    script:
    """
    samtools collate -@ ${task.cpus} -O ${sorted_bam} | \
        samtools fixmate -m -@ ${task.cpus} - - | \
        samtools sort -@ ${task.cpus} -o fixed.sorted.bam -
    samtools markdup -@ ${task.cpus} -s -f ${sample_id}.dup_metrics.txt \
        fixed.sorted.bam ${sample_id}.dedup.bam
    samtools index ${sample_id}.dedup.bam
    rm fixed.sorted.bam
    """
}
