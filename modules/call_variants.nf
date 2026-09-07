process CALL_VARIANTS {
    tag "$sample_id"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path ref
    path ref_fai

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi")

    script:
    """
    bcftools mpileup -f ${ref} -a AD,DP ${bam} | \
        bcftools call -mv -Oz -o ${sample_id}.raw.vcf.gz
    bcftools index -t ${sample_id}.raw.vcf.gz

    bcftools filter -s LOWQUAL -e 'QUAL<20 || INFO/DP<10' \
        -Oz -o ${sample_id}.vcf.gz ${sample_id}.raw.vcf.gz
    bcftools index -t ${sample_id}.vcf.gz
    """
}
