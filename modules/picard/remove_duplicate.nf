// Module to Remove Duplicate Reads from Bam File
process remove_duplicate {

    tag "${sample_type}_${sample_id}"
    publishDir "${params.outdir}/deduplicated_bams/${sample_type}", mode: 'copy'

    conda 'modules/picard/picard_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_type), val(sample_id), file("${sample_id}_deduplicated.bam")

    script:
    """
    picard MarkDuplicates \
    ${params.markduplicate} \
    -I=${bam} \
    -O=${sample_id}_deduplicated.bam \
    -M=${sample_id}_dedup_metrics.txt \
    ASSUME_SORTED=true

    bcftools index ${sample_id}_deduplicated.bam
    """
}