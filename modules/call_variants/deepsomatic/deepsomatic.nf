// Module to Call Variants Using deepsomatic
process deepsomatic_call {

    tag "${sample_type}_${sample_id}"

    publishDir "${params.outdir}/variants/deepsomatic/${sample_type}", mode: 'copy'

    container 'google/deepsomatic:1.10.0'

    input:
    tuple val(sample_type), val(sample_id), path(recalibrated_bam)
    path(reference)

    output:
    tuple val(sample_type), val(sample_id), file("${sample_id}_deepsomatics_variants_snp.vcf")

    script:
    """
    run_deepsomatic \
    ${params.deepsomatic} \
    --ref=${reference} \
    --reads_tumor=${recalibrated_bam} \
    --output_vcf=${sample_id}_deepsomatics_variants_snp.vcf.gz \
    --sample_name_tumor="${sample_id}"

    bcftools index ${sample_id}_deepsomatics_variants_snp.vcf.gz
    """
}