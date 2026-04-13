// Module to Call Variants Using GATK MuTect2
process gatk_mutect2 {

    tag "${sample_type}_${sample_id}"

    publishDir "${params.outdir}/variants/gatk/${sample_type}", mode: 'copy'

    conda 'modules/call_variants/gatk/gatk_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(recalibrated_bam)
    path(reference)

    output:
    tuple val(sample_type), val(sample_id), file("${sample_id}_gatk_variants_filter_snp.vcf")

    script:
    """
    gatk Mutect2 \
    -I ${recalibrated_bam} \
    -O ${sample_id}_gatk_variants.vcf \
    -R ${reference} \
    -tumor ${sample_id}

    gatk FilterMutectCalls \
    -R ${reference} \
    -V ${sample_id}_gatk_variants.vcf \
    -O ${sample_id}_gatk_variants_filter.vcf

    gatk SelectVariants \
    -R ${reference} \
    -V ${sample_id}_gatk_variants_filter.vcf \
    ${params.gatk} \
    -O ${sample_id}_gatk_variants_filter_snp.vcf

    bgzip ${sample_id}_gatk_variants_filter_snp.vcf

    bcftools index ${sample_id}_gatk_variants_filter_snp.vcf.gz
    """
}