// Module to Call Variants Using Bcftools mpileup
process bcftools_mpileup {

    tag "${sample_type}_${sample_id}"

    publishDir "${params.outdir}/variants/bcftools/${sample_type}", mode: 'copy'

    conda 'modules/call_variants/bcftools/bcftools_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(recalibrated_bam)
    path(reference)

    output:
    tuple val(sample_type), val(sample_id), file("${sample_id}_bcftools_variants_filtered_snp.vcf.gz")

    script:
    """
    bcftools mpileup ${params.bcftools} -f ${reference} ${recalibrated_bam} | bcftools call -mv -Ov -o ${sample_id}_bcftools_variants.vcf

    bcftools filter -i 'QUAL>=30' ${sample_id}_bcftools_variants.vcf -o ${sample_id}_bcftools_variants_filtered.vcf

    bcftools view -i 'TYPE="snp"' ${sample_id}_bcftools_variants_filtered.vcf -o ${sample_id}_bcftools_variants_filtered_snp.vcf

    bgzip ${sample_id}_bcftools_variants_filtered_snp.vcf

    bcftools index ${sample_id}_bcftools_variants_filtered_snp.vcf.gz
    """
}