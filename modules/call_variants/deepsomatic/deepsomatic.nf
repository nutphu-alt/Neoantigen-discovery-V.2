// Module to Call Variants Using deepsomatics
process deepsomatics {

    tag "${sample_type}_${sample_id}"

    publishDir "${params.outdir}/variants/deepsomatics/${sample_type}", mode: 'copy'

    container 'google/deepsomatic:latest'

    input:
    tuple val(sample_type), val(sample_id), path(recalibrated_bam)
    path(reference)

    output:
    tuple val(sample_type), val(sample_id), file("${sample_id}_deepsomatics_variants_snp.vcf")

    script:
    """
    run_deepsomatic \
    --model_type=WGS_TUMOR_ONLY \
    --ref=${reference} \
    --reads_tumor=${recalibrated_bam} \
    --output_vcf=${sample_id}_deepsomatics_variants_snp.vcf.gz \
    --sample_name_tumor="${sample_id}" \
    --num_shards=$(nproc) \
    --logging_dir=logs \
    --intermediate_results_dir=intermediate_results

    bcftools index ${sample_id}_deepsomatics_variants_snp.vcf.gz
    """
}