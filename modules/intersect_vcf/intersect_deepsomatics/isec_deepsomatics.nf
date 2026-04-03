// Intersect variant calls from Deepsomatics
process intersect_deepsomatics {

    tag "Intersecting Deepsomatics VCFs"

    publishDir "${params.outdir}/variants/deepsomatics/intersect", mode: 'copy'

    conda 'modules/intersect_vcf/intersect_somaticsniper/bcftools_env.yaml'

    input:
    tuple val(tumor_type), val(tumor_id), path(tumor_vcf)
    tuple val(normal_type), val(normal_id), path(normal_vcf)

    output:
    tuple val(tumor_id), file("all_tumor_deepsomatics.vcf.gz")

    script:
    """
    bcftools isec -n+4 -c none -w1 -O z -o all_tumor_deepsomatics.vcf.gz \
    ${tumor_vcf[0]} \
    ${tumor_vcf[1]} \
    ${tumor_vcf[2]} \
    ${tumor_vcf[3]} \

    bcftools index all_tumor_deepsomatics.vcf.gz

    bcftools isec -C -c none -w1 -O z -o tumor_specific_deepsomatics.vcf.gz \
    all_tumor_deepsomatics.vcf.gz \
    ${normal_vcf[0]} \
    ${normal_vcf[1]}
    """
}