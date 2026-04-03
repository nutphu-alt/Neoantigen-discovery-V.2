// Intersect VCFs from all variant callers
process intersec_all_tools {

    tag "Intersecting VCFs from all variant callers"

    publishDir "${params.outdir}/variants/all_callers", mode: 'copy'

    conda 'modules/intersect_vcf/intersec_all_tools/bcftools_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(vcf)

    output:
    tuple val(sample_id), file("tumor_specific_all_callers.vcf")

    script:
    """
    for f in ${vcf}; do
        bcftools index \$f
    done

    bcftools isec -n=3 -c none -w1 -o tumor_specific_all_callers.vcf \
    ${vcf[0]} \
    ${vcf[1]} \
    ${vcf[2]}
    """
}