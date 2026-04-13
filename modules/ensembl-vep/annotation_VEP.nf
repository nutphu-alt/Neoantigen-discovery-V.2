// Annotation variants using Ensembl VEP
process annotation_VEP {

    tag "Annotating VCFs with Ensembl VEP"

    publishDir "${params.outdir}/variants/annotated", mode: 'copy'
    
    container 'ensemblorg/ensembl-vep'

    input:
    tuple val(sample_id), path(vcf)
    path(vep_cache)

    output:
    tuple val(sample_id), file("${sample_id}_missense.txt")

    script:
    """
    vep \
    ${params.VEP} \
    -i ${vcf} \
    -o ${sample_id}_missense.txt
    """