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
    --cache \
    --dir_cache ${vep_cache} \
    --species ${params.species} \
    --tab \
    --offline \
    --force_overwrite \
    --symbol \
    --fork 4 \
    --filter "Consequence is missense_variant" \
    -i ${vcf} \
    -o ${sample_id}_missense.txt
    """