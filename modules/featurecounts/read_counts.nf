// Calculate read counts using featureCounts
process read_counts {

    tag "${sample_id}"
    publishDir "${params.outdir}/read_counts", mode: 'copy'

    conda 'modules/featurecounts/featurecounts_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(bam)
    path(gtf)

    output:
    tuple val(sample_type), val(sample_id), path("${sample_id}_read_counts.txt")

    script:
    """
    featureCounts ${params.featcounts} \
    -a ${gtf} \
    -o ${sample_id}_read_counts.txt \
    ${sample_id[0]}
    ${sample_id[1]}
    """
}