// Module to Build Reference Index for BWA Aligner
process bwa_build_index {  

    tag "bwa_build_index"
    publishDir "${params.outdir}/reference_index/bwa", mode: 'copy'

    conda 'modules/bwa/index_building/bwa_env.yaml'

    input:
    path(reference)
    

    output:
    path("bwa_index")

    script:
    """
    bwa index ${reference}
    """
}