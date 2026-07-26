// Module to Build Reference Index for STAR Aligner
process star_index {  

    tag "star_build_index"
    publishDir "${params.outdir}/reference_index/star", mode: 'copy'

    conda 'modules/star/index_building/star_env.yaml'

    input:
    path(reference)
    path(gtf)
    
    output:
    path("star_index")

    script:
    """
    STAR --runMode genomeGenerate \
    --genomeDir star_index \
    --genomeFastaFiles ${reference} \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang ${params.read_length - 1} \
    """
}