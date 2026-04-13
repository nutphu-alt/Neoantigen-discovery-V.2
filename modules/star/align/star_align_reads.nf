// Module to Align Reads to Reference Genome Using Star
process star_align_reads {

    tag "${sample_type}_${sample_id}"
    publishDir "${params.outdir}/aligned_bams/${sampe_type}", mode: 'copy'

    conda 'modules/star/align/star_env.yaml'

    input:
    tuple val(sampe_type), val(sample_id), path(reads)
    path(star_index_files)

    output:
    tuple val(sampe_type), val(sample_id), file("aligned_${sample_id}_sorted.bam")

    script:
    """
    STAR --runMode alignReads \
    --genomeDir ${star_index_files} \ 
    --readFilesIn ${reads[0]} ${reads[1]} \ 
    --outSAMtype BAM SortedByCoordinate aligned_${sample_id}.bam\ 
    ${params.star}

    samtools index aligned_${sample_id}_sorted.bam
    """
}


