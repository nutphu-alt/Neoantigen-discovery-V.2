// Module to Align Reads to Reference Genome Using Star
process star_align_reads {

    tag "${sample_type}_${sample_id}"
    publishDir "${params.outdir}/aligned_bams", mode: 'copy'
    
    conda 'modules/star/align/star_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(reads)
    path(star_index_files)

    output:
    tuple val(sample_type), val(sample_id), file("${sample_type}/aligned_${sample_id}_sorted.bam")

    script:
    """
    STAR --runMode alignReads \
    --genomeDir ${star_index_files} \
    --readFilesIn ${reads[0]} ${reads[1]} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix aligned_${sample_id} \
    ${params.star}

    samtools index aligned_${sample_id}_sorted.bam
    """
}


