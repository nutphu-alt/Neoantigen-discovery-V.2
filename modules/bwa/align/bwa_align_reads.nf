// Module to Align Reads to Reference Genome Using BWA
process bwa_align_reads {  

    tag "${sample_type}_${sample_id}"
    publishDir "${params.outdir}/aligned_bams/${sampe_type}", mode: 'copy'

    conda 'modules/bwa/align/bwa_env.yaml'

    input:
    tuple val(sampe_type), val(sample_id), path(reads)
    path(bwa_index_files)

    output:
    tuple val(sample_type), val(sample_id), file("aligned_${sample_id}_sorted.bam")

    script:
    """
    bwa mem ${bwa_index_files} ${reads[0]} ${reads[1]} > aligned_${sample_id}.sam
    
    samtools view -bS aligned_${sample_id}.sam -o aligned_${sample_id}.bam

    samtools sort aligned_${sample_id}.bam -o aligned_${sample_id}_sorted.bam

    samtools index aligned_${sample_id}_sorted.bam
    """
}