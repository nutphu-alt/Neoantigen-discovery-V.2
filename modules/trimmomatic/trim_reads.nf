// module to Trim Reads
process trim_reads {

    tag "${sample_type}_${sample_id}"
    publishDir "${params.outdir}/trimmed_reads/${sample_type}", mode: 'copy'

    conda 'modules/trimmomatic/trimmomatic_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(reads)
    path(adapters)

    output:
    tuple val(sample_type), val(sample_id), file("trimmed_${sample_id}_1_paired.fastq.gz"), file("trimmed_${sample_id}_2_paired.fastq.gz")

    script:
    """
    trimmomatic ${params.trimmomatic} \
        ${reads[0]} ${reads[1]} \
        trimmed_${sample_id}_1_paired.fastq.gz \
        trimmed_${sample_id}_1_unpaired.fastq.gz \
        trimmed_${sample_id}_2_paired.fastq.gz \
        trimmed_${sample_id}_2_unpaired.fastq.gz \
        ILLUMINACLIP:${adapters}:${params.trimmomatic_clip} \
        ${params.trimmomatic_qual}
    """
}
