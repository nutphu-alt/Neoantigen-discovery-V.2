// Module FastQC on Trimmed Reads
process fastqc_trimmed_reads {

    tag "${sample_type}_${sample_id}"
    
    publishDir "${params.outdir}/fastqc_reports/${sample_type}", mode: 'copy'

    conda 'modules/fastqc_trimmed_reads/fastqc_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(reads)

    output:
    tuple val(sample_type), val(sample_id), file("*_fastqc.html")
    
    script:
    """
    fastqc ${reads[0]} ${reads[1]} -o .
    """
}