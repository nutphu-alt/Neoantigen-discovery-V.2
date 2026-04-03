// Module to Calculate And Apply Base Quality Score Recalibration (BQSR)
process BQSR {  

    tag "${sample_type}_${sample_id}"
    publishDir "${params.outdir}/BQSR_table", mode: 'copy'

    conda 'modules/base_score_recalibration/gatk_env.yaml'

    input:
    tuple val(sample_type), val(sample_id), path(bam)
    path(reference)
    path(known_snps)
    path(known_indels)
    
    output:
    tuple val(sample_type), val(sample_id), file("${sample_id}_BQSR.table")
    tuple val(sample_type), val(sample_id), file("${sample_id}_recalibrated.bam"), emit: recalibrated_bam

    script:
    """
    # Add read group
    picard AddOrReplaceReadGroups \
    I=${bam} \
    O=${sample_id}_addreadgroup.bam \
    RGID=${sample_id} \
    ${params.addorreplacereadgroup} \
    RGSM=${sample_id}

    # Index known sites VCF files
    gatk IndexFeatureFile -I ${known_snps}
    gatk IndexFeatureFile -I ${known_indels}

    # Create the reference genome dictionary
    gatk CreateSequenceDictionary -R ${reference} -O ${reference}.dict

    # Create BQSR table
    gatk BaseRecalibrator \
    -I ${sample_id}_addreadgroup.bam \
    -R ${reference} \
    --known-sites ${known_snps} \
    --known-sites ${known_indels} \
    -O ${sample_id}_BQSR.table

    # Apply BQSR to BAM file (optional, uncomment if needed)
    gatk ApplyBQSR \
    -I ${sample_id}_addreadgroup.bam \
    -R ${reference} \
    --bqsr ${sample_id}_BQSR.table \
    -O ${sample_id}_recalibrated.bam

    # Index the recalibrated BAM file
    samtools index ${sample_id}_recalibrated.bam
    """
}