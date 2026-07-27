#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Parameters ---
params.mode                     = params.mode ?: 3
params.tumor_rna_reads          = params.tumor_rna_reads ?: 'data/raw_reads/tumor_rna/*.fastq.gz'
params.tumor_dna_reads          = params.tumor_dna_reads ?: 'data/raw_reads/tumor_dna/*.fastq.gz'
params.normal_dna_reads         = params.normal_dna_reads ?: 'data/raw_reads/normal_dna/*.fastq.gz'
params.outdir                   = params.outdir ?: 'results'
params.reference                = params.reference ?: 'data/reference/*.fa'
params.adapters                 = params.adapters ?: 'data/adapters/*.fa'
params.gtf                      = params.gtf ?: 'data/reference/*.gtf'
params.star_index               = params.star_index ?: "results/reference_index/star/*"
params.bwa_index                = params.bwa_index ?: "results/reference_index/bwa/*"
params.known_snps               = params.known_snps ?: 'data/reference/*.vcf.gz'
params.known_indels             = params.known_indels ?: 'data/reference/*.vcf.gz'
params.vep_cache                = params.vep_cache ?: 'data/vep_cache/*'


// --- Import modules ---
include { fastqc_raw_reads }                        from './modules/fastqc/fastqc_raw_reads/fastqc_raw_reads.nf'
include { trim_reads }                              from './modules/trimmomatic/trim_reads.nf'
include { fastqc_trimmed_reads }                    from './modules/fastqc/fastqc_trimmed_reads/fastqc_trimmed_reads.nf'
include { star_index }                              from './modules/star/index_building/star_build_index.nf'
include { star_align_reads }                        from './modules/star/align/star_align_reads.nf'
include { bwa_index }                               from './modules/bwa/index_building/bwa_build_index.nf'
include { bwa_align_reads as bwa_tumor
          bwa_align_reads as bwa_normal }           from './modules/bwa/align/bwa_align_reads.nf'
include { remove_duplicate }                        from './modules/picard/remove_duplicate.nf'
include { BQSR }                                    from './modules/base_score_recalibration/BQSR.nf'
include { bcftools_mpileup }                        from './modules/call_variants/bcftools/bcftools_mpileup.nf'
include { gatk_mutect2 }                            from './modules/call_variants/gatk/gatk_mutect2.nf'
include { deepsomatic_call }                        from './modules/call_variants/deepsomatic/deepsomatic.nf'
include { intersect_bcftools }                      from './modules/intersect_vcf/intersect_bcftools/isec_bcftools.nf'
include { intersect_gatk }                          from './modules/intersect_vcf/intersect_gatk/isec_gatk.nf'
include { intersect_deepsomatic }                   from './modules/intersect_vcf/intersect_deepsomatic/isec_deepsomatic.nf'
include { intersect_all_tools }                     from './modules/intersect_vcf/intersect_all_tools/isec_all.nf'
include { cache_download }                          from './modules/ensembl-vep/cache_download.nf'
include { annotation_VEP }                          from './modules/ensembl-vep/annotation_VEP.nf'
include { read_counts }                             from './modules/featurecounts/read_counts.nf'
include { tpm_calculation }                         from './modules/tpm_calculation/tpm_calculation.nf'
include { annotation_biomart }                      from './modules/biomart/annotation_biomart.nf'
include { merge_annotation }                            from './modules/merge_annotation/merge_annotation.nf'


//include { predict_neoantigens }   from './modules/predict_neoantigens.nf'



// --- Workflow definition ---
workflow {

    // Input channels
    tumor_rna_ch    = Channel.fromFilePairs(params.tumor_rna_reads)
                        .map { sample_id, reads -> tuple('tumor_rna', sample_id, reads) }
    tumor_dna_ch    = Channel.fromFilePairs(params.tumor_dna_reads)
                        .map { sample_id, reads -> tuple('tumor_dna', sample_id, reads) }
    normal_dna_ch   = Channel.fromFilePairs(params.normal_dna_reads)
                        .map { sample_id, reads -> tuple('normal_dna', sample_id, reads) }
    reads_ch        = tumor_rna_ch.concat(tumor_dna_ch, normal_dna_ch) // Combine all read channels

    reference       = file(params.reference)
    adapters        = file(params.adapters)
    gtf             = file(params.gtf)
    star_index      = file(params.star_index)
    bwa_index       = file(params.bwa_index)
    known_snps      = file(params.known_snps)
    known_indels    = file(params.known_indels)
    vep_cache       = file(params.vep_cache)


    // Workflow steps

    def executionMode = params.mode.toInteger()

    // Mode 1: FastQC on raw reads

    if(executionMode == 1){

        log.info "=========================================="
        log.info "=========================================="
        log.info "                  MODE 1                 "
        log.info "           FastQC on raw reads            "
        log.info "=========================================="
        log.info "=========================================="

        fastqc_raw = fastqc_raw_reads(reads_ch)

    }


    // Mode 2: Trim and FastQC on trimmed reads

    else if(executionMode == 2){

        log.info "=========================================="
        log.info "=========================================="
        log.info "                  MODE 2                  "
        log.info "           Trim reads + FastQC            "
        log.info "=========================================="
        log.info "=========================================="

        // Trim reads
        trim = trim_reads(reads_ch, adapters)

        // FastQC on trimmed reads
        fastqc_trimmed = fastqc_trimmed_reads(trim)

    }


    // Mode 3: Neoantigen analysis

    else if(executionMode == 3){

        log.info "=========================================="
        log.info "=========================================="
        log.info "                  MODE 3                  "
        log.info "           Neoantigen analysis            "
        log.info "=========================================="
        log.info "=========================================="


        // Trim reads
        trim = trim_reads(reads_ch, adapters)

        // Separate trimmed reads by sample type
        tumor_rna_trimmed_ch  = trim.filter { sample_type, sample_id, reads -> sample_type == 'tumor_rna' }
        tumor_dna_trimmed_ch  = trim.filter { sample_type, sample_id, reads -> sample_type == 'tumor_dna' }
        normal_dna_trimmed_ch = trim.filter { sample_type, sample_id, reads -> sample_type == 'normal_dna' }


        // Build or use existing STAR index
        if (params.star_build_new_index) {
        log.info "🧬 Building new STAR index from reference and GTF"
            (star_index_files) = star_index(reference, gtf)
        } else {
            log.info "✅ Using existing STAR index: ${params.star_index}"
            star_index_files = Channel.fromPath("${params.star_index}*").collect()
        }


        // Align RNA reads
        tumor_rna_aligned_ch  = star_align_reads(tumor_rna_trimmed_ch, star_index_files)


        // Build or use existing BWA index
        if (params.bwa_build_new_index) {
        log.info "🧬 Building new BWA index from reference"
            (bwa_index_files) = bwa_index(reference)
        } else {
            log.info "✅ Using existing BWA index: ${params.bwa_index}"
            bwa_index_files = Channel.fromPath("${params.bwa_index}*").collect()
        }


        // Align DNA reads
        tumor_dna_aligned_ch  = bwa_tumor(tumor_dna_trimmed_ch, bwa_index_files)
        normal_dna_aligned_ch = bwa_normal(normal_dna_trimmed_ch, bwa_index_files)


        // Combine all bam files
        all_bam_ch = tumor_rna_aligned_ch.concat(tumor_dna_aligned_ch, normal_dna_aligned_ch)


        // Remove Duplicates
        dedup = remove_duplicate(all_bam_ch)


        // Base Quality Score Recalibration
        base_recal = BQSR(dedup, reference, known_snps, known_indels)


        // Variant Calling
        bcftools = bcftools_mpileup(base_recal.recalibrated_bam, reference)
        gatk = gatk_mutect2(base_recal.recalibrated_bam, reference)
        deepsomatic = deepsomatic_call(base_recal.recalibrated_bam, reference)

            // Separate bcftools and GATK VCFs by sample type
        collapse_bcftools_vcf_ch = bcftools.map { type, id, reads -> 
            def new_type
            if (type in ['tumor_rna', 'tumor_dna']) {
                new_type = 'tumor'
            } else if (type == 'normal_dna') {
                new_type = 'normal'
            } else {
                new_type = type
            }
            tuple(new_type, id, reads)  
        }

        tumor_bcftools_vcf_ch  = collapse_bcftools_vcf_ch.filter { type, id, reads -> type == 'tumor' }
        normal_bcftools_vcf_ch = collapse_bcftools_vcf_ch.filter { type, id, reads -> type == 'normal' }

        collapse_gatk_vcf_ch = gatk.map { type, id, reads -> 
            def new_type
            if (type in ['tumor_rna', 'tumor_dna']) {
                new_type = 'tumor'
            } else if (type == 'normal_dna') {
                new_type = 'normal'
            } else {
                new_type = type
            }
            tuple(new_type, id, reads)  
        }

        tumor_gatk_vcf_ch  = collapse_gatk_vcf_ch.filter { type, id, reads -> type == 'tumor' }
        normal_gatk_vcf_ch = collapse_gatk_vcf_ch.filter { type, id, reads -> type == 'normal' }

        collapse_deepsomatic_vcf_ch = deepsomatic.map { type, id, reads -> 
            def new_type
            if (type in ['tumor_rna', 'tumor_dna']) {
                new_type = 'tumor'
            } else if (type == 'normal_dna') {
                new_type = 'normal'
            } else {
                new_type = type
            }
            tuple(new_type, id, reads)  
        }

        tumor_deepsomatic_vcf_ch  = collapse_deepsomatic_vcf_ch.filter { type, id, reads -> type == 'tumor' }
        normal_deepsomatic_vcf_ch = collapse_deepsomatic_vcf_ch.filter { type, id, reads -> type == 'normal' }

        isec_bcftools = intersect_bcftools(tumor_bcftools_vcf_ch, normal_bcftools_vcf_ch)
        isec_gatk = intersect_gatk(tumor_gatk_vcf_ch, normal_gatk_vcf_ch)
        isec_deepsomatic = intersect_deepsomatic(tumor_deepsomatic_vcf_ch, normal_deepsomatic_vcf_ch)
        isec_all_tools_ch = isec_bcftools.mix(isec_gatk, isec_deepsomatic)
        isec_all_tools = intersect_all_tools(isec_all_tools_ch)


        // Annotation
        if (params.download_vep_cache) {
            log.info "🧬 Downloading VEP cache for ${params.species} ${params.genome_version}"
            vep_cache = cache_download()
        } else {
            log.info "✅ Using existing VEP cache: ${params.vep_cache}"
            vep_cache = file(params.vep_cache)
        }

        VEP = annotation_VEP(isec_all_tools, vep_cache)

        // Calculate read counts
        featurecounts = read_counts(tumor_rna_aligned_ch, gtf)

        tpm = tpm_calculation(featurecounts)

        biomart = annotation_biomart(tpm)

        // Predict Neoantigens

        merge = merge_annotation(biomart, VEP)

    }



}