process tpm_calculation {

    tag "${sample_id}"
    publishDir "${params.outdir}/tpm_calculation", mode: 'copy'

    conda 'modules/tpm_calculation/r_env.yml'

    input:
    tuple val(sample_id), path(count_file)

    output:
    tuple val(sample_id), path("${sample_id}_tpm.csv")

    script:
    """
    Rscript -e "
    data <- read.csv("${count_file}", header=TRUE, skip=1)
    
    colnames(data) <- c("GeneID", "Chr", "Start", "End", "Strand", "Length", "RNA1", "RNA2")

    rpk1 <- data\$RNA1 / ((as.numeric(data\$Length) / 1000))
    rpk2 <- data\$RNA2 / ((as.numeric(data\$Length) / 1000))

    data\$TPM_RNA1 <- (rpk1 / sum(rpk1)) * 1e6
    data\$TPM_RNA2 <- (rpk2 / sum(rpk2)) * 1e6

    data\$mean_TPM <- apply(data[, c("TPM_RNA1", "TPM_RNA2")], 1, mean)

    write.csv(data, file = "${sample_id}_tpm.csv", row.names = FALSE)
    "
    """
}