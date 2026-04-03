process annotation_biomart {

    tag "${sample_id}"
    publishDir "${params.outdir}/annotation_biomart", mode: 'copy'

    conda 'modules/biomart/r_env_biomart.yml'

    input:
    tuple val(sample_id), path(tpm_file)

    output:
    tuple val(sample_id), path("${sample_id}_tpm_biomart_annotation.csv")

    script:
    """
    Rscript -e "
    library(biomartr)
    library(biomaRt)

    data <- read.csv("${tpm_file}", header=TRUE)
    attributes <- c("ensembl_gene_id_version", "mgi_symbol")

    if ${params.species} == 'Mus_musculus' {
        species <- 'mmusculus_gene_ensembl'
    } else if ${params.species} == 'Homo_sapiens' {
        species <- 'hsapiens_gene_ensembl'
    } else {
        stop("Unsupported species: ${params.species}")
    }

    mart = useMart("ensembl", dataset = species)

    annotations <- getBM(values = data\$Geneid, 
    mart = mart, 
    attributes = attributes, 
    filters = "ensembl_gene_id_version"
    )

    data_annotated <- merge(data, annotations[, c("ensembl_gene_id_version", attributes)], 
    by.x = "Geneid", by.y = "ensembl_gene_id_version", all.x = TRUE)

    matched <- data_annotated\$Geneid == data_annotated\$ensembl_gene_id_version.1

    # Check if there are NA values in the matched vector
    if (any(is.na(matched))) {
    message("There are NA values in the 'matched' vector. Investigate why.")
    na_rows <- which(is.na(matched))
    message(paste("Rows with NA values:", paste(na_rows, collapse = ", ")))
    }

    # Summarize the results
    if (all(matched, na.rm = TRUE)) {
    message("All data in column Geneid and ensembl_gene_id_version.1 are exactly the same.")
    } else {
    message("There are rows where data in column Geneid and ensembl_gene_id_version.1 differ.")
    # Optionally, find the rows with mismatches
    mismatched_rows <- which(!matched & !is.na(matched))
    message(paste("Rows with mismatches:", paste(mismatched_rows, collapse = ", ")))
    }

    # Select column to remove
    library(dplyr)

    data_annotated <- data_annotated %>%
    select(-Geneid, 
    -ensembl_gene_id_version.1)

    # Filter mean TPM >= 10
    data_annotated_TPM_filter <- data_annotated[data_annotated\$mean_TPM >= 10, ]

    write.csv(data_annotated_TPM_filter, "${sample_id}_tpm_biomart_annotation.csv", row.names =  FALSE)

    # Check if there are any duplicated rows
    any_duplicated_A <- anyDuplicated(data_annotated\$mgi_symbol)
    if (any_duplicated_A) {  
    print(paste("There are duplicates in column A starting at row", any_duplicated_A))
    } else {
    print("There are no duplicates in column A")
    }

    # Count and display duplicated rows
    count_A <- table(data_annotated\$mgi_symbol)
    duplicates_A <- subset(count_A, count_A > 1)
    print(duplicates_A)

    duplicate_rows <- data_annotated[duplicated(data_annotated\$mgi_symbol), ]
    print(duplicate_rows)
    )
    "
    """
}