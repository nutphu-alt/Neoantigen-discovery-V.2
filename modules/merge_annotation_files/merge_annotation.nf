process merge_annotation {
    
    tag "Merging Biomart and VEP annotations"

    publishDir "${params.outdir}/neoantigen_candidates", mode: 'copy'

    conda 'modules/merge_annotation_files/r_env_merge.yml'

    input:
    tuple val(sample_id), path(biomart_file)
    tuple val(sample_id), path(vep_file)

    output:
    tuple val(sample_id), file("${sample_id}_neoantigen_candidates.csv")

    script:
    """
    Rscript -e "
    library(tidyverse)

    biomart_data <- read.csv("${biomart_file}", header=TRUE)
    vep_data <- read.delim("${vep_file}", header=TRUE)

    # Rmove unwanted column
    biomart_data <- biomart_data %>%
    select(-1, -2, -3, -4, -5)

    symbol_vep <- vep_data\$SYMBOL
    symbol_biomart <- biomart_data\$mgi_symbol

    vep_filter <- vep_data[vep_data\$SYMBOL %in% biomart_data\$mgi_symbol,]
    biomart_filter <- biomart_data[biomart_data\$mgi_symbol %in% vep_data\$SYMBOL,]

    merged_data <- merge(vep_filter, biomart_filter, by.x = "SYMBOL", by.y = "mgi_symbol", all.x = TRUE)

    # Check for exact match between column A and C for each row
    matched <- merged_data\$SYMBOL == merged_data\$mgi_symbol

    # Summarize the results
    if (all(matched)) {
    message("All data in column Geneid and ensembl_gene_id_version.1 are exactly the same.")
    } else {
    message("There are rows where data in column Geneid and ensembl_gene_id_version.1 differ.")
    # Optionally, find the rows with mismatches
    mismatched_rows <- which(!matched)
    message(paste("Rows with mismatches:", mismatched_rows, sep = ", "))
    }

    # Find duplicated rows in mutation_filter\$SYMBOL
    duplicated_rows <- which(duplicated(vep_filter\$SYMBOL) | duplicated(vep_filter\$SYMBOL, fromLast = TRUE))

    # Output the indices of duplicated rows
    print(duplicated_rows)

    write.csv(merged_data, file = "${sample_id}_neoantigen_candidates.csv", row.names = FALSE)
    "
    """
}