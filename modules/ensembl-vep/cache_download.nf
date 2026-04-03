// Download cache
process cache_download {

    tag "Downloading VEP cache"

    publishDir "data/vep_cache", mode: 'copy'

    container 'ensemblorg/ensembl-vep'

    output:
    path 'vep_cache/'

    script:
    """
    vep_install -a cf -s ${params.species} -y ${params.genome_version} -c vep_cache/
    """
}