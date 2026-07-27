process TEST {

    input:
    tuple val(type), val(id), path(reads)

    script:
    """
    echo ${type}
    echo ${id}
    ls
    """
}