process design {
    // publishDir "${params.results}/${name}", pattern: '*.json', mode: 'copy', overwrite: true
    publishDir "${params.results}/${name}", pattern: '*.txt', mode: 'copy', overwrite: true
    // echo true
    errorStrategy 'ignore'
    memory '8 GB'
    cpus 8

    input:
        tuple(
            val(name),
            val(method),
            val(variant),
            path(config),
            val(data),
            )

    output:
        tuple(val(name), val(method), path("${name}.json"), path("${name}.tsv"), path("${name}.fna"), path('log.txt'))

    script:
    """
    design.py -n ${name} -p ${name} -q '${variant}' -m ${method} -c ${config} -d ${data} 2> log.txt
    # Note the ' around the variant; otherwise will fail bc/ eg
    # NM_032682.6:c.484C>T would write to file T
    """
}


process pcr {
    // publishDir "${params.results}", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 8
    container 'nanozoo/ispcr:33--2df9365'

    input:
        path(primers)

    output:
        path('ispcr.fna')

    """
    cat ${primers} > primers.tsv
    isPcr /Users/phi/tmp/primer/GRCh37_latest_genomic.fna primers.tsv -maxSize=4000 -minPerfect=12 -minGood=12 -out=fa ispcr.fna
    """
}


process blast {
    publishDir "${params.results}/${name}", mode: 'copy', overwrite: true
    cpus 8

    input:
        tuple(val(name), path(primers), val(blastdb))

    output:
        tuple(val(name), path('blast.tsv'))

    """
    blastn -task blastn -num_threads ${task.cpus} -query ${primers} -db ${blastdb} -evalue 100 -dust no -word_size 7 -outfmt 6 > blast.tsv
    """
}


process pseudo {
    input:
        tuple(val(name), val(method), path(candidates), path(annealing))

    output:
        path("${name}.tsv")

    """
    pseudo.py -m ${method} -c ${candidates} -a ${annealing} -o ${name}.tsv
    """
}


process cat {
    publishDir "${params.results}", mode: 'copy', overwrite: true

    input:
        path(results)

    output:
        path("summary.tsv")

    """
    cat ${results} > summary.tsv
    """
}