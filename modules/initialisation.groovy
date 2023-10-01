#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process help {
    label 'linux'
    input:
        val errorMessage
    when:
        params.help
    exec:
        "cat $errorMessage".execute().text.readLines().each{println it}
}

process initializeFolders {
    label 'linux'
    input:
        val results
    when:
        !params.help
    exec:
        file("$results/FASTQ/RAW").mkdirs()
        file("$results/FASTQ/FILTERED").mkdirs()
        file("$results/BAM").mkdirs()
        file("$results/REFERENCE_FILES").mkdirs()
        file("$results/COUNTING").mkdirs()
}

workflow initialisation {
    main:
        results = file(params.results)

        errorMessage = file("./bin/help.txt")
        help(errorMessage)
        initializeFolders(results)
}