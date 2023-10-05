#!/usr/bin/env nextflow

nextflow.enable.dsl=2  /* choice of nextflow version */


/* Process that print the help in the console  : */
process help {
    label 'linux'
    input:
    /* File containing the message to print */
        val errorMessage
    when:
    /* This process is executed only if there is the parameter 'help' used in the line command.
    By default (when not used), help = false (see nextflow.config) */
        params.help
    exec:
    /* Print the given file */
        "cat $errorMessage".execute().text.readLines().each{println it}
}




/* Process that create the folder for the results : */
process initializeFolders {
    label 'linux'
    input:
    /* Name of the result folder (by default to "results", see nextflow.config) */
        val results
    when:
    /* Executed if the parameter help is not used */
        !params.help
    exec:
    /* Creation of all the folder */
        file("$results/FASTQ/RAW").mkdirs()         /* for the downloaded fastq */
        file("$results/FASTQ/FILTERED").mkdirs()    /* for the trimmed fastq */
        file("$results/BAM").mkdirs()               /* for the bam (mapping) */
        file("$results/REFERENCE_FILES").mkdirs()   /* for the genome and annotation */
        file("$results/COUNTING").mkdirs()          /* for the reads counts */
}




workflow initialisation {
    main:
        /* Retrieval of the path for the result folder */
        results = file(params.results)

        /* File containing the error message */
        errorMessage = file("./bin/help.txt")

        /* Call of process */
        help(errorMessage)          /* if help message needed */
        initializeFolders(results)  /* if help message not needed */
}