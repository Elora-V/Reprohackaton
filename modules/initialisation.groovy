#!/usr/bin/env nextflow

nextflow.enable.dsl=2  /* choice of nextflow version */
/*
This Nextflow script provides help information and initializes the folder structure for storing results.
*/

/* Process for printing the help in the console  : */
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
    /* Print the content of the given file */
        "cat $errorMessage".execute().text.readLines().each{println it}
}




/* Process that create the folder for the results : */
process initializeFolders {
    label 'linux'
    output:val(params.results)
    when:
    /* Executed if the parameter help is not used */
        !params.help
    exec:
    /* Creation of all the folder */
        file("$params.results/FASTQ/RAW").mkdirs()         /* for the downloaded fastq */
        file("$params.results/FASTQ/FILTERED").mkdirs()    /* for the trimmed fastq */
        file("$params.results/BAM").mkdirs()               /* for the bam (mapping) */
        file("$params.results/REFERENCE_FILES").mkdirs()   /* for the genome and annotation */
        file("$params.results/COUNTING").mkdirs()          /* for the reads counts */
        file("$params.results/STATS").mkdirs()          /* for the stat analysis */
}


workflow initialisation {
    main:
        /* File containing the error message */
        errorMessage = file("./bin/help.txt")
        /* Call of process */
        help(errorMessage)          /* if help message needed */
        initializeFolders()         /* if help message not needed */
}