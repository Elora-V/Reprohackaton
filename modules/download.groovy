#!/usr/bin/env nextflow
/*
This pipeline downloads FASTQ files from SRA, as well as a reference genome and reference annotation.

*/

nextflow.enable.dsl=2 /* choice of nextflow version */

/* Process that download the fastq file */
process download_NCBI {
    label 'sratoolkit'
    input:
        /* run = list of the sra file */
        val(run)
    output:
        /* name of sra files and one downloaded sra file */
        tuple val(run), file("${run}.fastq.gz")
    publishDir path: "${params.results}/FASTQ/RAW", mode: 'copy'
    when:
        /* Executed when the fastq file aren't already downloaded (and not the help parameter) */
        !params.help && !mf.checkFile("$params.results/FASTQ/RAW", run, "q.gz")
    script:
        """
        prefetch $run # downloads all necessary files for the computer for fasterq-dump
        fasterq-dump $run --skip-technical --threads "${params.threads_download}" # extraction of fastq from SRA-accessions 
        #(with number of threads specified, default indicated in nextflow.config)
        # compressing fastq (pigz : parallel implementation of gzip)
        pigz ${run}.fastq --best
        """
}

/* Process that download the genome */
process download_reference_genome {
    label 'linux'
    output:
        /* genome file */
        file("reference.fasta")
    publishDir path: "${params.results}/REFERENCE_FILES", mode: 'copy'
    when:
        /* Executed when there is the genome parameter, and the genome isn't already downloaded (and not the help parameter) */
        !params.help && params.reference_genome != "none" && !mf.checkFile("$params.results/REFERENCE_FILES", "reference", "fasta")
    script:
        """
        # download the genome with wget (need the genome parameter from the user)
        wget -q -O reference.fasta "${params.reference_genome}"
        """
}


/* Process that download the annotation */
process download_reference_annotation {
    label 'linux'
    output:
        /* annotation file */
        file("reference.gff")
    publishDir path: "${params.results}/REFERENCE_FILES", mode: 'copy'
    when:
        /* Executed when there is the annotation parameter, and the annotation isn't already downloaded (and not the help parameter) */
        !params.help && params.reference_annotation != "none" && !mf.checkFile("${params.results}/REFERENCE_FILES", "reference", "gff")
    script:
        """
        # download the annotation with wget (need the annotation parameter from the user)
        wget -q -O reference.gff "${params.reference_annotation}"
        """
}


/* Create the class that will contain functions needed in the pipeline  (class define in the "lib" folder) */
mf = new functions()

workflow download {
    main:
        /* Create a channel that contains the name of the SRA :
         - get the file with the SRA list (by default in bin, named 'sra.txt', see nextflow.config) 
         - split the file by line 
         - suppress the space in each line
         - keep only lines that doesn't start by # (comment line) */
         
         /* Création d'un canal */
        SRA = Channel.from(file(params.sra)).splitText().map{it -> it.trim()}.filter(it -> it !=~ /^#/)

        /* Call the process */
        download_NCBI(SRA)                                 /* FASTQ file download process from SRA */
        download_reference_genome()  /* Reference genome download process */
        download_reference_annotation()                    /* Reference annotation download process */
    emit:
    /* Send by channel the downloaded files */
        fastq_files=download_NCBI.out
        ref=download_reference_genome.out
        annot=download_reference_annotation.out
}
