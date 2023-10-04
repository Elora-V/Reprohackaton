#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */


/* Process that download the fastq file */
process download_NCBI {
    label 'sratoolkit'
    input:
        /* run = list of the sra file */
        val(run)
        /* name of the result folder */
        val(results)
    output:
        /* name of sra file and the downloaded sra file */
        tuple val(run), file("${run}.fastq.gz")
    when:
        /* Executed when the fastq file aren't already downloaded (and not the help parameter) */
        !params.help && !mf.checkFile("$results/FASTQ/RAW", run, "q.gz")
    script:
        """

        prefetch $run     # downloads all necessary files for the computer for fasterq-dump
        fasterq-dump $run --skip-technical --threads "${params.threads_download}" # extraction of fastq from SRA-accessions 
        #(with number of threads specified, default indicated in nextflow.config)

        # compressing fastq (pigz : parallel implementation of gzip)
        pigz ${run}.fastq --best
    
        # suppressing an old potential fastq of the sra
        rm -rf $results/FASTQ/RAW/${run}*.fastq.gz
        # moving the compressed fastq into the wanted folder
        mv *.fastq.gz $results/FASTQ/RAW
        # create a symbolic link of this file at the base (here it's the current position : the 'base' of the repository)
        # (needed if the complete path is not describe in the output of the process)
        ln -s $results/FASTQ/RAW/$run* .
        # suppressing the non-compressed file
        rm -rf $run
        """
}



/* Process that download the genome */
process download_reference_genome {
    label 'linux'
    input:
        /* name of the result folder */
        val(results)
    output:
        /* genome file */
        file("reference.fasta")
    when:
        /* Executed when there is the genome parameter, and the genome isn't already downloaded (and not the help parameter) */
        !params.help && params.reference_genome != "none" && !mf.checkFile("$results/REFERENCE_FILES", "reference", "fasta")
    script:
        """
        # download the genome with wget (need the genome parameter from the user)
        wget -q -O $results/reference.fasta "${params.reference_genome}"
        # create a link at the base
        ln -s $results/reference.fasta reference.fasta 
        """
}



/* Process that download the annotation */
process download_reference_annotation {
    label 'linux'
    input:
        /* name of the result folder */
        val(results)
    output:
        /* annotation file */
        file("reference.gff")
    when:
        /* Executed when there is the annotation parameter, and the annotation isn't already downloaded (and not the help parameter) */
        !params.help && params.reference_annotation != "none" && !mf.checkFile("$results/REFERENCE_FILES", "reference", "gff")
    script:
        """
        # download the annotation with wget (need the annotation parameter from the user)
        wget -q -O $results/reference.gff "${params.reference_annotation}"
        # create a link at the base
        ln -s $results/reference.gff reference.gff
        """
}


/* Create the class that will contain functions needed in the pipeline  (class define in the "lib" folder) */
mf = new functions()



workflow download {

    main:
        /* Retrieval of the path for the result folder */
        results = file(params.results)

        /* Create a channel that contains the name of the SRA :
         - get the file with the SRA list (by default in bin, named 'sra.txt', see nextflow.config) 
         - split the file by line 
         - suppress the space in each line
         - keep only lines that doesn't start by # (comment line) */
        SRA = Channel.from(file(params.sra)).splitText().map{it -> it.trim()}.filter(it -> it !=~ /^#/)

        /* Call the process */
        download_NCBI(SRA, results)                                 /* download the fastq in results */
        download_reference_genome(results+"/REFERENCE_FILES")       /* download the genome in "/REFERENCE_FILES" in results */
        download_reference_annotation(results+"/REFERENCE_FILES")   /* download the annotation in "/REFERENCE_FILES" in results*/


    emit:
    /* Send by channel the downloaded files */
        fastq_files=download_NCBI.out
        ref=download_reference_genome.out
        annot=download_reference_annotation.out
}
