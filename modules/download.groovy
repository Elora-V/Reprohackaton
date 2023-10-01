#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process download_NCBI {
    label 'sratoolkit'
    input:
        val(run)
        val(results)
    output:
        tuple val(run), file("${run}.fastq.gz")
    when:
        !mf.checkFile("$results/FASTQ/RAW", run, "q.gz")
    script:
    """
        prefetch $run 
        fasterq-dump $run --skip-technical

        pigz ${run}.fastq --best
    
        rm -rf $results/FASTQ/RAW/${run}*.fastq.gz
        mv *.fastq.gz $results/FASTQ/RAW
        ln -s $results/FASTQ/RAW/$run* .
        rm -rf $run
        """
}

process download_reference_genome {
    input:
        val(results)
    output:
        file("reference.fasta")
    when:
        !params.help && params.reference_genome != "none" && !mf.checkFile("$results/REFERENCE_FILES", "reference", "fasta")
    script:
        """
        wget -q -O $results/reference.fasta "${params.reference_genome}"
        ln -s $results/reference.fasta reference.fasta
        """
}

process download_reference_annotation {
    input:
        val(results)
    output:
        file("reference.gff")
    when:
        !params.help && params.reference_annotation != "none" && !mf.checkFile("$results/REFERENCE_FILES", "reference", "gff")
    script:
        """
        wget -q -O $results/reference.gff "${params.reference_annotation}"
        ln -s $results/reference.gff reference.gff
        """
}

mf = new functions()

workflow download {
    results = file(params.results)
    main:
        //cherche dans le fichier conf sra les noms des run, supprime les espaces et ignore les commentaires en #
        SRA = Channel.from(file(params.sra)).splitText().map{it -> it.trim()}.filter(it -> it !=~ /^#/)
        download_NCBI(SRA, results)
        download_reference_genome(results+"/REFERENCE_FILES")
        download_reference_annotation(results+"/REFERENCE_FILES")
    emit:
        fastq_files=download_NCBI.out
        ref=download_reference_genome.out
}
