#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process download_NCBI {
    label 'sratoolkit'
    input:
        val(run)
        val(results)
    output:
        val(run)
        file("${run}.fastq.gz")
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

mf = new functions()

workflow download {
    results = file(params.results)
    main:
        //cherche dans le fichier conf sra les noms des run, supprime les espaces et ignore les commentaires en #
        SRA = Channel.from(file(params.sra)).splitText().map{it -> it.trim()}.filter(it -> it !=~ /^#/)
        download_NCBI(SRA, results)
    emit:
        run=fastq_files=download_NCBI.out[0]
        fastq_files=download_NCBI.out[1]
}
