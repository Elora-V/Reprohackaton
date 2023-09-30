#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process downloadNCBI {
    label 'sratoolkit'
    input:
        val(run)
        val(results)
    output:
        tuple val(run), file("*.fastq.gz")
    when:
        !mf.checkFile("$results/FASTQ/RAW", run, "q.gz")
    script:
        """
        prefetch $run -N 1000
    
        fasterq-dump $run --skip-technical
        pigz *.fastq --best
    
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
        SRA = Channel.from(file(params.sra)).splitText().map{it -> it.trim()}.filter(it -> it !=~ /^#/)
        downloadNCBI(SRA, results)
        downloadNCBI.out.branch{paired: it[1].size() == 2}.set{fastq}
    emit:
        all_paired_fastq = fastq.paired
}
