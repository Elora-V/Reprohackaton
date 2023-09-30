#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process trim_single_end {
    label 'cutadapt'
    input:
        val(run)
        file("${run}.fastq.gz")
        val(results)
    output:
        val(run)
        file("${run}_trimmed.fq.gz")
    when:
        !mf.checkFile("$results/FASTQ/FILTERED", run, "q.gz")
    script:
        """
        trim_galore -q 20 --phred33 --gzip --length 25 -o $results/FASTQ/FILTERED "${run}.fastq.gz"
        """
}

mf = new functions()

workflow clean_fastq {
    take: run
    take: fastq_files

    main:
        results = file(params.results)
        trim_single_end(run, fastq_files, results)
}