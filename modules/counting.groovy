#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process counting_reads {
    label 'subread'
    input:
        file("$results/BAM/*.bam")
        file(annotation_genome)
        val(results)
    output:
        file("counts.txt")
    when:
        !params.help && !mf.checkFile("$results/COUNTING", "counts", ".txt")
    script:
        """
        featureCounts -t gene -g ID -s 1 -F GTF -T "${params.threads_counting}" -a $annotation_genome -o $results/COUNTING/counts.txt "$results/BAM/*.bam"
        ln -s $results/COUNTING/counts.txt counts.txt
        """
}

mf = new functions()

workflow counting {
    take: fastq_files
    take: ref_genome
    take: annot_genome
    take: all_bam_files

    main:
        results = file(params.results)
        counting_reads(all_bam_files, annot_genome, results)
}