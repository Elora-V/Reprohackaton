#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */

process counting_reads {
    label 'subread'
    input:
        file(flag)
        file(annotation_genome)
        val(results)
    output:
        file("counts.txt")
        file("final_count_matrix.txt")
    when:
        !params.help && !mf.checkFile("$results/COUNTING", "counts", ".txt")
    script:
        """
        featureCounts -t gene -g ID -s 1 -F GTF -T "${params.threads_counting}" -a $annotation_genome -o $results/COUNTING/counts.txt $results/BAM/*.bam
        cut -f1,7,8,9,10,11,12 $results/COUNTING/counts.txt > $results/COUNTING/final_count_matrix.txt
        sed -i '1d' $results/COUNTING/final_count_matrix.txt
        ln -s $results/COUNTING/counts.txt counts.txt
        ln -s $results/COUNTING/final_count_matrix.txt final_count_matrix.txt
        """
}

/* Create the class that will contain functions needed in the pipeline  (class define in the "lib" folder) */
mf = new functions()

workflow counting {
    take: fastq_files
    take: ref_genome
    take: annot_genome
    take: all_bam_files

    main:
        /* Retrieval of the path for the result folder */
        results = file(params.results)

        counting_reads(all_bam_files, annot_genome, results)
}