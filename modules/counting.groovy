#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */


/* Process that count the number of reads for each gene */
process counting_reads {
    label 'subread'
    input:
        /* input needed to start the process but not used */
        file(flag)
        /* annotation */
        file(annotation_genome)
    output:
        /* raw file of featureCounts */
        file("counts.txt")
        /* modified file of counts */
        file("final_count_matrix.txt")
    publishDir path: "${params.results}/COUNTING", mode: 'copy'
    when:
        /* Executed when the counts.txt file doesn't exist (and not the help parameter) */
        !params.help && !mf.checkFile("$params.results/COUNTING", "counts", ".txt")
    script:
        """
        # Use of featureCounts to count reads :
        # - t : on gene
        # - g : get id of gene
        # - s : 1 for stranded reads
        # - F : format of annotation file
        # - T : number of threads (define by user or default define in nextflow.config)
        # - a : annotation file
        # - o : output file
        featureCounts -t gene -g ID -s 1 -F GTF -T "${params.threads_counting}" -a $annotation_genome -o counts.txt *.bam
        # Selection of column 1,7,8,9,10,11,12 (column with gene id and counts)
        cut -f1,7,8,9,10,11,12 counts.txt > final_count_matrix.txt
        # Suppression of the first line of this new file
        sed -i '1d' $results/COUNTING/final_count_matrix.txt
        """
}



/* Create the class that will contain functions needed in the pipeline  (class define in the "lib" folder) */
mf = new functions()


workflow counting {

    /* Retrieval of the emited files from the previous sub-workflow */
    take: fastq_files  
    take: ref_genome
    take: annot_genome
    take: all_bam_files

    main:
        /* Call the counting process */
        counting_reads(all_bam_files, annot_genome) /* executed only when all_bam_files is complete (previous step finished) */
}