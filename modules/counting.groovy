#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */

/* Process for counting reads using featureCounts */
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
        featureCounts -t gene -g ID -s 1 -F GTF -T "${params.threads_counting}" -a $annotation_genome -o "counts.txt" *.bam
        # Selection of column 1,7,8,9,10,11,12 (column with gene id and counts)
        cut -f1,7,8,9,10,11,12 "counts.txt" > "final_count_matrix.txt"
        # Suppression of the first line of this new file
        sed -i '1d' "final_count_matrix.txt"
        """
}

process analyse_stat {
    label 'r_stats'
    input:
        file("counts.txt")
        file("final_count_matrix.txt")
        val GSE139659_IPvsctrl
        val GeneSpecificInformation_NCTC8325
        val geneTranslation
        val script
    output:
        /* pdf figures for our data and the article ones */
        file("MA-plot.pdf")
        file("MA-plot_article.pdf")
        file("Comparaison_resultats.pdf")
    publishDir path: "${params.results}/STATS", mode: 'copy'
    when:
        /* Executed when the counts.txt file doesn't exist (and not the help parameter) */
        !params.help && !mf.checkFile("$params.results/STATS", "MA-plot", ".pdf")
    script:
        """
        Rscript $script "final_count_matrix.txt" $GSE139659_IPvsctrl $GeneSpecificInformation_NCTC8325 $geneTranslation
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
        counting_reads(all_bam_files, annot_genome) /* Executed only when all_bam_files is complete (previous step finished) */
        analyse_stat(counting_reads.out, file("./bin/GSE139659_IPvsctrl.complete.xls"), file("./bin/GeneSpecificInformation_NCTC8325.tsv"), file("./bin/geneTranslation.txt"), file("./bin/analysis_stat.r"))
}
