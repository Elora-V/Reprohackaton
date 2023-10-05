#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */

process trim_single_end {
    label 'cutadapt'
    input:
        /* list of the sra files and one fastq */
        tuple val(run), file("${run}.fastq.gz")
        /* name of the result folder */
        val(results)
    output:
        /* list of the sra files and one trimmed file */
        tuple val(run), file("${run}_trimmed.fq.gz")
    when:
        /* Executed when the trimmed file doesn't exist and the fastq exist (and not the help parameter) */
        !params.help && !mf.checkFile("$results/FASTQ/FILTERED", run, "q.gz") && mf.checkFile("$results/FASTQ/RAW", run, "q.gz")
    script:
        """
        # Conservation of reads with a minimal length of 25 with good quality
        trim_galore -q 20 --phred33 --gzip --length 25 -o $results/FASTQ/FILTERED "${run}.fastq.gz"
        # creation of link at the base
        ln -s $results/FASTQ/FILTERED/${run}_trimmed.fq.gz ${run}_trimmed.fq.gz
        """
}

process index_bowtie {
    label 'bowtie'
    input:
        /* annotation */
        file(ref_genome)
        /* name of the result folder */
        val results
        /* name of the index output */
        val index
    output:
        /* used to indicate the end of this process */
        val index
    when:
        /* Executed when the genome is downloaded (and not the help parameter) */
        !params.help && mf.checkFile("$results/REFERENCE_FILES", "reference", "fasta")
    script:
        """
        # building the index for the mapping with bowtie
        bowtie-build $ref_genome $results/BAM/$index
        # creation of link at the base
        ln -s $results/BAM/$index $index
        """
}

process mapping_bowtie {
    label 'bowtie'
    input:
        /* list of the sra files and one trimmed file */
        tuple val(run), file("${run}_trimmed.fq.gz")
        /* name of the result folder */
        val(results)
        /* output of previous process (used to link the two process) */
        val index
    output:
        /* one mapped file */
        file("${run}.bam")
    when:
        /* Executed when the mapped file doesn't exist, the genome is downloaded and the trimmed fastq exist (and not the help parameter) */
        !params.help && !mf.checkFile("$results/BAM", run, "bam") && mf.checkFile("$results/REFERENCE_FILES", "reference", "fasta") && mf.checkFile("$results/FASTQ/FILTERED", run, "q.gz")
    script:
        """
        # Mapping with bowtie (using the index) and sorting with samtools
        bowtie -p "${params.threads_mapping}" -S $results/BAM/$index <(gunzip -c "${run}_trimmed.fq.gz") | samtools sort -@ "${params.threads_sort_bam}" > "$results/BAM/${run}.bam"
        # Indexing the bam (mapped file) with samtools
        samtools index "$results/BAM/${run}.bam"
        # creation of link at the base
        ln -s $results/BAM/${run}.bam ${run}.bam
        """
}

/* Create the class that will contain functions needed in the pipeline  (class define in the "lib" folder) */
mf = new functions()

workflow process_fastq {

    /* Retrieval of the emited files from the previous sub-workflow */
    take: fastq_files
    take: ref_genome
    take: annot_genome

    main:
        /* Retrieval of the path for the result folder */
        results = file(params.results)

        /* Call of process */
        trim_single_end(fastq_files, results)                           /* Trimming of reads */
        index_bowtie(ref_genome, results, "index_mapping")              /* Index for mapping of reads */
        mapping_bowtie(trim_single_end.out, results, index_bowtie.out)   /* Mapping of reads */

    emit:
        /* Emition of a list of all the items emitted by the process mapping_bowtie */
        all_bam_files=mapping_bowtie.out.collect()
}