#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */

process trim_single_end {
    label 'cutadapt'
    input:
        /* list of the sra files and one fastq */
        tuple val(run), file("${run}.fastq.gz")
    output:
        /* list of the sra files and one trimmed file */
        tuple val(run), file("${run}_trimmed.fq.gz")
    publishDir path: "${params.results}/FASTQ/FILTERED", mode: 'copy'
    when:
        /* Executed when the trimmed file doesn't exist and the fastq exist (and not the help parameter) */
        !params.help && !mf.checkFile("$params.results/FASTQ/FILTERED", run, "q.gz") && mf.checkFile("$params.results/FASTQ/RAW", run, "q.gz")
    script:
        """
        # Conservation of reads with a minimal length of 25 with good quality
        trim_galore -q 20 --phred33 --gzip --length 25 -o . ${run}.fastq.gz
        """
}

process index_bowtie {
    label 'bowtie'
    input:
        /* annotation */
        file(ref_genome)
        /* name of the index output */
        val index
    output:
        /* used to indicate the end of this process */
        val index
        tuple file("${index}.1.ebwt"), file("${index}.2.ebwt"), file("${index}.3.ebwt"), file("${index}.4.ebwt"), file("${index}.rev.1.ebwt"), file("${index}.rev.2.ebwt")
    publishDir path: "${params.results}/BAM", mode: 'copy'
    when:
        /* Executed when the genome is downloaded (and not the help parameter) */
        !params.help && mf.checkFile("$params.results/REFERENCE_FILES", "reference", "fasta")
    script:
        """
        # building the index for the mapping with bowtie
        bowtie-build $ref_genome $index
        """
}

process mapping_bowtie {
    label 'bowtie'
    input:
        /* list of the sra files and one trimmed file */
        tuple val(run), file("${run}_trimmed.fq.gz")
        /* output of previous process (used to link the two process) */
        val index
        tuple file("${index}.1.ebwt"), file("${index}.2.ebwt"), file("${index}.3.ebwt"), file("${index}.4.ebwt"), file("${index}.rev.1.ebwt"), file("${index}.rev.2.ebwt")
    output:
        /* one mapped file */
        file("${run}.bam")
    publishDir path: "${params.results}/BAM", mode: 'copy'
    when:
        /* Executed when the mapped file doesn't exist, the genome is downloaded and the trimmed fastq exist (and not the help parameter) */
        !params.help && !mf.checkFile("$results/BAM", run, "bam")
    script:
        """
        # Mapping with bowtie (using the index) and sorting with samtools
        bowtie -p "${params.threads_mapping}" -S $index <(gunzip -c "${run}_trimmed.fq.gz") | samtools sort -@ "${params.threads_sort_bam}" > "${run}.bam"
        # Indexing the bam (mapped file) with samtools
        samtools index "${run}.bam"
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
        /* Call of process */
        trim_single_end(fastq_files)                           /* Trimming of reads */
        index_bowtie(ref_genome, "index_mapping")              /* Index for mapping of reads */
        mapping_bowtie(trim_single_end.out, index_bowtie.out)  /* Mapping of reads */

    emit:
        /* Emition of a list of all the items emitted by the process mapping_bowtie */
        all_bam_files=mapping_bowtie.out.collect()
}