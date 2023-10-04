#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */

process trim_single_end {
    label 'cutadapt'
    input:
        tuple val(run), file("${run}.fastq.gz")
        val(results)
    output:
        tuple val(run), file("${run}_trimmed.fq.gz")
    when:
        !params.help && !mf.checkFile("$results/FASTQ/FILTERED", run, "q.gz") && mf.checkFile("$results/FASTQ/RAW", run, "q.gz")
    script:
        """
        trim_galore -q 20 --phred33 --gzip --length 25 -o $results/FASTQ/FILTERED "${run}.fastq.gz"
        ln -s $results/FASTQ/FILTERED/${run}_trimmed.fq.gz ${run}_trimmed.fq.gz
        """
}

process index_bowtie {
    label 'bowtie'
    input:
        file(ref_genome)
        val results
        val index
    output:
        val index
    when:
        !params.help && mf.checkFile("$results/REFERENCE_FILES", "reference", "fasta")
    script:
        """
        bowtie-build $ref_genome $results/BAM/$index
        ln -s $results/BAM/$index $index
        """
}

process mapping_bowtie {
    label 'bowtie'
    input:
        tuple val(run), file("${run}_trimmed.fq.gz")
        val(results)
        val index
    output:
        file("${run}.bam")
    when:
        !params.help && !mf.checkFile("$results/BAM", run, "bam") && mf.checkFile("$results/REFERENCE_FILES", "reference", "fasta") && mf.checkFile("$results/FASTQ/FILTERED", run, "q.gz")
    script:
        """
        bowtie -p "${params.threads_mapping}" -S $results/BAM/$index <(gunzip -c "${run}_trimmed.fq.gz") | samtools sort -@ "${params.threads_sort_bam}" > "$results/BAM/${run}.bam"
        samtools index "$results/BAM/${run}.bam"
        ln -s $results/BAM/${run}.bam ${run}.bam
        """
}

/* Create the class that will contain functions needed in the pipeline  (class define in the "lib" folder) */
mf = new functions()

workflow process_fastq {
    take: fastq_files
    take: ref_genome
    take: annot_genome

    main:
        results = file(params.results)
        trim_single_end(fastq_files, results)
        index_bowtie(ref_genome, results, "index_mapping")
        mapping_bowtie(trim_single_end.out, results, "index_mapping")
    emit:
        all_bam_files=mapping_bowtie.out.collect()
}