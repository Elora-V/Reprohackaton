#!/usr/bin/env nextflow

nextflow.enable.dsl=2 /* choice of nextflow version */

/*
------------------------------------------------------------------------------------------------------------
			 S3P Pipeline - Reprohackathon
------------------------------------------------------------------------------------------------------------

	Staphylococcus Persistence Profiling Pipeline. M2 AMI2B

	#### Homepage / Documentation
	https://gitlab.dsi.universite-paris-saclay.fr/fiona.hak/reprohackathon

	#### Authors
	Celine Guo
	Shun Robert
	Elora Vigo
	Fiona Hak
	#### Version : 1.0.0

------------------------------------------------------------------------------------------------------------
*/

/* Importation of sub-workflows (in "modules") */
include { initialisation }	from		'./modules/initialisation.groovy'
include { download }	from		'./modules/download.groovy'
include { process_fastq }	from		'./modules/process_fastq.groovy'
include { counting }	from		'./modules/counting.groovy'

workflow {
    main:

		/* Call of the sub-workflows */
        initialisation() 							/* print the help or create folders */
        download()									/* download fastq, genome and annotation */
        process_fastq(download.out)					/* trim and map the fastq */
        counting(download.out, process_fastq.out)	/* count the reads */
}

