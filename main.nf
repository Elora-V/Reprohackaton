#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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

// sub-workflow import
include { initialisation }	from		'./modules/initialisation.groovy'
include { download }	from		'./modules/download.groovy'
include { process_fastq }	from		'./modules/process_fastq.groovy'
include { counting }	from		'./modules/counting.groovy'

workflow {
    main:
        initialisation()
        download()
        process_fastq(download.out)
        counting(download.out, process_fastq.out)
}

