# SP_cube : Staphylococcus Persistence Profiling Pipeline

__Celine Guo<sup>1</sup>, Shun Robert<sup>1</sup>, Elora Vigo<sup>1</sup> et Fiona Hak<sup>1</sup>__
<br>
<sub>1. Université Paris-Saclay

Ce dépôt contient un pipeline Nextflow d'analyse de données RNA-seq.
Les consignes du TP, le rapport et l'article associé aux données tests sont dans le dossier "docs".

## Installation
Le lancement du pipeline nécessite que le gestionnaire de worflow Nextflow soit installé.
Il peut être installé via conda sur une machine virtuelle de l'IFM (contenant conda et ubuntu) ou en local.

    conda create -n nextflow
    conda activate nextflow
    conda install -c bioconda nextflow

Une fois Nextflow installé, ce dépôt doit être cloné :

    git clone git@gitlab.dsi.universite-paris-saclay.fr:fiona.hak/reprohackathon.git

## Usage
Le pipeline peut être lancé via ses executables directement dans le dossier d'installation :

    cd reprohackathon
    nextflow main.nf <options>

Pour plus d'information sur les options disponible :

    nextflow main.nf --help

Lancement de l'executable pour le projet :

    nextflow main.nf --reference_genome "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000253.1&rettype=fasta" --reference_annotation "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=CP000253.1"

## Dépendances
Ce pipeline est basé sur plusieurs dépendances (les mêmes qu'utilisées dans l'article, sauf cutadapt), contenu dans des environnements docker :

- sratoolkit v3.0.7
- trim_galore v0.6.10
- FastQC v0.11.9
- cutadapt 4.4 with Python 3.10.12
- bowtie v0.12.7
- samtools v1.13
- featureCounts v1.4.6
