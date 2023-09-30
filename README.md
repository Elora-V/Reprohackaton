# S3P : Staphylococcus Persistence Profiling Pipeline

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

*EXÉCUTABLE A CRÉER*

## Dépendances
Ce pipeline est basé sur plusieurs dépendances (les mêmes qu'utilisées dans l'article), contenu dans des environnements docker :

- 


