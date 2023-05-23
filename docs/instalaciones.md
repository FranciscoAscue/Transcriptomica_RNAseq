---
layout: default
title: Instalaciones y archivos del curso
nav_order: 2
---


# Instalaciones y archivos del curso
{: .no_toc }

Antes de empezar con el curso instale los siguientes programas y paquetes, elija el SO operativo adecuado.

## Tabla de contenidos
{: .no_toc .text-delta }

1. TOC
{:toc}


## Instalacion de R y Rstudio para Windows.


[Descargar R ](https://cran.r-project.org/bin/windows/base/R-4.3.0-win.exe){: .btn } [Descargar Rstudio ](https://download1.rstudio.org/electron/windows/RStudio-2023.03.1-446.exe){: .btn } 


## Instalacion de R y Rstudio para Linux.

### Instalar R in Linux

Procure instalar la version 4 de R (cran-40), si tiene inconvenientes instalando esta version de R añada el repositorio de cran-40 con los siguientes pasos:

```bash
## GPG key ubuntu
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

#Output
#Executing: /tmp/apt-key-gpghome.cul0ddtmN1/gpg.1.sh --keyserver #!/usr/bin/env bashkeyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
#gpg: key 51716619E084DAB9: public key "Michael Rutter <marutter@gmail.com>" imported
#gpg: Total number processed: 1
#gpg:               imported: 1

## agrere el repositorio segun la version de ubuntu que tenga en este ejemplo es para ubuntu 20 

sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

#Output
#...
#Get:7 https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/ InRelease [3622 B]                  
#Get:8 https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/ Packages [15.6 kB]
#...

## finalmente actulizar los repositorios e instalar R
sudo apt update
sudo apt install r-base
```

Si tiene la version de ubuntu 18, 20 o 22 elija entre los siguientes instaladores:

### Instalar Rstudio en Linux

[Descargar Rstudio Ubuntu 18+](https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2022.07.2-576-amd64.deb){: .btn } 


[Descargar Rstudio Ubuntu 22](https://download1.rstudio.org/desktop/jammy/amd64/rstudio-2022.07.2-576-amd64.deb){: .btn } 


## Instalaciones de paquetes para el análisis de RNASeq

Para instalar los paquetes de R usaremos el gestor de paquetes de Bioconductor.

```r
## Instalar gestor de paquetes de Bioconductor
install.packages("BiocManager")

## Activar BiocManager
library(BiocManager)
```


```r
## Instalar los siguientes paquetes:

BiocManager::install("Rsubread")
BiocManager::install("dplyr")
BiocManager::install("org.Sc.sgd.db")
BiocManager::install("pheatmap")
BiocManager::install("GO.db")
BiocManager::install("DESeq2")
BiocManager::install("ggplot2")
BiocManager::install("pathview")
BiocManager::install("KEGGgraph")
BiocManager::install("ape")
BiocManager::install("ggplot2")
BiocManager::install("pheatmap")
BiocManager::install("pathview")
BiocManager::install("clusterProfiler")
BiocManager::install("biomaRt")
BiocManager::install("genefilter")

```  

## Instalacion de programas de linux para el curso

Para el analsis de datos en Linux se instalaran los siguiente programas:

```bash
sudo apt-get install ncbi-entrez-direct bowtie2 sortmerna subread samtools rna-star trim-galore fastqc seqtk
```

`ncbi-entrez-direct` : Descarga de secuencias desde el NCBI   
`bowtie2`   : Alineador de datos NGS fastq  
`sortmerna` : Alineador de ARN ribosomales   
`subread`   : Conteo de secuencias alineadas    
`samtools`  : Manejo de secuencias alineadas SAM/BAM  
`rna-star`  : Alineador de secuencias para RNA-Seq  
`trim-galore`: Filtrado de secuencias fastq de baja calidad  
`fastqc`    : Calculo de calidad de secuencias fastq  
`setqk`     : Manejo de secuencias fastq y fasta  s

## Descarga de archivos para el curso

Para la creacion del directorio de trabajo y descarga de los archivos debe ejecutar el siguiente comando en la terminal Linux:

```bash
wget -O - https://raw.githubusercontent.com/FranciscoAscue/Transcriptomica_RNAseq/master/downloadSeq.bash | bash

curl -s https://raw.githubusercontent.com/FranciscoAscue//Transcriptomica_RNAseq/master/downloadSeq.bash | bash
```

Este script debe crear el siguiente arbol de directorios para la organizacion de los archivos:

```

saccharomyces_rnaSeq
├── data
│   ├── ERR2929683_sample0.1.fastq.gz
│   ├── ERR2929684_sample0.1.fastq.gz
│   ├── ERR2929685_sample0.1.fastq.gz
│   ├── ERR2929686_sample0.1.fastq.gz
│   ├── ERR2929687_sample0.1.fastq.gz
│   ├── ERR2929688_sample0.1.fastq.gz
│   ├── GCF_000146045.2_R64_genomic.fna.gz
│   └── GCF_000146045.2_R64_genomic.gtf.gz
├── results
│   ├── counts
│   │   ├── final_counts.txt
│   │   └── final_counts.txt.summary
│   ├── filter
│   └── map
└── scripts

```
