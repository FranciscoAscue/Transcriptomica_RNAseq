---
layout: default
title: Instalaciones y archivos del curso
nav_order: 2
---


## Instalaciones del curso

Antes de empezar con el curso instale los siguientes programas y paquetes, elija el SO operativo adecuado.

## Tabla de contenidos
{: .no_toc .text-delta }


## Instalacion de R y Rstudio para Windows.


[Descargar R ](https://cran.r-project.org/bin/windows/base/R-4.2.1-win.exe){: .btn } [Descargar Rstudio ](https://download1.rstudio.org/desktop/windows/RStudio-2022.07.2-576.exe){: .btn } 


## Instalacion de R y Rstudio para Linux.

### Install R in Linux

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
BiocManager::install("org.At.tair.db")
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
