---
layout: default
title: Ejecucion de funciones y scripts en R
nav_order: 4
---

## Cargar paquetes de R. 
Rbowtie2 ( alineador )    
Rsamtools ( manejar archivos tipo SAM)  
ape (leer archivos gff)  
paletas de colores (viridis)

```r
library(Rbowtie2)
library(Rsamtools)
library(ape)
library(viridisLite)
library(viridis)
```

## Explorar nuestro genoma de A. thaliana a partir de gff (archivo de anotacion)

```r
## Cargar gff file

gff_file <- read.gff("sequence.gff3", na.strings = c(".","?"), GFF3 = TRUE)
```

![](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/counts-workflow.png)

## Algunas estadisticas

```r
tab <- as.matrix(table(gff_file$type))

rnames <- as.matrix(rownames(tab))

etiquetas <- paste0(rnames[c(4,7,10,16),],"=",round(100 * tab[c(4,7,10,16),]/sum(tab[c(4,7,10,16),]), 2), "%")                                                       
par(mfrow=c(1,2), adj = TRUE)

pie(tab[c(4,7,10,16),], labels = etiquetas, col = viridis(4))

#pie(tab[c(4,7,10,16),], col = viridis(4), labels = paste0(tab[c(4,7,10,16),], "%"))

barplot(tab[c(4,7,10,16),], col = viridis(4), width = 60)

```