---
layout: page
title: Análisis Diferencial y Anotacion de Genes
parent: Procesamiento de datos RNA-Seq
nav_order: 2
---

## 2.Análisis Diferencial de Genes 

## Introducción

En esta seccion vamos a obtener los recuentos de los 2 tratamientos que hemos visto para el experimento, antes de empezar se esta dejando a continuacion otra forma de realizar el conteo conjunto a travez de R, para ello se esta cargando el paquete `Rsubread`.
``` r
library(DESeq2)
library(ggplot2)
library(Rsubread)
```

Antes de realizar el conteo de reads se realiza un lista de todos los archivos BAM alineados anteriormente que pueden se reads SE (Single End) o PE (Paired End)

``` r
bamfilesSE <- list.files(path = "bamfile_SE", 
                       pattern = "*.bam$", 
                       full.names = TRUE )

bamfilesPE <- list.files(path = "bamfile_PE", 
                       pattern = "*.bam$", 
                       full.names = TRUE )
                       
bamfilesSE
```

Una vez listados los archivos que se usaran se procede al conteo de las secuencias, para ello se necesitara tambien el genoma de referencia para lo cual deben fijar la direccion donde se encuenta:

```r
conteo <- featureCounts( files = bamfilesSE, 
                      annot.ext = "GCF_000001735.4_TAIR10.1_genomic.gff", 
                      isGTFAnnotation = TRUE, 
                      GTF.featureType = "gene", 
                      GTF.attrType = "ID", 
                      isPairedEnd = FALSE, 
                      requireBothEndsMapped = FALSE, 
                      minMQS = 20, 
                      strandSpecific = 2 )
```
           
                             ==========     _____ _    _ ____  _____  ______          _____  
                             =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
                               =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
                                 ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
                                   ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
                             ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
                            Rsubread 2.6.4
                     
                     //========================== featureCounts setting ===========================\\
                     ||                                                                            ||
                     ||             Input files : 10 BAM files                                     ||
                     ||                                                                            ||
                     ||                           SRR390302_sortedByCoord.out.bam                  ||
                     ||                           SRR390303_sortedByCoord.out.bam                  ||
                     ||                           SRR390304_sortedByCoord.out.bam                  ||
                     ||                           SRR390306_sortedByCoord.out.bam                  ||
                     ||                           SRR390307_sortedByCoord.out.bam                  ||
                     ||                           SRR390308_sortedByCoord.out.bam                  ||
                     ||                           SRR390309_sortedByCoord.out.bam                  ||
                     ||                           SRR390311_sortedByCoord.out.bam                  ||
                     ||                           SRR390312_sortedByCoord.out.bam                  ||
                     ||                           SRR390313_sortedByCoord.out.bam                  ||
                     ||                                                                            ||
                     ||              Paired-end : no                                               ||
                     ||        Count read pairs : no                                               ||
                     ||              Annotation : GCF_000001735.4_TAIR10.1_genomic.gff (GTF)       ||
                     ||      Dir for temp files : .                                                ||
                     ||                 Threads : 1                                                ||
                     ||                   Level : meta-feature level                               ||
                     ||      Multimapping reads : counted                                          ||
                     || Multi-overlapping reads : not counted                                      ||
                     ||   Min overlapping bases : 1                                                ||
                     ||                                                                            ||
                     \\============================================================================//
   
El conteo de secuencias debe ser guardado para evitar perder el analisis.
```r
conteo$counts
write.table(conteo$counts , file = "rawSE_counts.txt", sep = "\t")
```

- Desde este punto en adelante los analisis que vamos a realizar siguen el mismo procedimiento:

### Lectura de tabla de conteos
```r
data <- read.table("results/counts/conteoCompleto/final_counts.txt", header = TRUE , row.names = 1)
head(data)
```
### Edición del nombre de las secuencias 
```r
colnames(data) <- gsub(".Aligned.sortedByCoord.out.bam", "", colnames(data), fixed = T)
colnames(data) <- gsub("..", "", colnames(data), fixed = T)
row.names(data) <- gsub("gene-", "", rownames(data), fixed = T)
head(data)

```

### Creacion de metadata 

```
### Metadata 
SampleID,Group,Replicate
ERR2929683,Test,Rep1
ERR2929684,Test,Rep2
ERR2929685,Test,Rep3
ERR2929686,Reference,Rep1
ERR2929687,Reference,Rep2
ERR2929688,Reference,Rep3
```

### Cargado de metadata (muestras.txt) y edición 
```r
metadata <- read.table("muestras.txt", sep = "\t", header = TRUE)
row.names(metadata) <- metadata$SampleID 
metadata <- metadata[match(colnames(data[,-c(1:5)]), metadata$SampleID),]
#### Replace Group
metadata$Group[metadata$Group == "Test"] <- 1
metadata$Group[metadata$Group == "Reference"] <- 0
metadata$Group <- factor(metadata$Group)
```

### Creacion de datos DESeq y Normalizacion de conteos 

```r
ddsMat <- DESeqDataSetFromMatrix(countData = data[,-c(1:5)], 
                                 colData = metadata, 
                                 design = ~Group)

ddsMat <- estimateSizeFactors(ddsMat)

sizeFactors(ddsMat)
normalized_counts <- counts(ddsMat, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)


ddsMat <- DESeq(ddsMat)
```

### Ajuste de resultados para posterior análisis 
```r
# Get results from testing with FDR adjust pvalues
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)

# Generate summary of testing. 
summary(results)
```
```r
mcols(results, use.names = T)
results
```
