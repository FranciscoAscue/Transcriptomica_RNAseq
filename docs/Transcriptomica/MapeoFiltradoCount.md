---
layout: page
title: Filtrado y Mapeo de RNA 
parent: Procesamiento de datos RNA-Seq
nav_order: 1
---

##  Filtrado de calidad 

## Filtrado de RNAr

Para esta seccion descargaremos las base de RNA ribosomales en el repositorio de github y instalaremos y generaremos los index para la base de datos

##### descarga de base de datos
```bash
wget -P sortmerna_db https://github.com/biocore/sortmerna/archive/2.1b.zip

unzip 2.1b.zip
./configure --prefix=$PWD
./build.sh


```

##### creacion del index de base de datos y filtrado de rRNA
```bash

sortmernaREF=sortmerna_db/rRNA_databases/silva-arc-16s-id95.fasta,sortmerna_db/index/silva-arc-16s-id95:\
sortmerna_db/rRNA_databases/silva-arc-23s-id98.fasta,sortmerna_db/index/silva-arc-23s-id98:\
sortmerna_db/rRNA_databases/silva-bac-16s-id90.fasta,sortmerna_db/index/silva-bac-16s-id95:\
sortmerna_db/rRNA_databases/silva-bac-23s-id98.fasta,sortmerna_db/index/silva-bac-23s-id98:\
sortmerna_db/rRNA_databases/silva-euk-18s-id95.fasta,sortmerna_db/index/silva-euk-18s-id95:\
sortmerna_db/rRNA_databases/silva-euk-28s-id98.fasta,sortmerna_db/index/silva-euk-28s-id98

### genear index

./sortmerna_db/indexdb_rna --ref $sortmernaREF

### loop para el filtrado de secuencias
for i in {4..8}
do
	sortmerna \
	--ref $sortmernaREF \
	--reads data/ERR292968${i}.fastq \
	--aligned results/filter/ERR292968${i}_sample_aligned.fq \
	--other results/filter/ERR292968${i}_sample_filtered.fq \
	--fastx \
	--log \
	-a 4 \
	-v

done
```
## Alinemiento de secuencias filtradas por STAR

Al igual que el alineador Bowtie2 se generara un indexado del genoma 

```bash

STAR --runMode genomeGenerate --genomeDir . --genomeSAindexNbases 7 --genomeFastaFiles GCF_000146045.2_R64_genomic.fna --sjdbGTFfile GCF_000146045.2_R64_genomic.gtf --runThreadN 40

```

#### script para el alineamiento de las secuencias

```bash

for i in {3..8}
do
		mkdir -p results/counts/ERR292968${i}
        STAR \
                --genomeDir data/reference/ \
                --readFilesIn results/filter/ERR292968${i}_sample_filtered.fq.fastq  \
                --runThreadN 4 \
                --outSAMtype BAM SortedByCoordinate \
                --limitBAMsortRAM 12000000000 \
                --quantMode GeneCounts

        mv Aligned.sortedByCoord.out.bam results/counts/ERR292968${i}
done

```


## Conteo de secuencias con featureCounts

Se calcula el conteo de las secuencias alineadas a los genes 

```bash
dir=( `find . -name "*bam"`)
featureCounts -a reference/GCF_000146045.2_R64_genomic.gtf -o trimmed/final_counts.txt -g 'gene_id' -T 10 $dir
```