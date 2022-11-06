---
layout: page
title: Filtrado y Mapeo de RNA 
parent: Procesamiento de datos RNA-Seq
nav_order: 1
---

##  


```bash
wget -P sortmerna_db https://github.com/biocore/sortmerna/archive/2.1b.zip
unzip 2.1b.zip


./install-sh
./configure --prefix=$PWD


./build.sh


```


```bash

sortmernaREF=sortmerna_db/rRNA_databases/silva-arc-16s-id95.fasta,sortmerna_db/index/silva-arc-16s-id95:\
sortmerna_db/rRNA_databases/silva-arc-23s-id98.fasta,sortmerna_db/index/silva-arc-23s-id98:\
sortmerna_db/rRNA_databases/silva-bac-16s-id90.fasta,sortmerna_db/index/silva-bac-16s-id95:\
sortmerna_db/rRNA_databases/silva-bac-23s-id98.fasta,sortmerna_db/index/silva-bac-23s-id98:\
sortmerna_db/rRNA_databases/silva-euk-18s-id95.fasta,sortmerna_db/index/silva-euk-18s-id95:\
sortmerna_db/rRNA_databases/silva-euk-28s-id98.fasta,sortmerna_db/index/silva-euk-28s-id98

for i in {4..8}
do
	sortmerna \
	--ref $sortmernaREF \
	--reads ERR292968${i}.fastq \
	--aligned trimmed/ERR292968${i}_sample_aligned.fq \
	--other trimmed/ERR292968${i}_sample_filtered.fq \
	--fastx \
	--log \
	-a 40 \
	-v

done
```


```bash
dir=( `find . -name "*bam"`)
featureCounts -a reference/GCF_000146045.2_R64_genomic.gtf -o trimmed/final_counts.txt -g 'gene_id' -T 10 $dir
```