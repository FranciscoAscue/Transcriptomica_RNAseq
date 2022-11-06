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



./sortmerna_db/indexdb_rna --ref $sortmernaREF


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


```bash

STAR --runMode genomeGenerate --genomeDir . --genomeSAindexNbases 7 --genomeFastaFiles GCF_000146045.2_R64_genomic.fna --sjdbGTFfile GCF_000146045.2_R64_genomic.gtf --runThreadN 40

```

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



```bash
dir=( `find . -name "*bam"`)
featureCounts -a reference/GCF_000146045.2_R64_genomic.gtf -o trimmed/final_counts.txt -g 'gene_id' -T 10 $dir
```