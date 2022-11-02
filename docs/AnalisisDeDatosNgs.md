---
layout: default
title: Análisis de datos NGS
nav_order: 4
---

Introducción a datos NGS
=========================

**Que son los datos NGS?**
Son las secuencias tal cual salen de la plataforma de secuenciación (Illumina, IonTorrent, PacBio, entre otros). Es decir los **reads** (lecturas).  


Los datos crudos (sin importar la plataforma y el método de laboratorio utilizado) conllevan cierto nivel de basura, es decir:

* *reads* de baja calidad que deben ser descartados por completo o en parte (*trimmed*) antes de proceder a los análisis biológicos de los datos. 

* Secuencias sobrerepresentadas, por ejemplo dímeros de adaptadores.

La longitud y distribución de la calidad de los *reads* varían de plataforma a plataforma, pero también el tipo de error más común:

* Illumina: errores puntuales en la asignación de un nucleótido
* IonTorrent y símiles: se les complica determinar el número correcto de nucleótidos en cadenas como AAAAA.

Dado que los datos crudos son muy pesados (espacio de disco) y que buena parte de los datos crudos son basura, en general los datos crudos a este nivel no se guardan en repositorios públicos (e.g. SRA) hasta que hayan pasado por los filtros del pre-procesamiento.

Como veremos más adelante, los filtros de pre-procesamiento ayudan a identificar las buenas secuencias de la basura a partir de su calidad asociada.


## Información en los archivos FASTQ 		
### Representación de secuencias

En cómputo las secuencias de ADN son una *string* (cadena) de caracteres. 

* Secuencias genómicas

{A,C,G,T}+

* Secuencias mRNA

{A,C,G,U}+



##### Secuencias simples: FASTA/Pearson Format

Línea 1: información de la secuencia

Línea 2: la secuencia.

```

>gi|365266830|gb|JF701598.1| Pinus densiflora var. densiflora voucher Pf0855 trnS-trnG intergenic spacer, partial sequence; chloroplast
GCAGAAAAAATCAGCAGTCATACAGTGCTTGACCTAATTTGATAGCTAGAATACCTGCTGTAAAAGCAAG
AAAAAAAGCTATCAAAAATTTAAGCTCTACCATATCTTCATTCCCTCCTCAATGAGTTTGATTAAATGCG
TTACATGGATTAGTCCATTTATTTCTCTCCAATATCAAATTTTATTATCTAGATATTGAAGGGTTCTCTA
TCTATTTTATTATTATTGTAACGCTATCAGTTGCTCAAGGCCATAGGTTCCTGATCGAAACTACACCAAT
GGGTAGGAGTCCGAAGAAGACAAAATAGAAGAAAAGTGATTGATCCCGACAACATTTTATTCATACATTC
AGTCGATGGAGGGTGAAAGAAAACCAAATGGATCTAGAAGTTATTGCGCAGCTCACTGTTCTGACTCTGA
TGGTTGTATCGGGCCCTTTAGTTATTGTTTTATCAGCAATTCGCAAAGGTAATCTATAATTACAATGAGC
CATCTCCGGAGATGGCTCATTGTAATGATGAAAACGAGGTAATGATTGATATAAACTTTCAATAGAGGTT
GATTGATAACTCCTCATCTTCCTATTGGTTGGACAAAAGATCGATCCA

```


##### Fastq: formato para secuencias de secuenciación de siguiente generación

Secuencia fasta + detalles calidad de la información (la Q es de Quality).

* Línea 1: Encabezado (*Header*): comienza con @. Sigue el identificador (*identifier*). Si son datos crudos contiene info del secuenciador que identifica a esta secuencia y el read pair (/1 o /2), si son datos ya procesados en SRA contiene una descripción de la secuencia.

* Línea 2: la secuencia. 

* Línea 3: Comienza con +. Puede ser sólo el símbolo + o repetir la info del Header. 

* Línea 4: Información de la calidad de secuenciación de cada base. Cada letra o símbolo representa a una base de la secuencia codificado en formato [ASCII](http://ascii.cl/). 


La info de calidad se codifica en ASCII porque esto permite en 1 sólo caracter codificar un valor de calidad. Por lo tanto la línea 2 y la 4 tienen el mismo número de caracteres. 
 

Ordenados de menor a mayor estos son los caracteres ASCII usados para representar calidad en FASTQ:

```
!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
```
Dependiendo del tipo y versión de plataforma de secuenciación se toman diferentes caracteres (pero sin desordenarlos):


![ascii_fastqplataformas.PNG](ascii_fastqplataformas.PNG)

(Tomé la imágen de [aquí](http://en.wikibooks.org/wiki/Next_Generation_Sequencing_(NGS)/Pre-processing#Sequence_Quality)
)

Pero en todos los casos el **valor máximo de calidad es = ~40** y los valores **< 20 se consideran bajos**.


**Ojo:**
 @ y + están dentro de los caracteres ASCII utilizados para codificar la calidad, lo que hace que usar `grep` en archivos FASTQ pueda complicarse, ojo.

**Ejemplos:**

Ejemplo de datos FASTQ recién salidos de Illumina:

```
@HWI-ST999:102:D1N6AACXX:1:1101:1235:1936 1:N:0:
ATGTCTCCTGGACCCCTCTGTGCCCAAGCTCCTCATGCATCCTCCTCAGCAACTTGTCCTGTAGCTGAGGCTCACTGACTACCAGCTGCAG
+
1:DAADDDF<B<AGF=FGIEHCCD9DG=1E9?D>CF@HHG??B<GEBGHCG;;CDB8==C@@>>GII@@5?A?@B>CEDCFCC:;?CCCAC
```
Y uno más:

```
@OBIWAN:24:D1KUMACXX:3:1112:9698:62774 1:N:0:
TAATATGGCTAATGCCCTAATCTTAGTGTGCCCAACCCACTTACTAACAAATAACTAACATTAAGATCGGAAGAGCACACGTCTGAACTCAGTCACTGACC
+
CCCFFFFFHHHHHIJJJJJJJJJJJJIIHHIJJJJJJJJJJJJJJJJJJJJIJJJJJJIJJJJIJJJJJJJHHHHFDFFEDEDDDDDDDDDDDDDDDDDDC

```
¿Quieres saber cuáles son las partes del Header? [Clic aquí](https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers). Y sí, el último ejemplo es real y por lo tanto hay un Illumina HiSeq2000 que se llama Obiwan :)

 Ejemplo de datos FASTQ del SRA:

```
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
```




Los datos FASTQ típicamente están comprimidos en formato 
`gzip` (.gz) o `tar` (.tar.gz o .tgz).

**Ejercicio**:

En la Unidad 1 vimos cómo descomprimir archivos .tar.gz ¿Cómo descomprimir un archivo .gz?


Análisis de datos NGS
=====================

<p align="center" width="100%">
    <img src="https://www.researchgate.net/profile/Victoria-Dominguez-Del-Angel/publication/322946559/figure/fig2/AS:590843312361473@1517879431449/General-steps-in-a-genome-assembly-workflow-Input-and-output-data-are-indicated-for-each.png">
</p>


## CONTROL CALIDAD

Antes de saltar a filtrar tus datos con filtros de calidad que la terminal ejecute muy obediente, lo mejor es ver algunos gráficos básicos que nos dicen mucho más que una serie semi-eterna de caracteres ASCII. 

[FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) es un programa que hace una serie de análisis básicos y estándar de calidad. La mayoría de las empresas de secuenciación efectúan este análisis y te mandan los resultados junto con tus datos crudos.

Los análisis de FASTQC son útiles para identificar problemas que pudieron surgir durante el laboratorio o durante la secuenciación. 

El análisis de FASTQC consiste en los siguientes campos:


* [Basic Statistics](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/1%20Basic%20Statistics.html)
* [Per Base Sequence Quality](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html)
* [Per Sequence Quality Scores](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html)
* [Per Base Sequence Content](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html)
* [Per Sequence GC Content](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html)
* [Per Base N Content](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html)
* [Sequence Length Distribution](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html)
* [Duplicate Sequences](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html)
* [Overrepresented Sequences](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html)
* [Adapter Content](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html)
* [Kmer Content](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html)
* [Per Tile Sequence Quality](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/12%20Per%20Tile%20Sequence%20Quality.html)

**Notas importantes:**

* FASTQ automáticamente dice si nuestra muestra "pasó" o "falló" la evaluación. Sin embargo debemos tomar esto dentro del **contexto de lo que esperamos de nuestra librería**, ya que FASTQ espera una distribución al diversa y al azar de nucleótidos, lo que puede no cumplirse en algunos protocolos.



Vamos a la página de dicho programa a ver ejemplos de:

* [Buenos datos Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
* [Malos datos Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)
* [Corrida Illumina contaminada con dímeros de adaptadores (*adapter dimers*)](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html)


**¿Qué que son los dímeros de adaptadores?**

Los adaptadores se ligan al ADN de nuestras muestras en un paso de ligación, sin embargo, también pueden ligarse entre sí y luego pegarse a la *flow cell* (que lo traduzca quién sepa cómo). Resultado: son secuenciados pero no proven datos útiles, simplemente la secuencia de los adaptadores repetida  muchas veces. Adelante veremos cómo lidiar con ellos bioinformáticamente, pero se recomienda intentar deshacerse de ellos desde el laboratorio (con pequeños, pequeños imanes como [Agencourt](https://www.beckmancoulter.com/wsrportal/wsrportal.portal;jsessionid=jhp2WT8G5B4zYXhzKCPnWk8J1n3TL1JRGLsDbp130t9VRWtbFrY4!-1744838288!-683502135?_nfpb=true&_windowLabel=UCM_RENDERER&_urlType=render&wlpUCM_RENDERER_path=%2Fwsr%2Fcountry-selector%2Findex.htm&_WRpath=%252Fwsr%252Fresearch-and-discovery%252Fproducts-and-services%252Fnucleic-acid-sample-preparation%252Fagencourt-ampure-xp-pcr-purification%252Findex.htm&intBp=true) o símiles de otras marcas). 


**¿Qué tanto importa el análisis FASTQC?**

Mucho, a partir del análisis FASTQC es que decidirás si tu secuenciación fue exitosa y qué parámetros de pre-procesamiento deberás utilizar para deshacerte del ruido y quedarte con **datos limpios**. 

Escoger los parámetros adecuados de pre-procesamiento es vital ya que todas las corridas de secuenciación son diferentes. Lo más seguro es que el default del programa o lo que Perenganos *et al* 2015 reportaron en su artículo magno no sea lo mejor para procesar **tus** datos.

Además entender bien tu FASTQC puede permitirte **rescatar** datos usables incluso dentro de una mala corrida. 

**Ejercicio**.

[Baja la versión de escritorio de FastQC](http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc). Corre en java y está disponible para Linux, Mac y lo crean o no, Windows. 

Abre el archivo `Unidad5/Prac_Uni5/FastqsEjemplo/datos/human_Illumina_dataset.fastq` y examina cada uno de los módulos de FastQC ¿Qué opinas del resultado? ¿Qué harías para limpiarlo?

Estos datos ejemplo vienen de [Galaxy Data Libraries](https://usegalaxy.org/library/list#folders/F5bee13e9f312df25/datasets/99fa250d93e003f7) y son de libre uso.

### Pre-procesamiento 

El pre-procesamiento se refiere al filtrado y edición de los datos crudos para quedarnos con los **datos limpios**, que son los que se analizarán para resolver preguntas biológicas.

El input son archivos .fastq y el output son también archivos .fastq (o más posiblemente sus versiones comprimidas).



```bash
fastqc -t 2 ngsexample/data/*fastq.gz -o ngsexample/results/quality/
```

## FILTRADO DE READS
El pre-procesamiento por lo común incluye los siguientes pasos:

##### Recortar secuencias por calidad (*Sequence Quality Trimming*)
Recortar (quitar) las bases de baja calidad de la secuencia. En Illumina por lo general se trata de las últimas bases (-3' end). Cuántas bases cortar puede determinarse tras la inspección visual de análisis FASTQC o automáticamente con un parámetro de calidad. 

**Factor a considerar**: algunos programas de ensamblado requieren que las secuencias a ensamblar sean del mismo largo, por lo que si ocuparás uno de esos programas es necesario recortar todos tus datos al mismo largo (incluso si se secuenciaron en lanes distintas y una tiene buena calidad).

##### Recortar secuencias (*Trimming*)
Recortar (quitar) x bases de la secuencia porque no nos interesa conservarlas para el análisis (por ejemplo barcodes o sitios de restricción). 

##### Filtrar secuencias por calidad
Remueve del conjunto de datos todas las secuencias que estén por debajo de un mínimo de calidad (número de bases con calidad <x, promedio de calidad <x y símiles).

##### Quitar adaptadores 
Busca la secuencia de los adaptadores y los recorta de las secuencias finales. También es posible limitar las secuencias finales a sólo aquellas con un adaptador especificado (en vez de otro que pudiera ser contaminación). 

##### Filtrar artefactos
Detecta primers de PCRs, quimeras y otros artefactos y los desecha de los datos finales.

##### Separar por barcodes "demultiplexear" (*demultiplexing*)
Identifica las secuencias que contienen uno o más *barcodes* (también llamado índices), que son secuencias cortas (4-8 bp por lo general) que se incluyen en los adaptadores y que son únicos por muestra. Esto permite identificar y separar todas las secuencias pertenecientes a una muestra de otra secuenciadas al mismo tiempo. 

Requiere que le demos la lista de barcodes utilizados y en qué extremo se localizan. Muchos programas tendrán como output un archivo llamado algo como GATCATGC.fastq.gz, donde se encuentran todas las secuencias del barcode GATCATGC. El nombre de tu muestra deberás ponerlo en un paso subsecuente.

**Ojo** Tu lista barcodes-nombremuestra es de la info más valiosa de tu proyecto, no la pierdas.

##### *Paired end merging*
Si se realizó secuenciación Illumina a ambos lados (*pair end*) es posible unir las lecturas si se detecta que coinciden (aunque sea parcialmente). Esto permite corregir errores de secuenciación al tomar la base de la lectura de mayor calidad.

##### Remover otras secuencias no deseadas
Busca secuencias no deseadas, como genoma de *E. coli*, restos de PhiX o partes del genoma que no son de interés (e.g. cpDNA).  

**Ejercicio**: busca un artículo con datos y análisis parecidos a los que tendrás en tu proyecto y determina con qué programas y parámetros realizaron el pre-procesamiento. 

#### Software para realizar el pre-procesamiento

##### Línea de comando con programa especializado

La mayoría de los pasos del pre-procesamiento puede hacerse con programas dedicados a esto como [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/), pero también algunos programas especializados (e.g. en datos RAD, como pyRAD o Stacks) tienen sus propios módulos para el pre-procesamiento. 

Es importante hacer el pre-procesamiento acorde a la naturaleza y calidad de nuestros datos. 


**Ejercicio**

Utiliza un contenedor con FastXtools y otro con FASTQC para limpiar las secuencias del archivo `Prac_Uni5/FastqsEjemplo/datos/human_Illumina_dataset.fastq`, y visualizar el resultado del clipeo antes y después con FastQC. Deberás utilizar docker en la forma `docker run -v [rutaabsoutaalosdatos] fastxtools [comandosfastxtools]`, es decir **no** entrando a los contenedores.

Primero conseguir las imágenes de los contenedores:

```
# docker pull biocontainers/fastxtools
# dokcer pull biocontainers/FASTQC
```

Como sabemos, necesitaremos montar un volumen a `Prac_Uni5/FastqsEjemplo/datos/human_Illumina_dataset.fastq` utilizando la respectiva ruta absoluta dentro de nuestros equipos. Para no andar copiando la ruta tan larga todo el tiempo, es buena idea crear una variable que la contenga.

```
pmisdatos=/Users/ticatla/hubiC/Science/Teaching/Mx/BioinfInvgRepro/BioinfInvRepro2017-II/Unidad5/Prac_Uni5/FastqsEjemplo/datos/
```

Ahora sí, ¿qué comandos correr de fastxtools y de fastqc? 


**Ejercicio**

Integra todo lo que hicimos antes en un script.

**Ejercicio**

Como mencionamos antes, además de programas dedicados sólo al pre-procesamiento como FASTX, hay programas especializados que incluyen una pipeline completa para pre-procesar y procesar los datos. Esto puede ser una buena idea para cierto tipo de datos que no necesariamente cumplen con las expectativas de herramientas como FASTX. 

El ejercicio es que formen equipos de 3 personas con un tipo de datos en común y que discutan entre sí con qué programas se analizan sus tipos de datos y si estos realizan también el pre-procesamiento. 

Algunos temas para formar los equipos:

* Mapeo a genomas de referencia
* Transcriptomas/genomas de novo
* RAD/GBS para genética de poblaciones 
* RAD/GBS para filogenética
* ChipsSeq para asociaciones fenotipo-genotipo
* Metabarcoding

```bash
#trim_galore 
trim_galore --quality 20 --fastqc --length 25 --output_dir ngsexample/results/trimmed/ ngsexample/data/*.gz
```


#### GTF y GFF
"Genetic transfer format" y "Genomic feature format" [Ref](http://www.ensembl.org/info/website/upload/gff.html)

Sirven para representar características genéticas con una función anotada (gen, mRNA, rRNA, tRNA, etc).

Contenido del formato GTF:
`seqname(#chr) program feature start end score strand frame attribute(gene_id; txpt_id, etc)`

Cada intervalo toma una línea, con info en difernetes columnas.
Columnas 1-9 separadas por tab y campos *dentro* de la col 9 separados por espacio.
La columna 9 es compuesta y puede tener varios atributos: (mínimo el identificador id del gen), separados por `;`.
Las coordenadas son 1-based

Ejemplo:

```
1 transcribed_unprocessed_pseudogene  gene        11869 14409 . + . gene_id "ENSG00000223972"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; 
1 processed_transcript                transcript  11869 14409 . + . gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; gene_name "DDX11L1"; gene_sourc e "havana"; gene_biotype "transcribed_unprocessed_pseudogene"; transcript_name "DDX11L1-002"; transcript_source "havana";
```
 

#### Formato SAM/BAM
"Sequence Alignment Map", su versión binaria (comprimida) BAM.
[Ref](https://samtools.github.io/hts-specs/SAMv1.pdf)

Formato para representar un alineamiento de NGS seq.

Programa asociado: [Samtools](https://github.com/samtools/samtools)

Alineamiento: Mapear las letras (bases) de dos o más secuencias, permitiendo algunos espaciadores (indels). 

Variaciones dentro de un alineamiento:
Indels
Substituciones

El ADN se alinea de forma continua al genoma

Pero el mRNA puede tener un alineamiento dividido (spliced aligment) -> <- y forman un continuo

Contenido del formato SAM:

Header: líneas que empiezan con `@` y dan información del alineamiento, la longitud de cada cromosoma, programa con el que se hizo, etc.

Alineamiento: una línea por cada alineamiento.

Contenido de las columnas: 
```
Read id
FLAG
Chr
Start
Mapping quality
CIGAR (aligment)
MateChr
MateStart
MateDist
QuerySeq
QueryBaseQuals
AlignmentScore
Edot-distance-to-reference
Number-of-hits
Strand
Hit-index
```

Ejemplo:

Un alineamiento así:

```
Coor     12345678901234  5678901234567890123456789012345
ref      AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
+r001/1        TTAGATAAAGGATA*CTG
+r002         aaaAGATAA*GGATA
+r003       gcctaAGCTAA
+r004                     ATAGCT..............TCAGC
-r003                            ttagctTAGGC
-r001/2                                        CAGCGGCAT
```

En SAM se codifica así:

```
@HD VN:1.5 SO:coordinate
@SQ SN:ref LN:45
r001   99 ref  7 30 8M2I4M1D3M = 37  39 TTAGATAAAGGATACTG *
r002    0 ref  9 30 3S6M1P1I4M *  0   0 AAAAGATAAGGATA    *
r003    0 ref  9 30 5S6M       *  0   0 GCCTAAGCTAA       * SA:Z:ref,29,-,6H5M,17,0;
r004    0 ref 16 30 6M14N5M    *  0   0 ATAGCTTCAGC       *
r003 2064 ref 29 17 6H5M       *  0   0 TAGGC             * SA:Z:ref,9,+,5S6M,30,1;
r001  147 ref 37 30 9M         =  7 -39 CAGCGGCAT         * NM:i:1
```


## ALINEAMIENTO

## 5.6. Mapeo a genoma de referencia 		
Para esta sección leeremos los artículos:

Otra vez:
[De Wit P, Pespeni MH, Ladner JT et al. (2012) The simple fool’s guide to population genomics via RNA-Seq: an introduction to high-throughput sequencing data analysis. Molecular Ecology Resources, 12, 1058–1067.](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12003/full)

Y este tutorial asociado [http://sfg.stanford.edu/mapping.html](http://sfg.stanford.edu/mapping.html)


## 5.7. Metabarcoding y símiles 

Para esta sección leeremos los artículos:

[Creer S, Deiner K, Frey S et al. (2016) The ecologist’s field guide to sequence-based identification of biodiversity. Methods in Ecology and Evolution, 7, 1008–1018.](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12574/full)

[Coissac E, Riaz T, Puillandre N (2012) Bioinformatic challenges for DNA metabarcoding of plants and animals. Molecular ecology, 21, 1834–1847.](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2012.05550.x/full)



## 5.8. Ensamblando datos RAD, GBS y símiles 		

Para esta sección leeremos el artículo:

[Andrews KR, Good JM, Miller MR, Luikart G, Hohenlohe PA (2016) Harnessing the power of RADseq for ecological and evolutionary genomics. Nature Reviews Genetics, 17, 81–92.](http://www.nature.com/nrg/journal/v17/n2/full/nrg.2015.28.html)

**Ejercicio**

* Divídanse en equipos por tipo de datos (transcriptoma, genoma, RAD, metabarcoding, etc) y tipo de análisis (ensamblado de novo, mapeo, filogenética, genética de poblaciones).
* Discutan los artículos relevantes que dejamos leer para las secciones anteriores y cualquier otro que hayas consultado como parte de la tarea opcional anterior.
* Contesta lo siguiente:

1. ¿Cuáles son las principales variantes del método de laboratorio para generar mis datos y cuándo es más útil cada una?
2. ¿Qué limitantes y posibles fuentes de error puede presentar este método (en el laboratorio o la bioinformática)? ¿Qué puede hacerse para amortiguarlos?
3. ¿El muestreo requiere algún diseño específico? Por ejemplo, si se quiere secuenciar un genoma *de novo* ¿qué individuo sería ideal? Si trabajo con trascriptomas, ¿cómo afecta el tejido, la edad, las condiciones, etc. mi muestreo?
4. Menciona al menos dos softwares principales que se utilicen para realizar la parte medular de los análisis bioinformáticos de este tipo de análisis (e.g. si es ensamblado *de novo* con qué se ensambla, no con qué se hace el pre-procesamiento) y cuáles son los pros y contras de cada uno. 

* Contesten las preguntas anteriores en un archivo markdown llamado `Resumen-EquipoX.md` (donde X será  cada uno de los tipos de datos). Si requieres poner imágenes ponles el prefijo `Resumen-EquipoX-FigN`. Pongan los integrantes del equipo al inicio del documento.
* Realicen un pull-request para subir dicho archivo(s) al master del repositorio de la clase en la ruta `Unidad5/Prac_Uni5`.
* Unx o varixs representantes por equipo pasen a exponer brevemente (5 mins) las conclusiones utilizando como apoyo los archivos anteriores.
		
## 5.9. Importancia de elección de parámetros en análisis bioinformáticos

Para esta sección leeremos el artículo:

[Mastretta-Yanes A, Arrigo N, Alvarez N et al. (2015) Restriction site-associated DNA sequencing, genotyping error estimation and de novo assembly optimization for population genetic inference. Molecular Ecology Resources, 15, 28–41.](http://onlinelibrary.wiley.com/wol1/doi/10.1111/1755-0998.12291/full)

Y lo discutiremos con el apoyo de esta [presentación](./S2_1640_Mastretta-Yanes_5_noanimation.pdf) y los [Materiales Suplementarios 1](./SupportingInformation_1_diagram_corrected.pdf) del artículo.


## 5.10. Ambiente computacional y reproducibilidad de análisis bioinformáticos 

Para esta sección leeremos el artículo:

[Beaulieu-Jones BK, Greene CS (2017) Reproducibility of computational workflows is automated using continuous analysis. Nature Biotechnology, advance online publication.](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3780.html)


### Preparacion del index 

```bash

## descarga de genoma de referencia

efetch -db nucleotide -id NC_001224.1 -format fasta > ngsexample/data/NC_001224.1.fasta

## creamos una carpeta para el indexado de genoma

mkdir -p ngsexample/data/index

## realizamos el indexado del genoma

bowtie2-build -threads 2 ngsexample/data/NC_001224.1.fasta ngsexample/data/index/mitocp

```
### Alineamiento de secuencias

```bash
#!/bin/bash

###bowtie2

###CONSTANTS

WD="~/Curso_transcriptomica/SARS"
REF="${WD}/data/reference/NC_000.fasta"
RES="${WD}/results"
READS="${WD}/data/reads/"
r1="${READS}/sars2_1.fq"
r2="${READS}/sars2_2.fq"
OD="${RES}/map"

###EXECUTION
echo "started at `date`"

echo "mkdir -p ${OD}"
mkdir -p ${OD}

bowtie -end-to-end -I 0 -X 1000 -p 30 -x ${REF} -1 $r1 -2$r2 -S ${OD}/cavtsc.sam

samtools view -u@ 8 ${OD}/cavtsc.sam | samtools sort -@ 40 -o ${OD}/cavtsc.sorted.bam -

samtools index ${OD}/cavtsc.sorted.bam

echo "Finished at `date`"
```

