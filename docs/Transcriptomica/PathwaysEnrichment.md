---
layout: page
title: Vias metabolicas de genes expresados diferencialmente
parent: Procesamiento de datos RNA-Seq
nav_order: 4
---

## Enrriquecimiendo de vias metabolicas

 
``` r
library(DESeq2)
library(org.Sc.sgd.db)
library(GO.db)
library(KEGGgraph)
library(clusterProfiler)
library(pathview)
```


### Anotacion de genes segun bases de datos
```r
columns(org.Sc.sgd.db)
keytypes(org.Sc.sgd.db)

# Add gene full name
results$description <-  AnnotationDbi::select(org.Sc.sgd.db, keys = row.names(results),
                                              columns = 'DESCRIPTION', 
                                              keytypes = 'ALIAS')

results$description

results$genename <-  AnnotationDbi::select(org.Sc.sgd.db, keys = row.names(results),
                                              columns = 'GENENAME', 
                                              keytypes = 'ALIAS')

# Add gene symbol
results$symbol <- row.names(results)

# Add ENTREZ ID
results$entrez <- AnnotationDbi::select(org.Sc.sgd.db, keys = row.names(results),
                                                                columns = 'ENTREZID', 
                                                                keytypes = 'ALIAS')
results$ORF <- AnnotationDbi::select(org.Sc.sgd.db, keys = row.names(results),'ORF','ALIAS', multiVals = "list")

head(results$ORF)
columns(org.Sc.sgd.db)
keytypes(org.Sc.sgd.db)
results$GO <- AnnotationDbi::select(org.Sc.sgd.db, keys = row.names(results),"GOALL","ALIAS", multiVals = "list")

### Add gene Ontology code

results$GO <- mapIds(org.Sc.sgd.db, keys = row.names(results),'GOALL','ORF')


```
### Enrriquecimiento de GO (gene ontology)

```r
columns(GO.db)
keytypes(GO.db)
enrich <- AnnotationDbi::select(GO.db, keys = results_sig$GO, columns = c('TERM', 'ONTOLOGY'), keytypes = 'GO')
enrich$symbol <- results_sig$symbol
enrich

enrich2 <- AnnotationDbi::select(GO.db, keys = results_sig$GO, columns = 'DEFINITION', keytypes = 'GO')
enrich2$symbol <- results_sig$symbol
enrich2

results_sig_entrez <- subset(results_sig, is.na(GO) == FALSE)

data <- enrich %>% group_by(ONTOLOGY) %>% summarise(N = n())
barplot(data$N)
# Create a matrix of gene log2 fold changes
# View the format of the gene matrix
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
results_sig_entrez

gene_matrix <- results_sig_entrez$log2FoldChange
rownames(results_sig_entrez)
names(gene_matrix) <- results_sig_entrez$entrez$ENTREZID

```
### Enrriquecimiento de KEGG

```r
enrichKEeeGG <- enrichKEGG(gene = rownames(results_sig_entrez),
           organism = 'sce',
           pvalueCutoff = 0.05, 
           qvalueCutoff = 0.10)


barplot(enrichKEeeGG, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)
```

### Plot de Pricipales vias metabolicas afectadas

```r
pathview(gene.data = gene_matrix,
         pathway.id = "00010", 
         species = "sce")
```

