install.packages("BiocManager")
install.packages("tidyverse")
install.packages("RSQLite")
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("DT")
install.packages("gsea")
install.packages("stats")
install.packages('conflicted')



BiocManager::install(version = "3.15")


BiocManager::install("recount")
BiocManager::install("SummarizedExperiment")
BiocManager::install("MatrixGenerics")
BiocManager::install("sparseMatrixStats")
BiocManager::install("GenomeInfoDbData")
BiocManager::install("DESeq2")
BiocManager::install("fgsea")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("org.Hs.eg.db", character.only = TRUE)
BiocManager::install("ReactomePA")
BiocManager::install("reactome.db")



library(tidyverse)
library(BiocManager)
library(recount)
library(RColorBrewer)
library(pheatmap)
library(DT)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(plotly)
library(styler)
library(kableExtra)
library(fgsea)
library(stats)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(org.Hs.eg.db)
library(ReactomePA)
library(reactome.db)
library(devtools)
library(conflicted)

library(fgsea)
library(data.table)
library(ggplot2)

data(examplePathways)
data(exampleRanks)

fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)

pathways <- reactomePathways(names(exampleRanks))
fgseaRes <- fgsea(pathways, exampleRanks, maxSize=500)
head(fgseaRes)




url <- download_study("SRP022028") 
url

load(file.path("SRP022028", 'rse_gene.Rdata'))
rse <- scale_counts(rse_gene)

colData(rse)

conditions <- c("patient","patient","patient","patient",
                "control","control","control","control")

rse$condition <- conditions

dds <- DESeqDataSet(rse, design = ~condition)

dds <-DESeq(dds)

result <- results(dds)
result

plotMA(result)

plotCounts(dds, gene=which.min(result$padj), intgroup="condition")


# This is to convert the Ensembl IDs into Entrez ID eventually to be read by fgsea. 
BiocManager::install("biomaRt")
library(biomaRt)

# This code should connect us to the most recent human dataset. 
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

listEnsembl()

datasets <- listDatasets(ensembl)
head(datasets)


EnsemblIDs <- rownames(result_dataframe)
EnsemblIDs <- gsub("\\.*","",EnsemblIDs)


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))


Test <- getBM(attributes = "ensembl_gene_id", "entrezgene_id",
              filters = "ensembl_gene_id",
              values = EnsemblIDs,
              mart = mart)













