install.packages("BiocManager")
install.packages("tidyverse")
install.packages("RSQLite")
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("DT")
install.packages("stats")
install.packages('conflicted')
install.packages("stringr")
install.packages("rebus")
install.packages('europepmc')
install.packages("kableExtra")
install.packages("plotly")
install.packages("devtools")
install.packages('naniar')


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
BiocManager::install('EnhancedVolcano')
BiocManager::install("recount3")

library("recount3")

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
library(conflicted)
library(stringr)
library(fgsea)
library(data.table)
library(ggplot2)
library(rebus)
library(europepmc)
library(devtools)
library(EnhancedVolcano)



BiocManager::install("clusterProfiler")











conflict_prefer("select","dplyr")


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



mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

mart2 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

sig_ensembl_stat <- result_dataframe %>%
  arrange(desc(stat)) %>%
  top_n(5,stat)
  
sig_ensembl_stat

sig_ensemblID <- rownames(sig_ensembl_stat)

sig_ensemblID <- str_replace(sig_ensemblID,
                          pattern = "\\..*", 
                          replacement = "")
sig_ensemblID



Test <- getBM(attributes = "ensembl_gene_id", "entrezgene_id",
              filter= "entrezgene_id",
              values = sig_ensemblID,
              mart = mart2)
head(Test)

mart2 = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

EnsemblIDs <- rownames(result_dataframe)
EnsemblIDs <- str_replace(EnsemblIDs,
                           pattern = "\\..*", 
                           replacement = "")











