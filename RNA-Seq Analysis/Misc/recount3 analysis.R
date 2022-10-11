library(tidyverse)
library(BiocManager)
library(recount)
library(RColorBrewer)
library(pheatmap)
library(DT)
library(DESeq2)
library(fgsea)
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
library(devtools)
library(conflicted)
library(stringr)
library(fgsea)
library(data.table)
library(biomaRt)
library(EnhancedVolcano)
library(gage)
library(gageData)

rse <- recount3::create_rse_manual(
       project = "SRP127016",
       project_home = "data_sources/sra",
       organism = "human",
       annotation = "gencode_v26",
       type = "gene")





























