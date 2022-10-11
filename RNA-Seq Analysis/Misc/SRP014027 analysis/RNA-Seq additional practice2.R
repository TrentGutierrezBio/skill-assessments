url <- download_study("SRP014027") 
load(file.path("SRP014027", 'rse_gene.Rdata'))

colData(rse_gene)

metadata <- colData(rse_gene)

metadata$conditions <- sapply(metadata$characteristics, function(x){x[3]}) %>%
                       grepl(pattern = "control") %>%
                       ifelse("control","treatment")

metadata$genotype <- grepl("Gata6Hopx",rse_gene$title) %>%
                     ifelse("Gata6Hopx","WT")

metadata %>% 
  as.data.frame() %>%
  dplyr::select(genotype,conditions) %>%
  kbl(caption = "Table 1: Metadata for SRP014027") %>%
  kable_styling(bootstrap_options = "striped", full_width = T, html_font = "Cambria")

rse_gene$genotype  <- metadata$genotype

rse_gene$condition <- metadata$conditions

dds <- DESeqDataSet(rse_gene, design = ~condition)
dds <- DESeq(dds)

result <- results(dds, contrast = c("condition","treatment","control"))

result_dataframe <- as.data.frame(result)

fullensemblID <- rownames(result_dataframe) %>%
  str_replace(pattern = "\\..*",replacement = "")

sampledist <- vst(dds) %>%
  assay() %>%
  t() %>%
  dist()


sampledistmatrix <- as.matrix(sampledist)

vsd_dds <- vst(dds)

rownames(sampledistmatrix) <- paste(vsd_dds$condition, vsd_dds$genotype, sep="-")
colnames(sampledistmatrix) <- NULL

colsBlue <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

distheatmap <- pheatmap(sampledistmatrix,
                        clustering_distance_rows=sampledist,
                        clustering_distance_cols=sampledist,
                        col=colsBlue)

plotPCA(vsd_dds, intgroup = c("condition","genotype"))


sigresults_df <- result %>%
  as.data.frame() %>%
  dplyr::mutate(ID = gsub(pattern = "\\..*", replacement = "", x= rownames(result))) %>%
  dplyr::relocate(ID,baseMean) %>%
  arrange(desc(stat)) %>%
  subset(padj < 0.05) 


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

sigresult_gene_id <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                           filters = "ensembl_gene_id",
                           values = sigresults_df$ID,
                           mart = ensembl) %>%
  dplyr::rename(ID = ensembl_gene_id, 
                Gene_name = hgnc_symbol)

sigresults_df <- left_join(sigresults_df,sigresult_gene_id, by= "ID") %>%
  dplyr::relocate(Gene_name,baseMean) %>%
  dplyr::select(ID,Gene_name,log2FoldChange,stat,pvalue,padj)

sigresults_df_filtered <- sigresults_df %>%
  drop_na(Gene_name) %>%
  dplyr::filter(Gene_name != "")

datatable(sigresults_df_filtered, class = 'cell-border stripe' , 
          caption = "Table 2: Patient vs Control Condition Significant Genes")


EnhancedVolcano(sigresults_df_filtered, 
                lab = sigresults_df_filtered$Gene_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                title ="volcano Slot of Signifcant DEGs",
                xlim = c(-10,10),
                ylim = c(0,20),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1.5,
                labSize = 2 )

topover <- dplyr::slice_max(sigresults_df_filtered,n = 10, order_by = stat)
topunder <- dplyr::slice_min(sigresults_df_filtered,n = 10, order_by = stat)
exprdata <- rbind(topover,topunder)

normalized_dds_counts <- counts(dds, normalized=TRUE)
row.names(normalized_dds_counts) <- str_replace(fullensemblID,
                                                pattern = "\\..*", 
                                                replacement = "")

sig_norm_dds_counts <- normalized_dds_counts[exprdata$ID,]
row.names(sig_norm_dds_counts) <- (exprdata$Gene_name)



heatmeta <- as.data.frame(metadata) %>%
  dplyr::select(conditions,genotype)


colsRdBu <- brewer.pal(11, "RdBu")

palette <- colorRampPalette(colsRdBu)

heatmap <- pheatmap(sig_norm_dds_counts,
                    color = palette(20),
                    cluster_rows = TRUE,
                    show_rownames = TRUE,
                    annotation = dplyr::select(heatmeta,conditions),
                    scale = "row") 

gene_list <- result_dataframe$stat %>%
  setNames(object = ,fullensemblID) %>%
  na.omit() %>%
  sort(decreasing = TRUE)

gse <- gseGO(geneList = gene_list , ont = "ALL",
             keyType = "ENSEMBL",
             pvalueCutoff = 0.05, verbose = TRUE,
             OrgDb= "org.Hs.eg.db", pAdjustMethod = "holm")

overGSEAplot <- gse %>%
  arrange(desc(rank)) %>%
  gseaplot2(geneSetID = 1:5, title = "Top 5 Overexpressed Pathways")
overGSEAplot


underGSEAplot <- gse %>%
  arrange(rank) %>%
  gseaplot2(geneSetID = 1:5, title = "Top 5 Underexpressed Pathways")
underGSEAplot











