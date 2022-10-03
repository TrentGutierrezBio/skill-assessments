#This additional practice is going to pull data from recount3 to see if it gives more information compared to recount2

library("recount3")

#Finds all available human projects. 
human_projects <- available_projects()


#We gotta find an interesting project and subset the info. 
proj_info_main <- subset(human_projects,
  project == "SRP068912" & project_type == "data_sources")

proj_info_test <- subset(human_projects,
  project == "SRP169501" & project_type == "data_sources")


rse_gene_SRPtest <- create_rse(proj_info_test)



#This creates a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_SRP068912 <- create_rse(proj_info_main)

dim(rse_gene_SRP068912)

MatrixGenerics::rowRanges(rse_gene_SRP068912)

recount3_cols <- colnames(colData(rse_gene_SRP068912))

# For studies from SRA, we can further extract the SRA attributes
rse_gene_SRP186857_expanded <-
  expand_sra_attributes(rse_gene_SRPSRP068912)
colData(rse_gene_SRP068912_expanded)[, ncol(colData(rse_gene_SRP068912)):ncol(colData(rse_gene_SRP068912_expanded))]

metadata <- colData(rse_gene_SRP186857_expanded)
metadata$conditions <- rse_gene_SRP186857_expanded$sra_attribute.source_name

metadata$conditions <- grepl("control",metadata$conditions) %>%
                       ifelse("Control","Treated") 

metadata$genotype <- metadata$sra_attribute.cell_line




metadata %>% 
  as.data.frame() %>%
  dplyr::select(genotype,conditions) %>%
  kbl(caption = "Metadata for SRP186857") %>%
  kable_styling(bootstrap_options = "striped", full_width = T, html_font = "Cambria")

rse_gene_SRP186857_expanded$condition <- metadata$conditions 
rse_gene_SRP186857_expanded$genotype  <- metadata$genotype

dds <- DESeqDataSet(rse_gene_SRP186857_expanded, design = ~condition)
dds <- DESeq(dds)

result <- results(dds, contrast = c("condition","Treated","Control"))

result_dataframe <- as.data.frame(result)


colData <- colData(rse_gene_SRP186857) %>% as.data.frame()
countData <- assays(rse_gene_SRP186857)$raw_counts %>% as.data.frame()

head(rownames(colData) == colnames(countData))

colData <- colData %>%
  dplyr::filter(sra.run_acc)



