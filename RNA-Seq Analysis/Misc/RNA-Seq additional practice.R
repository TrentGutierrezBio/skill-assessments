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

# Once you have your RSE object, you can transform the raw coverage
# base-pair coverage counts using transform_counts().
assay(rse_gene_SRP068912, "counts") <- transform_counts(rse_gene_SRP068912)

metadata(rse_gene_SRP068912)

recount3_cols <- colnames(colData(rse_gene_SRP068912))

# How many are there?
length(recount3_cols)

# View the first few ones
head(recount3_cols)

# Group them by source
recount3_cols_groups <- table(gsub("\\..*", "", recount3_cols))

# Common sources and number of columns per source
recount3_cols_groups[recount3_cols_groups > 1]

# Unique columns
recount3_cols_groups[recount3_cols_groups == 1]

# For studies from SRA, we can further extract the SRA attributes using expand_sra_attributes() as shown below.
rse_gene_SRP068912_expanded <-
  expand_sra_attributes(rse_gene_SRP068912)
colData(rse_gene_SRP068912_expanded)[, ncol(colData(rse_gene_SRP068912)):ncol(colData(rse_gene_SRP068912_expanded))]

recount3_cols


metadata_treatment <- colData(rse_gene_SRP068912_expanded)

metadata_treatment$genotype <- metadata_treatment$sra_attribute.cell_line

metadata_treatment$condition <- grepl("control",metadata_treatment$sra_attribute.source_name) %>%
                                ifelse("control","treated")

colData(metadata_treatment)

metadata_treatment %>% 
  as.data.frame() %>%
  dplyr::select(genotype,condition) %>%
  kbl(caption = "Table 1: Metadata for SRP068912") %>%
  kable_styling(bootstrap_options = "striped", full_width = T, html_font = "Cambria")



rse_gene_SRP068912_expanded$genotype  <- metadata_treatment$genotype
rse_gene_SRP068912_expanded$condition <- metadata_treatment$condition

dds <- DESeqDataSet(rse_gene_SRP068912_expanded, design = ~condition)
dds <- DESeq(dds)
 