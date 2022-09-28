#This additional practice is going to pull data from recount3 to see if it gives more information compared to recount2

library("recount3")

#Finds all available human projects. 
human_projects <- available_projects()


#We gotta find an interesting project and subset the info. 
proj_info <- subset(
  human_projects,
  project == "SRP186857" & project_type == "data_sources"
)

#This creates a RangedSummarizedExperiment (RSE) object at the gene level
rse_gene_SRP186857 <- create_rse(proj_info)

dim(rse_gene_SRP186857)

MatrixGenerics::rowRanges(rse_gene_SRP186857)

recount3_cols <- colnames(colData(rse_gene_SRP186857))

# For studies from SRA, we can further extract the SRA attributes
rse_gene_SRP186857_expanded <-
  expand_sra_attributes(rse_gene_SRP186857)
colData(rse_gene_SRP186857_expanded)[, ncol(colData(rse_gene_SRP186857)):ncol(colData(rse_gene_SRP186857_expanded))]

















