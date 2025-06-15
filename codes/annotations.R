library(MetaboAnalystR)
library(magrittr)
library(dplyr)
library(ggplot2)
library(biodb)

load("data/loadme.RData")

set1 <- setdiff(colnames(read.csv("results/all_data_used/3 groups comparison/data_filtered.csv", check.names = F)), c("samples", "response"))
set2 <- setdiff(colnames(read.csv("results/all_data_used/PD versus PR/data_filtered.csv", check.names = F)), c("samples", "response"))
set3 <- setdiff(colnames(read.csv("results/all_data_used/SD versus PD/data_filtered.csv", check.names = F)), c("samples", "response"))
set4 <- setdiff(colnames(read.csv("results/all_data_used/SD versus PR/data_filtered.csv", check.names = F)), c("samples", "response"))

metabolites <- setdiff(c(set1, set2, set3, set4), c("sample_id", "response")) %>% unique()
untargeted <- metabolites[grepl("m/z", metabolites)]
untargeted_df <- combined[, untargeted]

mybiodb <- biodb::newInst()
compUrl <- system.file("extdata", "chebi_extract.tsv", package='biodb')
compdb <- mybiodb$getFactory()$createConn('comp.csv.file', url=compUrl)

mz_values_pos <- read.csv("data/raw_data/Immunotherapy S LPOS_Profiling_2024-04-18/Immunotherapy S LPOS_Profiling_featureMetadata.csv",
                          row.names = 1) %>% dplyr::filter(Feature.Name %in% untargeted)
intensity <- read.csv("data/raw_data/Immunotherapy S LPOS_Profiling_2024-04-18/Immunotherapy S LPOS_Profiling_intensityData.csv", header = F) %>% dplyr::filter(Feature.Name %in% untargeted)
write.table(mz_values_pos[, 2:3], "~/Desktop/mz_values_pos.txt", sep = "\t", quote = F, col.names = F, row.names = F)

annotations <- compdb$annotateMzValues(mz_values_pos, mz.tol = 1e-3, ms.mode= "pos")


hmdb_results <- SearchHMDB(mz_values, mode = "pos", tol = 0.01)  # Adjust polarity/tolerance
