library(dplyr)
library(magrittr)

setwd("~/Desktop/metabolomics/")

load("data/loadme.RData")

# ## Add the sample column name
# metabolites <- lapply(file_names, function(x){
#   df <- metabolites[[x]]
#   df$samples <- rownames(df)
#   return(df)
# })
# names(metabolites) <- file_names


# merge all samples
# combined <- merge(metabolites$`LNEG_PeakPantheR_2024-04-22`, metabolites$`LNEG_Profiling_2024-04-18`, by = "samples")
combined <- merge(metabolites$`LNEG_PeakPantheR_2024-04-22`, metabolites$`LPOS_PeakPantheR_2024-04-22`, by = "samples")
# combined <- merge(combined, metabolites$`LPOS_PeakPantheR_2024-04-22`, by = "samples", all.x = TRUE)
# combined <- merge(combined, metabolites$`LPOS_Profiling_2024-04-18`, by = "samples", all.x = TRUE)
combined <- merge(combined, metabolites$`BI-LISA_2024-01-09`, by = "samples", all.x = TRUE)
combined <- merge(combined, metabolites$`BI-QUANT_2023-12-22`, by = "samples", all.x = TRUE)

# ## Add patient
# meta <- read.csv("data/metadata/lneg.csv") %>% dplyr::select(sample_id, subject_id)
# combined <- merge(meta, combined, by.x = "sample_id", by.y = "samples")

## metadata
meta <- readxl::read_xlsx("~/Desktop/metabolomics/data/metadata/2025_Jan_02_immunotherapy_lipidomics.xlsx") %>%
  dplyr::select(Sample.ID, Response_Rank_1)  %>%
  dplyr::filter(Response_Rank_1 %in% c("PR", "SD"))

## Immunotherapy
samples <- meta$Sample.ID

## combined
combined_subset <- combined %>% dplyr::filter(samples %in% samples)
combined_subset <- merge(meta, combined_subset, by.x = "Sample.ID", by.y = "samples")
combined_subset <- combined_subset %>% dplyr::rename(samples = Sample.ID, response = Response_Rank_1)
write.table(combined_subset, "~/Desktop/metabolomics/results/annotated/SD versus PR/data.txt",
            sep = "\t", row.names = F, col.names = T)


# ## transpose the dataframe
# combined <- as.data.frame(t(combined))
# colnames(combined) <- combined[1, ]
# combined <- combined[2:nrow(combined), ]
#
# ## Save file
# write.table(combined, "data/combined.txt", sep = "\t", quote = F, col.names = T, row.names = T)
#
# ## metadata
# meta <- readxl::read_xlsx("data/metadata/2025_Jan_02_immunotherapy_lipidomics.xlsx")


filtered_id <- read.csv("data_filtered.csv", row.names = 1, check.names = F) %>% colnames() %>% setdiff("response")
plot_df <- read.csv("data_norm.csv", row.names = 1, check.names = F)
plot_df$`K-EDTA` <- NULL

plot_df$sample <- rownames(plot_df)
plot_df <- merge(meta, plot_df, by.x = "Sample.ID", by.y = "sample")
plot_df$sample_response <- paste0(plot_df$Sample.ID, "_", plot_df$Response_Rank_1)
plot_df <- plot_df %>% tibble::column_to_rownames(var = "sample_response")
plot_df$Sample.ID <- NULL
plot_df$Response_Rank_1 <- NULL
plot_df <- plot_df %>% dplyr::select(all_of(filtered_id))

correlation <- cor(t(plot_df))

png("1_normalisation/after_normalisation_heatmap_correlation_filtered.png", width = 10, height =10, units = "in", res = 600)
pheatmap::pheatmap(correlation)
dev.off()




