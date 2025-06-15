library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggridges)


combinedData_path <- list.files(pattern = "_combinedData.csv", path = "data", full.names = TRUE, recursive = T)

## Function to clean up the data frame
clean_up <- function(df) {
  df_cleaned <- df[, colSums(is.na(df) | df == "") < nrow(df)]
  df_cleaned <- df_cleaned[rowSums(is.na(df_cleaned) | df_cleaned == "") <= 15, ]
  return(df_cleaned)
}

## Read in all dataframes
metabolites <- lapply(combinedData_path, function(x){
  df <- read.csv(x, skip = 1, check.names = F, row.names = 1) %>% clean_up
  df <- df[!(df[[3]] %in% c("Study Pool", "Study Pool Sample") |
               df[[4]] %in% c("Study Pool", "Study Pool Sample")), ]
  return(df)
})
file_names <- unlist(lapply(strsplit(combinedData_path, split = " "), function(x){x[[3]]}))
file_names <- unlist(lapply(strsplit(file_names, split = "/"), function(x){x[[1]]}))
names(metabolites) <- file_names

## Clean up
for(i in file_names[1:4]){
  rownames(metabolites[[i]]) <- NULL
  rownames(metabolites[[i]]) <- metabolites[[i]][[4]]
  metabolites[[i]] <- metabolites[[i]][18:ncol(metabolites[[i]])]
}

## Clean up 2
for(i in file_names[5:6]){
  rownames(metabolites[[i]]) <- NULL
  rownames(metabolites[[i]]) <- metabolites[[i]][[3]]
  metabolites[[i]] <- metabolites[[i]][17:ncol(metabolites[[i]])]
}


## plot
df_long <- bind_rows(lapply(file_names, function(name) {
  df <- apply(metabolites[[name]], 2, as.numeric)
  df_long <- tidyr::pivot_longer(as.data.frame(df), cols = everything(),
                          names_to = "Metabolite", values_to = "Value")
  df_long$Group <- name  # Add a group identifier
  return(df_long)
}))

quantiles <- quantile(df_long$Value)

df_long <- df_long %>% mutate(Quantile = case_when(Value <= quantiles[2] ~ "0-25%",
                                                   Value <= quantiles[3] ~ "25-50%",
                                                   Value <= quantiles[4] ~ "50-75%",
                                                   TRUE ~ "75-100%"))

p <- ggplot(df_long, aes(x = Value, y = Group, fill = Group)) +
  geom_density_ridges(alpha = 0.6) +
  facet_wrap(~ Quantile, scales = "free_x", nrow = 1) +
  theme_ridges() +
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Metabolite Intensity") +
  ylab("Datasets") +
  labs(title = "Ridge Plot of Metabolite Distributions by Quantile")
ggsave("plots/MetaboliteDistributionQuantileRidgePlot.png", p, width = 9, height = 4, dpi = 600, units = "in")




