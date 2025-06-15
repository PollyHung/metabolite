library(dplyr)
library(tidyr)
library(tidyverse)
library(survival)
library(survminer)
library(openxlsx)
library(readxl)
source("~/Desktop/metabolomics/scripts/functions.R")

## Load in Markers -------------------------------------------------------------
setwd("~/Desktop/metabolomics/results/")
set1 <- setdiff(colnames(read.csv("all_data_used/3 groups comparison/data_filtered.csv", check.names = F)), c("samples", "response"))
set2 <- setdiff(colnames(read.csv("all_data_used/PD versus PR/data_filtered.csv", check.names = F)), c("samples", "response"))
set3 <- setdiff(colnames(read.csv("all_data_used/SD versus PD/data_filtered.csv", check.names = F)), c("samples", "response"))
set4 <- setdiff(colnames(read.csv("all_data_used/SD versus PR/data_filtered.csv", check.names = F)), c("samples", "response"))

metabolites <- setdiff(c(set1, set2, set3, set4), c("sample_id", "response")) %>% unique()


sum(set4 %in% colnames(metabolites$`LNEG_PeakPantheR_2024-04-22`))
sum(set4 %in% colnames(metabolites$`LNEG_Profiling_2024-04-18`))
sum(set4 %in% colnames(metabolites$`LPOS_PeakPantheR_2024-04-22`))
sum(set4 %in% colnames(metabolites$`LPOS_Profiling_2024-04-18`))
sum(set4 %in% colnames(metabolites$`BI-LISA_2024-01-09`))
sum(set4 %in% colnames(metabolites$`BI-QUANT_2023-12-22`))


## Load in normalized data
norm_data <- read.csv("results/3 groups comparison/data_norm.csv", check.names = F)
colnames(norm_data)[[1]] <- "sample_id"
norm_data <- norm_data %>% dplyr::select(sample_id, all_of(metabolites))
metadata <- read.csv("~/Desktop/metabolomics/data/metadata/meta_for_survival_analysis.csv")


## merge them together
survival <- merge(metadata, norm_data, by = "sample_id")


## Univariate analysis to identify co-variates -------------------------
os_obj <- Surv(time = survival$OS_time, event = survival$OS_event)
pfs_obj <- Surv(time = survival$PFS_time, event = survival$PFS_event)

# For numeric variables: cox-ph
# note: stage was mutated from categorical BCLC stage to numeric
numerical <- c("age", "AFP", "MTD_cm")
x = length(numerical)
univariate <- data.frame(factors = numerical,
                         p.val_os = rep(1, x), HR_os = rep(1, x),
                         p.val_pfs = rep(1, x), HR_pfs = rep(1, x))
rownames(univariate) <- univariate$factors
for(factor in numerical){
  fit <- coxph(os_obj ~ unlist(survival[factor]), data = survival)
  fit <- summary(fit)
  univariate[factor, "p.val_os"] <- fit$coefficients[1, 5]
  univariate[factor, "HR_os"] <- fit$conf.int[1, 1]

  fit2 <- coxph(pfs_obj ~ unlist(survival[factor]), data = survival)
  fit2 <- summary(fit2)
  univariate[factor, "p.val_pfs"] <- fit2$coefficients[1, 5]
  univariate[factor, "HR_pfs"] <- fit2$conf.int[1, 1]
}

# > univariate
#                           factors    p.val_os        HR_os   p.val_pfs    HR_pfs
# age                           age 0.457353108 1.015093e+00 0.188238794 0.9733662
# AFP                           AFP 0.634786683 9.999847e-01 0.047394020 1.0001186
# MTD_cm                     MTD_cm 0.469942036 1.029681e+00 0.004103663 1.1138018


# For categorical variables: log-rank test
# note: sex was taken as categorical as it does not make sense numerically
categorical <- c("etiology", "etiology_rank", "etiology_viral", "No.of.Nodules", "BCLC")
x = length(categorical)
univariate_c <- data.frame(factors = categorical,
                           p.val_os = rep(1, x), Chisq_os = rep(1, x),
                           p.val_pfs = rep(1, x), Chisq_pfs = rep(1, x))
rownames(univariate_c) <- univariate_c$factors
for(factor in categorical){
  fit <- survdiff(os_obj ~ unlist(survival[factor]), data = survival)
  univariate_c[factor, "p.val_os"] <- fit$pvalue
  univariate_c[factor, "Chisq_os"] <- fit$chisq

  fit2 <- survdiff(pfs_obj ~ unlist(survival[factor]), data = survival)
  univariate_c[factor, "p.val_pfs"] <- fit2$pvalue
  univariate_c[factor, "Chisq_pfs"] <- fit2$chisq
}

# > univariate_c
#                       factors    p.val_os  Chisq_os  p.val_pfs  Chisq_pfs
# etiology             etiology 0.004288689 24.007995 0.04312992 17.3769545
# etiology_rank   etiology_rank 0.230955973  4.298784 0.05289074  7.6892229
# etiology_viral etiology_viral 0.062860378  3.460261 0.54600642  0.3645203
# No.of.Nodules   No.of.Nodules 0.289398186  4.979589 0.01341330 12.5985701 ## categorical or continuous????
# BCLC                     BCLC 0.097143350  4.663135 0.33865008  2.1655758

wb <- createWorkbook()
addWorksheet(wb, "numerical_factors")
addWorksheet(wb, "categorical_factors")
writeData(wb, "numerical_factors", univariate)
writeData(wb, "categorical_factors", univariate_c)
saveWorkbook(wb, "results/survival analysis/univariate.xlsx")


# covariates to adjust for OS: age, BCLC, etiology, etiology_viral
# covariates to adjust for PFS: age, AFP, MTD_cm, etiology, etiology_rank


## Overall Survival by Median --------------------------------------------------
os_HR <- createEmptyDF(metabolites)
setwd("results/")
for(metabolite in metabolites){
  # mutate the level
  survival <- survival %>%
    mutate(level = ifelse(!!sym(metabolite) > median(!!sym(metabolite), na.rm = TRUE), "High", "Low"))

  # record hazard values
  cox_model <- coxph(os_obj ~ level, data = survival) ## +age+BCLC+etiology+etiology_viral
  cox_model <- summary(cox_model)
  os_HR[metabolite, "p.value"] <- cox_model$coefficients[1, 5]
  os_HR[metabolite, "HR"] <- cox_model$conf.int[1, 1]
  os_HR[metabolite, "upper.95"] <- cox_model$conf.int[1, 4]
  os_HR[metabolite, "lower.95"] <- cox_model$conf.int[1, 3]

  # plot
  # fit <- survfit(os_obj ~ level, data = survival)
  # name <- gsub("m/z", "", metabolite)
  # # plot_survival(surv_data=survival, metabolite=name, title=metabolite, surv_type="os")
}
os_HR$p.adj <- p.adjust(os_HR$p.value, method = "fdr")
os_HR_sig <- os_HR %>% dplyr::filter(p.adj < 0.25)

wb <- createWorkbook()
addWorksheet(wb, "survival")
addWorksheet(wb, "survival_sig")
writeData(wb, "survival", os_HR)
writeData(wb, "survival_sig", os_HR_sig)
saveWorkbook(wb, "survival analysis/OS_no_adjustments.xlsx")



## Progression Free Survival by Median -----------------------------------------
pfs_HR <- createEmptyDF(metabolites)

for(metabolite in metabolites){
  # mutate the level
  survival <- survival %>%
    mutate(level = ifelse(!!sym(metabolite) > median(!!sym(metabolite), na.rm = TRUE), "High", "Low"))

  # record hazard values
  cox_model <- coxph(pfs_obj ~ level, data = survival) ## +age+AFP+MTD_cm+etiology
  cox_model <- summary(cox_model)
  pfs_HR[metabolite, "p.value"] <- cox_model$coefficients[1, 5]
  pfs_HR[metabolite, "HR"] <- cox_model$conf.int[1, 1]
  pfs_HR[metabolite, "upper.95"] <- cox_model$conf.int[1, 4]
  pfs_HR[metabolite, "lower.95"] <- cox_model$conf.int[1, 3]

  # # plot
  # fit <- survfit(pfs_obj ~ level, data = survival)
  # name <- gsub("m/z", "", metabolite)
  # plot_survival(surv_data=survival, metabolite=name, title=metabolite, surv_type="pfs")
}

pfs_HR$p.adj <- p.adjust(pfs_HR$p.value, method = "fdr")
pfs_HR_sig <- pfs_HR %>% dplyr::filter(p.adj < 0.25)

wb <- createWorkbook()
addWorksheet(wb, "survival")
addWorksheet(wb, "survival_sig")
writeData(wb, "survival", pfs_HR)
writeData(wb, "survival_sig", pfs_HR_sig)
saveWorkbook(wb, "survival analysis/PFS_no_adjustments.xlsx")





adj <- readxl::read_xlsx("survival analysis/OS.xlsx", sheet = 2)$metabolite
imp <- unique(c(os_HR_sig$metabolite, adj))

for(metabolite in imp){
  # mutate the level
  survival <- survival %>%
    mutate(level = ifelse(!!sym(metabolite) > median(!!sym(metabolite), na.rm = TRUE), "High", "Low"))

  # plot
  fit <- survfit(pfs_obj ~ level, data = survival)
  name <- gsub("m/z", "", metabolite)
  plot_survival(surv_data=survival, metabolite=name, title=metabolite, surv_type="os")
}
