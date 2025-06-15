# library packages
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(MetaboAnalystR)
library(reshape2)
library(tidyverse)
library(FSA)

# set working directory
setwd("/Users/polly_hung/Desktop/metabolomics/")
options(max.print = 1000000)
system("lsof -i -P -n | grep Rserve")
system("kill -9 5246")
HOME = "results/PD versus PR/"

# create directories
if(!dir.exists(HOME)){
  dir.create(HOME)
  dir.create(file.path(HOME, "1_normalisation"))
  dir.create(file.path(HOME, "2_univariate"))
  dir.create(file.path(HOME, "3_PCA"))
  dir.create(file.path(HOME, "4_PLS"))
  dir.create(file.path(HOME, "5_ROC"))
  setwd(HOME)
  print(paste0("directory created! We are now at: ", getwd()))
} else {
  setwd(HOME)
  print(paste0("directory exists! We are now at: ", getwd()))
}

# raw data
data <- read.delim("data.txt", check.names = FALSE)

## Preprocessing -----------
pre_mSet <- InitDataObjects("conc", "stat", FALSE)
pre_mSet <- Read.TextData(pre_mSet, "data.txt", "rowu", "disc")
pre_mSet <- SanityCheckData(pre_mSet)
pre_mSet <- ReplaceMin(pre_mSet)
pre_mSet <- PreparePrenormData(pre_mSet)
pre_mSet <- Normalization(pre_mSet,
                          rowNorm = "MedianNorm",
                          transNorm = "LogNorm",
                          scaleNorm = "MeanCenter",
                          ratio = FALSE,
                          ratioNum = 20)

# Data exploration
pre_mSet<-PlotNormSummary(pre_mSet, "1_normalisation/metabolites", format ="png", dpi=300, width=NA)
pre_mSet<-PlotSampleNormSummary(pre_mSet, "1_normalisation/samples", format = "png", dpi=300, width=NA)
data_norm <- pre_mSet$dataSet$norm
write.csv(data_norm, "data_norm.csv", row.names = TRUE, col.names = TRUE, quote = F)


## Univariate Analysis -----------
# Fold Change Analysis
pre_mSet <- FC.Anal(pre_mSet, 1.0, 0, FALSE)
pre_mSet <- PlotFC(pre_mSet, "2_univariate/fold_change", "png", 300, width=6)
fold_change <- pre_mSet$analSet$fc$sig.mat

# Perform T-test (parametric) and fold change test when comparing two groups
pre_mSet <- Ttests.Anal(pre_mSet, nonpar=F, threshp=0.15, paired=FALSE, equal.var=TRUE, "fdr", FALSE)
pre_mSet <- PlotTT(pre_mSet, imgName = "2_univariate/t_test", format = "png", dpi = 300, width=7)
t_test <- pre_mSet$analSet$tt$sig.mat

imp_var <- rownames(t_test)
data_filtered <- data %>% dplyr::select(samples, response, all_of(imp_var))
write.csv(data_filtered, "data_filtered.csv", row.names = F, col.names = T, quote = F)

## MetaboAnalystR Processing -------------
system("lsof -i -P -n | grep Rserve")
system("kill -9 5312")

mSet <- InitDataObjects("conc", "stat", FALSE)
mSet <- Read.TextData(mSet, "data_filtered.csv", "rowu", "disc")
mSet <- SanityCheckData(mSet) ## do a basic check on dataframe

mSet <- ReplaceMin(mSet) ## replace missing or zero values
mSet <- PreparePrenormData(mSet) ## prepare data for normalization
mSet <- Normalization(mSet,
                      rowNorm = "MedianNorm",
                      transNorm = "LogNorm",
                      scaleNorm = "MeanCenter",
                      ratio = FALSE,
                      ratioNum = 20)

mSet<-PlotNormSummary(mSet, "1_normalisation/filtered_norm", format ="png", dpi=300, width=NA)
mSet<-PlotSampleNormSummary(mSet, "1_normalisation/filtered_sample", format = "png", dpi=300, width=NA)


## Perform PCA analysis --------------------------------------------------------
mSet <- PCA.Anal(mSet)
w = 5
if(length(mSet$analSet$pca$cum.var) < 5){
  pc = length(mSet$analSet$pca$cum.var)
} else {
  pc = 5
}
# plots
mSet<-PlotPCAPairSummary(mSet, "3_PCA/pair_", format = "png", dpi = 300, width=w, pc)
mSet<-PlotPCAScree(mSet, "3_PCA/scree_", "png", dpi = 300, width=w, pc)
mSet<-PlotPCA2DScore(mSet, "3_PCA/score_", format = "png", dpi=300, width=w, 1, 2, 0.95, 1, 0)
mSet<-PlotPCALoading(mSet, "3_PCA/loading_", "png", 300, width=w, 1,2)
mSet<-PlotPCABiplot(mSet, "3_PCA/biplot_", format = "png", dpi = 300, width=w, 1, 2)


## Perform PLS analysis --------------------------------------------------------
# perform PLS analysis
set.seed(222)
mSet<-PLSR.Anal(mSet, reg=FALSE)
mSet<-PlotPLSPairSummary(mSet, "4_PLS/pair2_", "png", 300, width=w, pc)
mSet<-PlotPLS2DScore(mSet, "4_PLS/score2d_", "png", 300, width=w, 1,2,0.95,1,0)
mSet<-PlotPLS3DScoreImg(mSet, "4_PLS/score3d_", "png", 300, width=w, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "4_PLS/loading_", "png", 300, width=w, 1, 2)
mSet<-PLSDA.CV(mSet, choice = "Q2")
mSet$analSet$plsda$fit.info

# get permutation
mSet<-PlotPLS.Classification(mSet, "4_PLS/cv_", "png", 300, width=w)
mSet<-PlotPLS.Imp(mSet, "4_PLS/imp_", "png", 300, width=w, "vip", "Comp. 1", FALSE)
mSet<-PLSDA.Permut(mSet, 1000, "accu")
mSet<-PlotPLS.Permutation(mSet, "4_PLS/perm_", "png", 300, width=w)
pls_vip <- mSet$analSet$plsr$vip.mat %>% as.data.frame()

# get individual plot for the sample
performance_data <- mSet$analSet$plsda$fit.info %>% as.data.frame()
performance_data$Metrics <- rownames(performance_data)
performance_data <- performance_data %>%
  pivot_longer(cols = c("1 comps", "2 comps", "3 comps", "4 comps"),
               names_to = "Components", values_to = "Value")
performance_data$Components <- gsub(" comps", "", performance_data$Components)

p1 <- ggplot(performance_data, aes(x = factor(Components), y = Value, fill = Metrics)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(Value, 2)), position = position_dodge(0.9), vjust = -0.5) +
  labs(x = "Number of Components",
       y = "Performance Score") +
  ylim(0, 1.1) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB")) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "top",
        axis.title = element_text(size = 13),
        legend.title = element_text(size = 13))
ggsave("4_PLS/cv_new.png", p1, width = 4, height = 3.5, units = "in", dpi = 600)



## Perform ROC analysis --------------------------------------------------------
mSet <- SetAnalysisMode(mSet, "explore")
mSet <- PrepareROCData(mSet)
mSet <- PerformCV.explore(mSet, cls.method = "pls", rank.method = "pls", lvNum = 2)

fix(PlotROC) # w = 5, h = 5
fix(PlotProbView) # w = 6, h = 5
fix(PlotImpBiomarkers) # w = 5, h = 4
fix(PlotAccuracy) # w = 5, h = 4

mSet <- PlotROC(mSet, imgName = "5_ROC/all_models", format = "png", dpi = 600, mdl.inx = 0,
                avg.method = "threshold", show.conf = 1, show.holdout = 1)
mSet <- PlotProbView(mSet, imgName = "5_ROC/multi_roc_prob", format = "png", dpi = 600,
                     mdl.inx = 1, show = 1, showPred = 1)
mSet <- PlotAccuracy(mSet, imgName = "5_ROC/multi_roc_accuracy", format = "png", dpi = 600)
mSet <- PlotImpBiomarkers(mSet, imgName = "5_ROC/multi_roc_impvar",
                          format="png", dpi=600, mdl.inx = -1, measure="mean")
mSet<-PlotROC(mSet, imgName = "5_ROC/ROC_model",
              format = "png", dpi = 600,
              mdl.inx = 2, avg.method = "threshold",
              show.conf = 1, show.holdout = 1)

pls_vip2 <- mSet$analSet$imp.mat %>% as.data.frame()

save.image("image.RData")





















