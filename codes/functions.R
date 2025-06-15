## general purpose packages
library(dplyr)
library(tidyverse)
library(magrittr)
library(survival)
library(scales)
library(tsibble)
library(readxl)

## plot packages
library(ggplot2)
library(RColorBrewer)
library(plotly)
library(pheatmap)

## analysis specific packages
# library(peakPantheR)
# library(e1071)
# library(sva)
# library(pls)
# library(limma)
# library(glmnet)
# library(caret)
# library(ROCR)
# library(rms)
# library(MetaboAnalystR)
# library(stats)
# library(multcomp)
library(survminer)
# library(nnet)

## set seet
set.seed(100)

## plot_distribution_and_skewness
plot_distribution <- function(data,
                              plot_name_distribution){
  metabolites <- setdiff(colnames(data), metadata)
  p <- ggplot(pivot_longer(data, cols = metabolites,
                           names_to = "metabolite", values_to = "intensity"),
              aes(x = metabolite, y = intensity)) +
    geom_boxplot(outlier.size = 2, outlier.colour = "red", outlier.alpha = 0.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p_interactive <- ggplotly(p, width = 1800, height = 384)
  return(p_interactive)

  ggsave(paste0("/plots/normalisation/", variable, "/", plot_name_distribution),
         p, width = 18, height = 4, units = "in", dpi = 600)
}

plot_skewness <- function(data,
                          plot_name_skewness){
  metabolites <- setdiff(colnames(data), metadata)
  metabolites_data <- data
  skewness <- apply(metabolites_data, 2, skewness, na.rm = TRUE) %>% abs()
  skewness_value <- data.frame(metabolite = names(skewness), skewness = skewness)

  p <- ggplot(skewness_value, aes(x = skewness)) +
    geom_density() + labs(x = "Skewness", y = "Density") + theme_bw() +
    geom_vline(xintercept = c(1, 2), colour = "red")
  p_interactive <- ggplotly(p, width = 300, height = 200)
  return(p_interactive)

  ggsave(paste0("/plots/normalisation/", variable, "/", plot_name_skewness),
         p, width = 3, height = 2, units = "in", dpi = 600)
}


customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}



function (mSetObj = NA, imgName, format = "png", dpi = 72, mdl.inx,
          show, showPred)
{
  mSetObj <- .get.mSet(mSetObj)
  anal.mode <- mSetObj$analSet$mode
  smpl.nms <- rownames(mSetObj$dataSet$norm)
  prob.vec <- rep(0.5, length(smpl.nms))
  names(prob.vec) <- smpl.nms
  if (anal.mode == "explore") {
    if (mdl.inx == -1) {
      mdl.inx <- mSetObj$analSet$multiROC$best.model.inx
    }
    probs <- MergeDuplicates(unlist(mSetObj$analSet$multiROC$pred.cv[[mdl.inx]]))
  }
  else {
    probs <- MergeDuplicates(unlist(mSetObj$analSet$ROCtest$pred.cv))
  }
  prob.vec[names(probs)] <- probs
  nms <- names(prob.vec)
  ord.inx <- order(nms)
  prob.vec <- prob.vec[ord.inx]
  cls <- mSetObj$dataSet$cls[ord.inx]
  nms <- names(prob.vec)
  pred.out <- as.factor(ifelse(prob.vec > 0.5, 1, 0))
  act.cls <- as.numeric(cls) - 1
  prob.res <- data.frame(Probability = prob.vec, Predicted = pred.out,
                         Actual = act.cls)
  if (anal.mode == "explore") {
    write.table(prob.res, file = "roc_pred_prob.csv", sep = ",",
                col.names = TRUE)
  }
  else {
    write.table(prob.res, file = "roc_pred_prob1.csv", sep = ",",
                col.names = TRUE)
  }
  conf.res <- table(pred.out, act.cls)
  mSetObj$analSet$conf.table <- xtable::xtable(conf.res, caption = "Confusion Matrix (Cross-Validation)")
  mSetObj$analSet$conf.mat <- print(mSetObj$analSet$conf.table,
                                    type = "html", print.results = F, caption.placement = "top",
                                    html.table.attributes = "border=1 width=150")
  if (anal.mode == "test") {
    if (!is.null(mSetObj$dataSet$test.data)) {
      test.pred <- ifelse(mSetObj$analSet$ROCtest$test.res >
                            0.5, 1, 0)
      test.cls <- as.numeric(mSetObj$dataSet$test.cls) -
        1
      test.df <- data.frame(Prob_HoldOut = mSetObj$analSet$ROCtest$test.res,
                            Predicted_HoldOut = test.pred, Actual_HoldOut = test.cls)
      suppressMessages(write.table(test.df, file = "roc_pred_prob1.csv",
                                   sep = ",", append = TRUE, col.names = TRUE))
      test.res <- table(test.pred, test.cls)
      mSetObj$analSet$conf.mat.test <- print(xtable::xtable(test.res,
                                                            caption = "Confusion Matrix (Hold-out)"), type = "html",
                                             print.results = F, xtable.width = 120, caption.placement = "top",
                                             html.table.attributes = "border=1 width=150")
    }
  }
  imgName = paste(imgName, "dpi", dpi, ".", format, sep = "")
  w <- 9
  h <- 8
  if (anal.mode == "explore") {
    mSetObj$imgSet$roc.prob.plot <- imgName
    mSetObj$imgSet$roc.prob.name <- mdl.inx
  }
  else {
    mSetObj$imgSet$roc.testprob.plot <- imgName
    mSetObj$imgSet$roc.testprob.name <- mdl.inx
  }
  Cairo::Cairo(file = imgName, unit = "in", dpi = dpi, width = w,
               height = h, type = format, bg = "white")
  set.seed(123)
  y <- rnorm(length(prob.vec))
  max.y <- max(abs(y))
  ylim <- max.y * c(-1.05, 1.05)
  xlim <- c(0, 1)
  op <- par(mar = c(4, 4, 3, 6))
  pchs <- ifelse(as.numeric(cls) == 1, 1, 19)
  colors <- ifelse(show == 1, "darkgrey", "black")
  ROCR::plot(prob.vec, y, pch = pchs, col = colors, xlim = xlim,
             ylim = ylim, xlab = "Predicted Class Probabilities",
             ylab = "Samples")
  abline(h = 0, lty = 2, col = "grey")
  abline(v = 0.5, lty = 2, col = "grey")
  par(xpd = T)
  legend("right", inset = c(-0.11, 0), legend = unique(as.character(cls)),
         pch = unique(pchs))
  test.y <- test.x <- 0
  if (showPred) {
    if (anal.mode == "explore") {
      test.y <- rnorm(length(mSetObj$analSet$multiROC$test.res))
      test.x <- mSetObj$analSet$multiROC$test.res
    }
    else {
      test.y <- rnorm(length(mSetObj$analSet$ROCtest$test.res))
      test.x <- mSetObj$analSet$ROCtest$test.res
    }
    pchs <- ifelse(as.numeric(mSetObj$dataSet$test.cls) ==
                     1, 1, 19)
    points(test.x, test.y, pch = pchs, cex = 1.5, col = "red")
  }
  if (show == 1) {
    act.ones <- as.numeric(cls) - 1 == 1
    pred.vec <- ifelse(prob.vec > 0.5, 1, 0)
    wrong.inx <- (pred.vec != as.numeric(cls) - 1) & pred.vec ==
      1
    if (sum(wrong.inx) > 0) {
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx],
           pos = 4)
    }
    act.zeros <- as.numeric(cls) - 1 == 0
    pred.vec <- ifelse(prob.vec < 0.5, 0, 0.5)
    wrong.inx <- pred.vec != as.numeric(cls) - 1 & pred.vec ==
      0
    if (sum(wrong.inx) > 0) {
      text(prob.vec[wrong.inx], y[wrong.inx], nms[wrong.inx],
           pos = 2)
    }
    if (showPred) {
      nms <- rownames(mSetObj$dataSet$test.data)
      act.ones <- as.numeric(mSetObj$dataSet$test.cls) -
        1 == 1
      act.zeros <- as.numeric(mSetObj$dataSet$test.cls) -
        1 == 0
      pred.vec <- ifelse(test.x > 0.5, 1, 0.5)
      wrong.inx <- (pred.vec != as.numeric(mSetObj$dataSet$test.cls) -
                      1) & act.ones
      if (sum(wrong.inx) > 0) {
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx],
             pos = 4, cex = 0.9)
      }
      pred.vec <- ifelse(test.x < 0.5, 0, 0.5)
      wrong.inx <- pred.vec != as.numeric(mSetObj$dataSet$test.cls) -
        1 & act.zeros
      if (sum(wrong.inx) > 0) {
        text(test.x[wrong.inx], test.y[wrong.inx], nms[wrong.inx],
             pos = 2, cex = 0.9)
      }
    }
  }
  par(op)
  dev.off()
  return(.set.mSet(mSetObj))
}


plot_survival <- function(surv_data, metabolite, title,
                          surv_type){

  p <- ggsurvplot(fit, data = surv_data, title = title,
                  risk.table = TRUE, risk.table.y.text.col = TRUE, risk.table.height = 0.2,
                  risk.table.y.text = FALSE,
                  pval = TRUE, pval.method = TRUE, pval.size = 3,
                  conf.int = TRUE, conf.int.style = "ribbon",
                  xlab = "Time (months)", ylab = "Survival Probability",
                  break.time.by = 2, test.for.trend = FALSE, ggtheme = theme_light(),
                  ncensor.plot = TRUE, ncensor.plot.height = 0.20,
                  surv.median.line = "hv", surv.scale = "percent",
                  legend = "top", legend.labs = c("High", "Low"),
                  fontsize = 3, size = 0.5, font.family = "Arial", font.legend = c(9),
                  censor.size = 1, censor.shape = 124,
                  palette = c("#FF407D", "#1B3C73"))
  p <- customize_labels(p, font.title = c(9, "bold"),
                        font.subtitle = c(9, "italic", "darkgrey"),
                        font.x = c(9),
                        font.y = c(9),
                        font.xtickslab = c(8))
  metabolite <- gsub("/", "_", metabolite)

  png(paste0("survival analysis/", surv_type, "/", metabolite, ".png"), width = 3, height = 5.8, units = "in", res = 600)
  print(p)
  dev.off()
}


createEmptyDF <- function(metabolites){
  df <- data.frame(metabolite = metabolites,
                      p.value = rep(1, length(metabolites)),
                      HR = rep(1, length(metabolites)),
                      upper.95 = rep(1, length(metabolites)),
                      lower.95 = rep(1, length(metabolites)))
  rownames(df) <- df$metabolite
  return(df)
}







