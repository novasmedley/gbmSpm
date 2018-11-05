#---------------------------------------- Description---------------------------------------- 
# Assumes binary classification
# ROC and PR curves

#---------------------------------------- binary classification---------------------------------------- 

#' depends on ROCR
#' 
#' @description Get the ROC and PR curve from model probabilies.
#' 
#' @param prob An object of class prediction, predicted values and ground truth, input to \code{\link[ROCR]{performance}}
#' @param round logical, whether to round roc and pr axis values to 3 decimals.
#' 
#' @return list of dataframes \code{roc} and \code{pr} x and y data points, and the area under the curves.
#' @export
getROC_PR <- function(prob, round=TRUE) {
  r <- ROCR::performance(prob, "tpr", "fpr")
  p <- ROCR::performance(prob, "prec","rec")
  
  roc <- data.frame(x=r@x.values[[1]], 
                    y=r@y.values[[1]])
  pr <- data.frame(x=p@x.values[[1]], 
                   y=p@y.values[[1]])
  
  aucROC <- unlist(slot(ROCR::performance(prob, 'auc'),"y.values"))
  aucPR <- MLmetrics::PRAUC(y_pred = prob@predictions[[1]], 
                             y_true = prob@labels[[1]])
  if (round) {
    aucPR <- round(aucPR, digits=3)
    aucROC <- round(aucROC, digits=3)
  }
  results <- list(roc=roc, 
                  aucROC=aucROC, 
                  pr=pr, 
                  aucPR=aucPR)
  return(results)
}

#' Select the roc threshold with the max F1 score and calculate classification metrics
#' 
#' @param trainPred An object of class prediction, predicted values and ground truth, input to \code{\link[ROCR]{performance}}
#' @return dataframe of various metrics when using threshold to do classification, 
#' including sensitivity, specificity, TP, FP, TN, FN, precision, recall, f1score
#' 
#' @export
selectROCThreshold <- function(trainPred) {
  pr <- ROCR::performance(trainPred, "prec","rec")
  acc <- ROCR::performance(trainPred,'acc')
  ss <- ROCR::performance(trainPred,'sens','spec')
  f1Scores <- 2*( (pr@y.values[[1]]*pr@x.values[[1]]) / (pr@y.values[[1]]+pr@x.values[[1]]) )
  
  # specificy a single threshold (instead of multiple, as in the ROC and PR curves)
  # and extract performance stats
  cuttoffInd <- which.max(f1Scores)
  cutData <- list()
  cutData$cutoff <- trainPred@cutoffs[[1]][cuttoffInd] 
  cutData$fp <- trainPred@fp[[1]][cuttoffInd]
  cutData$tp <- trainPred@tp[[1]][cuttoffInd]
  cutData$tn <- trainPred@tn[[1]][cuttoffInd]
  cutData$fn <- trainPred@fn[[1]][cuttoffInd]
  cutData$spec <- ss@x.values[[1]][cuttoffInd]
  cutData$sens <- ss@y.values[[1]][cuttoffInd]
  cutData$prec <- pr@y.values[[1]][cuttoffInd]
  cutData$recall <- pr@x.values[[1]][cuttoffInd]
  cutData$f1score <- f1Scores[cuttoffInd]
  cutData$accuracy <- acc@y.values[[1]][cuttoffInd]
  cutData$cutoffInd <- cuttoffInd
  cutData <- unlist(cutData, recursive=F)
  return(cutData)
}

#---------------------------------------- ROC and PR  --------------------------------------- 

#' Plot multiple performance curves in one graph
#' 
#' @param ptitle plot title
#' @param fn plot filename
#' @param rocs list of ROCs, each containing dataframe with columns:
#' \itemize{
#' \item \code{x} false positive rate
#' \item \code{y} true positive rate
#' \item \code{L1} names of different ROCs to differeniate x and y (important if there are multiple ROCs)
#' }
#' @param prs list of precision-recall curves, each containing dataframe with columns:
#' \itemize{
#' \item \code{x} recall
#' \item \code{y} precision
#' \item \code{L1} names of different curvers to differeniate x and y (important if there are multiple)
#' }
#' @param aucs area under the ROC values matching the order of \code{names}
#' @param names names of the different ROCs in \code{L1}, to be included in legend name
#' @param outDir directory to put \code{fn}, the plot
#' @param perfTheme ggplot2 theme() object to include plot formats
#' @param scaleX scale_x_continuous, format of x scale
#' @param scaleY scale_y_continuous, format of y scale
#' @param scaleColor scale_color_manual, specify colors
#' 
#' @export
plotROC <- function(ptitle, fn, rocs, aucs, names, outDir, perfTheme, scaleX, scaleY, scaleColor) {
  names(rocs) <- paste0(names, ' (area=', aucs, ')')
  rocs <- melt(rocs, id.vars=c('x','y'))
  rocP <- ggplot(rocs, aes(x=x, y=y, color=factor(L1) )) + 
    ggtitle(ptitle) +
    geom_line(size=.5) + 
    geom_abline(intercept = 0, slope = 1, colour = "gray") +
    ylab("True positive rate") + 
    xlab("False positive rate") + 
    theme_bw() + perfTheme + scaleX + scaleY + scaleColor
    
  print(rocP)
  cat(file.path(outDir, paste0(fn,'.png')))
  dev.print(png, file.path(outDir, paste0(fn,'.png')), res=300,height=3, width=4.1, units="in")
  rocP
}

#' @rdname plotROC
#' @export
plotPR <- function(ptitle, fn, prs, aucs, names, outDir, perfTheme, scaleX, scaleY, scaleColor) {
  names(prs) <- paste0(names, ' (area=', aucs, ')')
  prs <- melt(prs, id.vars=c('x','y'))
  prP <- ggplot(prs, aes(x=x, y=y, color=factor(L1))) + 
    ggtitle(ptitle) +
    geom_line(size=.5) + 
    ylab("Precision") + 
    xlab("Recall") +
    theme_bw() + perfTheme + scaleX + scaleY + scaleColor
  print(prP)
  out <- file.path(outDir, paste0(fn,'.png'))
  cat(out)
  dev.print(png, out, res=300,height=3, width=4.1, units="in")
  prP
}



#' Multiple ROC plots in one and also faceted
#' 
#' @description Plots multiple ROC or PR curves in one plot (i.e., different models), 
#' and faceted by another variable (i.e., different classification task)
#' 
#' @inheritParams plotROC
#' 
#' @param fn plot filename
#' @param rocsList named list of ROCs, each containing dataframe with columns:
#' \itemize{
#' \item \code{x} false positive rate
#' \item \code{y} true positive rate
#' \item \code{L1} names of different ROCs to differeniate x and y (important if there are multiple ROCs, i.e., training vs testing)
#' }
#' @param prsList  named list of  precision-recall curves, each containing dataframe with columns:
#' \itemize{
#' \item \code{x} recall
#' \item \code{y} precision
#' \item \code{L1} names of different ROCs to differeniate x and y (important if there are multiple ROCs, i.e., training vs testing)
#' }
#' @param aucsList named list of area under the ROC values
#' @param listOrder order of plotting for legend labeling of multipe ROC or PR curves
#' 
#' @export
plotMultiROC <- function(fn, rocsList, aucsList, listOrder, legendTitle,printIt, outDir, perfTheme, scaleX, scaleY, perfGuide) {
  rocs <- melt(rocsList, id.vars=c('x','y'))
  names <- names(rocsList)
  rocs$L1<- as.factor(rocs$L1)
  rocs$L1 <- factor(rocs$L1, levels=levels(rocs$L1)[listOrder])

  aucs <- melt(aucsList)
  aucs$label <- paste0('area=', aucs$value)
  aucs$x <- 0.55 # align at this position on x axis
  aucs$y <- 0.1*seq(0,length(unique(aucs$L2))-1) # text align on y axis
  aucs$L1<- as.factor(aucs$L1)
  aucs$L1 <- factor(aucs$L1, levels=levels(aucs$L1)[listOrder])
  
  p <- ggplot(rocs, aes(x=x, y=y, color=factor(L2) )) + 
    geom_abline(intercept = 0, slope = 1, colour = "gray") +
    geom_line(size=.35) + 
    geom_text(data=aucs, aes(x,y, label=label), show.legend = F, hjust=0, vjust=.05,size=2.5) + # aligned text to the left starting
    ylab("True positive rate") + 
    xlab("False positive rate") +
    theme_linedraw() + perfTheme + scaleX + scaleY + perfGuide +
    theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank()) +
    labs(color=legendTitle)
  p <- p + facet_grid(. ~ L1)
  
  out <- file.path(outDir, paste0(fn,'.png') )
  if (printIt) {
    print(p)
    dev.print(png, out, res=800,height=3, width=8.75, units="in")
  }  
  p
}
#' @rdname plotMultiROC
#' 
#' @export
plotMultiPR <- function(fn, prsList, aucsList, listOrder, legendTitle, printIt, outDir, perfTheme, scaleX, scaleY, perfGuide) {
  prs <- melt(prsList, id.vars=c('x','y'))
  names <- names(prsList)
  prs$L1<- as.factor(prs$L1)
  prs$L1 <- factor(prs$L1, levels=levels(prs$L1)[listOrder])
  
  aucs <- melt(aucsList)
  aucs$label <- paste0('area=', aucs$value)
  aucs$x <- .55
  aucs$y <- rev(.95-(0.1*seq(0,length(unique(aucs$L2))-1)))
  aucs$L1<- as.factor(aucs$L1)
  aucs$L1 <- factor(aucs$L1, levels=levels(aucs$L1)[listOrder])
  
  # aucs$y[aucs$L2=='train'] <- 0
  p <- ggplot(prs, aes(x=x, y=y, color=factor(L2) )) + 
    geom_line(size=.35) + 
    geom_text(data=aucs, aes(x,y, label=label), show.legend = F, hjust=0, size=2.5) +
    ylab("Precision") + 
    xlab("Recall") +
    theme_linedraw() + perfTheme + scaleX + scaleY + perfGuide +
    theme(panel.grid.minor=element_blank(),panel.grid.major = element_blank()) +
    labs(color=legendTitle)
  p <- p + facet_grid(. ~ L1)
  out <- file.path(outDir, paste0(fn,'.png') )
  if (printIt) {
    print(p)
    dev.print(png, out, res=800,height=3, width=8.75, units="in")
  }
  p
}