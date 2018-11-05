# ----------- Description -------

# plot ROC curves and compare for different approaches

library(gbmSpm)

library(gridExtra) # graph
library(grid) # graph
library(scales) # graph
library(wesanderson) # graph color
# library(viridis) # graph color
# library(RColorBrewer) # graph color

# install.packages('devtools') # if not installed, run below two lines
# devtools::install_github("baptiste/egg") 
library(egg) # arrange graphs
library(cowplot) # arrange graphs
# constants ------------------------------------------------------------

dataDir <- '~/Data/Agile'
spmLogitDir <- file.path(dataDir, 'analysis_spm') # temporal pattern approach
volLogitDir <- file.path(dataDir, 'analysis_vol', 'logits') # tumor vol approaches
volThresDir <- file.path(dataDir, 'analysis_vol', 'threshold') # tumor vol approaches
  
outDir <- file.path(dataDir, 'compare_approaches')# saved performance plots and coefficients
figDir <- file.path(dataDir, 'figures')# manuscript figures
finalDir <- file.path(dataDir, 'final_model') # saved performance plots and coefficients

months <- c('2-month', '6-month', '9-month', '12-month')
canUse <- c('continuous','percentChange', 'rateChange', 'baselineSize')
canUseNames <- c('mm^3', 'percent change', 'rate change', 'baseline')

# performance plot constants
perfTheme <- theme(legend.text=element_text(size=7, color="grey30"),
                   legend.title=element_text(size=7), 
                   legend.text.align =1,
                   legend.position='right', 
                   legend.key.size=unit(.8,'line'),
                   legend.margin=margin(l = -0.15, unit='cm'),
                   legend.background = element_rect(fill = "transparent", colour = "transparent"),
                   axis.text=element_text(size=7, color="#666666"),
                   axis.title.y = element_text(margin=margin(0,5,0,0),size=7, angle=90),
                   axis.title.x = element_text(margin=margin(5,0,-5,0),size=7),
                   panel.spacing = unit(.35, "lines"),
                   strip.text.x = element_text(size=9,margin = margin(.05,0,.05,0, "cm")))
scaleX <- scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2))
scaleY <- scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2))
perfGuide <- guides(colour = guide_legend(override.aes = list(size=4)))

# 
# ---- get SPM logit CV results -----

# cv logit results
topCVS <- readRDS(file.path(spmLogitDir,'spm_logit_top_cv_results.rds'))
bestLogitSpm <- lapply(months, function(m) {
  z <- topCVS[which(topCVS$month==m),]
  b <- z[which.max(z$rocAUC),]
  b <- b[,colnames(b) %in% c('rocAUC','rocAUCSD','tType','month','model')]
  rownames(b) <- b$model
  b
})
bestLogitSpm <- do.call('rbind', bestLogitSpm)
colnames(bestLogitSpm) <- c('rocAUC','rocAUCSD','month','name','model')

# ---- get volume logit performance -----

volLogitResults <- readRDS(file.path(volLogitDir,'tumorVol_logits_CV_results.rds'))

# # single plots
# lapply(months, function(i) plotCVROC(results, i))
# lapply(months, function(i) plotCVPR(results, i))

# group plots
volLogitResults$line <- as.factor(volLogitResults$line)
levels(volLogitResults$line) <- c('continuous','discretized', 'continuous with\ncovariates', 'discretized with\ncovariates' )
volLogitResults$line <- factor(volLogitResults$line, levels=levels(volLogitResults$line)[c(1,3,2,4)])
volLogitResults$month <- as.factor(volLogitResults$month)
volLogitResults$month <- factor(volLogitResults$month, levels=levels(volLogitResults$month)[c(2,3,4,1)])

crossValTheme <- theme(axis.text=element_text(size=7),
                       axis.text.x = element_text(margin=margin(10,0,-5,0), size=6),
                       axis.title.y = element_text(angle=0, vjust=.5,size=7, margin=margin(0,5,0,0)),
                       axis.title.x = element_blank(),
                       legend.text=element_text(size=7, color='grey30'), 
                       legend.title=element_text(size=7),
                       legend.position="right",
                       legend.margin=margin(l = -0.15, unit='cm'),
                       legend.key.height=unit(1,'line'),
                       legend.key.width = unit(.85,'line'),
                       panel.grid.minor.x = element_blank(),
                       panel.spacing = unit(.2, "lines"),
                       panel.grid.major.x = element_blank(),
                       strip.text.x = element_text(size=9,margin = margin(.05,0,.05,0, "cm")))
crossValGuide <- guides(fill = guide_legend(override.aes = list(size=1)))
crossValLabels <- labs(fill='Variable(s)',x='Tumor Volume Type')
crossValFillColor <-  scale_fill_manual(values = c(wes_palette("Darjeeling2")[c(4,2)], wes_palette('Darjeeling1')[c(3,4)]) )

volLogitResults$type <- factor(volLogitResults$type, levels=levels(volLogitResults$type)[c(3,2,4,5,1)]) # must match names below
crossValXLabels <- scale_x_discrete(labels=c('rate\nchange', parse(text='mm^3'), 'baseline', 'percent\nchange', 'response\ncriteria'))

dev.new()

gROC <- ggplot(volLogitResults, aes(x=type, y=rocAUC, fill=line) ) +
  geom_bar(position=position_dodge(), stat='identity') +
  scale_y_continuous('Average Area\n Under the\n ROC Curve', limits=c(0,1), breaks=seq(0,1,0.2)) +
  theme_linedraw() + 
  crossValTheme + crossValLabels + crossValFillColor  + crossValXLabels + crossValGuide +
  geom_errorbar(aes(ymin=rocAUC-rocAUCSD, ymax=rocAUC+rocAUCSD), width=0.1, size=.2, position=position_dodge(.9))
gROC <- gROC + facet_grid(month~.) + theme(axis.title.x=element_blank())
print(gROC)
# dev.print(png, file.path(subDir, paste0('tumorVol_logit_all_CV_aucroc','.png')), res=600,height=4, width=6, units="in")

gPR <- ggplot(volLogitResults, aes(x=type, y=prAUC, fill=line) ) +
  geom_bar(position=position_dodge(), stat='identity') +
  scale_y_continuous('Average Area\n Under the\n PR Curve', limits=c(0,1), breaks=seq(0,1,0.2)) +
  theme_linedraw() + 
  crossValTheme + crossValLabels + crossValFillColor + crossValXLabels + crossValGuide +
  geom_errorbar(aes(ymin=prAUC-prAUCSD, ymax=prAUC+prAUCSD), width=0.1, size=.2, position=position_dodge(.9))
gPR <- gPR + facet_grid(month~.) #+ guides(fill=F) + theme(strip.background = element_blank(),
                                  #                        strip.text.x = element_blank())
print(gPR)
# dev.print(png, file.path(subDir, paste0('tumorVol_logit_all_CV_aucpr','.png')), res=600,height=4, width=6, units="in")

tvLogitPlots <- ggarrange(gROC, gPR, heights=c(1.33,.91) , ncol=1)

setEPS(); postscript(file.path(figDir, paste0('s6_tumorVol_logit_all_CV_aucPR_aucROC.eps')), height=9.31, width = 7.5)
pp <- plot_grid(gROC, gPR, labels = c("A", "B"), ncol=1)
print(pp)
dev.off()


#
#
# ---- get volume threshold performance -----
#
# prep
volResults <- readRDS(file.path(volThresDir,'tumorVol_threshold_results.rds'))

vol.ROCs <- lapply(volResults, '[[',2)
vol.auc.roc <- lapply(volResults, '[[',3)
vol.PR <- lapply(volResults, '[[', 4)
volauc.pr <- lapply(volResults, '[[',5)
vol.cutoffs <- lapply(volResults, '[', c(2,6))
vol.cutoffs <- lapply(vol.cutoffs, function(i) {
  r <- cbind(i$roc,i$cutoffs)
  colnames(r) <- c('x','y','cut')
  r
})
rocsList <- list(vol.ROCs[1:4],vol.ROCs[5:8], vol.ROCs[9:12], vol.ROCs[13:16])
names(rocsList) <- months
rocsList <- lapply(rocsList, function(i) {
  names(i) <- canUseNames
  i
})
aucsList <- list(vol.auc.roc[1:4],vol.auc.roc[5:8], vol.auc.roc[9:12], vol.auc.roc[13:16])
names(aucsList) <- months
aucsList <- lapply(aucsList, function(i) {
  names(i) <- canUseNames
  i
})
prsList <- list(vol.PR[1:4],vol.PR[5:8], vol.PR[9:12], vol.PR[13:16])
names(prsList) <- months
prsList <- lapply(prsList, function(i) {
  names(i) <- canUseNames
  i
})
praucsList <- list(volauc.pr[1:4],volauc.pr[5:8], volauc.pr[9:12], volauc.pr[13:16])
names(praucsList) <- months
praucsList <- lapply(praucsList, function(i) {
  names(i) <- canUseNames
  i
})

# plot
scaleColor <- scale_color_brewer(type='qual',palette='Set1')

# single plots
plotROC('2-month', 'volThres_roc_2-month', vol.ROCs[1:4], vol.auc.roc[1:4], canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)
plotROC('6-month', 'volThres_roc_6-month', vol.ROCs[5:8], vol.auc.roc[5:8], canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)
plotROC('9-month', 'volThres_roc_9-month', vol.ROCs[9:12], vol.auc.roc[9:12],canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)
plotROC('12-month', 'volThres_roc_12-month', vol.ROCs[13:16], vol.auc.roc[13:16],canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)

plotPR('2-month', 'volThres_pr_2-month', vol.PR[1:4], volauc.pr[1:4], canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)
plotPR('6-month', 'volThres_pr_6-month', vol.PR[5:8], volauc.pr[5:8], canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)
plotPR('9-month', 'volThres_pr_9-month', vol.PR[9:12], volauc.pr[9:12], canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)
plotPR('12-month', 'volThres_pr_12-month', vol.PR[13:16], volauc.pr[13:16], canUseNames, volThresDir, perfTheme, scaleX, scaleY, scaleColor)

# group plots
tvColors <- wes_palette("Darjeeling1")[rev(c(1,4,2,5))]
tvLabels <- c('baseline',parse(text='mm^3'),'percent\nchange', 'rate\nchange')
tvLegendKey <- theme(legend.key.height=unit(1.25,'line'))
scaleCol <- scale_color_manual(values=tvColors,labels=tvLabels)
lo <- c(2,3,4,1)
tvROC <- plotMultiROC('volThres_allROCs', rocsList, aucsList, lo,
                      legendTitle = 'Tumor\nVolume Type', printIt = T, volThresDir,
                      perfTheme, scaleX, scaleY, perfGuide) + scaleCol + tvLegendKey

tvPR <- plotMultiPR('volThres_allPRs', prsList, praucsList, lo,
                    legendTitle = 'Tumor\nVolume Type', printIt = T, volThresDir,
                    perfTheme, scaleX, scaleY, perfGuide) + scaleCol + tvLegendKey


# adjust group plots
tvROC <- tvROC + theme(legend.text.align=0) ; print(tvROC)
tvPR <- tvPR + guides(color=F) + 
  theme(legend.text.align=0,strip.background = element_blank(), strip.text.x = element_blank()); print(tvPR)

tvPlots <- ggarrange(tvROC, tvPR, heights = c(1,.92), ncol=1)

setEPS(); postscript(file.path(figDir, paste0('s4_volThres_aucPR_aucROC.eps')), height=5, width = 9.3)
ggarrange(tvROC, tvPR, heights = c(1, .92), ncol=1)
dev.off()

# cuttoff plots
# plot thresholds for tumor volume (no logit
plotROC_threshold <- function(fn, vol_cutoffs, ind, aucs, listOrder=lo) {
  names(vol_cutoffs) <- paste0(months,' (area=',aucs,')')
  vol_cutoffs <- melt(vol_cutoffs,id.vars=c('x','y','cut'))
  
  vol_cutoffs$L1<- as.factor(vol_cutoffs$L1)
  vol_cutoffs$L1 <- factor(vol_cutoffs$L1, levels=levels(vol_cutoffs$L1)[listOrder])

  p <- ggplot(vol_cutoffs, aes(x=x, y=y, color=cut)) + 
    geom_abline(intercept = 0, slope = 1, colour = "gray") +
    geom_line(size=1) + 
    # geom_text(data=aucs, aes(x,y, label=label), show.legend = F, hjust=0, vjust=.05,size=2.5) + # aligned text to the left starting
    theme_linedraw() +
    theme(legend.text=element_text(size=5, color="#666666"),
          legend.title=element_text(size=7), 
          legend.position = 'right', 
          legend.key.size=unit(.6,'line'),
          legend.key.height = unit(1.75,'line'),
          legend.background = element_rect(fill = "transparent", colour = "transparent"),
          legend.margin=margin(l = -0.15, unit='cm')) +
    theme(axis.text=element_text(size=7, color="#666666"),
          axis.title.y = element_text(margin=margin(0,10,0,0),size=7),
          axis.title.x = element_text(margin=margin(10,0,0,0),size=7) ) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing = unit(.35, "lines"),
          strip.text.x = element_text(size=9,margin = margin(.05,0,.05,0, "cm"))) +
    scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.2)) +
    ylab("True positive rate") + 
    xlab("False positive rate")
  
  # scale_color_gradientn(colours = rainbow(10))
  mm3 <- c(0,1000,2000,4000,6000,8000,10000,20000,30000) #,250000)
  pc <- c(seq(-100,200,by=20))
  ratec <- c(seq(-100,1000, by=100))
  # myColors <- c('#5e4fa2','#3288bd','#66c2a5','#abdda4','#e6f598','#ffffbf','#fee08b','#fdae61','#f46d43','#d53e4f','#9e0142')
  if (ind=='mm3') p <- p + facet_grid(. ~ L1) + scale_color_distiller(palette='Spectral', name=parse(text='mm^3'), breaks=mm3,limits=c(0,30000), labels=c(mm3[1:8],'>=30000'), oob=squish)
  if (ind=='p') p <- p + facet_grid(. ~ L1) + scale_color_distiller(palette='Spectral', name='percent change', breaks=pc, limits=c(-100,max(pc)), labels=c(pc[1:length(pc)-1],'>=200'), oob=squish)
  if (ind=='r') p <- p + facet_grid(. ~ L1) + scale_color_distiller(palette = 'Spectral', name='rate change', breaks=ratec, limits=c(-100,1000), labels=c(ratec[1:length(ratec)-1],'>=1000'), oob=squish) +
    theme(legend.key.height = unit(1.5,'line'))
  if (ind=='b') p <- p + facet_grid(. ~ L1) + scale_color_distiller(palette='Spectral', name='baseline', breaks=mm3, limits=c(0,max(mm3)), labels=c(mm3[1:length(mm3)-1],'>=30000'), oob=squish)
  p
}


pmm3 <- plotROC_threshold(fn='threshold_roc_mm3',vol_cutoffs=vol.cutoffs[c(1,5,9,13)], ind='mm3', aucs=vol.auc.roc[c(1,5,9,13)])
ppercent <- plotROC_threshold('threshold_roc_percent',vol.cutoffs[c(2,6,10,14)], 'p', vol.auc.roc[c(2,6,10,14)])
prate <- plotROC_threshold('threshold_roc_rate',vol.cutoffs[c(3,7,11,15)], 'r', vol.auc.roc[c(3,7,11,15)])
pbaseline <- plotROC_threshold('threshold_roc_baseline',vol.cutoffs[c(4,8,12,16)], 'b', vol.auc.roc[c(4,8,12,16)])

ggarrange(prate, pmm3, pbaseline, ppercent, heights=c(2,2,2,2) , ncol=1)

setEPS(); postscript(file.path(figDir, paste0('s5_volThres_cutoffs.eps')), height=8, width = 8.66)
pp <- plot_grid(prate, pmm3, pbaseline, ppercent, labels = c("A", "B",'C','D'), ncol=1)
print(pp)
dev.off()

#
#
#-----------compare all three approaches to select best model --------------------
# approach 1: tumor vol without modeling results
vol.auc.roc <- lapply(volResults, '[[',3)
bestVol <- vol.auc.roc[c(names(which.max(vol.auc.roc[1:4]  )),
                         names(which.max(vol.auc.roc[5:8]  )),
                         names(which.max(vol.auc.roc[9:12] )),
                         names(which.max(vol.auc.roc[13:16])))]
bestVol <- data.frame(rocAUC=unlist(bestVol), 
                      rocAUCSD=rep(0,length(bestVol)),
                      lower = rep(NA,length(bestVol)),
                      higher = rep(NA,length(bestVol)),
                      name=unlist(lapply(strsplit(names(bestVol), '\\.'), '[[',2)),
                      month=months,
                      approach=1)

# approach 2: tumor vol model CV results
bestVolLogit <- lapply(months, function(m) {
  z <- volLogitResults[which(volLogitResults$month==m),]
  b <- z[which.max(z$rocAUC),]
  b <- b[,colnames(b) %in% c('rocAUC','rocAUCSD','month')]
  b$lower <- b$rocAUC - b$rocAUCSD
  b$higher <- b$rocAUC + b$rocAUCSD
  b$name <- rownames(b)
  b
})
bestVolLogit <- do.call('rbind', bestVolLogit)
bestVolLogit$approach <- 2

# approach 3: logits with spm patterns
bestLogitSpm <- bestLogitSpm[,c('rocAUC','rocAUCSD','month','name')]
bestLogitSpm$lower <- bestLogitSpm$rocAUC - bestLogitSpm$rocAUCSD
bestLogitSpm$higher <- bestLogitSpm$rocAUC + bestLogitSpm$rocAUCSD
bestLogitSpm$approach <- 3

# plot comparision across results
comp <- do.call('rbind', list(bestVol,bestVolLogit,bestLogitSpm))
comp$approach <- as.factor(comp$approach)
comp$approach <- factor(comp$approach, levels=levels(comp$approach), 
                        labels=c('tumor volume','tumor volume\nwith covariates', 'temporal patterns\nwith covariates'))
pd <- position_dodge(0.6) # move them .05 to the left and right

# include 12-month

# don't include
comparedPlot <- ggplot(comp[!(comp$month=='12-month'),], aes(x=month, y=rocAUC, fill=factor(approach)) ) +
  geom_bar(stat='identity',position=position_dodge()) +
  scale_y_continuous(limits = c(.5,1), oob=rescale_none) +
  theme_linedraw() +
  # ggtitle('') +
  perfTheme +
  geom_errorbar(aes(ymin=lower, ymax=higher), width=.2, position=position_dodge(.9)) +
  labs(fill='Best predictor\nin approach:') + xlab('Task') + ylab('AUC') +
  scale_fill_manual(values = wes_palette("Zissou1")[3:5] ) + 
  theme(legend.text.align=0, legend.key.height=unit(1.5,'line')) +
  guides(fill = guide_legend(override.aes = list(size=4))) 
print(comparedPlot)
dev.print(png, file.path(outDir,paste0('compare_approaches.','png')), res=600,height=5, width=7, units="in")




# ---- get final model performancs------

# model performance results, including cutoffs
retrainResults <- readRDS(file.path(finalDir, 'retrain_test_results.rds')) 

# group plot best model
# extract performance results
rocsList <- lapply(retrainResults, function(fit) {
  rocs <- list(train=fit$train$performance$roc,
               test =fit$test$performance$roc )
  rocs
})
aucsList <- lapply(retrainResults, function(fit) {
  aucs <- list(train=fit$train$performance$aucROC,
               test =fit$test$performance$aucROC )
  aucs
})
prsList <- lapply(retrainResults, function(fit) {
  prs <- list(train=fit$train$performance$pr,
               test =fit$test$performance$pr )
  prs
})
praucsList <- lapply(retrainResults, function(fit) {
  auc.pr <- list(train=fit$train$performance$aucPR,
               test =fit$test$performance$aucPR)
  auc.pr
})
names(rocsList) <- months # assumes months and order of retrainResults are the same
names(aucsList) <- months
names(prsList) <- months
names(praucsList) <- months

#
# plot ROC and PR of final model's test and train results
#
scaleColor <- scale_color_manual(values = c(wes_palette("Darjeeling1")[1], wes_palette('FantasticFox1')[3]) )

finalModelROC <- plotMultiROC(fn = '',
                              rocsList = rocsList[1:3],
                              aucsList = aucsList[1:3],
                              # listOrder = lo,
                              outDir = NULL,
                              legendTitle = 'Data partition', 
                              printIt=F, 
                              perfTheme = perfTheme,
                              scaleX = scaleX,
                              scaleY = scaleY,
                              perfGuide = perfGuide) + 
  theme(legend.text.align=0) + 
  scaleColor 
finalModelPR <- plotMultiPR(fn = '',
                            prsList = prsList[1:3], 
                            aucsList = praucsList[1:3], 
                            legendTitle = 'data\npartition',
                            outDir = NULL,
                            # listOrder = lo,
                            printIt=F, 
                            perfTheme = perfTheme, 
                            scaleX = scaleX,
                            scaleY = scaleY, 
                            perfGuide = perfGuide) + 
  guides(colour=F) + 
  theme(legend.text.align=0) + 
  scaleColor + 
  theme(legend.text.align=0,strip.background = element_blank(), strip.text.x = element_blank())
finalModelPlots <- ggarrange(finalModelROC, finalModelPR, heights=c(1,.92),ncol=1, padding = unit(0, "line"))
# dev.print(tiff, file.path(spmDir, paste0('all_ROC_PR','.tiff') ), res=600,height=5, width=7, units="in")

#
# plot prediction/probability distribution by label
#

testLabels <- lapply(retrainResults, function(fit) train=fit$test$preds@labels )
names(testLabels) <- months
testLabels <- unlist(testLabels, recursive = F)
testLabels <- melt(testLabels)
testLabels$value <- as.character(testLabels$value)
testLabels$value[testLabels$value=='dead'] <- '<= x-months'
testLabels$value[testLabels$value=='alive'] <- '> x-months'
testLabels$value <- as.factor(testLabels$value)

testPredictions <- lapply(retrainResults, function(fit) train=fit$test$preds@predictions )
names(testPredictions) <- months
testPredictions <- unlist(testPredictions, recursive = F)
testPredictions <- melt(testPredictions)

testPredictions$labels <- testLabels$value
names(testPredictions) <- c('prediction','month', 'Truth')
p <- ggplot(testPredictions, aes(x=prediction, color=Truth, fill=Truth, linetype=Truth) ) + 
  geom_density(alpha=0.1) +
  theme_bw() +
  theme(legend.text=element_text(size=7, color="grey30"),
        legend.title=element_text(size=7),
        legend.text.align =1,
        legend.position='right',
        legend.key.size=unit(.8,'line'),
        legend.margin=margin(l = -0.15, unit='cm'),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text.y=element_text(size=7, color="#666666"),
        axis.text.x=element_text(size=7,color='#666666'),
        axis.title.y = element_text(margin=margin(0,5,0,0),size=7, angle=90),
        axis.title.x = element_text(margin=margin(5,0,0,0),size=7),
        strip.text.x = element_text(color='white',size=9,margin = margin(.05,0,.05,0, "cm")),
        strip.background =  element_rect(fill='black'))+
  xlab('Model estimate (probability of residual survival \u2264 x-months)') + ylab('Density') +
  scale_color_manual(values = wes_palette("Darjeeling1")[c(1,2)], labels=c('\u2264 9 months','> 9 months') ) +
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(1,2)], labels=c('\u2264 9 months','> 9 months') ) +
  scale_linetype_manual(values=c('solid','twodash'), labels=c('\u2264 9 months','> 9 months')) +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25))
p <- p + facet_wrap( ~ month, scales='free') + theme(panel.spacing = unit(.5, "lines"), legend.text.align=0)
print(p)
dev.print(png, file.path(finalDir,paste0('testPrediction_density.','png')), res=600,height=3, width=8, units="in")


p <- ggplot(testPredictions[!(testPredictions$month=='12-month'),], 
            aes(x=prediction, color=Truth, fill=Truth, linetype=Truth) ) + 
  geom_density(alpha=0.1) +
  theme_bw() +
  theme(legend.text=element_text(size=7, color="grey30"),
        legend.title=element_text(size=7),
        legend.text.align =1,
        legend.position='right',
        legend.key.size=unit(.8,'line'),
        legend.margin=margin(l = -0.15, unit='cm'),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text.y=element_text(size=7, color="#666666"),
        axis.text.x=element_text(size=7,color='#666666'),
        axis.title.y = element_text(margin=margin(0,5,0,0),size=7, angle=90),
        axis.title.x = element_text(margin=margin(5,0,0,0),size=7),
        strip.text.x = element_text(color='white',size=9,margin = margin(.05,0,.05,0, "cm")),
        strip.background =  element_rect(fill='black'))+
  xlab('Model estimate (probability of residual survival \u2264 x-months)') + ylab('Density') +
  scale_color_manual(values = wes_palette("Darjeeling1")[c(1,2)] , labels=c('\u2264 9 months','> 9 months')) +
  scale_fill_manual(values = wes_palette("Darjeeling1")[c(1,2)] , labels=c('\u2264 9 months','> 9 months')) +
  scale_linetype_manual(values=c('solid','twodash'), labels=c('\u2264 9 months','> 9 months')) +
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,0.25))
p <- p + facet_wrap( ~ month, scales='free') + theme(panel.spacing = unit(.5, "lines"), legend.text.align=0)
# cairo_ps(file.path(spmDir, paste0('all_ROC_PR_density.eps')), height=7, width = 7)
# finalModelPlots <- ggarrange(finalModelROC, finalModelPR, p, heights=c(1,.92,.92),ncol=1)
# dev.off()
#



#
# ---- combine plots ------
# unicode fails, use cairo_ps

# setEPS(); postscript(file.path(figDir, paste0('Fig2_performance.eps')), height=7, width = 5)
# pp <- plot_grid(comparedPlot, finalModelPlots, p, labels = c("A", "B",'C','D'), ncol=1,
#           rel_heights = c(.8,1,.5))
# print(pp)
# dev.off()

cairo_ps(file.path(figDir, paste0('Fig2_performance.eps')), height=7, width = 5)
pp <- plot_grid(comparedPlot, finalModelPlots, p, labels = c("A", "B",'C','D'), ncol=1,
                rel_heights = c(.8,1,.5))
print(pp)
dev.off()
