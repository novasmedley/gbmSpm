#---------------------------------------- Description---------------------------------------- 
# model survival prediction using tumor measurements

# Load Sources and Libs
library(gbmSpm)
library(ggplot2)
library(scales)
library(ggrepel)
library(egg)
library(wesanderson)


dataDir <- '~/Data/Agile/analysis_vol/logits'
outDir <- dataDir

date <- format(Sys.Date(), format="%Y.%m.%d")

demographics <- c('+ageDecade+gender+ethnicity+MSP+laterality+location')
demoPreds <- c('ageDecade','gender','ethnicity','MSP','laterality:  right', 'laterality:  left',
               'location: frontal', 'location: temporal', 'location: parietal', 'location: other')
# set up what features should be used in logit modeling
# no lasso, use all features named

# tumor volume types, each individually and should be in the same order
# continuous types
contPreds <- c('continuous','percentChange', 'rateChange', 'baselineSize') 
# discretizied types
discPreds <- c('discreteContinuous','discretePercent','discreteRateChange',
               'responseCriteria', 'discreteBaseline')
# tumor volume types, each individually + covariaties
adjContPreds <- lapply(contPreds, function(i) c(i,demoPreds))
adjDiscPreds <- lapply(discPreds, function(i) c(i,demoPreds))

# list of different combinations of predictors to use as features in logit
predCombos <- c(contPreds, discPreds, adjContPreds, adjDiscPreds)
names(predCombos) <- c(contPreds, 
                       discPreds, 
                       paste0('adj.',contPreds), 
                       paste0('adj.d.',discPreds))
# graph names
comboNames <- c('mm^3','percent\nchange', 'rate\nchange', 'baseline', # continuous tyes
                'mm^3','percent\nchange', 'rate\nchange', 'response\ncritieria', 'baseline') # will be discreized types




# -------- logit CV plot functions -----------
pd <- position_dodge(0.2)

plotCVROC <- function(results, month) {
  filename <- file.path(subDir,paste0('volume_logit_CVresults_aucroc_',month) )
  g <- ggplot(results[results$month==month,], aes(x=type, y=ROC, line=type, linetype=covariates, shapes=discretized,color=discretized, group=line) ) +
    geom_point(size=1, position=pd) + geom_line(size=0.4, position=pd) +
    scale_y_continuous('average\nAUCROC', limits=c(0,1), breaks=seq(0,1,0.1), minor_breaks=seq(0,1,0.05)) +
    theme_linedraw() + 
    # ggtitle(paste0('Using only Tumor Volume: 10-fold CV Results for ',month, 'Survival Prediction') )+
    theme(axis.text=element_text(size=10),
          axis.title.y = element_text(angle = 0, vjust = .5, size=10, margin=margin(0,5,0,0)),
          axis.title.x = element_text(size=10,margin=margin(10,0,0,0)) ,
          legend.text=element_text(size=10, color='grey30'), 
          legend.title=element_text(size=10),
          # legend.position=c(0.8,0.8),
          legend.key.size=unit(.7,'line') ) +   geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD), width=0.01, position=pd) +
    labs(color='discretized',x='Tumor Volume Types')  +
    scale_color_manual(values = rev(c('#f1a340','#998ec3')))
  
  print(g)
  dev.print(png, paste0(filename,'.png'), res=600,height=3, width=6, units="in")
}
plotCVPR <- function(results, month) {
  filename <- file.path(subDir,paste0('volume_logit_CVresults_aucpr_',month) )
  g <- ggplot(results[results$month==month,], aes(x=type, y=PR, line=type, linetype=covariates, shapes=discretized,color=discretized, group=line) ) +
    geom_point(size=1) + geom_line(size=0.4) +
    scale_y_continuous('average\nAUCPR', limits=c(0,0.8), breaks=seq(0,0.8,0.1), minor_breaks=seq(0,0.8,0.05)) +
    theme_minimal() + 
    # ggtitle(paste0('Using only Tumor Volume: 10-fold CV Results for ',month, 'Survival Prediction') )+
    theme(axis.text=element_text(size=10),
          axis.title.y = element_text(angle = 0, vjust = .5, size=10, margin=margin(0,10,0,0)),
          axis.title.x = element_text(size=10,margin=margin(10,0,0,0)) ,
          legend.text=element_text(size=15, color='grey30'), 
          legend.title=element_text(size=15),
          # legend.position=c(0.8,0.8),
          legend.key.size=unit(.7,'line') ) +   # geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD), width=0.01, position=pd) +
    labs(color='discretized',x='Tumor Volume Types')  +
    scale_color_manual(values = rev(c('#f1a340','#998ec3')))
  
  print(g)
  dev.print(png, paste0(filename,'.png'), res=600,height=3, width=6, units="in")
}




#---------------------------------------- tumor vol logit results ----------------------------------------
# process cross validation training data
months <- c('2-month','6-month','9-month','12-month')

models2m <- readRDS(file.path(dataDir, 'tumorVol_logits_cvResults_survivalIn60.rds'))
models6m <- readRDS(file.path(dataDir, 'tumorVol_logits_cvResults_survivalIn180.rds'))
models9m <- readRDS(file.path(dataDir, 'tumorVol_logits_cvResults_survivalIn270.rds'))
models12m <- readRDS(file.path(dataDir, 'tumorVol_logits_cvResults_survivalIn360.rds'))

# extract averaged CV results
results <- lapply(list(models2m, models6m, models9m, models12m), function(month) {
  lapply(month, function(model) model$avg.CVPerformance )
})
names(results) <- months
results <- unlist(results, recursive = F, use.names = T)
results <- do.call('rbind', results)

# metadata
combos <- length(predCombos)
nModels <- length(months)
trainMonths <- c(rep( '2-month',combos), rep('6-month',combos), 
                 rep( '9-month',combos), rep('12-month',combos))
trainTypes <- as.factor(rep( comboNames, 2*nModels))
trainTypes <- factor(trainTypes, levels=levels(trainTypes)[c(5,2,4,1,3,6)] )
trainDiscrete <- as.factor(  rep(  c( rep('no', length(contPreds)  ), # plus 1 for the percentchange + baseline models
                                            rep('yes', length(discPreds) )   ), 2*nModels ))# 2x for adj and not adj, and 4x model months    
trainDiscrete <- factor(trainDiscrete, levels=levels(trainDiscrete)[c(2,1)])
trainCovars <- as.factor(rep( c(rep('no',length(contPreds)+length(discPreds) ) ,
                                       rep('yes',length(contPreds)+length(discPreds)) )   ,nModels) ) # 3x model months
trainCovars <- factor(trainCovars, levels=levels(trainCovars)[c(2,1)])
trainLines <- paste0(trainCovars, trainDiscrete)
trainGroups <- paste0(trainTypes,trainLines)

# attach metadata
results$month <- trainMonths
results$type <- trainTypes
results$discretized <- trainDiscrete
results$covariates <- trainCovars
results$line <- trainLines
results$group <- trainGroups

# save info
write.table(results, file=file.path(outDir,'tumorVol_logits_CV_results.txt'), sep='\t', row.names=F,col.names=T)
saveRDS(results, file=file.path(outDir, 'tumorVol_logits_CV_results.rds'))

# get just the best results
bestVolLogit <- lapply(months, function(m) {
  z <- results[which(results$month==m),]
  b <- z[which.max(z$rocAUC),]
  b <- b[,colnames(b) %in% c('rocAUC','rocAUCSD','month')]
  b$name <- rownames(b)
  b
})
bestVolLogit <- do.call('rbind', bestVolLogit)
write.table(bestVolLogit, file=file.path(outDir,'topPerformers_tumorVol_logits_CV.txt'), sep='\t', row.names=F,col.names=T)




