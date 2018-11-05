library(gbmSpm)

library(gridExtra) # graph
library(grid) # graph
library(scales) # graph
library(wesanderson) # graph color


#------ local functions to format results after reading in -------------------------------------------------------------------------- 

getCVResults <- function(dir, fn, name) {
  cvResultsFile <- file.path(dir,fn)
  cvResults <- read.table(cvResultsFile, header=T, sep=',')
  
  # identify the cspade parameters the logit's features were based on
  temp <- do.call(rbind, lapply(cvResults$model, function(n) regmatches(n, gregexpr('[0-9]+',n))[[1]][2:5] ) )
  colnames(temp) <- c('support', 'gap', 'length', 'size')
  temp <- as.data.frame(temp)
  temp$support <- paste0( '.', as.character(temp$support))
  cvResults <- cbind(cvResults, temp)
  
  cvResults <- cvResults[order(cvResults$model),]
  cvResults$name <- rep(name, nrow(cvResults))
  
  ind <- length(strsplit(as.character(cvResults$model[1]), '\\/')[[1]])
  temp <- strsplit(as.character(cvResults$model),'\\/')
  cvResults$model <- unlist(lapply(temp, '[[', ind))
  cat('duplicates?: ', which(duplicated((cvResults))) )
  return(cvResults)
}

getCVByMonth <- function(cvResults, month, metricColName) {
  cv <- cvResults[cvResults$month==month,]
  cv <- cv[order(-cv[,metricColName]),]
  cv$model <- factor(cv$model, levels=unique(cv$model)) # order levels by decreasing ROC
  cv$name <- as.factor(cv$name)
  cv$group <- paste0(cv$tType, cv$name)
  cv$tType <- as.factor(cv$tType)
  cv$tType <- factor(cv$tType, levels=levels(cv$tType)[c(3,1,2,4)])
  cv
}

rankModels <- function(cvsResults, months, topN, metricColName) {
  r <- lapply(months, function(m) {
    models <- unique(cvsResults$model)
    ranked <- lapply(models, function(m) {
      results <- max(cvsResults[cvsResults$model==m,metricColName])
    })
    names(ranked) <- models
    ranked <- unlist(ranked)
    ranked <- ranked[order(-ranked)]
    names(ranked)[1:topN]
  })
  r <- unlist(r)
  r <- r[!duplicated(r)]
  return(r)
}

#------plotting functions-------------------------------------------------------------------------- 
# https://learnr.wordpress.com/2009/04/29/ggplot2-labelling-data-series-and-adding-a-data-table/
vplayout <- function(...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}
subplot <- function(x, y) viewport(layout.pos.row = x,
                                   layout.pos.col = y)
mmplot <- function(a, b) {
  vplayout()
  print(a, vp = subplot(1, 1))
  print(b, vp = subplot(2, 1))
}
# end
Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2,0.5), c("null", "null")))
pd <- position_dodge(0.6) # move them .05 to the left and right

# ROC
plotCVResults <- function(cv, month, filename) {
  g <- ggplot(cv, aes(x=model, y=rocAUC, colour=tType, line=tType, group=group) ) +
    geom_point() + geom_line() +
    scale_y_continuous('Average\nArea Under the ROC Curve\n',limits=c(0.4,0.7),breaks=seq(0.4,0.7,.05), minor_breaks=seq(0.4,0.7,0.01)) +
    theme_linedraw() + 
    ggtitle(month) +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.text=element_text(size=15),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,10,0,0)),
          legend.text=element_text(size=12, color='grey30'), 
          legend.title=element_text(size=12),
          # legend.position=c('none'),
          legend.position = c(0.7,0.8),
          legend.key.size=unit(.8,'line'),
          plot.title = element_text(size=15, hjust=0.5)) +  #geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD), width=0.001, position=pd) +
    # geom_line(position=pd) + geom_point(position=pd) + 
    labs(color='tumor volume type') +
    scale_color_manual(values = wes_palette("Darjeeling1")[-4] )
  
  table <- cv[!(duplicated(cv$model)),]
  table <- table[, c('model','size','length','gap','support')]
  # colnames(table) <- c('model', 'events', 'visits', 'days', 'patients')
  table <- melt(table, id='model')
  
  t <- ggplot(table, aes(x=model,y=variable, label=format(value))) +
    geom_text(size = 4.2, color='grey30', hjust=0.6) + theme_bw() + scale_y_discrete('cSPADE\nparameters') +  scale_x_discrete(expand=c(0.04,0)) +
    theme(panel.grid.major = element_blank(), legend.position = "none",
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text=element_text(size=12),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,11,0,5)) ) +
    theme(plot.margin = unit(c(0,.19, 0, 0), "lines")) + xlab(NULL) + ylab(NULL)
  
  dev.new()
  mmplot(g, t)
}

# PR
plotCVResults_PR <- function(cv, month, filename, lp) {
  g <- ggplot(cv, aes(x=model, y=prAUC, colour=tType, line=tType, group=group) ) +
    geom_point() + geom_line() +
    scale_y_continuous('Average\nArea Under the PR Curve\n',limits=c(0,1),breaks=seq(0,1,.05), minor_breaks=seq(0,1,0.01)) +
    theme_linedraw() + 
    ggtitle(month) +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.text=element_text(size=15),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,13,0,0)),
          legend.text=element_text(size=12, color='grey30'), 
          legend.title=element_text(size=12),
          legend.position=lp,
          # legend.position=c('none'),
          legend.key.size=unit(.8,'line'),
          plot.title = element_text(size=15, hjust=0.5)) +  #geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD), width=0.001, position=pd) +
    # geom_line(position=pd) + geom_point(position=pd) + 
    labs(color='tumor volume type') +
    scale_color_manual(values = wes_palette("Darjeeling1")[-4] )
  
  table <- cv[!(duplicated(cv$model)),]
  table <- table[, c('model','size','length','gap','support')]
  # colnames(table) <- c('model', 'events', 'visits', 'days', '% patients')
  table <- melt(table, id='model')
  
  t <- ggplot(table, aes(x=model,y=variable, label=format(value))) +
    geom_text(size = 4.2, color='grey30', hjust=0.6) + theme_bw() + scale_y_discrete('cSPADE\nparameters') +  scale_x_discrete(expand=c(0.04,0)) +
    theme(panel.grid.major = element_blank(), legend.position = "none",
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text=element_text(size=12),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,11,0,5)) ) +
    theme(plot.margin = unit(c(0,.19, 0, 0), "lines")) + xlab(NULL) + ylab(NULL)
  
  mmplot(g, t)
}

# id best model parameters ------------------------------------------------------------

# Parse cross-validation log files
dataDir <- '~/gbm_spm_example_multi/logits'
theseDirs <- list.dirs(path=dataDir, full.names=F, recursive=F)
theseDirs

# types and typenames are matched
types <- c('percent','rateChange') # tumor volume types expected in log filename
typenames <- c('percent change','rate change') # corresponding pretty names of tumor volume types
months <- c('2-month','6-month','9-month', '12-month') # month names for classification labels

cvs <- lapply(theseDirs, function(f) {
  cat(f, '\n')
  logFiles <- list.files(path=file.path(dataDir,f), pattern='\\.log$') # locate log files
  
  # extract results for each tumor volume type
  results <- lapply(seq(1,length(types)), function(i) {
    n <- logFiles[grepl(types[i],logFiles)] # read in a log file
    if (length(n) > 0) {
      log <- file.path(dataDir,f,n)
      con <- file(log)
      open(con)
      resultsCols <- read.table(con, skip=114, nrow=1, stringsAsFactors = F)  # expected results colnames
      r2 <- read.table(con, skip=0, nrow=1) # expected location of CV results for each classification
      r6 <- read.table(con, skip=25, nrow=1)
      r9 <- read.table(con, skip=25, nrow=1)
      r12 <- read.table(con, skip=25, nrow=1)
      close(con)
      
      r <- do.call('rbind', list(r2,r6,r9,r12)) # combine and format
      r <- r[,-1]
      colnames(r) <- c(resultsCols[1,])
      
      r$month <- months # attach additional info
      r$tType <- rep(typenames[i], nrow(r))
      r$model <- rep(f, nrow(r))
      temp <- do.call(rbind, lapply(r$model, function(n) regmatches(n, gregexpr('[0-9]+',n))[[1]][2:5] ) )
      colnames(temp) <- c('support', 'gap', 'length', 'size')
      temp <- as.data.frame(temp)
      temp$support <- paste0( '.', as.character(temp$support))
      r <- cbind(r, temp)
      r$name <- 'without'
      r
    } else {
      cat('\t ... file not found for ', typenames[i] , '\n')
    }
  })
  results <- do.call('rbind',results)
})
cvs <- do.call('rbind',cvs) 
cat('duplicates?: ', which(duplicated((cvs))) )
cvs


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

# individual plots
spmROC <- lapply(months, function(m) plotCVResults(getCVByMonth(cvResults=cvs, m, 'rocAUC'), 
                                                   m, 
                                                   '' ) )
spmPR <- lapply(months, function(m) plotCVResults_PR(getCVByMonth(cvResults=cvs, m, 'prAUC'),
                                                     m,
                                                     '',
                                                     c(0.6,0.8)) )


# cv logit results
bestLogitSpm <- lapply(months, function(m) {
  z <- cvs[which(cvs$month==m),]
  b <- z[which.max(z$rocAUC),]
  b <- b[,colnames(b) %in% c('rocAUC','rocAUCSD','tType','month','model', 'gap','length','name')]
  rownames(b) <- b$model
  b
})
bestLogitSpm <- do.call('rbind', bestLogitSpm)
colnames(bestLogitSpm) <- c('rocAUC','rocAUCSD','month','name','model','maxgap','maxlength')
bestLogitSpm$label <- getClassLabels() # assumes row order already matches
bestLogitSpm$fitname <- c('fit2m','fit6m','fit9m','fit12m')
bestLogitSpm

# approach 3: logits with spm patterns
comp <- bestLogitSpm[,c('rocAUC','rocAUCSD','month','name')]
comp$lower <- comp$rocAUC - comp$rocAUCSD
comp$higher <- comp$rocAUC + comp$rocAUCSD
comp$approach <- 3

# plot comparision across results
comp <- do.call('rbind', list(comp))
comp$approach <- as.factor(comp$approach)
comp$approach <- factor(comp$approach, levels=levels(comp$approach), 
                        labels=c( 'temporal patterns\nwith covariates'))
pd <- position_dodge(0.6) # move them .05 to the left and right

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

# ---- retrain and apply to test set ------
seed <- 9
metric <- 'rocAUC'
expDir <- dirname(dataDir) 
expDir
dataPartitions <- file.path(expDir, 'train_test_partitions.rds') # expected ids of train and test, if doesn't exist, will create it
dataPartitions
partitions <- readRDS(file=dataPartitions) # data partitions

# metadata, unwanted features, and 
# unwanted labels should not be in training data for glmnet
# only one label will be kept in this set, depending on which model is being fit
#
needToRemove <- c('id','iois','eventID', # remove ids
                  bestLogitSpm$label, # remove labels
                  'IDH1') # not interested

# applies to all rows in 'bestLogitSpm':
selectedTType <- 'rate'
selectedTTypeFN_spm <- 'featureVectors_rateChange.rds'
selectedTTypeFN_cvLogit <- 'logits_CVresults_rateChange.rds'

# model performance results, including cutoffs
retrainResults <- lapply(seq(1, nrow(bestLogitSpm)), function(i) {
  selected <- bestLogitSpm[i,]
  
  cat('\n....', selected$label ,' ...\n')
  
  # cv results
  cvLogitPath <- file.path(expDir,'logits', selected$model, selectedTTypeFN_cvLogit)
  cvResults <- readRDS(cvLogitPath)
  cvResults <- cvResults[[selected$fitname]]
  
  # get lambda from best cv results
  cvResults <- cvResults$avg.CVPerformance[which.max(cvResults$avg.CVPerformance$rocAUC),]
  cat('\n\n','Summary to repeated CV results (lambdas for retraining) ...\n\n')
  print(cvResults)
  
  
  # get data
  dataPath <- file.path(expDir,'spm', selected$model, selectedTTypeFN_spm)
  data <- readRDS(dataPath) # features and labels for each clinical visit
  names <- colnames(data)
  data <- data.frame(id=row.names(data),data)
  colnames(data) <- c('id',names)
  cat('... \n overall samples: ', nrow(data), '\n')
  
  data <- prepLaterality(data) 
  data <- prepLocation(data)
  
  cat('... \n removing clinical visits with no history of length maxgap:', selected$maxgap, '\n')
  data <- removeVisits(data = data, 
                       maxgap = as.numeric(as.character(selected$maxgap)), 
                       maxlength = as.numeric(as.character(selected$maxlength)), 
                       tType = selectedTType, 
                       save=F, 
                       outDir=outputDir)
  start.local <- Sys.time()
  cat('\n....', cvResults$lambda ,' ...\n')
  label <- selected$label
  
  fit <- retrainLogit(data = data[data$id %in% partitions$train,],  #training data
                      formula = as.formula(paste0(selected$label,'~.')), 
                      labelName = selected$label, 
                      lambda = cvResults$lambda,
                      lasso=TRUE, # using lamba
                      needToRemove=needToRemove, 
                      createModelMatrix=FALSE,
                      metric=metric, 
                      verbose=TRUE,
                      seed=seed)
  cat('\n'); print(Sys.time() - start.local)
  
  trainResults <-  getTrainResults(fit) 
  testResults <- getTestResults(fit = fit,
                                data = data[data$id %in% partitions$test,], #testing data 
                                labelName = selected$label, 
                                formula = as.formula(paste0(selected$label,'~.')), 
                                needToRemove=needToRemove, 
                                createModelMatrix=FALSE)
  
  testResults$selectedThreshold <-  trainResults$selectedThreshold # use the threshold selected from train
  
  list(fit=fit, train=trainResults, test=testResults)
})
names(retrainResults) <- bestLogitSpm$fitname


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
