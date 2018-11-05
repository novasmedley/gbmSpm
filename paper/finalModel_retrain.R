#------ description ----

# extract data from final model for survival predictions
# save data (de-identified) for further plotting 

dataDir <- '~/Data/Agile/analysis_spm' # all data used for logits

library(gbmSpm)
library(RMySQL)

library(gridExtra)
library(grid)
library(wesanderson)
library(viridis)
library(scales) # to truncate geom_bars

# install.packages('devtools')
# devtools::install_github("baptiste/egg")
library(egg)

#- functions ----
getPredictions <- function(result, info, month, patSurvival) {
  predictions <- data.frame(pred=result$pred@predictions, labels=result$pred@labels)
  colnames(predictions) <- c('prediction', 'label')
  rownames(predictions) <- info$id
  
  predictions$pat <- as.factor(unlist(lapply(info$id, function(i) strsplit(i, '\\.')[[1]][1])))
  predictions$eventID <- as.integer(unlist(lapply(info$id, function(i) strsplit(i, '\\.')[[1]][3])))
  predictions$label_pred <- as.integer(ifelse(predictions$prediction>result$selectedThreshold['cutoff',], 2,1))
  predictions$correct <- ifelse(as.integer(predictions$label)==predictions$label_pred, 'yes','no')
  predictions$patSurv <- patSurvival$survival[match(predictions$pat,patSurvival$iois)]
  predictions$patStat <- patSurvival$status[match(predictions$pat,patSurvival$iois)]
  
  if(month=='2-month') { predictions$classLine <- predictions$patSurv - 60
  } else if (month=='6-month') {predictions$classLine <- predictions$patSurv - 180
  } else if (month=='9-month'){ predictions$classLine <- predictions$patSurv - 270
  } else { predictions$classLine <- predictions$patSurv - 360
  }
  
  ind <- which(predictions$patStat=='ALIVE')
  predictions$patSurv[ind] <- NA
  predictions$classLine[ind] <- NA
  
  predictions$name <- paste0(predictions$pat,' (',predictions$patSurv,')')
  row.names(predictions) <- seq(1,nrow(predictions))
  return(predictions)
}
getID <- function(data) {
  ids <- row.names(data)
  patients <- unique(unlist(lapply(ids, function(i) strsplit(i, '\\.')[[1]][1])))
  return(list(id=ids, pat=patients))
}
getPatID <- function(data) {
  ids <- row.names(data)
  patients <- unlist(lapply(ids, function(i) strsplit(i, '\\.')[[1]][1]))
  return(patients)
}
getDayID <- function(data) {
  ids <- row.names(data)
  days <- unlist(lapply(ids, function(i) strsplit(i, '\\.')[[1]][3]))
  return(days)
}



# constants and data ----------
dataDir <- '~/Data/Agile/analysis_spm'
outputDir <- '~/Data/Agile/final_model'
months <- c('2-month','6-month','9-month', '12-month')
metric <- 'aucROC'
seed <- 9

# database connection
group <- 'agile' # group name in MySQL configuration file
patTabName <- 'patient'
eventTabName <- 'agent6'

dataPartitions <- file.path(dataDir, 'train_test_partitions.rds') # expected ids of train and test, if doesn't exist, will create it


# cv logit results
topCVS <- readRDS(file.path(dataDir,'spm_logit_top_cv_results.rds'))
bestLogitSpm <- lapply(months, function(m) {
  z <- topCVS[which(topCVS$month==m),]
  b <- z[which.max(z$rocAUC),]
  b <- b[,colnames(b) %in% c('rocAUC','rocAUCSD','tType','month','model', 'gap','length')]
  rownames(b) <- b$model
  b
})
bestLogitSpm <- do.call('rbind', bestLogitSpm)
colnames(bestLogitSpm) <- c('rocAUC','rocAUCSD','month','name','model','maxgap','maxlength')

# labels
bestLogitSpm$label <- getClassLabels() # assumes row order already matches
bestLogitSpm$fitname <- c('fit2m','fit6m','fit9m','fit12m')


# applies to all rows in 'bestLogitSpm':
selectedTType <- 'rate'
selectedTTypeFN_spm <- 'featureVectors_rateChange.rds'
selectedTTypeFN_cvLogit <- 'logits_CVresults_rateChange.rds'

#---------------------------------------- log data ---------------------------------------- 

logfn <- file(file.path(outputDir,
                        paste0(format(Sys.Date(), format="%Y.%m.%d"),
                               '_finalModel_retrain','.log')), 
              open='wt')

sink(logfn, type='output', split=T)

cat(format(Sys.Date(), format="%Y.%m.%d"), '\n')
cat('selected pattern dataset:', selectedModel, '\n')
cat('selected tumor vol variable type: ', getVolTypeName(selectedTType) , '\n')


#---------------------------------------- get data  partitions ---------------------------------------- 
partition <- getTrainTestPartition(dataPartitionFile = dataPartitions,
                                   createParts = F,
                                   data = data,
                                   verbose = T,
                                   seed = NULL,
                                   database = NULL,
                                   personTable = NULL)

#---------------------------------------- get patient survival data   ---------------------------------------- 
patSurvival <- getSurvData(database=group, personTable=patTabName) # includes MGMT biomarker

#---------------------------------------- get event data   ---------------------------------------- 





# metadata, unwanted features, and 
# unwanted labels should not be in training data for glmnet
# only one label will be kept in this set, depending on which model is being fit
#
needToRemove <- c('id','iois','eventID', # remove ids
                  bestLogitSpm$label, # remove labels
                  'IDH1') # not interested


#---------------------------------------- retrain data  ---------------------------------------- 
start.total <- Sys.time()
cat('\n\n','Retrain logit with selected lambda from repeated CV...\n\n')

partitions <- readRDS(file=dataPartitions) # data partitions

results <- lapply(seq(1, nrow(bestLogitSpm)), function(i) {
  selected <- bestLogitSpm[i,]
  
  cat('\n....', selected$label ,' ...\n')
  
  # cv results
  cvLogitPath <- file.path(dataDir,'logits', selected$model, selectedTTypeFN_cvLogit)
  cvResults <- readRDS(cvLogitPath)
  cvResults <- cvResults[[selected$fitname]]
  
  # get lambda from best cv results
  cvResults <- cvResults$avg.CVPerformance[which.max(cvResults$avg.CVPerformance$rocAUC),]
  cat('\n\n','Summary to repeated CV results (lambdas for retraining) ...\n\n')
  print(cvResults)
  
  
  # get data
  dataPath <- file.path(dataDir,'spm', selected$model, selectedTTypeFN_spm)
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
names(results) <- bestLogitSpm$fitname
saveRDS(results, file.path(outputDir, 'retrain_test_results.rds'))

close(logfn)
# ------------------------ get data used for training and test in each model------------------------
trainInfo <- lapply(results, function(model) getID(model$fit$trainingData) ) 
testInfo <- lapply(results, function(model) getID(model$test$data) )

# ------------------------  get plotable performance of best model at the patient level -------------------------------------------------------------------------------- 

# note that labels are 1 and 2 (2 being the positive class)

testPredictions <- lapply(seq(1,length(months)), function(i) getPredictions(result = results[[i]]$test, 
                                                                            info = testInfo[[i]], 
                                                                            month = months[i],
                                                                            patSurvival = patSurvival) )

trainPredictions <- lapply(seq(1,length(months)), function(i) getPredictions(result = results[[i]]$train, 
                                                                             info = trainInfo[[i]], 
                                                                             month = months[i],
                                                                             patSurvival = patSurvival) )
testPredictions <- lapply(testPredictions, function(tp) {
  tp$label <- as.integer(tp$label)
  tp$label[tp$label==1] <- 'alive'
  tp$label[tp$label==2] <- 'dead'
  tp
})

saveRDS(trainPredictions, file=file.path(outputDir,'plotReady_trainPredictions.rds'))
saveRDS(testPredictions, file=file.path(outputDir,'plotReady_testPredictions.rds'))


# ------------------------  proportion of clinical events per patient----
spmFreq <- lapply(seq(1,length(months)), function(i) {
  r <- lapply(list(trainPredictions[[i]], testPredictions[[i]]), function(d) {
    nPat <- length(unique(d$pat))
    nVisits <- length(unique(d$eventID))
    freq <- lapply(split(d,d$pat), function(p) length(unique(p$eventID)) )
    names(freq) <- seq(1,length(freq))
    freq <- unlist(freq)
    list(nPat=nPat,nVisits=nVisits,freq=freq)
  })
  
  
  return(list(train=r[[1]], test=r[[2]]))
})
names(spmFreq) <- months
saveRDS(spmFreq, file=file.path(outputDir,'spm_sampleFreq.rds'))



# ------------------------   get coefficients of best model -------------------------------------------------------------------------------- 
getModelStruct <- function(fit) {
  print(table(fit$trainingData$.outcome)) # class label sizes
  
  cfs <- coef(fit$finalModel, fit$bestTune$lambda)
  cfs2 <- summary(coef(fit$finalModel, fit$bestTune$lambda))
  cfs2$var <- cfs@Dimnames[[1]][cfs2$i]
  cfs2 <- cfs2[,c('x','var')]
  names(cfs2) <- c('coefficient','variable')
  cat('number of variables in model:',nrow(cfs2),',','\n')
  cfs2
}

allCoefs <- lapply(results, function(result) getModelStruct(result$fit) )  # coefs
allCoefs <- Reduce(function(a,b) merge(a,b, by='variable',all=T), allCoefs)
allCoefs[,2:5] <- lapply(seq(2,5), function(i) round(exp(allCoefs[,i]), digits=2) ) # convert to odds ratios
colnames(allCoefs) <- c('variable', paste0(c('2m','6m','9m','12m'),'.','rateChange'))
write.table(allCoefs, file=file.path(outputDir,'allCoefficients.csv'), sep=',', row.names=F, col.names=T)

# model specific coefs
vars <- allCoefs$variable
vars <- gsub('\\},\\{',' --> ',vars)
vars <- gsub('<\\{','',vars)
vars <- gsub('\\}>','',vars)
vars <- gsub('`','',vars)
allCoefs$variable <- vars
write.table(na.omit(allCoefs[,c(1,2)]), file=file.path(outputDir,'coefs_2m.txt'), sep='\t',row.names=F)
write.table(na.omit(allCoefs[,c(1,3)]), file=file.path(outputDir,'coefs_6m.txt'), sep='\t',row.names=F)
write.table(na.omit(allCoefs[,c(1,4)]), file=file.path(outputDir,'coefs_9m.txt'), sep='\t',row.names=F)
write.table(na.omit(allCoefs[,c(1,5)]), file=file.path(outputDir,'coefs_12m.txt'), sep='\t',row.names=F)
