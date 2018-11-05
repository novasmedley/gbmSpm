#---------------------------------------- Description---------------------------------------- 
# Predict survival using tumor measurements. Tumor volume is used a single predictor. 
# There is no modeling involved
# Should have very similar results to logit models using only tumor volume.

library(gbmSpm)

outDir <- '~/Data/Agile/analysis_vol/threshold'
dataPartFile <- '~/Data/Agile/analysis_spm/train_test_partitions.rds'
seed <- 9
group <- 'agile' # group name in MySQL configuration file
patTabName <- 'patient'
eventTabName <- 'agent6'

#---------------------------------------- log data ---------------------------------------- 

logfn <- file(file.path(outDir,
                        paste0(format(Sys.Date(), format="%Y.%m.%d"),
                               '_volThreshold','.log')), 
              open='wt')

sink(logfn, type='output', split=T)

cat(format(Sys.Date(), format="%Y.%m.%d"), '\n')

#---------------------------------------- thresholding method --------------------------------------- 
# types of tumor volumes that can be used for thresholding
canUse <- c('continuous','percentChange', 'rateChange', 'baselineSize') # names
canUseNames <- c('mm^3', 'percent change', 'rate change', 'baseline') # pretty names

getResults <- function(trainData,labelName) {
  # no model needs to be created, straight to prediction using tumor vol as the reponse
  # ordered factors lower < upper = lower is negative class, upper is positive class
  # here, alive is the negative class, dead is the positive class
  
  d <- lapply(canUse, function(tVolType) {

    # predict survival with tumor volume
    trainProb <- ROCR::prediction(trainData[,tVolType], trainData[,labelName]) 
    
    # get performance stats
    perf <- getROC_PR(trainProb) # overall performance (all thresholds)
    
    cutData <- selectROCThreshold(trainProb) # single threshold performance
    cutData[6:11] <- round(cutData[6:11], digits=3)
    
    return(list(cutData=as.data.frame(cutData), 
                roc=perf$roc, 
                auc=perf$aucROC, 
                pr=perf$pr, 
                auc.pr=perf$aucPR,
                cutoffs=trainProb@cutoffs[[1]],
                trainProb=trainProb,
                trainLabels=trainData[,labelName]))
  })
  names(d) <- paste0(labelName,'.',canUse)
  return(d)
}

#---------------------------------------- get data --------------------------------------- 

# get the same data used for tumor logits
# follows the steps as runSPM.R, but will then only keep tumor volume data
# patient 5902 is a good example of how different tumor vol types can be

# extract the same data as the ones used for SPM
data <- getData(database = group, 
                personTable = patTabName, 
                eventTable = eventTabName)

# by default eventName is the volumetric response critera
# will get all types eventually
data$event <- cleanData(data$event) 

# keep tumor volume only
data$event <- data$event[data$event$agent=='measurement',] 
data$event <- prepTumorVol(data$event) # get all tumor vol types

data <- merge(data$person, data$event, by='iois')
data <- prepDemographics(data, getDemographics(database=group,personTable = patTabName))# need gender for data partitioning

data <- prepSurvivalLabels(data) # class labels for each event
data$id <- paste0(data$iois,'.',data$eventID) # clinical ids

cat('...overall samples: ', nrow(data), '\n')
cat('...overall patients: ', length(unique(data$iois)), '\n')


#----------------------------------------  get results ---------------
#  tumor vol AS A PREDICTOR - no modeling

# there is no parameter selection, therefore no need for cross-validation
# use training partition to evaluate the accuracy of using tumor vol as a prediction

labels <- getClassLabels()

# data partition file used in spm logits
partition <- getTrainTestPartition(dataPartitionFile = dataPartFile ,
                                   createParts = F, # used one created already
                                   data = data,
                                   verbose = T,
                                   seed = seed,
                                   database = group,
                                   personTable = patTabName)
trainData <- data[data$id %in% partition$train  , ]
testData <- data[!(data$id %in% partition$train) , ]
rm(data)

# get results
cat('\n... getting predictions ...\n')
volResults <- lapply(labels, function(ln) getResults(trainData = trainData,
                                                     labelName = ln))
volResults <- unlist(volResults, recursive=F)
saveRDS(volResults, file=file.path(outDir,'tumorVol_threshold_results.rds'))

# print top results
cat('\n... best predictors in each task: ...\n')

vol.auc.roc <- lapply(volResults, '[[',3)
bestVol <- vol.auc.roc[c(names(which.max(vol.auc.roc[1:4]  )),
                         names(which.max(vol.auc.roc[5:8]  )),
                         names(which.max(vol.auc.roc[9:12] )),
                         names(which.max(vol.auc.roc[13:16])))]
print(bestVol)

#----------------------------------------  get meta data ---------------

# proportion of clinical events per patient
tumorVolFreq <- lapply(list(trainData, testData), function(i) {
    nPat <- length(unique(i$iois))
    nVisits <- length(unique(i$startdate))
    freq <- lapply(split(i,i$iois), function(p) length(unique(p$startdate)) )
    names(freq) <- seq(1,length(freq))
    freq <- unlist(freq)
    list(nPat=nPat,nVisits=nVisits,freq=freq)
})
names(tumorVolFreq) <- c('trainPart','testPart')
saveRDS(tumorVolFreq, file=file.path(outDir,'tumorVol_threshold_visitFreq.rds'))

#----------------------------------------  get cut off table ---------------

cutData <- Reduce(function(a,b) cbind(a,b), lapply(volResults, '[[',1) )
colnames(cutData) <- names(volResults)
cutData <- t(cutData)
cutData <- as.data.frame(cutData)
# check rowname match the following attached columns:
cutData$month <- c( rep('2-month', 4), 
                    rep('6-month', 4), 
                    rep('9-month', 4),
                    rep('12-month', 4))
cutData$volumeType <- canUseNames
cutData$samples <- unlist(lapply(canUse, function(i) length(na.omit(trainData[,i])))) # number of samples in train
cutData <- cutData[,c(13,14,seq(1,12))]
write.table(cutData, file=file.path(outDir,'tumorVol_selectedThreshold_performances.txt'), sep='\t')

close(logfn)