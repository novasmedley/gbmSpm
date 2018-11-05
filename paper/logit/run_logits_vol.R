#---------------------------------------- Description---------------------------------------- 
# model survival prediction using tumor measurements
# Lasso is not used as predictors are predfined.


#---------------------------------------- cmd line arguments---------------------------------------- 
# e.g., Rscript run_logits_vol.R --outDir ~/Data/analysis_vol --saveLogit yes --labelName survivalIn60

library(gbmSpm)
library(optparse)

option_list = list(make_option(c('--labelName'), type='character', default=NULL, help='logit label name', metavar='character'),
                   make_option(c('--dataPartFile'), type='character', default=NULL, help='data partition rds file, full path', metavar='character'),
                   make_option(c('--outDir'), type='character', default=NULL, help='output directory', metavar='character')
                   )
opt_parser = OptionParser(option_list = option_list)
opt = parse_args((opt_parser))

# testing
# opt <- list()
# opt$outDir <- '~/test'
# opt$labelName <- 'survivalIn60'
# opt$dataPartFile <- '~/test/train_test_partitions.rds'
# end testing

# library(doMC)
#---------------------------------------- constants ---------------------------------------- 
# registerDoMC(cores = 3)
date <- format(Sys.Date(), format="%Y.%m.%d")

# same as runLogits_spm.R so we may get the same data partitions and training procedure
seed <- 9
nFolds <- 10
cvReps <- 3
metric <- 'rocAUC' 
# no lasso parameters since no need for feature selection, they're set in this analysis

group <- 'agile' # group name in MySQL configuration file
patTabName <- 'patient'
eventTabName <- 'agent6'


#---------------------------------------- setup experiment ---------------------------------------- 
# outputDir is the main folder holding logit experiments
# subDir the name of the logit experiment folder within outputDir/logits
logitFolder <- file.path(opt$outDir, 'logits')
ifelse(!dir.exists(logitFolder), dir.create(logitFolder, recursive = T), F) 


logFile <- file.path(logitFolder, 
                     paste0('log_tumorVol_logits_',opt$labelName, '_data.txt'))
#---------------------------------------- log data ---------------------------------------- 

logfn <- file(logFile, open='wt')

sink(logfn, type='output', split=T)

cat(format(Sys.Date(), format="%Y.%m.%d"), '\n')
cat('label:', opt$labelName, '\n')
cat('seed: ', seed, '\n')
cat('nFolds: ', nFolds, '\n')
cat('cvReps: ', cvReps, '\n')
cat('lasso: FALSE','\n\n')

#---------------------------------------- get tumor data from database ------------------------
# covariates to consider
demographics <- c('+ageDecade+gender+ethnicity+MSP+laterality+location')
demoPreds <- c('ageDecade','gender','ethnicity','MSP',
               'laterality:  right', 'laterality:  left',
               'location: frontal', 'location: temporal', 'location: parietal', 'location: other')

# follows the steps as runSPM.R, but will then only keep tumor volume data
# patient 5902 is a good example of how different tumor vol types can be

# extract the same data as the ones used for SPM
data <- getData(database = group, 
                personTable = patTabName, 
                eventTable = eventTabName)
demo <- getDemographics(database=group, 
                        personTable = patTabName)
tumorInfo <- data$event # save tumor location and laterility strings before event cleaning

# by default eventName is the volumetric response critera
# will get all types eventually
data$event <- cleanData(data$event) 

# keep tumor volume only
data$event <- data$event[data$event$agent=='measurement',] 
data$event <- prepTumorVol(data$event) # get all tumor vol types

data <- merge(data$person, data$event, by='iois')
data <- prepDemographics(data,demo) # merge ethnicity data
data <- prepSurvivalLabels(data)
data <- getTumorLocation(data, tumorInfo) # get first tumor location

data <- prepLaterality(data)
data <- prepLocation(data)
data$id <- paste0(data$iois,'.',data$eventID) # clinical ids

cat('...overall samples: ', nrow(data), '\n')
cat('...covariate in training data?: ', demoPreds %in% colnames(data))

# general info
# pcRange <- data$percentChange
# pcRange <- pcRange[!(is.na(pcRange))]
# pcRange <- range(pcRange)
# 
# cRange <- range(data$continuous)

# data partition file used in spm logits
partition <- readRDS(file=opt$dataPartFile)

close(logfn)

#---------------------------------------- set predictors ---------------------------------------- 

# set up what features should be used in logit modeling
# no lasso, use all features named

# tumor volume types, each individually
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


  
#---------------------------------------- fit logit models ---------------------------------------- 

# classification task
f <- as.formula(paste0(opt$labelName,'~.'))

# logit models using different predictor lists

cvResults <- lapply(seq(1,length(predCombos)), function(i)  {
  # log model pipeline
  subLogFile <- file.path(logitFolder, 
                          paste0('log_tumorVol_logits_',opt$labelName, '_', names(predCombos)[i],'.txt'))
  if (file.exists(subLogFile)) file.remove(subLogFile)
  logfn_2 <- file(subLogFile, open='wt')
  sink(logfn_2, type='output', split=T)
  
  cat(format(Sys.Date(), format="%Y.%m.%d"), '\n')
  cat('label:', opt$labelName, '\n')
  cat('seed: ', seed, '\n')
  cat('nFolds: ', nFolds, '\n')
  cat('cvReps: ', cvReps, '\n')
  cat('lasso: FALSE','\n\n')
  cat('predictors in logit: ', names(predCombos)[i])
  
  # get cross validation results
  results <- runCV(data = data[data$id %in% partition$train,], 
                   formula= f, 
                   labelName=opt$labelName,
                   lasso=FALSE,
                   llength = NULL,
                   lmax = NULL,
                   predictors=predCombos[[i]], 
                   metric='rocAUC', # just quiets warning from caret, not used since no lasso parameter selection
                   seed=seed, 
                   folds=nFolds, 
                   cvReps=cvReps,
                   createModelMatrix = TRUE)
  
  # to compare a single ROC curve against other approaches
  trainResults <- retrainLogit(data = data[data$id %in% partition$train,],
                          formula = f, 
                          labelName = opt$labelName,
                          lambda = NULL, 
                          lasso = FALSE,
                          predictors=predCombos[[i]], 
                          createModelMatrix = TRUE,
                          metric='rocAUC', 
                          seed=seed, 
                          verbose=TRUE)
  saveRDS(trainResults, 
          file=file.path(logitFolder, paste0('tumorVol_logits_trainResults_',opt$labelName, '_', names(predCombos)[i],'.rds')))
  
  close(logfn_2)
  
  results
})
names(cvResults) <- names(predCombos)

saveRDS(cvResults, file=file.path(logitFolder,  
                                  paste0('tumorVol_logits_cvResults_',opt$labelName,'.rds')))

