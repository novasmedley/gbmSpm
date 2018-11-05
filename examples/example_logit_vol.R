#---------------------------------------- Description---------------------------------------- 
# model survival prediction using tumor measurements
# Lasso is not used as predictors are predfined.


#---------------------------------------- set up ---------------------------------------- 
library(gbmSpm)
library(doMC)

# constants
opt <- list()
opt$outDir <- '~/test'
opt$saveLogit <- 'no'
opt$labelName <- 'survivalIn60'

registerDoMC(cores = 3)
date <- format(Sys.Date(), format="%Y.%m.%d")

# same as runLogits_vol.R so we may get the same data partitions and training procedure
seed <- 99
nFolds <- 3
cvReps <- 1
metric <- 'rocAUC'

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

#---------------------------------------- get tumor data ------------------------
# covariates to consider
demographics <- c('+ageDecade+gender+ethnicity+MSP+laterality+location')
demoPreds <- c('ageDecade','gender','ethnicity','MSP',
               'laterality:  right', 'laterality:  left',
               'location: frontal', 'location: temporal', 'location: parietal', 'location: other')

# follows the steps as runSPM.R, but will then only keep tumor volume data

# extract the same data as the ones used for SPM
data("fake_data")
fake_demo <- fake_data$demo
fake_tumorInfo <- fake_data$events
survData <- fake_data$person
  
# by default eventName is the volumetric response critera
# will get all types eventually
fake_data$events <- cleanData(fake_data$events) 

# keep tumor volume only
fake_data$events <- fake_data$events[fake_data$events$agent=='measurement',] 
fake_data$events <- prepTumorVol(fake_data$events) # get all tumor vol types

fake_data <- merge(fake_data$person, fake_data$events, by='iois')
fake_data <- prepDemographics(fake_data, fake_demo) # merge ethnicity data
fake_data <- prepSurvivalLabels(fake_data)
fake_data <- getTumorLocation(fake_data, fake_tumorInfo) # get first tumor location
fake_data <- prepLaterality(fake_data)
fake_data <- prepLocation(fake_data)
fake_data$id <- seq(1,nrow(fake_data))

cat('...overall samples: ', nrow(fake_data), '\n')
cat('...covariate in training data?: ', demoPreds %in% colnames(fake_data))

# general info
# pcRange <- data$percentChange
# pcRange <- pcRange[!(is.na(pcRange))]
# pcRange <- range(pcRange)
# 
# cRange <- range(data$continuous)

# training and testing partitions
train.ids <- getTrainTestPartition(data=fake_data,
                                   database=NULL, 
                                   personTable=NULL,
                                   survData = survData,
                                   seed=1,
                                   verbose=T)
fake_data <- fake_data[fake_data$id %in% train.ids,] #training data

close(logfn)
#---------------------------------------- classification task ---------------------------------------- 
# forumla
f <- as.formula(paste0(opt$labelName,'~.'))


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
# logit models using different predictor lists

# log model pipeline
i <- 1 # just do one

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

# fit model
fit <- runCV(data = fake_data, 
             formula= f, 
             labelName=opt$labelName,
             lasso=FALSE,
             llength = NULL,
             lmax = NULL,
             predictors=predCombos[[i]], 
             metric='rocAUC', 
             seed=seed, 
             folds=nFolds, 
             cvReps=cvReps,
             createModelMatrix = TRUE)
close(logfn_2)

if (opt$saveLogit=='yes') {
  mfn <- file.path(logitFolder, 
                   paste0('tumorVol_logits_',opt$labelName, '_', names(predCombos)[i],'.rds'))
  saveRDS(fit, file=mfn)
}
