# Load spm pattern data,perform logit analysis using spm events + demographics information 

# saves fit data for each survival prediction

# get help: in terminal, Rscript runSpmLogits.R -h

#-------------------------------------------------------------------------------- #
library(gbmSpm)
library(optparse)

option_list = list(make_option(c('--tType'), type='character', default='c', help='tumor volume type: 
                               p (percent change), 
                               c (continuous, original volume units), 
                               r (volumetric response criteria), 
                               rate (rate change);
                               see "getVolTypeName"', metavar='character'),
                   make_option(c('--maxgap'), type='numeric', default=60, help='cSPADE parameter: max gap between elements', metavar='character'),
                   make_option(c('--maxlength'), type='numeric', default=2, help='cSPADE parameter: max length elements', metavar='character'),
                   make_option(c('--lmax'), type='numeric', default=0.2, help='max lamba for lasso', metavar='character'),
                   make_option(c('--llength'), type='numeric', default=1000, help='number of lambas to consider for lasso', metavar='character'),
                   make_option(c('--dataFolder'), type='character', default=NULL, help='relative directory of spm pattern data file name', metavar='character'),
                   make_option(c('--prefix'), type='character', default=NULL, help='where logit data will be stored', metavar='character'),
                   make_option(c('--dir'), type='character', default=NULL, help='input and output directory', metavar='character'),
                   make_option(c('--saveLogit'), type='character', default='no', help='save logit fit CV results, yes or no', metavar='character')
                   )
opt_parser = OptionParser(option_list = option_list)
opt = parse_args((opt_parser))

# testing
# opt <- list()
# opt$tType <- 'rate'
# opt$dataFolder <- 'sup0.4g60l2z2'
# opt$prefix <- 'logits'
# opt$dir <- '~/test'
# opt$maxgap <- 60
# opt$maxlength <- 2
# opt$saveLogit <- 'yes'
# opt$llength <- 100
# opt$lmax <- 0.2
# # end testing

#---------------------------------------- setup experiment ---------------------------------------- 

# Note:
#' Assumes feature have been created in dataFolder
#' Will not overwrite partition files it they already exist (if experiment changes, directories need to change too)

# constants
seed <- 9
nFolds <- 10
cvReps <- 3
metric <- 'rocAUC'
lasso <- TRUE

# database connection
group <- 'agile' # group name in MySQL configuration file
patTabName <- 'patient'
eventTabName <- 'agent6'

# data output
opt$dataFile <- file.path(opt$dataFolder,
                          switch(opt$tType, 
                                 p = 'featureVectors_percentChange.rds', 
                                 c ='featureVectors_continuous.rds', 
                                 r ='featureVectors_responseCriteria.rds', 
                                 rate='featureVectors_rateChange.rds'))

# outputDir is the main folder holding spm data and logit experiments
# subDir the name of the logit experiment folder within outputDir/logits
logitFolder <- file.path(opt$dir, opt$prefix, opt$dataFolder)
ifelse(!dir.exists(logitFolder), dir.create(logitFolder, recursive = T), F) 

dataPath <- file.path(opt$dir,'spm',opt$dataFile) # spm data folder

dataPartitions <- file.path(opt$dir, 'train_test_partitions.rds') # expected ids of train and test, if doesn't exist, will create it

#---------------------------------------- log data ---------------------------------------- 

logfn <- file(file.path(logitFolder,
                        paste0(format(Sys.Date(), format="%Y.%m.%d"),
                               '_logit_',getVolTypeName(opt$tType),'.log')), 
              open='wt')

sink(logfn, type='output', split=T)

cat(format(Sys.Date(), format="%Y.%m.%d"), '\n')
cat('pattern dataset:', opt$dataFile, '\n')
cat('tumor vol variable type: ', getVolTypeName(opt$tType), '\n')
cat('seed: ', seed, '\n')
cat('nFolds: ', nFolds, '\n')
cat('cvReps: ', cvReps, '\n')

if (lasso) {
  cat('lambda length: ', opt$llength, '\n\n')
  cat('lambda max: ', opt$lmax, '\n\n')
}


#---------------------------------------- get data  partitions ---------------------------------------- 

#' create data partitions once
#' 
#' get partitions based on original dataset used for SPM generation,
#' so that train and test clinical visits remain the exact same across different SPM parameters

if (file.exists(dataPartitions)) {
  cat('\n\n ... using data partition file found at: ', dataPartitions, '\n\n')
} else {
  cat('\n\n ... data partition file not found.. will create now: ')
  data <- getData(database=group, eventTable=eventTabName, personTable=patTabName) # event data
  demo <- getDemographics(database=group, personTable=patTabName) # includes MGMT biomarker
  
  # 'c'  tType does not remove initial imaging volumes 
  # (e.g., rate change needs two visits, the very first one is ignores)
  data$event <- cleanData(data$event, tType = 'c') # pool of visits considered
  
  data <- merge(data$event, data$person, by='iois', all.x=T)
  data <- prepDemographics(data, demo) # need gender info
  data <- prepSurvivalLabels(data) # get survival labels, these also change
  data$id <- paste0(data$iois,'.',data$eventID) # clinical ids
  
  train.ids <- getTrainTestPartition(data=data,
                                     database=group, 
                                     personTable=patTabName,
                                     seed=seed,
                                     verbose=T)
  saveRDS(train.ids, file=dataPartitions)
}

#---------------------------------------- get event data   ---------------------------------------- 

data <- readRDS(dataPath) # features and labels for each clinical visit for a SPM parameterization
names <- colnames(data)
data <- data.frame(id=row.names(data),data)
colnames(data) <- c('id',names)
cat('...overall samples: ', nrow(data), '\n')

data <- prepLaterality(data) 
data <- prepLocation(data)

cat('...removing clinical visits with no history of length maxgap:', opt$maxgap, '\n')
data <- removeVisits(data, 
                     maxgap=opt$maxgap, 
                     maxlength=opt$maxlength, 
                     tType=opt$tType, 
                     save=F, 
                     outDir=logitFolder)

labels <- getClassLabels()

# metadata, unwanted features, and 
# unwanted labels should not be in training data for glmnet
# only one label will be kept in this set, depending on which model is being fit
#
needToRemove <- c('id','iois','eventID', # remove ids
                  labels, # remove labels
                  'IDH1') # not interested

if (opt$saveLogit=='yes') {
  save = T
} else {
  save = F
}



#---------------------------------------- fit data  ---------------------------------------- 

# data partitions
partitions <- readRDS(file=dataPartitions)
data <- data[data$id %in% partitions$train,] #training data

start.total <- Sys.time()

cat('\n\n','Running lasso logit...\n\n')

results <- lapply(labels, function(ln) {
  cat('\n\n\n....', ln ,' ...\n\n'); start.local <- Sys.time()
  form <- as.formula(paste0(ln,'~.'))
  aucs <- runCV(data=data,
               formula=form,
               lasso = TRUE,
               llength=opt$llength,
               lmax=opt$lmax,
               labelName=ln,
               needToRemove=needToRemove,
               metric=metric,
               seed=seed,
               folds=nFolds,
               cvReps=cvReps,
               createModelMatrix = FALSE)
  cat('\n\n'); print(Sys.time() - start.local)
  aucs
})
names(results) <- c('fit2m','fit6m','fit9m','fit12m')

if (opt$saveLogit=='yes') {
  cat('saving logit results...')
  saveRDS(results, 
          file=file.path(logitFolder,paste0('logits_CVresults_',getVolTypeName(opt$tType),'.rds')))
}

cat('\n\n', 'Total CV training time: '); print(Sys.time() - start.total)
