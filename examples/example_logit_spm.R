# Load spm pattern data,perform logit analysis using spm events + demographics information 

# saves fit data for each survival prediction

# get help: in terminal, Rscript runSpmLogits.R -h

#-------------------------------------------------------------------------------- #
library(gbmSpm)
library(optparse)


#-------------------------------------------------------------------------------- #

# get cmd line arguments
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

# OR use run in RStudio using :

# opt <- list()
# opt$dir <- '~/gbm_spm_example'
# opt$tType <- 'rate'
# opt$dataFolder <- 'sup0.4g60l2z2' # previously created spm patterns
# opt$dataFile <- file.path(opt$dataFolder, 'featureVectors_rateChange.rds') # data input file
# opt$maxgap <- 60
# opt$maxlength <- 2
# 
# opt$prefix <- 'logits'
# opt$saveLogit <- 'yes'
# opt$llength <- 100
# opt$lmax <- 0.2


# constants
seed <- 9
nFolds <- 2
cvReps <- 2
metric <- 'rocAUC'
lasso <- TRUE

# output dirs
opt$dataFile <- file.path(opt$dataFolder,
                          switch(opt$tType, 
                                 p = 'featureVectors_percentChange.rds', 
                                 c ='featureVectors_continuous.rds', 
                                 r ='featureVectors_responseCriteria.rds', 
                                 rate='featureVectors_rateChange.rds'))
logitFolder <- file.path(opt$dir, opt$prefix, opt$dataFolder)
ifelse(!dir.exists(logitFolder), dir.create(logitFolder, recursive = T), F) 
dataPartitions <- file.path(opt$dir, 'train_test_partitions.rds') # expected ids of train and test, if doesn't exist, will create it

# input path
dataPath <- file.path(opt$dir,'spm',opt$dataFile) # spm data folder
cat('... expecting data file: ', dataPath, '\n\n')

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
  cat('lambda length: ', opt$llength, '\n')
  cat('lambda max: ', opt$lmax, '\n\n')
}


#---------------------------------------- get data  partitions ---------------------------------------- 
# data partitioning depends on the fully cleaned dataset
# and ensures the exact same clinical visits are maintained in training or testing, regardless of 
# feature engineering methods.

data("fake_data")

cat('\n\n ... data partition .. will create now: ')
# 'c'  tType does not remove initial imaging volumes 
# (e.g., rate change needs two visits, the very first one is ignores)
fake_demo <- fake_data$demo
survData <- fake_data$person
fake_data$events <- cleanData(fake_data$events, tType = 'c') # pool of visits considered

fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
fake_data <- prepDemographics(fake_data, fake_demo) # need gender info
fake_data <- prepSurvivalLabels(fake_data) # get survival labels, these also change
fake_data$id <- paste0(fake_data$iois,'.',fake_data$eventID) # clinical ids

part.ids <- getTrainTestPartition(data=fake_data,
                                   database=NULL, 
                                   personTable=NULL,
                                   survData=survData,
                                   seed=seed,
                                   verbose=T)

saveRDS(part.ids, file=dataPartitions)

#---------------------------------------- get event data   ---------------------------------------- 

# read in features generated from spm

data <- readRDS(dataPath) # features and labels for each clinical visit
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


# data partitions
data$id <- as.character(data$id)
data <- data[data$id %in% part.ids$train,] #training data

#---------------------------------------- fit data  ---------------------------------------- 

start.total <- Sys.time()

cat('\n\n','Running lasso logit...\n\n')

results <- lapply(labels, function(ln) {
  cat('\n\n\n....', ln ,' ...\n\n'); start.local <- Sys.time()
  form <- as.formula(paste0(ln,'~.'))
  aucs <- runCV(data = data,
                formula = form,
                lasso = TRUE,
                llength = opt$llength,
                lmax = opt$lmax,
                labelName = ln,
                needToRemove = needToRemove,
                metric = metric,
                seed = seed,
                folds = nFolds,
                cvReps = cvReps,
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