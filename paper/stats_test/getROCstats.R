library(pROC)

# apply retrained models onto the training partition for comparison

months <- c('2-month','6-month','9-month', '12-month')
dataDir <- '~/Data/Agile'
spmLogitDir <- file.path(dataDir, 'final_model') # temporal pattern approach
volLogitDir <- file.path(dataDir, 'analysis_vol', 'logits') # tumor vol approaches
volThresDir <- file.path(dataDir, 'analysis_vol', 'threshold') # tumor vol approaches

#------------- log ------------- 

logfn <- file(file.path(spmLogitDir,
                        paste0(format(Sys.Date(), format="%Y.%m.%d"),
                               '_trainingROC_pValues','.log')), 
              open='wt')

sink(logfn, type='output', split=T)

cat(format(Sys.Date(), format="%Y.%m.%d"), '\n')
cat('retrained spm logit file: ', file.path(spmLogitDir, 'retrain_test_results.rds'), '\n')
cat('retrained vol logit file: ', file.path(volLogitDir,'tumorVol_logits_CV_results.rds') , '\n')
cat('vol thresholidng file: ', file.path(volThresDir,'tumorVol_threshold_results.rds') , '\n')

#------------- spm logit results ------------- 
retrainResults <- readRDS(file.path(spmLogitDir, 'retrain_test_results.rds')) 

# are glmnet models
length(unlist(retrainResults$fit2m$train$preds@labels))
length(unlist(retrainResults$fit6m$train$preds@labels))
length(unlist(retrainResults$fit9m$train$preds@labels))
length(unlist(retrainResults$fit12m$train$preds@labels))

length(unlist(retrainResults$fit2m$test$preds@labels))
length(unlist(retrainResults$fit6m$test$preds@labels))
length(unlist(retrainResults$fit9m$test$preds@labels))
length(unlist(retrainResults$fit12m$test$preds@labels))

spm2mroc <- roc(response = unlist(retrainResults$fit2m$train$preds@labels),  predictor = unlist(retrainResults$fit2m$train$preds@predictions) )
spm6mroc <- roc(response = unlist(retrainResults$fit6m$train$preds@labels),  predictor = unlist(retrainResults$fit6m$train$preds@predictions) )
spm9mroc <- roc(response = unlist(retrainResults$fit9m$train$preds@labels),  predictor = unlist(retrainResults$fit9m$train$preds@predictions) )
spm12mroc<- roc(response = unlist(retrainResults$fit12m$train$preds@labels), predictor = unlist(retrainResults$fit12m$train$preds@predictions) )

spm2mroc_test <- roc(response = unlist(retrainResults$fit2m$test$preds@labels),  predictor = unlist(retrainResults$fit2m$test$preds@predictions) )
spm6mroc_test <- roc(response = unlist(retrainResults$fit6m$test$preds@labels),  predictor = unlist(retrainResults$fit6m$test$preds@predictions) )
spm9mroc_test <- roc(response = unlist(retrainResults$fit9m$test$preds@labels),  predictor = unlist(retrainResults$fit9m$test$preds@predictions) )
spm12mroc_test<- roc(response = unlist(retrainResults$fit12m$test$preds@labels), predictor = unlist(retrainResults$fit12m$test$preds@predictions) )

#------------- read in tumor logit ------------- 
volLogitResults <- readRDS(file.path(volLogitDir,'tumorVol_logits_CV_results.rds'))

#' best vol logit
#' 2m: continuous (no discretization, no covariates)
#' 6m: continuous (no discretization, no covariates)
#' 9m: discreteRateChange (discretized, no covariates)
#' 12m: adj.d.discreteRateChange (discretized, with covariates)

model2m <- readRDS(file.path(volLogitDir, paste0('tumorVol_logits_trainResults_','survivalIn60', '_', 'continuous','.rds')))
model6m <- readRDS(file.path(volLogitDir, paste0('tumorVol_logits_trainResults_','survivalIn180', '_', 'continuous','.rds')))
model9m <- readRDS(file.path(volLogitDir, paste0('tumorVol_logits_trainResults_','survivalIn270', '_', 'discreteRateChange','.rds')))
model12m <- readRDS(file.path(volLogitDir, paste0('tumorVol_logits_trainResults_','survivalIn360', '_', 'adj.d.discreteRateChange','.rds')))

fixVarNames <- function(trainingData) {
  colnames(trainingData) <- gsub('discrete', '`discrete', colnames(trainingData))
  colnames(trainingData) <- gsub('age', '`age', colnames(trainingData))
  colnames(trainingData) <- gsub(')', ')`', colnames(trainingData)) 
  colnames(trainingData) <- gsub(']', ']`', colnames(trainingData))
  
  colnames(trainingData) <- gsub('`location', '`\\\\`location', colnames(trainingData)) 
  colnames(trainingData) <- gsub('`laterality', '`\\\\`laterality', colnames(trainingData)) 
  colnames(trainingData) <- gsub('`TRUE', '\\\\`TRUE`', colnames(trainingData)) 
  
  colnames(trainingData) <- gsub('(Intercept)`', '(Intercept)', colnames(trainingData)) # convert back intercept 
  return(trainingData)
}

# these are glm models
tvlogit2mroc  <- roc(response = model2m$trainingData$.outcome,  predictor = predict.glm(model2m$finalModel,  newdata=fixVarNames(model2m$trainingData), type='response')) 
tvlogit6mroc  <- roc(response = model6m$trainingData$.outcome,  predictor = predict.glm(model6m$finalModel,  newdata=fixVarNames(model6m$trainingData), type='response'))
tvlogit9mroc  <- roc(response = model9m$trainingData$.outcome,  predictor = predict.glm(model9m$finalModel,  newdata=fixVarNames(model9m$trainingData), type='response'))
tvlogit12mroc <- roc(response = model12m$trainingData$.outcome, predictor = predict.glm(model12m$finalModel, newdata=fixVarNames(model12m$trainingData), type='response'))



#-------------  get tumor vol no model ------------- 
volResults <- readRDS(file.path(volThresDir,'tumorVol_threshold_results.rds'))

#' best vol only: all via continuous
#' 
tv2mroc <- roc(response = unlist(volResults$survivalIn60.continuous$trainProb@labels),  predictor = unlist(volResults$survivalIn60.continuous$trainProb@predictions))
tv6mroc <- roc(response = unlist(volResults$survivalIn180.continuous$trainProb@labels), predictor = unlist(volResults$survivalIn180.continuous$trainProb@predictions))
tv9mroc <- roc(response = unlist(volResults$survivalIn270.continuous$trainProb@labels), predictor = unlist(volResults$survivalIn270.continuous$trainProb@predictions))
tv12mroc <-roc(response = unlist(volResults$survivalIn360.continuous$trainProb@labels), predictor = unlist(volResults$survivalIn360.continuous$trainProb@predictions))


#-------------  get tests ------------- 

# 2-month
cat('\n\n ... 2-month results: ')
print('2-month: spm vs tv logit'); roc.test(spm2mroc,tvlogit2mroc, method='delong'); cat('\n')
print('2-month: spm vs tv only'); roc.test(spm2mroc, tv2mroc, method='delong') ; cat('\n')

# 6-month
cat('\n\n ... 6-month results: ')
print('6-month: spm vs tv logit'); roc.test(spm6mroc,tvlogit6mroc, method='delong'); cat('\n')
print('6-month: spm vs tv only'); roc.test(spm6mroc, tv6mroc, method='delong') ; cat('\n')

# 9-month
cat('\n\n ... 9-month results: ')
print('9-month: spm vs tv logit'); roc.test(spm9mroc,tvlogit9mroc, method='delong'); cat('\n')
print('9-month: spm vs tv only'); roc.test(spm9mroc, tv9mroc, method='delong') ; cat('\n')

# 12-month
cat('\n\n ... 12-month results: ')
print('12-month: spm vs tv logit'); roc.test(spm12mroc,tvlogit12mroc, method='delong'); cat('\n')
print('12-month: spm vs tv only'); roc.test(spm12mroc, tv12mroc, method='delong') ; cat('\n')


# training confidence intervals
cat('\n\n ... training partition confidence intervals: ')
cat('\n... spm logit 2m \n'); ci.auc(spm2mroc, method='bootstrap', progress = 'none') ; cat('\n')
cat('\n... spm logit 6m \n'); ci.auc(spm6mroc, method='bootstrap', progress = 'none') ; cat('\n')
cat('\n... spm logit 9m \n'); ci.auc(spm9mroc, method='bootstrap', progress = 'none') ; cat('\n')
cat('\n... spm logit 12m \n');ci.auc(spm12mroc, method='bootstrap', progress = 'none') ; cat('\n')

cat('\n... vol logit 2m \n'); ci.auc(tvlogit2mroc, method='bootstrap', progress = 'none'); cat('\n')
cat('\n... vol logit 6m \n'); ci.auc(tvlogit6mroc, method='bootstrap', progress = 'none'); cat('\n')
cat('\n... vol logit 9m \n'); ci.auc(tvlogit9mroc, method='bootstrap', progress = 'none'); cat('\n')
cat('\n... vol logit 12m \n');ci.auc(tvlogit12mroc, method='bootstrap', progress = 'none'); cat('\n')

cat('\n... vol 2m \n'); ci.auc(tv2mroc, method='bootstrap', progress = 'none'); cat('\n')
cat('\n... vol 6m \n'); ci.auc(tv6mroc, method='bootstrap', progress = 'none'); cat('\n')
cat('\n... vol 9m \n'); ci.auc(tv9mroc, method='bootstrap', progress = 'none'); cat('\n')
cat('\n... vol 12m \n'); ci.auc(tv12mroc, method='bootstrap', progress = 'none'); cat('\n')

# testing confidence intervals of selected model
cat('\n\n ... testing partition confidence intervals from selected model (spm logit): ')
cat('\n... spm logit 2m \n'); ci.auc(spm2mroc_test, method='bootstrap', progress = 'none') ; cat('\n')
cat('\n... spm logit 6m \n'); ci.auc(spm6mroc_test, method='bootstrap', progress = 'none') ; cat('\n')
cat('\n... spm logit 9m \n'); ci.auc(spm9mroc_test, method='bootstrap', progress = 'none') ; cat('\n')
cat('\n... spm logit 12m \n');ci.auc(spm12mroc_test, method='bootstrap', progress = 'none') ; cat('\n')

close(logfn)
