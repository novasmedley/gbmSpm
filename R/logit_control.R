# logit modeling, prediction, and plotting functions

# includes extracting glmnet model for information (e.g., coefs, avg ROCs)

#---------------------------------------- general logit modeling ----------------------------------------
# modification of caret's defaultSummary for summaryFunction


#' Sets the performance metric to be area under the AUC and PR curve
#' @details Modifies \code{summaryFunction} used in \code{\link[caret]{trainControl}}. 
#' Creates additional summary functions for \code{summaryFunction}. 
#' 
#' Possible methods for evaluating cross-validation performance in \code{\link[caret]{trainControl}} 
#' can then be \code{prAUC} or \code{rocAUC} for area under the precision-recall curve 
#' and area under the ROC curve, respectively. 
#' 
#' Their performance names are found in \code{\link{getPerformanceNames}}.
#' 
#' @return named vector of two functions: \code{prAUC} or \code{rocAUC}
#' @export 
#' @examples 
#' # not run
#' trainControl(method='cv', classProbs=T, index = index, summaryFunction=aucs)

aucs <- function(data, lev = NULL, model = NULL) {
  # modification of caret's defaultSummary for summaryFunction
  
  # print(data)
  # print(lev)
  # print(lev[1])
  # print(lev[2])
  prAUC <- MLmetrics::PRAUC(y_pred = data$dead, y_true = data$obs)
  lvls <- levels(data$obs)
  if (length(lvls) > 2) 
    stop(paste("Your outcome has", length(lvls), "levels. The twoClassSummary() function isn't appropriate."))
  if (!all(levels(data[, "pred"]) == lvls)) 
    stop("levels of observed and predicted data do not match")
  
  data$y = as.numeric(data$obs == lvls[2])
  rocAUC <- ModelMetrics::auc(ifelse(data$obs == lev[2], 0, 1), data[, lvls[1]])
  out <- c(prAUC,rocAUC)
  names(out) <- c("prAUC","rocAUC")
  
  out
}


#' Get colnames that will be kept for averaging across cross validation.
#' 
#' @details A results dataframe is returned from caret. These are the column
#' names that are numeric and can be averaged across multiple results.
#' @seealso \code{\link{aucs}}
#' @return char
#' 
#' 
#' @export
getPerformanceNames <- function(){
  c("prAUC","rocAUC", "prAUCSD","rocAUCSD")
}

#'Fit logit model using cross-validation
#'
#'@details Trains logits using \code{\link{glmnet}} and \code{\link{caret}}
#'
#'@param trainData dataframe, rows are samples to be classified, columns are features including sample ids
#'@param formula char or formula object
#'@param lasso logical, whether to use lasso regularization
#'@param llength num, number of lambdas to consider up to \code{lmax}
#'@param lmax num, maximum lambda to consider, cannot be NULL if lambda is NULL
#'@param lambda num, consider one lambda, \code{lmax} should be NULL if this has a value
#'@param metric char, see \code{aucs}
#'@param seed int, seed for split
#'@param folds int, number of folds for cross-validation
#'@param index list of ints, each element is an array, see \code{\link[caret:trainControl]{index}}
#'@param verbose 
#'@export
#' @examples
#' # use training partition to create folds for CV
#' data("features_ratechange_sup0.4g60l2z2") # features and labels for each clinical visit
#' t <- 'rate'
#' maxgap <- 60
#' maxlen <- 2
#' 
#' # format
#' names <- colnames(feats)
#' feats <- data.frame(id=row.names(feats),feats)
#' colnames(feats) <- c('id',names)
#' feats <- prepLaterality(feats)
#' feats <- prepLocation(feats)
#' feats <- removeVisits(feats,
#'                      maxgap=maxgap,
#'                      maxlength=maxlen,
#'                      tType=t,
#'                      save=F,
#'                      outDir=NA)
#' labels <- getClassLabels()
#' needToRemove <- c('id','iois','eventID', # remove ids
#'                   labels, # remove labels
#'                   'IDH1') # not interested
#' 
#' # data partitions
#' train.ids <- sample(feats$id, size=floor(0.80*nrow(feats)), replace = F) # random
#' feats <- feats[feats$id %in% train.ids,] #training data
#' ind <- getTrainingFolds(trainEvents=feats,
#'                         folds=3,
#'                         seed=1,
#'                         verbose=T)
#' feats <- prepLogitData(data = feats,
#'                       formula = 'survivalIn60 ~ .',
#'                       labelName = 'survivalIn60',
#'                       needToRemove=needToRemove)
#' 
#' # logit
#' fit <- getFit(formula = 'survivalIn60 ~ .',
#'               trainData = feats,
#'               index = ind,
#'               lasso = TRUE,
#'               llength = 100,
#'               lmax = 0.2,
#'               metric = 'rocAUC',
#'               seed = 1,
#'               verbose = T)
#'               
getFit <- function(formula, trainData, index, metric=c("prAUC","rocAUC"), cv=TRUE,
                   lasso=F, llength=NULL, lmax=NULL, seed, lambda=NULL, verbose=TRUE) {
  if (cv) {
    ctrl <- trainControl(method='cv',classProbs=T, index = index,
                         sampling='down', # TO DO: use weights
                         summaryFunction=aucs, savePredictions='final',
                         verboseIter=verbose, returnData = T, returnResamp = 'final')
    set.seed(seed) 
    ctrl$seeds <- sample.int(1e4,(length(index)*1)+1)
  } else {
    ctrl <- trainControl(method='none',classProbs=T,
                         sampling='down',
                         summaryFunction=aucs, savePredictions='final',
                         verboseIter=verbose, returnData = T, returnResamp = 'final')
    ctrl$seeds <- seed
  }
    
  if (lasso) {
    # lasso regularization (feature selection)
    if (!is.null(lmax)) {
      # lasso with grid search for lambda
      lambda.seq <- exp(seq(log(.0001*lmax),log(lmax),length.out=llength))
      
       train(as.formula(formula), 
              data=trainData,
              method ="glmnet", 
              family="binomial", 
              tuneGrid=expand.grid(alpha=1,lambda=lambda.seq),
              na.action=na.omit, 
              trControl=ctrl, 
              metric=metric)
      
    } else if (!is.null(lambda)) {
      # lasso without grid search, lambda is set
      train(as.formula(formula), 
            data=trainData,
            method ="glmnet", 
            family="binomial", 
            tuneGrid=expand.grid(alpha=1,lambda=lambda),
            na.action=na.omit, 
            trControl=ctrl, 
            metric=metric)
    }
  } else {
    # no lasso (use all features)
    train(as.formula(formula), 
          data=trainData,
          method ="glm", 
          family="binomial",
          na.action=na.omit,
          trControl=ctrl,
          metric=metric)
  }
}


#' Remove unwanted or keep interested columns from logit data
#' 
#' @details 
#' 
#' @param data dataframe, rows are samples, cols are features plus some metadata
#' not meant for modeling and will be removed
#' @param predictors char, names of columns in \code{data} that should be in logit fit data
#' @param needToRemove char, names of columns in \code{data} that should not be in logit fit data
#' @param labelName char, column name of binary label
#' 
#' @return \code{data} with removed columns, ready for logit modeling
#' @export
#' @examples
#' data("features_ratechange_sup0.4g60l2z2") # features and labels for each clinical visit
#' t <- 'rate'
#' maxgap <- 60
#' maxlen <- 2
#' 
#' # format
#' names <- colnames(feats)
#' feats <- data.frame(id=row.names(feats),feats)
#' colnames(feats) <- c('id',names)
#' feats <- prepLaterality(feats)
#' feats <- prepLocation(feats)
#' feats <- removeVisits(feats,
#'                      maxgap=maxgap,
#'                      maxlength=maxlen,
#'                      tType=t,
#'                      save=F,
#'                      outDir=NA)
#' labels <- getClassLabels()
#' needToRemove <- c('id','iois','eventID', # remove ids
#'                   labels, # remove labels
#'                   'IDH1') # not interested
#' 
#' # data partitions
#' train.ids <- sample(feats$id, size=floor(0.80*nrow(feats)), replace = F) # random
#' feats <- feats[feats$id %in% train.ids,] #training data
#' ind <- getTrainingFolds(trainEvents=feats,
#'                         folds=3,
#'                         seed=1,
#'                         verbose=T)
#' feats <- setDataCols(data=feats, 
#'                     labelName='survivalIn60', 
#'                     needToRemove=needToRemove)
#'                     
setDataCols <- function(data, labelName, predictors=NULL, needToRemove=NULL) {
  if (!(is.null(predictors))){
    data = data[ , c(predictors,labelName)]
  } else if (!(is.null(needToRemove))) {
    remove <- needToRemove[-which(needToRemove %in% labelName)]
    data = data[ , !(colnames(data) %in% remove)]
  } else {
    stop('both predictors and needToRemove cannot be NULL')
  }
}

#' Create data for logit model.
#' 
#' 
#' @details Removes non-features and non-labels. 
#' Creates dataframe ready to serve as input to logit fit, 
#' see \code{model.matrix}. Removes rows with any 'NA' values.
#' 
#' @inheritParams setDataCols
#' @param formula char or formula object
#' @param createModelMatrix logical, call \code{\link[stats]{model.matrix}}
#' 
#' @return dataframe ready for logit fitting
#' @export
#' 
#' @examples
#' 
#' # use training partition to create folds for CV
#' data("features_ratechange_sup0.4g60l2z2") # features and labels for each clinical visit
#' t <- 'rate'
#' maxgap <- 60
#' maxlen <- 2
#' 
#' # format
#' names <- colnames(feats)
#' feats <- data.frame(id=row.names(feats),feats)
#' colnames(feats) <- c('id',names)
#' feats <- prepLaterality(feats)
#' feats <- prepLocation(feats)
#' feats <- removeVisits(feats,
#'                      maxgap=maxgap,
#'                      maxlength=maxlen,
#'                      tType=t,
#'                      save=F,
#'                      outDir=NA)
#' labels <- getClassLabels()
#' needToRemove <- c('id','iois','eventID', # remove ids
#'                   labels, # remove labels
#'                   'IDH1') # not interested
#' 
#' # data partitions
#' train.ids <- sample(feats$id, size=floor(0.80*nrow(feats)), replace = F) # random
#' feats <- feats[feats$id %in% train.ids,] #training data
#' ind <- getTrainingFolds(trainEvents=feats,
#'                         folds=3,
#'                         seed=1,
#'                         verbose=T)
#' feats <- prepLogitData(data = feats,
#'                       formula = 'survivalIn60 ~ .',
#'                       labelName = 'survivalIn60',
#'                       needToRemove=needToRemove)
#'                       
prepLogitData <- function(data, formula, labelName, 
                          predictors=NULL, needToRemove=NULL,
                          createModelMatrix=FALSE) {
  data <- setDataCols(data=data, 
                      labelName=labelName, 
                      predictors=predictors,
                      needToRemove=needToRemove)
  if (createModelMatrix==T) {
    labels <- data[,labelName]
    data <- model.matrix(formula, data) # dummy vars and intercept
    data <- as.data.frame(data)
    data <- cbind(labels, data) # tack labels back on for modeling
    colnames(data)[1] <- labelName # ensure label is named correctly for modeling
  }
  return(data)
}

#' Perform cross-validation on training partition.
#' 
#' @details Removes rows for which a column is NA.
#' 
#' @inherit getTrainingFolds
#' @inherit prepLogitData
#' @inherit getFit
#' 
#' @export
getCV <- function(data, formula, labelName, 
                  lasso=FALSE, llength=NULL, lmax=NULL, 
                  predictors=NULL, needToRemove=NULL, createModelMatrix=FALSE,
                  metric=c("prAUC","rocAUC"),
                  seed, folds, verbose=TRUE) {
  data <- na.omit(data)
  ind <- getTrainingFolds(trainEvents=data,
                          folds=folds,
                          seed=seed,
                          verbose=verbose)
  data <- prepLogitData(data = data,
                        formula = formula,
                        labelName =labelName,
                        predictors = predictors,
                        needToRemove = needToRemove,
                        createModelMatrix = createModelMatrix)
  fit <- getFit(formula = formula,
                trainData = data,
                index = ind,
                lasso = lasso,
                llength =llength,
                lmax = lmax,
                metric = metric,
                seed = seed,
                verbose = verbose)
  return(fit)
}

#' Run logit modeling
#' 
#' @description Patients are first split into training and testing partitions. 
#' Next, samples with NA features will be removed. 
#' Then, training partition is split for cross-validation so that
#' no patient has events in both validation and training. 
#' Both partitions are prepped for logit fitting and cross-validation is run.
#' 
#' @inheritParams prepLogitData 
#' @inheritParams getCV

#' @seealso 
#' \code{\link{prepLogitData}}, 
#' \code{\link{getCV}}, 
#' \code{\link{getPerformanceNames}}
#' 
#' @return list of logit fits, averaged cross-validation results (if \code{cvReps} > 1), and test data partition 
#' (one that is formated like training data and is ready to be called by subsequent trained model)
#' 
#' @export
runCV <- function(data, cvReps, 
                  formula, labelName,
                  lasso=FALSE, llength=NULL, lmax=NULL, 
                  predictors=NULL, needToRemove=NULL, createModelMatrix=FALSE,
                  metric=c("prAUC","rocAUC"), seed, folds) {

  if (cvReps==1) {
    cat('\n\n....',labelName, 'CV with seed: ',seed,'\n\n')
    fit <- getCV(data = data, 
                 formula= formula, 
                 labelName=labelName, 
                 lasso = lasso, 
                 llength=llength,  
                 lmax=lmax, 
                 predictors=predictors, 
                 needToRemove=needToRemove, 
                 createModelMatrix=createModelMatrix,
                 metric=metric, 
                 folds=folds,
                 seed=seed,
                 verbose=TRUE)
    print(fit$results[ which.max(fit$results[, metric]),])
    
    r <- (list(fit=fit)) # only 1 fitting

  } else if (cvReps > 1) {
    # get repeated cross validation for training partition
    cat('\n\n....run ',cvReps,'CV repeats:\n\n')
    seeds <- seq(1,cvReps)*seed 
    results <- lapply(seeds, function(s) {
      cat('\n\n....',labelName, 'CV with seed: ',s,'\n\n')
      fit <- getCV(data = data, 
                   formula= formula, 
                   labelName=labelName, 
                   lasso = lasso, 
                   llength=llength,  
                   lmax=lmax, 
                   predictors=predictors, 
                   needToRemove=needToRemove, 
                   createModelMatrix=createModelMatrix,
                   metric=metric, 
                   folds=folds,
                   seed=s,
                   verbose=FALSE)
      fit$results
    })
    
    # unravel repeated CV results to get a sense of the averaged results
    avgPerf <- lapply(results, function(res) res[,getPerformanceNames()])
    avgPerf <- Reduce("+", avgPerf) / length(avgPerf) # average over repeated CV
    avgPerf$lambda <- results[[1]]$lambda # attach lambdas

    cat('\n\n....',labelName, 'highest averaged results over repeated CV: \n\n')
    print(avgPerf[which.max(avgPerf[,metric]),])

     r<- list(avg.CVPerformance=avgPerf)
  }
  
  return(r)
}

#' Retrain logit after hyperparameter selection.
#' 
#' @description Used for fitting the full training partition.
#' 
#' @inheritParams prepLogitData 
#' @inheritParams getCV
#' @param lambda numberic, lambda value
#' 
#' @export
retrainLogit <- function(data, formula, labelName,
                         lambda=NULL, lasso=TRUE,
                         predictors=NULL, needToRemove=NULL, createModelMatrix=FALSE,
                         metric=c("prAUC","rocAUC"), seed, verbose=TRUE) {
  
  data <- na.omit(data)
  data <- prepLogitData(data = data,
                        formula = formula,
                        labelName =labelName,
                        predictors = predictors,
                        needToRemove = needToRemove,
                        createModelMatrix = createModelMatrix)
  fit <- getFit(formula = formula,
                trainData = data,
                metric = metric,
                lambda = lambda,
                cv=FALSE,
                index = NULL,
                lasso = lasso,
                llength = NULL,
                lmax = NULL,
                seed = seed,
                verbose = verbose)
  return(fit)
}


#---------------------------------------- extracting information from glmnet ---------------------------------------- 
#' Get logit coefficients
#' 
#' @description Assumes \code{fit} contains the following (returned from \code{\link[caret]{trainControl}}):
#' \itemize{
#' \item \code{traininData} data used for training, with \code{.outcome} as the classification label column name
#' \item \code{fintalMode} the retrained, final model using selected \code{lamba} from cross-validation
#' \item \code{bestTune} values of selected hyperparameter during cross-validation (i.e, \code{labmda} value)
#' }
#' 
#' @param fit logit fit
#' @export
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


#---------------------------------------- extracting single model results ---------------------------------------- 
#' get training results of fitted model retrained on entire training partition
#' @description  asdf
#' @param model glmnet fitted
#' @export
getTrainResults <- function(fit) {
  trainProb <- predict(fit, newData=fit$trainingData, type='prob')
  trainPred <- ROCR::prediction(trainProb$dead, fit$trainingData$.outcome)
  overall_perf <- getROC_PR(trainPred)
  threshold_perf <- as.data.frame(selectROCThreshold(trainPred))
  return(list(preds=trainPred, performance=overall_perf, selectedThreshold=threshold_perf))
}

#' Apply 
#' 
#' @param fit glmnet fitted model
#' @param data testing data, which will under the same preprocessing as training, see \code{prepLogitData}
#' @inheritParams prepLogitData
#' 
#' @return named list of test predictions, ROC and PR curve performances, and the original data
#' 
#' @export
getTestResults <- function(fit, data, formula, labelName,
                           predictors=NULL, needToRemove=NULL, createModelMatrix=FALSE) {
  
  data <- na.omit(data)
  data <- prepLogitData(data = data,
                        formula = formula,
                        labelName =labelName,
                        predictors = predictors,
                        needToRemove = needToRemove,
                        createModelMatrix = createModelMatrix)
  
  labelInd <- which(colnames(data)==labelName)
  colnames(data)[labelInd] <- '.outcome'
  
  testProb <- predict(fit, newdata=data, type='prob')
  testPred <- ROCR::prediction(testProb$dead, data$.outcome )
  overall_perf <- getROC_PR(testPred)

  return(list(preds=testPred, performance=overall_perf, data=data))
}
getBestResults <- function(results, models, month, labelName) {
  bestModel <- getBestModel(month,results)
  train <- getTrainResults(models[[bestModel]])
  test <- getTestResults(models[[bestModel]], labelName)
  
  return(list(bestModel=bestModel, train=train, test=test))
}