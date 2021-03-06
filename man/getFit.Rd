% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logit_control.R
\name{getFit}
\alias{getFit}
\title{Fit logit model using cross-validation}
\usage{
getFit(formula, trainData, index, metric = c("prAUC", "rocAUC"), cv = TRUE,
  lasso = F, llength = NULL, lmax = NULL, seed, lambda = NULL,
  verbose = TRUE)
}
\arguments{
\item{formula}{char or formula object}

\item{trainData}{dataframe, rows are samples to be classified, columns are features including sample ids}

\item{index}{list of ints, each element is an array, see \code{\link[caret:trainControl]{index}}}

\item{metric}{char, see \code{aucs}}

\item{lasso}{logical, whether to use lasso regularization}

\item{llength}{num, number of lambdas to consider up to \code{lmax}}

\item{lmax}{num, maximum lambda to consider, cannot be NULL if lambda is NULL}

\item{seed}{int, seed for split}

\item{lambda}{num, consider one lambda, \code{lmax} should be NULL if this has a value}

\item{verbose}{}

\item{folds}{int, number of folds for cross-validation}
}
\description{
Fit logit model using cross-validation
}
\details{
Trains logits using \code{\link{glmnet}} and \code{\link{caret}}
}
\examples{
# use training partition to create folds for CV
data("features_ratechange_sup0.4g60l2z2") # features and labels for each clinical visit
t <- 'rate'
maxgap <- 60
maxlen <- 2

# format
names <- colnames(feats)
feats <- data.frame(id=row.names(feats),feats)
colnames(feats) <- c('id',names)
feats <- prepLaterality(feats)
feats <- prepLocation(feats)
feats <- removeVisits(feats,
                     maxgap=maxgap,
                     maxlength=maxlen,
                     tType=t,
                     save=F,
                     outDir=NA)
labels <- getClassLabels()
needToRemove <- c('id','iois','eventID', # remove ids
                  labels, # remove labels
                  'IDH1') # not interested

# data partitions
train.ids <- sample(feats$id, size=floor(0.80*nrow(feats)), replace = F) # random
feats <- feats[feats$id \%in\% train.ids,] #training data
ind <- getTrainingFolds(trainEvents=feats,
                        folds=3,
                        seed=1,
                        verbose=T)
feats <- prepLogitData(data = feats,
                      formula = 'survivalIn60 ~ .',
                      labelName = 'survivalIn60',
                      needToRemove=needToRemove)

# logit
fit <- getFit(formula = 'survivalIn60 ~ .',
              trainData = feats,
              index = ind,
              lasso = TRUE,
              llength = 100,
              lmax = 0.2,
              metric = 'rocAUC',
              seed = 1,
              verbose = T)
              
}
