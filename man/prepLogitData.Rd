% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logit_control.R
\name{prepLogitData}
\alias{prepLogitData}
\title{Create data for logit model.}
\usage{
prepLogitData(data, formula, labelName, predictors = NULL,
  needToRemove = NULL, createModelMatrix = FALSE)
}
\arguments{
\item{data}{dataframe, rows are samples, cols are features plus some metadata
not meant for modeling and will be removed}

\item{formula}{char or formula object}

\item{labelName}{char, column name of binary label}

\item{predictors}{char, names of columns in \code{data} that should be in logit fit data}

\item{needToRemove}{char, names of columns in \code{data} that should not be in logit fit data}

\item{createModelMatrix}{logical, call \code{\link[stats]{model.matrix}}}
}
\value{
dataframe ready for logit fitting
}
\description{
Create data for logit model.
}
\details{
Removes non-features and non-labels. 
Creates dataframe ready to serve as input to logit fit, 
see \code{model.matrix}. Removes rows with any 'NA' values.
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
                      
}
