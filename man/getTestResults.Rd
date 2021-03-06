% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logit_control.R
\name{getTestResults}
\alias{getTestResults}
\title{Apply}
\usage{
getTestResults(fit, data, formula, labelName, predictors = NULL,
  needToRemove = NULL, createModelMatrix = FALSE)
}
\arguments{
\item{fit}{glmnet fitted model}

\item{data}{testing data, which will under the same preprocessing as training, see \code{prepLogitData}}

\item{formula}{char or formula object}

\item{labelName}{char, column name of binary label}

\item{predictors}{char, names of columns in \code{data} that should be in logit fit data}

\item{needToRemove}{char, names of columns in \code{data} that should not be in logit fit data}

\item{createModelMatrix}{logical, call \code{\link[stats]{model.matrix}}}
}
\value{
named list of test predictions, ROC and PR curve performances, and the original data
}
\description{
Apply
}
