% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spm_getPatterns.R
\name{vectorizePatterns}
\alias{vectorizePatterns}
\title{Create a binary vector of patterns observed at a clinical visit.}
\usage{
vectorizePatterns(patternDays)
}
\arguments{
\item{patternDays}{list of patterns found by cSPADE, 
for each pattern there is a list of patientIDs (see \code{\link{getPatIDs}}), 
for each patientID there is a list of eventIDs when the pattern was observed}
}
\value{
dataframe, 
rows of clinical visits, 
cols of binary features (0: observed, 1: not observed)
}
\description{
Create a binary vector of patterns observed at a clinical visit.
}
\details{
For every clinical visit, each pattern has a value of 
1 (present) or 0 (absent).
}
\examples{
data("patterns_ratechange_sup0.4g60l2z2") # days when pattern occurred
features <- vectorizePatterns(patternDays)
}
