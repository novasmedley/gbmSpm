% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_discretize.R
\name{setContinDiscrete}
\alias{setContinDiscrete}
\alias{setPercentDiscrete}
\alias{setRateDiscrete}
\alias{setBaselineDiscrete}
\title{Discretize tumor volumes into intervals, where each interval is a category (an event).}
\usage{
setContinDiscrete(volumes, categories = 10)

setPercentDiscrete(volumes, categories = 10)

setRateDiscrete(volumes, categories = 10)

setBaselineDiscrete(volumes, categories = 10)
}
\arguments{
\item{volumes}{a dataframe, rows are tumor volumes}

\item{categories}{ number of categories or a vector with boundaries 
    (all values outside the boundaries will be set to NA). }
}
\value{
\code{volumes} with \code{eventName} as an additional column
}
\description{
\code{setContinDiscrete}: convert volumes, measured in original units, to intervals.

\code{setPercentDiscrete}: convert volumes, measured in percent change, to intervals.

\code{setRateDiscrete}: convert volumes, measured in rate change, to intervals.

\code{setBaselineDiscrete}: convert baseline volumes, measured in original units, to intervals.
}
\details{
Volumes are discretized into 10 intervals, where each is an event. Event names are named after the intervals.
Note: \itemize{
\item intervals are based on distribution of values in \code{volumes}
\item each interval or bin will be used as different categories of tumor volume
\item also see \code{prepTumorVol}
}

\code{volumes} contains tumor measurement events as rows and at least the following columns: \itemize{
\item agent_detail = event attribute, i.e., measured tumor volume in mm3 char
}
}
\examples{
data("fake_data")
vols <- fake_data$events[fake_data$events$agent=='measurement', ]

# format data
vols$startdate <- as.Date(as.character(vols$startdate), format='\%Y-\%m-\%d')
vols$enddate <- as.Date(as.character(vols$enddate), format='\%Y-\%m-\%d')
vols <- resetEventIDsToDays(vols)
vols <- removeData(vols) # get rid of excess details and empty or no data

vols_dist <- setContinDiscrete(vols)

data("fake_data")
vols <- fake_data$events[fake_data$events$agent=='measurement', ]

# format data
vols$startdate <- as.Date(as.character(vols$startdate), format='\%Y-\%m-\%d')
vols$enddate <- as.Date(as.character(vols$enddate), format='\%Y-\%m-\%d')
vols <- resetEventIDsToDays(vols)
vols <- removeData(vols) # get rid of excess details and empty or no data

# convert vol to correct type for all patients
vols <- split(vols, vols$iois) # list of dataframes, where each is a patient's set of tumor volumes
vols <- lapply(vols, function(v) setPercentChange(v))
vols <- do.call('rbind', vols)

# discretize
vols_dist <- setPercentDiscrete(vols)

}
\seealso{
\code{\link[arules]{discretize}}
}
