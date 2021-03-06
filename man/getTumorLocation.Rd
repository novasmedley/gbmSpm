% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_preprocess.R
\name{getTumorLocation}
\alias{getTumorLocation}
\title{Get the tumor location and laterality.}
\usage{
getTumorLocation(data, tumorInfo)
}
\arguments{
\item{data}{dataframe, 
rows are clinical visits to be classified, 
columns are features including clinical visits ids 
and at least the following columns:
\itemize{
\item \code{id} char, clinical visit id in the format "patientID.eventID"
\item \code{iois} char, patient id
\item \code{agent_detail} char, event attribute
}}

\item{tumorInfo}{a dataframe of pre-cleaned events and searches for tumor location}
}
\value{
\code{data} with additional columns containing tumor location dummy variables and tumor laterality
}
\description{
Get the tumor location and laterality.
}
\details{
Takes the first event with tumor location information in \code{tumorInfo} and
convert location into dummy variables for logit since one person can have multiple locations.
\code{tumorInfo} contains pre-cleaned data (\code{\link{cleanData}}
is used for cleaning events for SPM) \code{agent_details}.
Also retrieves tumor lateriality.

Assumes location and laterility is in \code{agent_details} column.
Assumes location is after "Location :" in string.
Assumes location and lateriality is unchanged for the remaining events (is static)
}
\examples{

data("fake_data")
t <- 'rate'
fake_tumorInfo <- fake_data$events # save tumor location and laterility strings before event cleaning
fake_demo <- fake_data$demo
fake_data$events <- cleanData(fake_data$events, tType = t)

# collect patient info for each event
# note that only eventName, iois, and eventID columns are used for SPM
fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
fake_data <- prepDemographics(fake_data, fake_demo) # prep for each event, since age does change
fake_data <- prepSurvivalLabels(fake_data) # get survival labels, these also change

fake_data <- getTumorLocation(fake_data, fake_tumorInfo) # get first tumor location

}
