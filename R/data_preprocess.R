#
# Cleans and preps (e.g., discretizes continuous variables) data for SPM.
# Function expects input data from a database. See 'getData.R' for data types and format.
#

#---------------------------------- clean events for SPM -----------------------------------

#' Clean data to desired format ready for 'arulesSequences'.
#'
#' @details
#'
#' @param event: dataframe, events as rows and at least the following columns:
#' \itemize{
#' \item \code{iois} patient id, int
#' \item \code{agent} event type, char
#' \item \code{agent_detail} event attribute, char
#' \item \code{startdate} event start date, char will be formatted 'YYYY-MM-DD' Date
#' \item \code{enddate} event end date, char will be formatted 'YYYY-MM-DD' Date
#' }
#' @param tType: char, type of tumor data to use: see 'getVolTypeName' for valid values
#' @param removeChemo: bool, remove chemo event types
#'
#' @return event: dataframe, cleaned and formated information with added columns:
#' \itemize{
#' \item \code{eventID} event id, num, see \code{\link{resetEventIDsToDays}}
#' \item \code{endEventID} end event id, num, see \code{\link{resetEventIDsToDays}}
#' \item \code{eventName} unique event name, char, see \code{\link{createEventNames}}
#' }
#' @export
#' 
#' @examples
#' data("fake_data")
#' cleaned <- cleanData(fake_data$events, tType='r',  removeChemo=T)
#' 
cleanData <- function(event, tType='r', removeChemo=T) {
  cat('\n\n','Cleaning data...\n\n')
  
  event$startdate <- as.Date(as.character(event$startdate), format='%Y-%m-%d')
  event$enddate <- as.Date(as.character(event$enddate), format='%Y-%m-%d')
  
  # note order of operations here matters
  
  event <- removeData(event) # get rid of excess details and empty or no data
  event <- resetEventIDsToDays(event)# Set eventIDs relative to each person. EventIDs will be in days after the first observation date.
  
  event <- removePreRadiation(event) # removes events prior to the FIRST radiation treatment
  event$agent_detail <- unlist(event$agent_detail)
  
  event <- createEventNames(event, tType, removeChemo) # Event names serve as an event in SPM, e.g., "tumorINCR", tumor volume is the event, and INCR is the attribute.)  
  
  event <- event[ !(is.na(event$eventName)) , ] # remove NA events, e.g., from percentChange
  event$eventName <- as.factor(event$eventName)
  
  return(event)
}

#' Formats event strings and removes incorrect or missing information.
#'
#' @details
#' \code{event} contains events as rows and at least the following columns:
#' \itemize{
#' \item agent = event type, char
#' \item agent_detail = event attribute, char
#' }
#' 
#' @param event a dataframe containing events, see Details
#' @return the \code{event} dataframe with removed events
#'
#' @export
#' @examples
#' data("fake_data")
#' events <- fake_data$events
#' events$startdate <- as.Date(as.character(events$startdate), format='%Y-%m-%d')
#' events$enddate <- as.Date(as.character(events$enddate), format='%Y-%m-%d')
#' 
#' cleaned <- removeData(events)
#' 
removeData <- function(event) {
  event$agent_detail <- lapply(event$agent_detail, function(x) gsub("\\].*","]", x) ) # remove extra details
  
  #### Remove tumor volume measurements which are empty (no value inserted) or  if agent_detail has '-' as value
  cat('...', length(which(event$agent_detail =='')), 
      ' events were removed due to empty agent_details \n')
  event <- event[which(event$agent_detail !=''),]
  
  cat('...', length(which(event$agent_detail =='-')), 
      ' events were removed due to agent_detail = "-" \n')
  event <- event[which(event$agent_detail !='-'),]
  
  cat('...', length(which((event$agent=='mentalStatus' & event$agent_detail=='9'))), 
      ' events were removed due to mental status = 9 \n')
  event <- event[!(event$agent=='mentalStatus' & event$agent_detail=='9'),]
  
  cat('...', length(which((event$agent=='overallNeuro' & event$agent_detail=='3'))), 
      ' events were removed due to overallNeuro = 3 \n')
  event <- event[!(event$agent=='overallNeuro' & event$agent_detail=='3'),]
  
  return(event)
}

#' @title Convert event dates into relative days.
#' 
#' @description Set \code{eventID} relative to each person. 
#' \code{eventID} will be in days after the first observation date.
#' 
#' @param event a dataframe containing all events across all patients, see Details
#' @return the \code{event} dataframe with additional columns: and added columns: \itemize{
#' \item \code{eventID}
#' \item \code{endEventID}
#' }
#' 
#' @details
#' \code{event} contains events as rows and at least the following columns:
#' \itemize{
#' \item iois = patient id, char
#' \item startdate = event start date in the format 'YYYY-MM-DD', char
#' \item enddate = event end date in the format 'YYYY-MM-DD', char
#' }
#' Relative days refer to the number of days since baseline.
#' Baseline date is the earliest event observed in a patient.
#' 
#' \code{eventID}: number of days since first event (grouped by patient id) to the start date of event
#' \code{endEventID}: number of days since first event (grouped by patient id) to the end date of event
#' 
#' @export
#' 
#' @examples
#' data("fake_data")
#' events <- fake_data$events
#' events$startdate <- as.Date(as.character(events$startdate), format='%Y-%m-%d')
#' events$enddate <- as.Date(as.character(events$enddate), format='%Y-%m-%d')
#' 
#' cleaned <- removeData(events)
#' cleaned <- resetEventIDsToDays(cleaned)
#' 
resetEventIDsToDays <- function(event) {
  ####  Set eventIDs relative to each person. EventIDs will be in days after the first observation date.
  # This format will make calculating days until death easier during KM plotting.
  event <- split(event, event$iois)
  event <- lapply(event, function(x) {
    earliestDate <- min(x$startdate)
    x$eventID <- as.integer(x$startdate - earliestDate) + 1 # Add 1 is to ensure eventID day substraction will show day differences
    x$endEventID <- as.integer(x$enddate - earliestDate) + 1
    x
  })
  event <- do.call('rbind', event) # collapse lists into a dataframe
  
  return(event)
}

#' Remove events which occurred prior to a patient's first radiation
#'
#' @param event a dataframe containing all events across all patients (assumes more than 1 patient), see Details
#' @return the \code{event} dataframe with removed 'Radiation' events 
#'
#' @details
#' \code{event} contains events as rows and at least the following columns:
#' \itemize{
#' \item iois = patient id, char
#' \item agent = event type, char
#' \item eventID = event id, num
#' }
#' 
#' Assumes radiation is labelled as 'Radiation' in \code{agent}.
#' @export
#' 
#' @examples
#' data("fake_data")
#' events <- fake_data$events
#' events$startdate <- as.Date(as.character(events$startdate), format='%Y-%m-%d')
#' events$enddate <- as.Date(as.character(events$enddate), format='%Y-%m-%d')
#' 
#' cleaned <- removeData(events)
#' cleaned <- resetEventIDsToDays(cleaned)
#' cleaned <- removePreRadiation(cleaned)
#' 
removePreRadiation <- function(event) {
  removedE <- nrow(event)
  ids <- unique(event$iois)
  removedP <- length(ids)
  event <- split(event, event$iois)
  
  event <- lapply(event, function(x) {
    x <- x[order(x$eventID),]
    if('Radiation' %in% x$agent) {
      index <- min(which(x$agent=='Radiation'))
      if( index < nrow(x) ) {x <- x[(index+1):nrow(x),]
      } else { x <- rep(NA,ncol(x))}
    } else {
      cat('...patient without radiation events: ', x$iois[1], 'with ', nrow(x), ' events \n')
      x <- NA
    }
    x
  })
  event <- do.call('rbind', event) # collapse lists into a dataframe
  event <- event[which(!is.na(event$iois)),] # remove the rows with NA (patients with empty dataframe because never had radiation)
  
  removedE <- removedE - nrow(event)
  cat('...', removedE, 'events were removed due to event being prior to first radiation or no radiation events in patient \n')
  
  
  n <- removedP - length(unique(event$iois))
  cat('...', n, 'patient(s) were removed \n')
  event$agent_detail <- unlist(event$agent_detail)
  return(event)
}

#' Create unique event names for each attribute of an event type.
#' 
#' Event names serve as an unique event in SPM. This will merge agent and agent_detail into a new column, eventName
#' Here, each attribute of an event type is stated in the event name. Events with continuous values will first discretized.
#' (e.g., "tumorINCR", tumor volume is the event from agent column, and INCR is the attribute converted from agent_detail column.)
#'
#' @details
#' \code{event} contains events as rows and at least the following columns:
#' \itemize{
#' \item \code{iois} patient id, char
#' \item \code{agent} event type, char
#' \item \code{agent_detail} event attribute, char
#' \item \code{startdate} event start date in the format 'YYYY-MM-DD', char
#' \item \code{enddate} event end date in the format 'YYYY-MM-DD', char
#' \item \code{eventID} event id
#' \item \code{eventName} distinct names for an event (e.g., mentalStatus0: agent is mental status, agent_detail is 0)
#' }
#'
#' @param event: dataframe, all events across all patients, each row is an event, returned from 'getData' in 'getData.R'
#' @inheritParams getVolTypeName
#' @param removeChemo: bool, remove chemo event types
#'
#' @return event: dataframe, cleaned and formated information that is ammendable to SPM
#' @export
#' 
#' @examples
#' data("fake_data")
#' events <- fake_data$events
#' events$startdate <- as.Date(as.character(events$startdate), format='%Y-%m-%d')
#' events$enddate <- as.Date(as.character(events$enddate), format='%Y-%m-%d')
#' 
#' cleaned <- removeData(events)
#' cleaned <- resetEventIDsToDays(cleaned)
#' cleaned <- removePreRadiation(cleaned)
#' 
#' cleaned <- createEventNames(event=cleaned,
#'                             tType='r',
#'                             removeChemo=T)
#'                             
createEventNames <- function (event, tType, removeChemo) {
  removedE <- nrow(event)
  ids <- unique(event$iois)
  removedP <- length(ids)
  
  event <- split(event, event$agent)
  
  # ignore certain agents
  cat('...', 'keep event agents: \n')
  if (removeChemo) {
    keepNames <- c('KPS', 'measurement', 'mentalStatus','neuroFunc', 'overallNeuro','Radiation', 'Surgery')
    concatThese <- c('neuroFunc','overallNeuro','mentalStatus')
  } else {
    keepNames <- c('Chemotherapy','KPS', 'measurement', 'mentalStatus','neuroFunc', 'overallNeuro','Radiation', 'Surgery')
    concatThese <- c('neuroFunc','overallNeuro','mentalStatus','Chemotherapy') # concat agent and agent detail
  }
  eNames <- names(event)
  ignored <- eNames[which(!eNames %in% keepNames)]
  cat(keepNames, sep='\n')
  event <- event[which(names(event) %in% keepNames)]
  
  
  # concat event names (agent) and values (agent_detail) for event types which are alrady discrete
  for (i in 1:length(concatThese)) { 
    temp <- event[[concatThese[i]]]
    
    if (length(nrow(temp)) == 0 ) next
    
    temp$eventName <- unlist(lapply(seq(1,nrow(temp)), function(x) paste0(temp[x, ]$agent,temp[x,]$agent_detail) ))
    temp$eventName <- unlist(lapply(temp$eventName, function(x) gsub('-','Neg', x)))
    event[[concatThese[i]]] <- temp # set changes
  }
  
  # concat event names (agent) and values (agent_detail) for event types which are continuous and will be discretized
  # - discretize tumor volumes
  temp <- discretizeTumVols(event$measurement, tType)
  event$measurement <- temp[, colnames(temp) %in% c(colnames(event$measurement),'eventName')]
  
  ####  Discretize KPS
  
  temp <- event$KPS
  
  temp <- split(temp, temp$iois) # a list where each item in list is a person's measurements
  for (i in 1:length(temp)) {
    if (length(nrow(temp[[i]])) == 0 ) next
    temp[[i]] <- discretizeKPS(temp[[i]])
  }
  t <- do.call('rbind', temp) # set changes
  event$KPS <- do.call('rbind', temp) # set changes
  
  # copy other events (agent) over to eventNames, since agent_detail will not be used in event name
  if(removeChemo) {
    temp <- lapply(c('Radiation','Surgery'), function(x) event[[x]]$eventName <<- event[[x]]$agent) #remove spaces
  } else {
    temp <- lapply(c('Radiation','Surgery','Chemotherapy'), function(x) event[[x]]$eventName <<- event[[x]]$agent) #remove spaces
  }
  
  # events for pattern mining
  event <- do.call('rbind', event) # collapse lists into a dataframe
  event$eventName <- unlist(event$eventName)
  
  removedE <- removedE - nrow(event) ; cat('...', removedE, 'events were removed due to no events in agents of interest \n')
  
  event <- event[which(!is.na(event$iois)),]
  n <- removedP - length(unique(event$iois)) ; cat('...', n, 'patients were removed due to no events in agents of interest \n')
  
  ids <- ids[-which(ids %in% unique(event$iois))] ; cat('...these iois were removed:'); print(ids); cat('\n')
  
  return(event)
}


#---------------------------------- parse data  -----------------------------------

#' Get the full name of the tumor volume type.
#'
#' Map shortcut name to fullname for type of tumor volume information.
#' Acceptable types: \itemize{
#'         \item 'c' for continuous, left as is, e.g., volume was measured in cubic mm.
#'         \item 'p' for percent change
#'         \item 'r' for volumetric response criteria
#'         \item 'rate' for rate change
#'         \item 'rc' for both 'c' and 'rate'
#' }
#' See \code{prepTumorVol} for full details.
#'
#' @param tType char, shortcut for type of tumor volume information
#' @return full name of tumor volume type
#' @export
#' @examples
#' getVolTypeName('p')
#' getVolTypeName('r')
#' getVolTypeName('rate')
#' 
getVolTypeName <- function(tType) {switch(tType, 
                                          p = 'percentChange', 
                                          c='continuous', 
                                          r='responseCriteria', 
                                          rate='rateChange', 
                                          rc='rateChange_continuous') }


#' Calculate numeric age.
#'
#' @param birthDate: date of birth in format 'YYYY-MM-DD', date
#' @param date: current date to evaluate age in format 'YYYY-MM-DD', date
#' @return \code{age}, numeric age at \code{date}
#' @export
#' @examples
#' data("fake_data")
#' birthdates <- fake_data$demo$DOB
#' today <- as.Date('2000-10-30')
#' ages <- getAge(birthdates, date=today)
#' 
getAge <- function(birthDate, date) {
  age <- as.numeric(difftime(date,birthDate, unit='weeks'))/52.25
}

#' Parse the tumor location string to get laterality.
#'
#' @param str: tumor location string, char
#' @return tumor laterality, char
#' @export
#' @examples
#' data("fake_data")
#' 
#' # do for all patients
#' lat_strs <- fake_data$events[fake_data$events$agent=='Surgery',] # event has tumor location information
#' lat_strs$laterality <- unlist(lapply(lat_strs$agent_detail, function(i) getLaterality(i) ))
#' 
#' # single patient
#' lat_str <- fake_data$events[fake_data$events$agent=='Surgery',][1,]
#' lat <- getLaterality(lat_str$agent_detail)
#' 
getLaterality <- function(str) {
  str <- tolower(str)
  l <- c(grepl('left', str), grepl('right',str), grepl('bi',str), grepl('both',str))
  names(l) <- c('left','right','both','both')
  n <- names(which(l==T))
  
  ifelse(is.null(n), NA, n)
}

#' Parse the tumor location string to get location
#'
#' @param str: tumor location string, char
#' @return tumor location, char
#'
#' @details Possible locatios are  c('frontal','temporal', 'occipital','parietal','cerebellar', 'thalamic', 'sellar',
#' 'corpus callosum', 'pineal', 'midbrain').
#' 
#' @export
#' @examples
#' data("fake_data")
#' 
#' # single patient
#' loc_str <- fake_data$events[fake_data$events$agent=='Surgery',][1,]
#' loc <- getLocation(loc_str$agent_detail)
#' 
getLocation <- function(str) {
  str <- tolower(str)
  n <- c('frontal','temporal', 'occipital','parietal','cerebellar', 'thalamic', 'sellar', 
         'corpus callosum', 'pineal', 'midbrain')
  l <- if (is.na(str)) { 
    rep(NA, length(n))
  } else {
    unlist(lapply(n, function(i) grepl(i, str)))
  }
  
  
  names(l) <- paste0('location: ', n)
  
  
  return(l)
}

#---------------------------------- prep clinical visits for logit modeling -----------------------------------

#' Remove clinical visits with no history of length \code{maxgap}.
#'
#' @details clinical visits that do not have a history of maxgap time cannot be evaluated
#' (such vists would have indeterminate values for temporal patterns item1 --> item2 when the maxgap exisits between these items)
#'   e.g., at the first observed visit for a patient, I search for item2's existence, and then look backwards in time to see
#'   if item1 happend within the maxgap. I wouldn't be able to determien this pattern's value because no visits are before the first visit.
#' Therefore, such clinical visits should be removed.
#' 
#' \code{data} dataframe, 
#' rows are clinical visits to be classified, 
#' columns are features including clinical visits ids 
#' and at least the following columns:
#' \itemize{
#' \item \code{id} char, clinical visit id in the format "patientID.eventID"
#' \item \code{iois} char, patient id
#' \item \code{eventID} char, event id
#' }
#' 
#' @inheritParams arulesSequences::cspade
#' @param data dataframe, see Details
#' @param tType char, tumor volume type nickname, see \code{\link{getVolTypeName}}
#' @param save bool, save removed data?
#' @param outDir char, full path to output directory
#' 
#' @return 
#'
#' @export
#' 
#' @examples
#' 
#' data("features_ratechange_sup0.4g60l2z2") # features and labels for each clinical visit
#' maxgap <- 60 # constants used to create features
#' maxlen <- 2
#' t <- 'rate'
#' 
#' # format
#' names <- colnames(feats)
#' feats <- data.frame(id=row.names(feats),feats)
#' colnames(feats) <- c('id',names)
#' cat('...overall samples: ', nrow(feats), '\n')
#' feats <- prepLaterality(feats) 
#' feats <- prepLocation(feats)
#' 
#' cat('...removing clinical visits with no history of length maxgap:', maxgap, '\n')
#' feats <- removeVisits(feats, 
#'                      maxgap=maxgap, 
#'                      maxlength=maxlen, 
#'                      tType=t, 
#'                      save=F, 
#'                      outDir=NA)
#' 
removeVisits <- function(data, maxgap, maxlength, tType, save=TRUE, outDir) {
  ids <- strsplit(as.character(data$id), '\\.')
  data$iois <- as.integer(unlist(lapply(ids, '[[',1)))
  data$eventID <- as.integer(unlist(lapply(ids, '[[',2)))
  data <- split(data, data$iois)
  
  removed <- list()
  data <- lapply(data, function(patient) {
    patient <- patient[order(patient$eventID),] # events in chronological order
    daysSinceFirstVist <- patient$eventID - patient$eventID[1]
    withinMaxGap <- which(daysSinceFirstVist <= (maxgap*(maxlength-1)) )
    cutoff <- max(withinMaxGap)
    
    # save removed visits for inspection later
    if (save) {
      removedVisits <- patient[1:cutoff,]
      removedVisits$daysSinceFirstVist <- daysSinceFirstVist[1:cutoff]
      removed <<- c(removed,list(removedVisits))
    }
    patient <- patient[-c(1:cutoff),]
  })
  
  if (save) {
    removed <- do.call('rbind', removed)
    cat('...saved at ', outDir, '\n')
    saveRDS(removed, file=file.path(outDir, paste0('removedLogitData_maxgap_',maxgap,'_',getVolTypeName(tType),'.rds')))
  }
  data <- do.call('rbind',data)
  cat('...visits left for training: ',nrow(data), '\n')
  return(data)
}

#' Get the tumor location and laterality.
#'
#' @details
#' Takes the first event with tumor location information in \code{tumorInfo} and
#' convert location into dummy variables for logit since one person can have multiple locations.
#' \code{tumorInfo} contains pre-cleaned data (\code{\link{cleanData}}
#' is used for cleaning events for SPM) \code{agent_details}.
#' Also retrieves tumor lateriality.
#'
#' Assumes location and laterility is in \code{agent_details} column.
#' Assumes location is after "Location :" in string.
#' Assumes location and lateriality is unchanged for the remaining events (is static)
#'
#' @param data dataframe, 
#' rows are clinical visits to be classified, 
#' columns are features including clinical visits ids 
#' and at least the following columns:
#' \itemize{
#' \item \code{id} char, clinical visit id in the format "patientID.eventID"
#' \item \code{iois} char, patient id
#' \item \code{agent_detail} char, event attribute
#' }
#' @param tumorInfo a dataframe of pre-cleaned events and searches for tumor location
#' @return \code{data} with additional columns containing tumor location dummy variables and tumor laterality
#' @export
#' 
#' @examples
#' 
#' data("fake_data")
#' t <- 'rate'
#' fake_tumorInfo <- fake_data$events # save tumor location and laterility strings before event cleaning
#' fake_demo <- fake_data$demo
#' fake_data$events <- cleanData(fake_data$events, tType = t)
#' 
#' # collect patient info for each event
#' # note that only eventName, iois, and eventID columns are used for SPM
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
#' fake_data <- prepDemographics(fake_data, fake_demo) # prep for each event, since age does change
#' fake_data <- prepSurvivalLabels(fake_data) # get survival labels, these also change
#' 
#' fake_data <- getTumorLocation(fake_data, fake_tumorInfo) # get first tumor location
#' 
getTumorLocation <- function(data, tumorInfo) {
  cat('parsing tumor locations...\n')
  
  tumorInfo <- split(tumorInfo, tumorInfo$iois)
  tumorLocation <- lapply(tumorInfo, function(p) {
    p <- p[order(p$startdate),]
    t <- grep('Location :*', p$agent_detail)
    ifelse(length(t), gsub('^.*Location :','',p$agent_detail[t])[1], NA)
  })
  cat('...patients without tumor location : ', names(which(is.na(tumorLocation))), '\n' )
  
  tumorLocation <- data.frame(iois=names(tumorLocation), location=unlist(tumorLocation))
  
  tumorLocation$tumorLaterality <- lapply(tumorLocation$location, function(s) getLaterality(s))
  
  dummyVars <- lapply(tumorLocation$location, function(s) getLocation(s))
  dummyVars <- do.call('rbind', dummyVars)
  tumorLocation <- cbind(tumorLocation, dummyVars)
  
  tumorLocation <- tumorLocation[-2]
  tumorLocation$tumorLaterality <- as.factor(unlist(tumorLocation$tumorLaterality))
  
  data <- merge(data, tumorLocation, by='iois')
  
  return(data)
}

#' Get all tumor volume types.
#'
#' @details
#' \code{volumes} contains tumor measurement events as rows and at least the following columns: \itemize{
#' \item iois = patient id, char
#' \item agent_detail = event attribute, i.e., measured tumor volume in mm3 char
#' \item startdate = event start date in the format 'YYYY-MM-DD', char
#' \item eventName = event name, char
#' }
#'
#' Tumor volume types are: \itemize{
#' \item volume = volume in original measured units (e.g., mm^3)
#' \item baseline size = earliest tumor volume observed in a patient
#' \item percent change = \eqn{(Volume_{current} - Volume_{previous})/Volume_{previous} }
#' \item rate change = \eqn{(Volume_{current} - Volume_{previous})/(DaysSinceBaseline_{current} - DaysSinceBaseline_{previous}) }
#' }
#'
#' @param volumes a dataframe containing  tumor measurement events acrosss patients, see Details
#' @return \code{volumes} dataframe, with additional columns for each type of tumor volume information, see Details
#' @export
#' 
#' @examples
#' data("fake_data")
#' fake_demo <- fake_data$demo
#' fake_tumorInfo <- fake_data$events
#' 
#' # by default eventName is the volumetric response critera
#' # will get all types eventually
#' fake_data$events <- cleanData(fake_data$events) 
#' 
#' # keep tumor volume only
#' fake_data$events <- fake_data$events[fake_data$events$agent=='measurement',] 
#' fake_data$events <- prepTumorVol(fake_data$events) # get all tumor vol types
#' 
prepTumorVol <- function(volumes) {
  
  # discretize tumor vols
  
  volumes$discreteContinuous <- discretizeTumVols(volumes, tType='c')$eventName
  
  temp <- discretizeTumVols(volumes, tType='b')[,c('baselineSize','eventName')] 
  volumes$baselineSize <- temp$baselineSize
  volumes$discreteBaseline <- temp$eventName
  
  temp <- discretizeTumVols(volumes, tType='p')[,c('percentChange','eventName')] 
  volumes$percentChange <- temp$percentChange
  volumes$discretePercent <- temp$eventName
  
  temp <- discretizeTumVols(volumes, tType='rate')[,c('rateChange','eventName')]
  volumes$rateChange <- temp$rateChange
  volumes$discreteRateChange <- temp$eventName
  
  # double checking response criteria function
  volumes <- subset(volumes, select=-c(eventName))
  volumes$responseCriteria <- discretizeTumVols(volumes, tType='r')$eventName 
  
  # now can format agent_detail
  colnames(volumes)[which(colnames(volumes)=='agent_detail')] <- 'continuous' # change name of tumor vol measurement
  volumes$continuous <- as.numeric(unlist(volumes$continuous,recursive = F))
  return(volumes)
}

#' Prepare tumor laterality feature for logit modeling.
#' 
#' @details Converts tumor laterality feature into two boolean features: "left" and "right".
#' Searches for "left", "right", "both" in \code{tumorLaterality} 
#' and sets "both" as "left" and "right."
#' 
#' @param data dataframe, 
#' rows are clinical visits to be classified, 
#' columns are features including clinical visits ids,
#' and at least the following column:
#' \itemize{
#' \item \code{tumorLaterality} char, tumor laterality
#' }
#' @return \code{data} dataframe with replaced laterality columns
#' @export
#' 
#' @examples
#' 
#' data("features_ratechange_sup0.4g60l2z2")  # features and labels for each clinical visit
#' 
#' # format
#' names <- colnames(feats)
#' feats <- data.frame(id=row.names(feats),feats)
#' colnames(feats) <- c('id',names)
#' cat('...overall samples: ', nrow(feats), '\n')
#' 
#' feats <- prepLaterality(feats) 
#' 
prepLaterality <- function(data) {
  # convert laterality to left or right, and absorb bilateral
  laterality <- data$tumorLaterality
  ind <- which(laterality=='both')
  laterality <- data.frame(right=laterality=='right', left=laterality=='left')
  laterality$right[ind] <- TRUE
  laterality$left[ind] <- TRUE
  colnames(laterality) <- paste('laterality: ',colnames(laterality))
  
  replaceInd <- which(colnames(data)=='tumorLaterality')
  data <- cbind(data[1:(replaceInd-1)], laterality, 
                data[(replaceInd+1):ncol(data)]) # insert new laterality variables
  return(data)
}

#' Prepare tumor location feature for logit modeling.
#' 
#' Locations with few observations are merged into 'other'. 
#' Categories include: frontal, temporal, parietal, and other.
#' 
#' @param data dataframe, 
#' rows are clinical visits to be classified, 
#' columns are features including clinical visits ids,
#' and at least the columns returned from \code{\link{getLocation}}
#' @return \code{data} dataframe with replaced location columns 
#' @export
#' 
#' @examples
#' 
#' data("features_ratechange_sup0.4g60l2z2")  # features and labels for each clinical visit
#' 
#' # format
#' names <- colnames(feats)
#' feats <- data.frame(id=row.names(feats),feats)
#' colnames(feats) <- c('id',names)
#' cat('...overall samples: ', nrow(feats), '\n')
#' 
#' feats <- prepLaterality(feats) 
# feats <- prepLocation(feats)
#' 
prepLocation <- function(data) {
  # merge tumor locations
  removeThis <- which(colnames(data) %in% c('location: sellar'))
  data <- data[, -removeThis]
  mergeThese <- c('location: occipital', 
                  'location: thalamic', 
                  'location: corpus callosum', 
                  'location: cerebellar',
                  'location: pineal', 
                  'location: midbrain')
  replaceInd <- which(colnames(data) %in% mergeThese)
  otherLoc <- data[, replaceInd]
  otherLoc$'location: other' <- apply(otherLoc, 1, function(i) sum(i))
  otherLoc$'location: other' <- as.logical(otherLoc$'location: other')
  data <- data[,-replaceInd]
  
  # insert new laterality variables
  splitCol <- replaceInd[1]
  if (splitCol == ncol(data)) {
    data <- cbind(data, 'location: other'=otherLoc$'location: other') 
  } else {
    nextCol <- splitCol+1
    data <- cbind(data[1:splitCol], 
                  'location: other'=otherLoc$'location: other', 
                  data[nextCol:ncol(data)])
  }

  return(data)
}

#' Add and prep demographic information (static 'events') to events.
#'
#' Ethnicity is either 'white' or 'non-white'.
#' Age is broken down to decades, with 70-90 decades merged.
#'
#' @param data dataframe, 
#' rows are clinical visits to be classified, 
#' columns are features including clinical visits ids 
#' and at least the following columns:
#' \itemize{
#' \item \code{iois} char, patient id
#' \item \code{startdate} = event start date in the format 'YYYY-MM-DD', char
#' }
#' 
#' @param demo dataframe, patient demograpics and biomarker information, containing columns:
#' \itemize{
#' \item \code{iois} char, patient id
#' \item \code{ethnicity} char, patient ethnicity, uses value 'WHITE' to differentiate into three categories: WHITE, non-WHITE, and unknown
#' \item \code{DOB} char, date of birth in the format 'YYYY-MM-DD'
#' \item \code{MSP} char, MGMT promoter methytation status, 'U' = unmethylated, 'M' = methylated
#' }
#' returns
#' - data: dataframe, with formated demographics
#' @export
#' 
#' @examples
#' data("fake_data")
#' t <- 'rate'
#' fake_tumorInfo <- fake_data$events # save tumor location and laterility strings before event cleaning
#' fake_demo <- fake_data$demo
#' 
#' fake_data$events <- cleanData(fake_data$events, tType = t)
#' 
#' # collect patient info for each event
#' # note that only eventName, iois, and eventID columns are used for SPM
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
#' fake_data <- prepDemographics(fake_data, fake_demo) # prep for each event, since age does change
#' 
prepDemographics <- function(data, demo, ageCategories = c(seq(20,70,10),90)) {
  data <- merge(data, demo, by='iois')
  
  data$gender <- as.factor(data$gender)
  
  h <- which(data$ethnicity=='', arr.ind=T)
  data$ethnicity[-(which(data$ethnicity=='WHITE',arr.ind=T))] <- 'non-WHITE' # merge ethnicity data
  data$ethnicity[h] <- 'unknown'
  data$ethnicity <- as.factor(data$ethnicity)
  
  data$DOB <- as.Date(as.character(data$DOB), format='%Y-%m-%d')
  data$ageDecade <- getAgeDecade(data$DOB, data$startdate, ageCategories)
  data$age <- getAge(data$DOB, data$startdate)
  
  # gene info
  data$MSP <- as.factor(data$MSP) # mgmt methyl
  
  return(data)
}

#' Prep survival information (labels for all classifictation tasks) to events.
#'
#' @param data dataframe, 
#' rows are clinical visits to be classified, 
#' columns are features including clinical visits ids 
#' and at least the following columns:
#' \itemize{
#' \item \code{eventID} num, event id
#' \item \code{survival} num, days since baseline
#' \item \code{status} char, patient status at the end of data collection, 'DECEASED' or 'ALIVE'
#' }
#' 
#' @return data dataframe, 
#' rows are clinical visits to be classified, 
#' columns are features including clinical visits ids 
#' and the following added columns:
#' \itemize{
#' \item \code{daysToDeath} num, event id
#' \item \code{survivalIn60} char, 'dead' or 'alive' in 60 days or ~2 months
#' \item \code{survivalIn180} char, 'dead' or 'alive' in 180 days or ~6 months
#' \item \code{survivalIn270} char, 'dead' or 'alive' in 270 days or ~9 months
#' \item \code{survivalIn360} char, 'dead' or 'alive' in 360 days or ~12 months
#' }
#' @export
#' 
#' @examples
#' data("fake_data")
#' t <- 'rate'
#' fake_tumorInfo <- fake_data$events # save tumor location and laterility strings before event cleaning
#' fake_demo <- fake_data$demo
#' 
#' fake_data$events <- cleanData(fake_data$events, tType = t)
#' 
#' # collect patient info for each event
#' # note that only eventName, iois, and eventID columns are used for SPM
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
#' fake_data <- prepDemographics(fake_data, fake_demo) # prep for each event, since age does change
#' fake_data <- prepSurvivalLabels(fake_data) # get survival labels, these also change
#' 
prepSurvivalLabels <- function(data) {
  data$daysToDeath <- data$survival-data$eventID
  data$survivalIn60 <- as.factor(ifelse(data$daysToDeath <= 60 & data$status=='DECEASED','dead','alive' ))
  data$survivalIn180 <- as.factor(ifelse(data$daysToDeath <= 180 & data$status=='DECEASED','dead','alive' ))
  data$survivalIn270 <- as.factor(ifelse(data$daysToDeath <= 270 & data$status=='DECEASED','dead','alive' ))
  data$survivalIn360 <- as.factor(ifelse(data$daysToDeath <= 360 & data$status=='DECEASED','dead','alive' ))
  return(data)
}






#' Get events that overlap with chemotherapy, specifically just the top 5 most used drug.
#'
#' args
#' - event: dataframe of all events called from Agile database (getData())
#'
#' returns
#' - overlaps: dataframe, events that occurred during chemotherapy
#' @export
getChemoOverlapEvents <- function(event, tType) {
  event <- cleanData(event,tType, removeChemo=F) 
  
  chemo <- event[event$agent=='Chemotherapy',] # 1476 chemo events, 1299 after data cleaning
  event <- event[!event$agent=='Chemotherapy',] # all other events, to use for spm
  
  drug <- as.data.frame(table(chemo$agent_detail))
  colnames(drug) <- c('name','freq')
  drug <- drug[order(-drug$freq),]
  topDrug <- drug[1:5,]
  
  allChemo <- chemo
  chemo <- chemo[chemo$agent_detail %in% topDrug$name,] # 923 with top drugs
  
  #-------------- get event overlaps -----------#
  
  # overlap data ranges in original chemo data with pattern data 
  
  event <- split(event, event$iois) # all non-chemo agent events
  chemo <- split(chemo,chemo$iois) # take only the top chemo types
  
  overlaps <- lapply(names(chemo), function(n) { # iterate through patients with chemo
    c <- chemo[[n]]
    c <- c[order(c$eventID),] # order chemo events
    e <- event[[n]]
    
    o <- lapply(seq(1,nrow(c)), function(i) {
      b <- which(e$eventID %between% c(c$eventID[i], c$endEventID[i])) # dates that overlap
      m <- e[b,] # get patterns
    })
    o <- do.call('rbind',o)
    o <- o[!duplicated(o), ] # remove events that may have occurred over two chemo ranges
  })
  overlaps <- do.call('rbind',overlaps)
  
  return(overlaps)
}