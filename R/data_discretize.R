# Helper functions for 'cleanData.R'
# Discretize events with continuous attributes.
# Discretized events are placed into 'eventName' column

#' Get age decade.
#' 
#' @param birthDate date of birth, in the format 'YYYY-MM-DD', chars will be converted to Date
#' @param date date to evaluate age decade, in the format 'YYYY-MM-DD', chars will be converted to Date
#' @inheritParams arules::discretize
#' 
#' @return age decade
#' 
#' @seealso \code{\link[arules]{discretize}}
#'
#' @export
#'
#' @examples
#'
# data(fake_data)
# ageCategories <- seq(0,100,10)
# birthdates <- fake_data$demo$DOB
# today <- as.Date("2000-10-30")
# 
# age_decades <- getAgeDecade(birthdates, today, ageCategories)

getAgeDecade <- function(birthDate, date, categories) {
  age <- getAge(birthDate, date)
  age <- floor(age)
  decade <- arules::discretize(age, method='fixed', categories=categories, ordered=F)
}

#' Discretize KPS according to expert opinion
#'
#' KPS is categorized into UNCH, SigDECR, SigINCR, DECR, and INCR.
#'
#' KPS categories are \itemize{
#' \item initial: first KPS in sequence
#' \item UNCH: unchanged, KPS is the same as the previous KPS
#' \item SigDECR: significant decrease, KPS drops to <= 60 OR KPS was already <= 60 and drops further
#' \item SigINCR: significant increase, KPS is > 60 AND previous KPS was <=60
#' \item DECR: decrease, KPS drops but remains above 60 AND psrevious KPS was above 60
#' \item INCR: increase, KPS improves AND is above 60 AND previous KPS was above 60
#' }
#' @param kps a dataframe of KPS scores, with columns: \itemize{
#' \item \code{startdate} date of KPS evaluation)
#' \item \code{agent_detail} KPS}
#' 
#' @return A dataframe with original data \code{kps} 
#' plus discretized KPS categories in \code{eventName} column.
#' 
#' @export
#' 
#' @examples
#' data("fake_data")
#' 
#' # get one person's set of KPS exams
#' kps_set <- fake_data$events[fake_data$events$agent=='KPS' & fake_data$events$iois=='2' , ]
#' 
#' # convert to categorical
#' kps_set_discretized <- discretizeKPS(kps_set)

discretizeKPS <- function(kps) {
  kps <- kps[order(kps$startdate, decreasing=F), ]
  kps$eventName <- lapply(seq(1,nrow(kps)), function(x) # for each person's measurements, find changes that occur
    if (x==1) name <- 'kpsInitial'
    else {
      previous <- as.numeric(kps[(x-1),]$agent_detail)
      current <- as.numeric(kps[x,]$agent_detail)
      if (current==previous) { name <- 'kpsUNCH'  
      } else if (current <= 60 || previous <= 60) {
        if (current < previous) {name <- 'kpsSigDECR'
        } else name <- 'kpsSigINCR'
      } else {
        if (current < previous) {name <- 'kpsDECR'
        } else name <- 'kpsINCR'
      }
    })
  kps$eventName <- unlist(kps$eventName)
  
  return(kps)
}


#' Discretize tumor volume.
#'
#' @details \code{tumorVols} contains tumor measurement events as rows and at least the following columns: \itemize{
#' \item iois = patient id, char
#' \item agent_detail = event attribute, i.e., measured tumor volume in mm3 char
#' \item startdate = event start datein the format 'YYYY-MM-DD', char
#' \item eventName = event name, char
#' }
#'
#' @param tumorVols a dataframe of tumor volume, see Details
#' @param tType char, abbrev. type of tumor data to use, see \code{\link{getVolTypeName}} for valid values
#'
#' @return A dataframe with original data in \code{tumorVols} plus discretized volumes in \code{eventName} column.
#'
#' @export
#' 
#' @examples
#' data("fake_data")
#' vols <- fake_data$events[fake_data$events$agent=='measurement', ]
#' 
#' # format data
#' vols$startdate <- as.Date(as.character(vols$startdate), format='%Y-%m-%d')
#' vols$enddate <- as.Date(as.character(vols$enddate), format='%Y-%m-%d')
#' vols <- resetEventIDsToDays(vols)
#' vols <- removeData(vols) # get rid of excess details and empty or no data
#' 
#' vols_dist <- discretizeTumVols(vols, tType = 'rate')

discretizeTumVols <- function(tumorVols, tType) {
  
  if (tType=='c') { 
    tumorVols <- setContinDiscrete(tumorVols) # continuous values are already in agent_detail so not repeated here.
  } else {
    if (tType=='p') {
      tumorVols <- split(tumorVols, tumorVols$iois) # list of dataframes, where each is a patient's set of tumor volumes
      tumorVols <- lapply(tumorVols, function(i) setPercentChange(i))
      tumorVols <- do.call('rbind', tumorVols)
      tumorVols <- setPercentDiscrete(tumorVols)
      
    } else if (tType=='rc') {
      # continuous values are already in agent_detail so not repeated here.
      continuousVols <- setContinDiscrete(tumorVols) 
      tumorVols <- split(tumorVols, tumorVols$iois)
      tumorVols <- lapply(tumorVols, function(i) setRatechange(i))
      tumorVols <- do.call('rbind', tumorVols)
      tumorVols <- setRateDiscrete(tumorVols)
      tumorVols <- tumorVols[,!(colnames(tumorVols) %in% c('rateChange'))]
      
      tumorVols <- do.call('rbind', list(tumorVols,continuousVols))
      
    } else if (tType=='b') {
      tumorVols <- split(tumorVols, tumorVols$iois) 
      for (i in 1:length(tumorVols)) {
        tumorVols[[i]]$baselineSize <- tumorVols[[i]][1,]$agent_detail # add baseline
      }
      tumorVols <- do.call('rbind', tumorVols)
      tumorVols <- setBaselineDiscrete(tumorVols)
      
    } else if (tType=='r') {
      tumorVols <- split(tumorVols, tumorVols$iois) 
      tumorVols <- lapply(tumorVols, function(i) setResponse(i))
      tumorVols <- do.call('rbind', tumorVols)
      
    } else if (tType=='rate') {
      tumorVols <- split(tumorVols, tumorVols$iois) 
      tumorVols <- lapply(tumorVols, function(i) setRatechange(i))
      tumorVols <- do.call('rbind', tumorVols)
      tumorVols <- setRateDiscrete(tumorVols)
    }
  }
  
  return(tumorVols)
}


#' Calculate different tumor volume information types.
#' 
#' \code{setPercentChange}: get volume percent changes from a sequence of volume measurements.
#'
#' @details  Percent changes, rate changes, and volumetric response criteria (see \code{prepTumorVol}) 
#' are based on volumes in \code{volumes}, where rows are ordered (e.g,. by time).
#'
#' \code{volumes} contains tumor measurement events as rows and at least the following column: \itemize{
#' \item agent_detail = tumor volume in original units (e.g, cubic mm)
#' \item eventID = baseline day of occurence, accepts only one volume per day
#' }
#'
#' @param volumes dataframe, a single patient's tumor volume events
#' @return \code{volumes} with \code{percentChange} as an additional column
#' @export
#' @examples
#' data("fake_data")
#' vols <- fake_data$events[fake_data$events$agent=='measurement', ]
#' 
#' # format data
#' vols$startdate <- as.Date(as.character(vols$startdate), format='%Y-%m-%d')
#' vols$enddate <- as.Date(as.character(vols$enddate), format='%Y-%m-%d')
#' vols <- resetEventIDsToDays(vols)
#' vols <- removeData(vols) # get rid of excess details and empty or no data
#' 
#' # get vols for one patient
#' vols_set <- vols[vols$iois=='1',] 
#' p <- setPercentChange(vols_set)
#' 
setPercentChange <- function(volumes) {
  volumes$agent_detail <- as.numeric(volumes$agent_detail)
  volumes <- volumes[order(volumes$eventID, decreasing=F), ] # order tumor volume by date
  
  # get percent changes and,
  # set Inf and NaN percent change to NA
  # NaN: V1=0 and V2=V1, V1=0 and V2=0
  # Inf: V2>0 and V1=0
  
  volumes$percentChange <-volumes$agent_detail # create another colume to write over
  
  volumes$percentChange <- lapply(seq(1,nrow(volumes)), function(x) {
    if (x==1) { 
      percent <- NA
    } else if (volumes$agent_detail[x-1]==0) { 
      percent <- NA
    } else {
      percent <- 100*(volumes$agent_detail[x] - volumes$agent_detail[x-1])/volumes$agent_detail[x-1] # percent change = 100*(V2-V1)/V1
    }
  })
  
  volumes$percentChange <- unlist(volumes$percentChange)
  return(volumes)
}
#' 
#' \code{setRatechange}: get volume rate changes from a sequence of volume measurements.
#' @rdname setPercentChange
#' @export
#' @examples
#' data("fake_data")
#' vols <- fake_data$events[fake_data$events$agent=='measurement', ]
#' 
#' # format data
#' vols$startdate <- as.Date(as.character(vols$startdate), format='%Y-%m-%d')
#' vols$enddate <- as.Date(as.character(vols$enddate), format='%Y-%m-%d')
#' vols <- resetEventIDsToDays(vols)
#' vols <- removeData(vols) # get rid of excess details and empty or no data
#' 
#' vols_set <- vols[vols$iois=='55',]
#' r <- setRatechange(vols_set)
#' 
setRatechange <- function(volumes) {
  volumes$agent_detail <- as.numeric(volumes$agent_detail)
  volumes <- volumes[order(volumes$eventID, decreasing=F), ] # order tumor volume by date
  volumes$rateChange <- volumes$agent_detail # create another colume to write over
  
  volumes$rateChange <- lapply(seq(1,nrow(volumes)), function(x) {
    if (x==1) { 
      rate <- NA
    } else {
      rate <- (volumes$agent_detail[x] - volumes$agent_detail[x-1])/(volumes$eventID[x] - volumes$eventID[x-1])
    }
  })
  
  volumes$rateChange <- unlist(volumes$rateChange)
  return(volumes)
}


#' \code{setResponse}: get volumetric response criteria by mapping continuous volumes to response criteria categories
#' @rdname setPercentChange
#' @examples
#' data("fake_data")
#' vols <- fake_data$events[fake_data$events$agent=='measurement', ]
#' 
#' # format data
#' vols$startdate <- as.Date(as.character(vols$startdate), format='%Y-%m-%d')
#' vols$enddate <- as.Date(as.character(vols$enddate), format='%Y-%m-%d')
#' vols <- resetEventIDsToDays(vols)
#' vols <- removeData(vols) # get rid of excess details and empty or no data
#' 
#' vols_set <- vols[vols$iois=='55',]
#' r <- setResponse(vols_set)
#' 
setResponse <- function(volumes)  {
  volumes$agent_detail <- as.numeric(volumes$agent_detail)
  volumes <- volumes[order(volumes$eventID, decreasing=F), ] # order tumor volume by date
  
  volumes$eventName <- lapply(seq(1,nrow(volumes)), function(x)
    if (x==1) {
      name <- 'tumorInitial'
    } else {
      if ( volumes$agent_detail[x]==0 ) name <- 'tumorCompleteRes'
      else if ( volumes$agent_detail[x] >= volumes$agent_detail[x-1]*1.40 ) name <- 'tumorProgression'
      else if ( volumes$agent_detail[x] <= volumes$agent_detail[x-1]*0.65 ) name <- 'tumorPartialRes'
      else name <- 'tumorStableDis'
    })
  
  volumes$eventName <- unlist(volumes$eventName)
  return(volumes)
}

#' Discretize tumor volumes into intervals, where each interval is a category (an event).
#' 
#' \code{setContinDiscrete}: convert volumes, measured in original units, to intervals.
#'
#' @details  Volumes are discretized into 10 intervals, where each is an event. Event names are named after the intervals.
#' Note: \itemize{
#' \item intervals are based on distribution of values in \code{volumes}
#' \item each interval or bin will be used as different categories of tumor volume
#' \item also see \code{prepTumorVol}
#' }
#'
#' \code{volumes} contains tumor measurement events as rows and at least the following columns: \itemize{
#' \item agent_detail = event attribute, i.e., measured tumor volume in mm3 char
#' }
#'
#' @param volumes a dataframe, rows are tumor volumes
#' @inheritParams arules::discretize
#' @seealso 
#' \code{\link[arules]{discretize}}
#' @return \code{volumes} with \code{eventName} as an additional column
#' @export
#' @examples
#' data("fake_data")
#' vols <- fake_data$events[fake_data$events$agent=='measurement', ]
#' 
#' # format data
#' vols$startdate <- as.Date(as.character(vols$startdate), format='%Y-%m-%d')
#' vols$enddate <- as.Date(as.character(vols$enddate), format='%Y-%m-%d')
#' vols <- resetEventIDsToDays(vols)
#' vols <- removeData(vols) # get rid of excess details and empty or no data
#' 
#' vols_dist <- setContinDiscrete(vols)
#' 
setContinDiscrete <- function(volumes, categories=10) {
  volumes$agent_detail <- as.numeric(volumes$agent_detail)
  # TO-DO: do this for all other discretize methods since arules has been updated to stop if duplicated breaks.
  # there are duplicate breaks in continouous tumor volume because there are many 0
  # this will just take the unique breaks
  breaks <- quantile(volumes$agent_detail, probs = seq(0,1, length.out = categories+1), na.rm = TRUE)
  volumes$eventName <- cut(volumes$agent_detail, breaks = unique(breaks),include.lowest = T, right = F, ordered_result = F)
  volumes$eventName <- gsub(',', ' to ', volumes$eventName)
  volumes$eventName<- paste0('tumorVol:', volumes$eventName)
  return(volumes)
}

#' \code{setPercentDiscrete}: convert volumes, measured in percent change, to intervals.
#' @rdname setContinDiscrete
#' @export
#' 
#' @examples
#' data("fake_data")
#' vols <- fake_data$events[fake_data$events$agent=='measurement', ]
#' 
#' # format data
#' vols$startdate <- as.Date(as.character(vols$startdate), format='%Y-%m-%d')
#' vols$enddate <- as.Date(as.character(vols$enddate), format='%Y-%m-%d')
#' vols <- resetEventIDsToDays(vols)
#' vols <- removeData(vols) # get rid of excess details and empty or no data
#' 
#' # convert vol to correct type for all patients
#' vols <- split(vols, vols$iois) # list of dataframes, where each is a patient's set of tumor volumes
#' vols <- lapply(vols, function(v) setPercentChange(v))
#' vols <- do.call('rbind', vols)
#' 
#' # discretize
#' vols_dist <- setPercentDiscrete(vols)
#' 
setPercentDiscrete <- function(volumes, categories=10) {
  volumes$eventName <- discretize(volumes$percentChange, method="frequency", categories = 10, ordered=F)
  volumes$eventName <- gsub(',', ' to ', volumes$eventName)
  volumes$eventName[!(is.na(volumes$eventName))]  <- paste0('tumorVol%Change:', volumes$eventName[!(is.na(volumes$eventName))] )
  return(volumes)
}
#' \code{setRateDiscrete}: convert volumes, measured in rate change, to intervals.
#' @rdname setContinDiscrete
#' @export
setRateDiscrete <- function(volumes, categories=10) {
  volumes$eventName <- discretize(volumes$rateChange, method="frequency", categories = 10, ordered=F)
  volumes$eventName <- gsub(',', ' to ', volumes$eventName)
  volumes$eventName[!(is.na(volumes$eventName))]  <- paste0('tumorVolRateChange:', volumes$eventName[!(is.na(volumes$eventName))] )
  return(volumes)
}
#' \code{setBaselineDiscrete}: convert baseline volumes, measured in original units, to intervals.
#' @rdname setContinDiscrete
#' @export
setBaselineDiscrete <- function(volumes, categories=10) {
  volumes$baselineSize <- as.double(volumes$baselineSize)
  volumes$eventName <- discretize(volumes$baselineSize, method="frequency", categories = 10, ordered=F)
  volumes$eventName <- gsub(',', ' to ', volumes$eventName)
  volumes$eventName[!(is.na(volumes$eventName))]  <- paste0('tumorVolBaselineSize:', volumes$eventName[!(is.na(volumes$eventName))] )
  return(volumes)
}
