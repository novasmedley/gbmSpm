#
# Calls to pull information from database, where authentication and connection details
# are in a local configuation file.
#

#' Get patient event and survival data from a database
#' 
#' @details Note: tumor lateriality and location is searched for in agent_detail regardless of agent/event type.
#'
#' @param database a MySQL database, where connection is
#' authenticated through a configuration file with group named \code{database})
#' @param eventTable  an event-based table in \code{database}, containing at least:
#' \itemize{
#' \item iois = patient id, char
#' \item agent = event type, char
#' \item agent_detail = event attribute, char
#' \item startdate = event start datein the format 'YYYY-MM-DD', char
#' \item enddate = event end date in the format 'YYYY-MM-DD', char
#' \item recurrenceNum = patient's recurrence count (since baseline) during event, int
#' \item daysSinceBaseline = days since baseline, int
#' }
#' @param personTable a patient-based table in \code{database}, containing at least:
#' \itemize{
#' \item iois = patient id, char
#' \item survival = total survival days at time of data collection, int
#' \item status = patient status at the end of data collection , 'DECEASED' or 'ALIVE', char
#' \item DOB = date of birth in the format 'YYYY-MM-DD', char
#' \item gender = gender, char
#' \item ethnicity = ethnicity, char
#' \item MSP = MGMT promoter methytation status, char ('U' = unmethylated, 'M' = methylated)
#' \item IDH1 = status of IDH1 gene, char ('WT' = wild type, 'MUT' = mutated)
#' }
#'
#' @return \code{getData} returns a named list containing two dataframes, 'person' for patient data 
#' (returned by \code{getSurvData}) and 'event' for event data (returned by \code{getEventData})
#' @export
getData <- function(database, eventTable, personTable){
  event <- getEventData(database, eventTable)
  person <- getSurvData(database, personTable)
  return(list(person=person, event=event))
}



#' \code{getEventData}: get event data from \code{eventTable} in database
#' @rdname getData
#' @export
getEventData <- function(database, eventTable){
  conn <- dbConnect(MySQL(), group=database)
  
  query <- paste("select iois, agent, agent_detail, startdate, enddate, recurrenceNum from", eventTable,";")
  cat('...', 'ran query: "', query, '"\n')
  RS <- dbSendQuery(conn, query)
  event <- fetch(RS, n=-1) 
  event$agent_detail <- unlist(event$agent_detail)
  
  z <- dbClearResult(dbListResults(conn)[[1]])
  z <- dbDisconnect(conn)
  
  cat('...', nrow(event), ' total intial events called from database')
  return(event)
}


#' \code{getSurvData}: get patient survival days and survival status from \code{personTable} in database
#' @rdname getData
#' @export
getSurvData <- function(database, personTable){
  conn <- dbConnect(MySQL(), group=database) 
  
  query <- paste("select iois, survival, status from", personTable,";")
  cat('...', 'ran query: "', query, '"\n')
  RS <- dbSendQuery(conn, query)
  person <- fetch(RS, n=-1)
  
  z <- dbClearResult(dbListResults(conn)[[1]])
  z <- dbDisconnect(conn)
  
  cat('...', nrow(person), ' total people called from database')
  return(person)
}


#' \code{getDemographics}: get patient demographics and two biomarkers from \code{personTable} in database
#' @rdname getData
#' @export
getDemographics <- function(database, personTable) {
  conn <- dbConnect(MySQL(), group=database) # authentication made through configuration file
  
  query <- paste('SELECT iois, DOB, gender, ethnicity, MSP, IDH1 FROM', personTable,';')
  cat('...', 'ran query: "', query, '"\n')
  
  RS <- dbSendQuery(conn, query)
  info <- fetch(RS, n=-1) 
  
  z <- dbClearResult(dbListResults(conn)[[1]])
  z <- dbDisconnect(conn) # prevent an output when disconnecting
  
  return(info)
}
