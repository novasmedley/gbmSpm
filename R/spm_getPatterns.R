#
# 1. obtain sequential patterns from SPM and find the days when pattern occurred for each patient
# 2. match static data with SPM patterns to get feature vectors for modeling
#


#' Get patient IDs.
#' 
#' @details Given a string in the form "patientID.eventID", return just the patientID.
#' @param ids an array of ids in the form "patientID.eventID"
#' @return an array of patientIDs
#' @export
#' @examples
#' data("features_ratechange_sup0.4g60l2z2")
#' feats$patID <- getPatIDs(rownames(feats))

getPatIDs <- function(ids){
  patIDs <- lapply(as.character(ids), function(n) unlist(strsplit(n,'\\.'))[1])
  return(unlist(patIDs))
}


#' Convert events into transactions.
#' 
#' @details Creates a transaction file with the following 'n' columns: \itemize{
#' \item column 1: sequenceID
#' \item column 2: timeStamp/eventID
#' \item column 3 = number of items
#' \item columns 4-n = items in event
#' }
#' 
#' Transaction file is to be read by \code{\link[arulesSequences]{read_baskets}}.
#' 
#' @param event dataframe of events,
#' where rows are single events, and
#' least the following columns:
#' \itemize{
#' \item \code{iois} patient id
#' \item \code{eventID} event id, see \code{\link{resetEventIDsToDays}}
#' \item \code{eventName} unique event name, see \code{\link{createEventNames}}
#' }
#' @param filename path to save transactions file
#'  
#' @export
#' @examples
#' data("fake_data")
#' outputDir <- '~/test'
#' 
#' # save tumor location and laterility strings before event cleaning
#' fake_tumorInfo <- fake_data$events 
#' fake_demo <- fake_data$demo
#' 
#' # clean data
#' fake_data$events <- cleanData(fake_data$events, tType = 'rate')
#' cat('...',nrow(fake_data$events), " events left for SPM after cleaning", '\n')
#' 
#' # collect patient info for each event
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T) 
#' 
#' # prep for each event, since age does change
#' # get survival labels, these also change
#' fake_data <- prepDemographics(fake_data, fake_demo) 
#' fake_data <- prepSurvivalLabels(fake_data) 
#' 
#' # get first tumor location
#' fake_data <- getTumorLocation(fake_data, fake_tumorInfo) 
#' 
#' # transactions
#' trans <- createTransactions(event = fake_data, 
#'                             filename = file.path(outputDir, 'example_transactions.txt'))

createTransactions <- function(event, filename) {
  cat('\n\n','Creating transactions file...\n\n')
  
  event$iois <- as.integer(event$iois)
  event <- event[order(event$iois,event$eventID),] # integer identifies must be ordered
  
  tempa <- event[, c('iois','eventID','eventName')] # treat each person as a sequence
  tempa <- dlply(tempa, .(iois,eventID), c)
  
  if (file.create(filename)) file.remove(filename) # delete old file
  
  fileConn <- file(filename, "w")
  for (i in 1:length(tempa)) {
    names <- as.character(tempa[[i]]$eventName)
    x <- c(tempa[[i]]$iois[1], tempa[[i]]$eventID[1], length(names),lapply(gsub(';','\t', names), function(x) x)) # format transactions
    write.table(as.data.frame(x), file=filename, sep='\t', append=T, quote=F, col.names=F, row.names=F) 
  }
  close(fileConn)
  
  cat('...', 'file written to: ', filename)
}

#' Get sequential patterns
#'
#' @inheritParams createTransactions
#' @param transFilename asdf
#' @param createT logical, create or load previously created transactions
#' @param support see slot in \code{\link[arulesSequences]{SPparameter-class}}
#' @param maxgap see slot in \code{\link[arulesSequences]{SPparameter-class}}
#' @param maxlen see slot in \code{\link[arulesSequences]{SPparameter-class}}
#' @param maxsize see slot in \code{\link[arulesSequences]{SPparameter-class}}
#' 
#' @return a named list containing: \itemize{
#' \item freqseq, the output from \code{\link{arulesSequences::cspade}}
#' \item data, the input \code{data}
#' }
#' 
#' @seealso \code{\link{createTransactions}}
#' 
#' @export
#' @examples
#' data("fake_data")
#' outputDir <- '~/test'
#' 
#' # save tumor location and laterility strings before event cleaning
#' fake_tumorInfo <- fake_data$events 
#' fake_demo <- fake_data$demo
#' 
#' # clean data
#' fake_data$events <- cleanData(fake_data$events, tType = 'rate')
#' cat('...',nrow(fake_data$events), " events left for SPM after cleaning", '\n')
#' 
#' # collect patient info for each event
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T) 
#' 
#' # prep for each event, since age does change
#' # get survival labels, these also change
#' fake_data <- prepDemographics(fake_data, fake_demo) 
#' fake_data <- prepSurvivalLabels(fake_data) 
#' 
#' # get first tumor location
#' fake_data <- getTumorLocation(fake_data, fake_tumorInfo) 
#' 
#' # spm
#' pSPM <- getSeqPatterns(event = fake_data,
#'                        transFilename = file.path(outputDir, 'example_transactions.txt'),
#'                        createT = T,
#'                        support = 0.4,
#'                        maxgap = 60,
#'                        maxlen = 2,
#'                        maxsize = 2)

getSeqPatterns <- function(event, transFilename, createT, 
                           support, maxgap, maxlen, maxsize) {
  
  if (createT) {
    # takes a few seconds - only eventName, iois, and eventID columns are used for SPM
    createTransactions(event, filename=transFilename) 
  } else {
    cat('\n\n','Reading from transactions file:  ', transFilename)
  }
  cat('\n\n','Performing SPM...\n\n')
  trans <- read_baskets(con=transFilename, 
                        sep="[\t]+", 
                        info=c("sequenceID","eventID","SIZE")) # persons as sequences
  freqseq <- cspade(trans, 
                    parameter = list(support=support, 
                                     maxgap=maxgap, 
                                     maxlen=maxlen, 
                                     maxsize=maxsize), 
                    control = list(verbose=TRUE, tidLists=T)) # get frequent sequences
  return(list(freqseq=freqseq, data=event))
}

#' Create a binary vector of patterns observed at a clinical visit.
#' 
#' @details For every clinical visit, each pattern has a value of 
#' 1 (present) or 0 (absent).
#' 
#' @param patternDays list of patterns found by cSPADE, 
#' for each pattern there is a list of patientIDs (see \code{\link{getPatIDs}}), 
#' for each patientID there is a list of eventIDs when the pattern was observed
#' 
#' @return dataframe, 
#' rows of clinical visits, 
#' cols of binary features (0: observed, 1: not observed)
#' 
#' @export
#' 
#' @examples
#' data("patterns_ratechange_sup0.4g60l2z2") # days when pattern occurred
#' features <- vectorizePatterns(patternDays)

vectorizePatterns <- function(patternDays) {
  # convert to list of patterns,
  # where each pattern has a list of clinical visit ids (patientID.eventID),
  # still shows when the pattern occurred
  patternDays <- lapply(patternDays, function(t) { 
    t <- t[!(is.na(t))]
    p <- lapply(seq(1,length(t)), function(i) paste0(names(t)[i],'.',t[[i]]))
    unique(unlist(p))
  })
  
  # for every clinical visit, set 1 or 0 for when pattern occurs
  g <- unique(unlist(patternDays)) # get unique clinical visits
  z <- rep(0,length(g))
  patternDays <- lapply(patternDays, function(p) {
    b <- g %in% p
    z[b] <- 1
    z
  })
  patternDays <- do.call('cbind', patternDays)
  rownames(patternDays) <- g
  return(patternDays)
}


#' Get column names of labels and non-spm features to potentially include in logit fitting data.
#' 
#' @return char, names of features and labels to look for in fitting data
#' @export
#' @examples
#' nonPatternColNames <- c('iois','eventID', # identification
#'                         getDemoNames(),
#'                         getBiomarkerNames(), 
#'                         getTumorNames(),
#'                         getClassLabels())
getTumorNames <- function() {
  n <- c('frontal','temporal', 'occipital','parietal','cerebellar', 'thalamic', 'sellar', 
         'corpus callosum', 'pineal', 'midbrain')
  n <- paste0('location: ', n)
  n <- c('tumorLaterality',n)
}

#' @rdname getTumorNames
#' 
#' 
#' @export
getDemoNames <- function() {
  d <- c('gender','ethnicity','ageDecade')
}

#' @rdname getTumorNames
#' 
#' 
#' @export
getBiomarkerNames <- function() {
  b <- c('MSP')
}

#' @rdname getTumorNames
#' 
#' 
#' @export
getClassLabels <- function() {
  l <- c('survivalIn60', 'survivalIn180', 'survivalIn270', 'survivalIn360')
}


#' Create feature vectors that can be used for logit modeling
#' 
#' @inheritParams vectorizePatterns
#' @seealso \code{\link{vectorizePatterns}}
#' @param event dataframe, 
#' rows are single events (used as input to cSPADE), 
#' columns are event details plus patient demogrpahics, 
#' tumor laterality, survival labels and MGMT biomarker
#' 
#' @return dataframe, where rows are clinical visits, 
#' and columns are features of the visit that can be used for logit modeling:
#' binary temporal (cSPADE patterns) events, demographic data, biomarkers, and tumor lateriality.
#' 
#' @export
#' @examples
#' data("fake_data")
#' outputDir <- '~/test'
#' 
#' # save tumor location and laterility strings before event cleaning
#' fake_tumorInfo <- fake_data$events 
#' fake_demo <- fake_data$demo
#' 
#' # clean data
#' fake_data$events <- cleanData(fake_data$events, tType = 'rate')
#' cat('...',nrow(fake_data$events), " events left for SPM after cleaning", '\n')
#' 
#' # collect patient info for each event
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T) 
#' 
#' # prep for each event, since age does change
#' # get survival labels, these also change
#' fake_data <- prepDemographics(fake_data, fake_demo) 
#' fake_data <- prepSurvivalLabels(fake_data) 
#' 
#' # get first tumor location
#' fake_data <- getTumorLocation(fake_data, fake_tumorInfo) 
#' 
#' # spm
#' pSPM <- getSeqPatterns(event = fake_data,
#'                        transFilename = file.path(outputDir, 'example_transactions.txt'),
#'                        createT = T,
#'                        support = 0.4,
#'                        maxgap = 60,
#'                        maxlen = 2,
#'                        maxsize = 2)
#' pSPM$patterns <- as(pSPM$freqseq, "data.frame")
#' pSPM$patterns$sequence <- as.character(pSPM$patterns$sequence)
#' 
#' # days when pattern occur
#' patternDays <- findPatternDays(pSPM$patterns, pSPM$data, maxgap=60)
#' 
#' # feature vectors to supply to logits
#' feat_vecs <- getFeatureVectors(patternDays, events=pSPM$data)

getFeatureVectors <- function(patternDays, events) {
  # prep static features
  # get their names
  nonPatternColNames <- c('iois','eventID', # identification
                          getDemoNames(),
                          getBiomarkerNames(), 
                          getTumorNames(),
                          getClassLabels())
  
  nonPatternCols <- events[,nonPatternColNames]
  
  features <- vectorizePatterns(patternDays)
  
  #-------- create label and demographics info for input data
  # get clinical visits
  info <- row.names(features)
  info <- lapply(info, function(n) unlist(strsplit(n,'\\.')))
  info <- lapply(info,function(i) as.integer(i))
  info <- do.call('rbind',info)
  info <- as.data.frame(info)
  colnames(info) <- c('iois','eventID')
  info <- as.data.frame(info) 
  

  #-------- attach features to clinical visits
  # get overall list of visits 
  nonPatternCols[,1:2] <- lapply(nonPatternCols[,1:2],function(i) as.integer(i))
  nonPatternCols <- nonPatternCols[!duplicated(nonPatternCols[,1:2]),]
  
  # attach non-pattern features
  info <- merge(info,nonPatternCols,by=c('iois','eventID'),all=T) 
  
  # # see events which did not pass support threshold and therefore not seen in patternDays; cSPADE filtered it out
  # info[!(info$eventID %in% t$eventID),]
  # t[!(t$eventID %in% info$eventID),]
  # see <- pSPM$data[pSPM$data$eventID %in% t$eventID[!(t$eventID %in% info$eventID)],]
  
  # attach pattern features
  features <- as.data.frame(features)
  id <- strsplit(rownames(features), '\\.') # get id info
  id <- lapply(id, function(i) as.integer(i))
  features <- cbind(iois = unlist(lapply(id,'[[',1)), 
                    eventID = unlist(lapply(id,'[[',2)), features)
  
  # include visits when patterns were not observed!!
  features <- merge(info,features, by=c('iois','eventID'), all=T)
  
  # set 0 for visits where no patterns from spm were observed 
  # (e.g., if pattern was below support)
  cat('\n... setting NA values to 0 -- will return warning if non-spm features are NA beforehand')
  features[is.na(features)] <- 0 
  
  # sort so data splits will be the same across models
  features <- features[with(features, order(iois,eventID)),] 
  
  row.names(features) <- lapply(seq(1,nrow(features)), function(i) {
    paste0(features$iois[i],'.',features$eventID[i])
  })
  features <- features[,!(colnames(features) %in% c('iois','eventID'))]
  
  return(features)
}

#' Run sequential pattern mining
#' 
#' @details This will automatically create three types of files: \itemize{
#' \item transactions - returned from \code{\link{createTransactions}}, saved in .txt
#' \item patterns - returned from \code{\link{findPatternDays}}, saved in .rds 
#' \item features - returned from \code{\link{getFeatureVectors}}, saved in .rds 
#' }
#' 
#' @inheritParams getVolTypeName
#' @inheritParams getSeqPatterns
#' @inheritParams createTransactions
#' 
#' @param suppList one or more supports to be used as a parameter, see slot in
#' code{\link[arulesSequences]{SPparameter-class}}
#' @param outputDir full path to output folder, e.g., where to save spm results
#' 
#' @seealso \code{\link{getSeqPatterns}}, \code{\link{getVolTypeName}}
#' 
#' @export
#' @examples
#' data("fake_data")
#' outputDir <- '~/test'
#' 
#' # save tumor location and laterility strings before event cleaning
#' fake_tumorInfo <- fake_data$events 
#' fake_demo <- fake_data$demo
#' 
#' # clean data
#' fake_data$events <- cleanData(fake_data$events, tType = 'rate')
#' cat('...',nrow(fake_data$events), " events left for SPM after cleaning", '\n')
#' 
#' # collect patient info for each event
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T) 
#' 
#' # prep for each event, since age does change
#' # get survival labels, these also change
#' fake_data <- prepDemographics(fake_data, fake_demo) 
#' fake_data <- prepSurvivalLabels(fake_data) 
#' 
#' # get first tumor location
#' fake_data <- getTumorLocation(fake_data, fake_tumorInfo) 
#' 
#' runSPM(event=fake_data,
#'        suppList=suppList, 
#'        maxgap=60, 
#'        maxlen=2, 
#'        maxsize=2, 
#'        tType='rate', 
#'        outputDir='~/test')
#' 
runSPM <- function(event, suppList, maxgap, maxlen, maxsize, tType, outputDir,createT ) {
  tp <- getVolTypeName(tType)
  cat('\n\n','Using spm patterns with', tp,'...\n')
  
  #----- find patterns with the lowest min supp #
  minSupp <- min(suppList)
  subDir <- paste0('sup',minSupp,'g',maxgap,'l',maxlen,'z',maxsize)
  spmDir <- file.path(outputDir,subDir)
  ifelse(!dir.exists(spmDir), dir.create(spmDir), FALSE)
  
  # note that seq. patterns already use discretized
  fn <- file.path(outputDir,paste0('transactions_',tp,'.txt'))
  if (file.exists(fn)) {
    cat('...using file that already exists: ', fn, '\n')
    ct <- F
  } else {
    cat('...creating: ', fn, '\n')
    ct <- T
  }
  pSPM <- getSeqPatterns(event = event,
                         transFilename = fn, 
                         createT = ct,
                         support = minSupp, 
                         maxgap = maxgap, 
                         maxlen = maxlen, 
                         maxsize = maxsize)
  pSPM$patterns <- as(pSPM$freqseq, "data.frame")
  pSPM$patterns$sequence <- as.character(pSPM$patterns$sequence)
  
  cat('\n\n','Getting spm events data ready for logit...\n')
  cat('...creating pattern and feature files in: ', spmDir, '\n')
  # Find person with the Event
  # slow because repeating search of items for each itemset in pattern
  patternDays <- findPatternDays(pSPM$patterns, pSPM$data, maxgap) 
  saveRDS(patternDays, file=file.path(spmDir, paste0('patterns_',getVolTypeName(tType) ,'.rds') ) )
  
  features <- getFeatureVectors(patternDays, events=pSPM$data)
  saveRDS(features, file=file.path(spmDir, paste0('featureVectors_',getVolTypeName(tType) ,'.rds') ) )
  
  #----- find patterns with other values of min supp in suppList #
  if (length(suppList >1)) {
    cat('\n\n','...subset data based on different support thresholds...\n')
    temp <- lapply(suppList[-1], function(supp) {
      cat('......supp: ', supp, '\n')
      p <- pSPM$patterns$sequence[pSPM$patterns$support >= supp] # subset of patterns
      days <- patternDays[p] # days when only subset patterns occurs
      feats <- getFeatureVectors(days, pSPM$data) # convert to feature vector representation
      
      subDir <- paste0('sup',supp,'g',maxgap,'l',maxlen,'z',maxsize)
      spmDir <- file.path(outputDir,subDir)
      ifelse(!dir.exists(file.path(outputDir,subDir)), dir.create(file.path(outputDir,subDir)), FALSE)
      ifelse(!dir.exists(spmDir), dir.create(spmDir), FALSE)
      saveRDS(days, file=file.path(spmDir, paste0('patterns_', tp ,'.rds') ) )
      saveRDS(feats, file=file.path(spmDir, paste0('featureVectors_',tp,'.rds') ) )
    })
  }
}
