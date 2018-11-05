# These functions provide additional parsing of data coming out of arules and arulessequences packages.


# Takes in a formal class sequences, and an index to identify which sequence to find IDs for. This function taks advantage of
# the information within the arulessequences package, tiDLists.
# Returns a list indicating which patient has/doesn't have the sequence.
getSequenceID <- function(sequences, index) {
  # get the tidList information and corresponding sequenceID order
  tIDs <- as(sequences@tidLists[index], 'list')
  itemNames <- names(tIDs)
  seqIDs <- sequences@tidLists@transactionInfo$sequenceID 
  
  # convert tIDs to sequenceIDs
  tIDs <- lapply(tIDs, function(x) seqIDs[x])
  
  return(tIDs)
}

# Takes in a frequent sequence output by cspade and converts the sequence into readable itemsets 
# that can be read by the function getEventDayOfPattern.
# returns itemsets as list information
parseSequence <- function (sequence) {
  str <- substring(sequence,2,nchar(sequence)-1)
  itemsets <- regmatches(str, gregexpr('\\{.*?\\}', str))
  itemsets <- lapply(itemsets, function(x) substring(x,2,nchar(x)-1))
  itemsets <- itemsets[[1]] # remove upper level list
  
  return(itemsets)
}

# Takes in a single patient's set of events, and a parsed sequence containing itemset information.
# returns day of first element each time pattern is found in patient sequence, or NA if not found
getLastDayofEachPattern <- function(patEvents, itemsets, maxgap) {
  
  # get a list of dates where events match each itemset in pattern/sequence
  hasItemsets <- lapply(itemsets, function(x) {
    s <- unlist(strsplit(x,','))
    
    # get a list of start dates with event for each event in s
    l <- lapply(s, function(y) patEvents$eventID[grepl(y, patEvents$eventName, fixed=T)==T] )
    b <- Reduce(intersect, l) # check if items are on the same day for each itemset, return dates where this is true
  })
  
  # check if all itemset dates are within maxgap of each other
  
  if (length(hasItemsets[[1]])==0) { return (NA) # pattern not found in patient events
  } else if (length(hasItemsets)==1) { return (hasItemsets[[1]]) # only one itemset in sequenece, and has days when it occured (ie, found)
  } else { 
    matches <- lapply(hasItemsets[[1]], function(currentDays){  # for each possible start date of this pattern
      z <- 2
      while (z <= length(hasItemsets) ) { # check if the rest of the itemset dates match pattern
        nextDays <- hasItemsets[[z]]
        matchDays <- c()
        
        temp <- lapply(currentDays, function(d) { # currentDays are matches from the prior itemset days, and are days to check the next days against
          nextDays <- nextDays[nextDays > d] # only check events that occur AFTER a current day
          if (!(length(nextDays)==0)) {
            matched <- nextDays[nextDays-d <= maxgap]
            if (!(length(matched)==0)) matchDays <<- c(matchDays,matched) # a running list of days that match pattern for the days in currentDays
          }
        })
        rm(temp)
        
        if ( length(matchDays)>0 ) {
          if (z==length(hasItemsets)) {
            return(matchDays)  # reached the end of itemsets and found match(es)
          } else {
            currentDays <- matchDays # found match and still have more to check, so set the next days to check as currentDays
            z <- z+1
          }
        } else { 
          return(NA) # itemset not found in current chain of days
        }
      }
    })
    
    patternEndDates <- unlist(matches)
    patternEndDates <- patternEndDates[!is.na(patternEndDates)]
    if (length(patternEndDates)==0) patternEndDates <- NA
    
    return(patternEndDates)
  }
  
}

#' Find People with a sequence pattern in their events
#' @return returns a list of patterns, where each pattern has a list of patients, and each patient has a list of days when event occurred (eventID)
#' @export
#' @examples
#' data("fake_data")
#' fake_tumorInfo <- fake_data$events # save tumor location and laterility strings before event cleaning
#' fake_demo <- fake_data$demo
#' fake_data$events <- cleanData(fake_data$events, tType = 'rate')
#' cat('...',nrow(fake_data$events), " events left for SPM after cleaning", '\n')
#' 
#' # collect patient info for each event
#' # note that only eventName, iois, and eventID columns are used for SPM
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
#' fake_data <- prepDemographics(fake_data, fake_demo) # prep for each event, since age does change
#' fake_data <- prepSurvivalLabels(fake_data) # get survival labels, these also change
#' fake_data <- getTumorLocation(fake_data, fake_tumorInfo) # get first tumor location
#' 
#' # pattern
#' pSPM <- getSeqPatterns(data = fake_data,
#'                        transFilename = 'example_transactions.txt',
#'                        createT = T,
#'                        supp = 0.2,
#'                        maxgap = 60,
#'                        maxlen = 2,
#'                        maxsize = 2)
#' pSPM$patterns <- as(pSPM$freqseq, "data.frame")
#' pSPM$patterns$sequence <- as.character(pSPM$patterns$sequence)
#' 
#' # days when pattern occur
#' patternDays <- findPatternDays(pSPM$patterns, pSPM$data, maxgap=60) 
#' 
findPatternDays <- function(patterns, events, maxgap) {
  test1 <- T # automatic inspection
  test2 <- F # manual inspection
  
  # Get the day when the pattern first occurred, the day should be the last eventID of the item(set) in pattern
  cat('\n\n','Searching for sequence patterns in events...\n\n')
  
  numOverallSequences <- length(unique(events$iois)) # the total set of sequences (here, each person is a sequence)
  events <- split(events, events$iois)
  
  # TO-DO: build hash table instead of looking up each itemset in a patient every time
  start.total <- Sys.time()
  count <- 1
  patternsDays <- lapply(patterns$sequence, function(seq) {
    itemsets <- parseSequence(seq)
    eventDay <- lapply(events, function(pe) getLastDayofEachPattern(pe, itemsets, maxgap))
    count <<- count + 1    
    if(count%%10==0) cat(count,'...',(Sys.time() - start.total),'\n')
       
    eventDay
  })
  names(patternsDays) <- patterns$sequence
  cat('... finished search...\n')
  print(Sys.time() - start.total)
  
  # TO-DO: fix testing to take in a subset of frequent sequences if they are too long.
  
  #  begin --- testing, TESTS FOR getLastDayOfFirstPattern function
  # Testing: check if support values match patients found with itemset pattern of frequent sequences
  
  if (test1==T) {
    cat('... testing...\n')
    
    # my results
    myResults <- lapply(patternsDays, function(pat) { pBools <- ifelse(is.na(pat), 0,1 ) })
    myResults <- as.data.frame(do.call(rbind, myResults))
    mySupports <- unlist(apply(myResults, 1, function(c) sum(c)))/numOverallSequences
    
    compareResults <- lapply(seq(1:nrow(patterns)), function(i) { ifelse (mySupports[i]==patterns$support[i],'yes','no') })
    unmatchedSeqs <- patterns[which(compareResults=='no'),]
    cat('...','number of unmatched frequent sequences support FROM TIDLIST (my results vs. arulesSequences, should match): ', nrow(unmatchedSeqs))
  }
  
# test 2 has not been adjusted
  if (test2==T) {
    cat('... testing...\n')
    # ind <- 17
    # e <- freqseq_df$sequence[ind] # testing pattern
    # e
    # eTIDs <- getSequenceID(freqseq, ind) # trans ID list containing pattern
    # notE <- patient$iois[!(patient$iois %in% eTIDs[[1]])] # patients without this pattern
    # # e <- '<{Radiation}>'
    # temp <- events[events$iois==1687,]
    # getEventDayOfPattern(temp, parseSequence(e))
    # patEvents <- temp
    # itemsets <- parseSequence(e)
    
    temp <- split(events, events$iois)
    freqseq_df$found <- lapply(seq(1,nrow(freqseq_df)), function(i) {
      e <- freqseq_df$sequence[i]
      r <- unlist(lapply(temp, function(p) getEventDayOfPattern(p, parseSequence(e)))) # returns if patients have the sequence or not
      r_pos <- names(r[r!=-1])
      
      # testing
      r_neg <- names(r[r==-1])
      pos <- getSequenceID(freqseq, i) # trans ID list containing pattern
      neg <- patient$iois[!(patient$iois %in% pos[[1]])] # patients without this pattern
      
      falseNegative <- r_neg[!(r_neg %in% neg)]  # is actually positive but was found negative
      falsePositive <- neg[!(neg %in% r_neg)] # is actually negative but was found positive
      # testing
      
      supp <- length(r_pos)/numOverallSequences
    })
    
    a <- lapply(seq(1:nrow(freqseq_df)), function(i) {
      if (freqseq_df$found[i]==freqseq_df$support[i]) {'yes'
      } else 'no'
    })
    
    unmatchedSeqs <- freqseq_df[which(a=='no'),]
    cat('number of unmatched frequent sequences support found MANUALLY: ', nrow(unmatchedSeqs))
  }
  # end --- testing
  
  return(patternsDays)
}

