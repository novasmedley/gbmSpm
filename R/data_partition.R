#
# Prepate data splits in data partitions (training and testing partitions) 
# and cross-validation (training and validation folds).
#

#' Partition patients into testing and training datasets based on mean survival and gender.
#'
#' @details Patients are separated in two groups, patients with lower than mean overall survival
#' and patients with higher or equal to mean overall survival. Each group is
#' partitioned into 75\% training and 25\% testing using gender stratification. The training
#' and testing partitions have similar gender distributions within the two groups.
#' 
#' The training partitions from the low survival group is combined with the high survival group.
#' The testing partitions are also combined in this way.
#' Patients in training are not in testing.
#' 
#' Note: this function calls \code{\link{getSurvData}} to obtain overall survival days since the parameter \code{data}
#' only contains labels and features for classification.
#' 
#' Ratios between folds work better when number of folds are higher (e.g., 10 versus 2).
#'
#' @param data a dataframe holding the event data ready for logit modeling, where each row is an event/clinical visit and the columns
#' contain features of the event and the labels. It must at least contain a sample \code{id} column in the form "patientID.eventID".
#' @param verbose boolean, True for print, False for silence
#' @param seed int, seed for split
#' 
#' @param database a MySQL database, see \code{\link{getData}}
#' @param personTable a patient-based table in \code{database}, see \code{\link{getData}}
#' @param survData if not calling database, provide the survival data, see \code{\link{getSurvData}}
#' 
#' @return The eventIDs of the samples to put in training partition.
#' @export
#'
#' @examples
#' data("fake_data")
#' seed <- 1
#' fake_demo <- fake_data$demo
#' survData <- fake_data$person
#' fake_data$events <- cleanData(fake_data$events, tType = 'c') # pool of visits considered
#' 
#' fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
#' fake_data <- prepDemographics(fake_data, fake_demo) # need gender info
#' fake_data <- prepSurvivalLabels(fake_data) # get survival labels, these also change
#' fake_data$id <- paste0(fake_data$iois,'.',fake_data$eventID) # clinical ids
#' 
#' train.ids <- getTrainTestPartition(data=fake_data,
#'                                    database=NULL, 
#'                                    personTable=NULL,
#'                                    survData=survData,
#'                                    seed=seed,
#'                                    verbose=T)
#'                                    
getTrainTestPartition <- function(data, database=NULL, personTable=NULL,verbose=TRUE, seed, survData=NULL) {
  # extract information from data
  splitInfo <- unique(data[, c('iois','gender')])

  # get overall survival days from database to calculate mean survival days of patients in data
  if (is.null(database)){
    person <- survData
  } else {
    person <- getSurvData(database, personTable)
  }
  cat('\n...','unique iois in data (passed in): ', length(unique(data$iois)))
  
  person <- person[person$iois %in% splitInfo$iois,]
  meanSurv <- mean(person$survival)
  person$meanSurv <- ifelse(person$survival<meanSurv, yes=1, no=0)

  # use gender and mean survival to split data
  person <- merge(person,splitInfo)
  person$meanSurv_gender <- paste(person$meanSurv, person$gender)
  lowSurv <- person[person$meanSurv==1,]
  hiSurv <- person[person$meanSurv==0,]

  set.seed(seed)
  lowSurv_parts <- createDataPartition(lowSurv$gender, p=0.75, times=1, list=F)
  lowSurv_train <- lowSurv[lowSurv_parts,]
  lowSurv_test <- lowSurv[-lowSurv_parts,]
  cat('\n...','gender distribution of  people (with OS < mean OS) in training partition: ')
  print(table(lowSurv_train$meanSurv_gender))
  cat('\n...','gender distribution of  people (with OS < mean OS) in testing partition: ')
  print(table(lowSurv_test$meanSurv_gender))

  set.seed(seed)
  hiSurv_parts <- createDataPartition(hiSurv$gender, p=0.75, times=1, list=F)
  hiSurv_train <- hiSurv[hiSurv_parts,]
  hiSurv_test <- hiSurv[-hiSurv_parts,]
  cat('\n...','gender distribution of people (with OS > mean OS) in training partition: ')
  print(table(hiSurv_train$meanSurv_gender))
  cat('\n...','gender distribution of  people (with OS > mean OS) in testing partition: ')
  print(table(hiSurv_test$meanSurv_gender))

  trainPeople <- rbind(lowSurv_train,hiSurv_train)
  trainIDs <- data[data$iois %in% trainPeople$iois,]$id
  testIDs <- data[!(data$iois %in% trainPeople$iois),]$id
  
  # print some basic info to log
  if (verbose) {
    testPeople <- rbind(lowSurv_test,hiSurv_test)
    trainEvents <- data[data$iois %in% trainPeople$iois,]
    testEvents <- data[data$iois %in% testPeople$iois,]
    cat('\n... num people in training partition: ', nrow(trainPeople))
    cat('\n... num people in testing partition: ', nrow(testPeople))
    cat('\n... num events in training partition: ', nrow(trainEvents))
    cat('\n... num events in testing partition: ', nrow(testEvents))
    cat('\n')
    cat('\n... Which patients are in testing and training partitions?: ')
    which(testPeople$iois %in% trainPeople$iois==TRUE)
    cat('\n')
    
    cat('\n... Event label distribution in:')
    labelNames <- getClassLabels()
    z <- lapply(labelNames, function(label) {
      cat('\n...... ',label,' training and testing partition, respectively:')
      print(table(trainEvents[ , label]))
      print(table(testEvents[ , label]))
    })
  }
  return(list(train=trainIDs, test=testIDs))
}

#' Create patient-independent folds in cross-validation in training partition.
#'
#' @details Patients in validation folds will not have events in training folds 
#' and vice versa.
#'
#' @param trainEvents a dataframe holding the event data ready for logit modeling, 
#' where each row is an event/clinical visit and the columns
#' contain features of the event and the labels. It must contain an \code{id} column 
#' in the form "patientID.eventID".
#' @param folds number of folds
#' @param seed int, seed for split
#' @param verbose logical, True for print, False for silence
#'  
#' @return \code{folds} a list of arrays indicating row 
#' position integers corresponding to fold split
#' 
#' @seealso \code{\link[caret]{groupKFold}},
#' \code{\link{getPatIDs}},
#' \code{\link{getClassLabels}},
#' \code{\link{getData}}
#' 
#' @export
#' 
#' @examples
#' 
#' # use training partition to create folds for CV
#' data("features_ratechange_sup0.4g60l2z2") # features and labels for each clinical visit
#' t <- 'rate'
#' maxgap <- 60
#' maxlen <- 2
#' names <- colnames(feats)
#' feats <- data.frame(id=row.names(feats),feats)
#' colnames(feats) <- c('id',names)
#' feats <- prepLaterality(feats)
#' feats <- prepLocation(feats)
#' feats <- removeVisits(feats,
#'                      maxgap=maxgap,
#'                      maxlength=maxlen,
#'                      tType=t,
#'                      save=F,
#'                      outDir=NA)
#' labels <- getClassLabels()
#' needToRemove <- c('id','iois','eventID', # remove ids
#'                   labels, # remove labels
#'                   'IDH1') # not interested
#' 
#' # data partitions
#' train.ids <- sample(feats$id, size=floor(0.80*nrow(feats)), replace = F) # random
#' feats <- feats[feats$id %in% train.ids,] #training data
#' ind <- getTrainingFolds(trainEvents=feats,
#'                         folds=5,
#'                         seed=1,
#'                         verbose=T)
#'                         
getTrainingFolds <- function(trainEvents, folds, seed, verbose) {
  iois <- getPatIDs(trainEvents$id)
  set.seed(seed)
  folds <- groupKFold(iois, k=folds)
  
  if (verbose) {  
    cat('\n...ratio of alive:dead (assumed order in label column) events in CV folds: \n')
    labelNames <- getClassLabels()
    z <- lapply(labelNames, function(label) {
      ratio <- lapply(folds, function(fold) { 
        t <- table(trainEvents[, label ][fold]); as.numeric(t[1]/t[2])
      })
      cat('......', label, ':',  round(unlist(ratio),digits=2), '\n')
    })

  }
  return(folds)
}

