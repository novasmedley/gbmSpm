#' Command line script to run spm by first calling Agile database and creates the feature vectors for logit modeling.
#' The feature vectors include temporal patterns identified by SPM and static variables such as ethnicity.
#' 
#' returns nothing, but saves the patterns in .RData object 
#' - note that the patterns are the feature vectors used for modeling.
#' - the patterns are formated so that it is ready for a logit model


# get help: in terminal, "Rscript /PATHTO/example_spm.R -h", where PATHTO is the path to the example_spm.R

#-------------------------------------------------------------------------------- #
library(gbmSpm)
library(optparse)

#-------------------------------------------------------------------------------- #

# get cmd line arguments
option_list = list(make_option(c('--tType'), type='character', default=NULL, help='tumor volume type: 
                               p (percent change), 
                               c (continuous, original volume units), 
                               r (volumetric response criteria), or
                               rate (rate change),
                               see "getVolTypeName" in "cleanData.R"', metavar='character'),
                   make_option(c('--minSupp'), type='numeric', default=NULL, help='cSPADE, the min support to consider', metavar='character'),
                   make_option(c('--minSuppList'), type='character', default='yes', help='explore other min support values to prevent redundant work, see "runSPM"', metavar='character'),
                   make_option(c('--maxgap'), type='numeric', default=NULL, help='cSPADE max time between CONSECUTIVE elements in sequence', metavar='character'),
                   make_option(c('--maxlength'), type='numeric', default=NULL, help='cSPADE max length of elements in sequence', metavar='character'),
                   make_option(c('--maxsize'), type='numeric', default=NULL, help='cSPADE max number of items in an element of a sequence', metavar='character'),
                   make_option(c('--outDir'), type='character', default=NULL, help='full path to output directory', metavar='character')
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args((opt_parser))

# OR use run in RStudio using :

# opt<- list()
# opt$tType <- 'rate'
# opt$minSupp <- 0.4
# opt$minSuppList <- 'no'
# opt$maxgap <- 60
# opt$maxlength <- 2
# opt$maxsize <- 2
# opt$outDir <- '~/gbm_spm_example'

#--------------------------------------------------------------------------------  #

# constants

outputDir <- opt$outDir
spmDir <- file.path(outputDir,'spm')
logDir <- file.path(outputDir,'spm_logs')
lapply(c(outputDir, spmDir, logDir), function(i) ifelse(!dir.exists(i), dir.create(i, showWarnings = F, recursive = T), F) )


chemoOverlap <- F

# itemsets come from items in one transaction (so here, they occur on the same day), so items in a visit
# sequences are made of itemsets (also called elments)

suppList <- opt$minSupp
if (opt$minSuppList=='yes') suppList <- seq(opt$minSupp, 0.4, .05)  # percentage of patients/sequences with pattern

# MAX GAP:
# max time between CONSECUTIVE elements in sequence, i.e., max days between visits
# note, patients be be intermitten if considered 'well'
maxgap <- opt$maxgap

# MAX LEN:
# max length of elements in sequence
# --> read as a sequence of 3 visits, with maxgap days between the visits
maxlen <- opt$maxlen  

# MAX SIZE
# max number of items in an element of a sequence 
# --> choosen because usually only 5 items in neuro evaluations, and they occur on diff days than chemo and radiation
# BUT not really necessary, 5 neuro evaluations + imaging + radiation + chemo could be up to 8 items in one visit
maxsize <- opt$maxsize # could change up to 8...


#---------------------------------------- log data ---------------------------------------- 
logfn <- file(file.path(logDir,paste0(paste0('g',maxgap,'l',maxlen,'z',maxsize),'_',
                                      getVolTypeName(opt$tType),'.log')), open='wt')

sink(logfn, type='output', split=T)

cat(format(Sys.Date(), format="%Y.%m.%d"), '\n')
cat('maxgap: ', maxgap, '\n')
cat('maxlen: ', maxlen, '\n')
cat('maxsize: ', maxsize, '\n')
cat('tumor vol variable type: ', getVolTypeName(opt$tType), '\n\n')


#---------------------------------------- load and prep data ---------------------------------------- 

# load dummy data (no database to call)
data("fake_data")
fake_tumorInfo <- fake_data$events # save tumor location and laterility strings before event cleaning
fake_demo <- fake_data$demo

fake_data$events <- cleanData(fake_data$events, opt$tType)

cat('...',nrow(fake_data$events), " events left for SPM after cleaning", '\n')


# collect patient info for each event
# note that only eventName, iois, and eventID columns are used for SPM
fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
fake_data <- prepDemographics(fake_data, fake_demo) # prep for each event, since age does change
fake_data <- prepSurvivalLabels(fake_data) # get survival labels, these also change
fake_data <- getTumorLocation(fake_data, fake_tumorInfo) # get first tumor location


#---------------------------------------- run spm  ---------------------------------------- 

runSPM(event=fake_data,
       suppList=suppList, 
       maxgap=maxgap, 
       maxlen=maxlen, 
       maxsize=maxsize, 
       tType=opt$tType, 
       outputDir=spmDir)
