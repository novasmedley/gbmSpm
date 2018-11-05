# Plot training and testing predictions for a patient's set of visits
# Overalp this with clinical information about the patient


# libs
library(gbmSpm)
library(ggplot2)
library(grid)
library(gridExtra)
library(lattice)
library(ggrepel)
library(wesanderson)
library(ROCR)
# 
# install.packages('devtools')
# devtools::install_github("baptiste/egg")
library(egg) # if not installed, run previous two lines

source('renameCoefs.R')

# dirs
dataDir <- '~/Data/Agile/final_model' # all data used for logits
outDir <- '~/Data/Agile/figures'

# model data
results <- readRDS(file.path(dataDir, 'retrain_test_results.rds'))

# trainResults <- readRDS(file.path(dataDir, 'trainResults.rds')) # model performance results, including cutoffs
# testResults <- readRDS(file.path(dataDir, 'testResults.rds'))
testPredictions <- readRDS(file.path(dataDir,'plotReady_testPredictions.rds')) # model predictions for each visit
trainPredictions <- readRDS(file.path(dataDir,'plotReady_trainPredictions.rds'))

allCoefs <- read.table(file.path(dataDir, 'allCoefficients.csv'), sep=',', header=T, stringsAsFactors = F)

# spm data
events <- readRDS(file.path(dataDir,'cleaned_events_rateChange.rds')) # events data used for SPM

# chemo data
chemo <- readRDS(file.path(dataDir,'chemo_events.rds'))

#----- get data for specific model ----------

# prep coefficients for comparision against model variables
# (strings must match)

coefs <- na.omit(allCoefs[,c('variable','X9m.rateChange')]) # 6-month coefs
colnames(coefs) <- c('variable', 'oddsRatio')
coefs$variable <- gsub('`','', coefs$variable)
# printcoef9m <- coef6m[order(-coef6m$oddsRatio),]
# printcoef9m$variable <- renameCoefs(printcoef9m$variable)
# printcoef9m$variable <- gsub('\\}, \\{',' --> ',printcoef9m$variable)
# printcoef9m$variable <- gsub('<\\{','',printcoef9m$variable)
# printcoef9m$variable <- gsub('\\}>','',printcoef9m$variable)

# coef9m <- read.table(file.path(dataDir, 'renamed_coefs_9m.txt'), header=T, stringsAsFactors = F)

# model predictions for each visit
i <- 3 # 9-month
temp_test <- testPredictions[[i]]
temp_train <- trainPredictions[[i]]
temp_train$partition <- 'train'
temp_test$partition <- 'test'
allPreds <- rbind(temp_train,temp_test)
# rm(temp_test,temp_train)

# best cutoff, using test (there's one for train too)
bestPredThreshold <- results$fit9m$train$selectedThreshold['cutoff',1]

# merge data
modelData <- do.call('rbind', list(results$fit9m$fit$trainingData,results$fit9m$test$data))
colnames(modelData)[which(colnames(modelData)=='.outcome')] <- 'survivalIn270'
modelData$iois <- getPatIDs(row.names(modelData))
modelData$eventID <- as.integer(unlist(lapply(row.names(modelData), function(i) strsplit(i, '\\.')[[1]][3])))

# fix chemo details
chemo$agent_detail <- gsub(".*\\[(.*)\\].*", "\\1", chemo$agent_detail)
chemo$agent_detail <- gsub('TEMODAR Schering-Plough', 'Temodar', chemo$agent_detail, fixed=T)
chemo$agent_detail <- gsub('Avastin (Bevacizumab) Genetech', 'Avastin', chemo$agent_detail, fixed=T)
chemo$agent_detail <- gsub('\\((.*)\\)', '', chemo$agent_detail)

# fix level order of eventName for tumor volume
tumorVol <- events[events$agent=='measurement',]
tumorVol$eventName <- as.character(tumorVol$eventName)
tumorVol$eventName <- renameCoefs_truncated(tumorVol$eventName)
tumorVol$eventName <- gsub('  ',' ', tumorVol$eventName)
tumorVol$eventName <- gsub('tumor vol. rate change:', '', tumorVol$eventName)
tumorVol$eventName <- as.factor(tumorVol$eventName); levels(tumorVol$eventName)
tumorVol$eventName <- factor(tumorVol$eventName, levels=levels(tumorVol$eventName)[c(3,2,5,1,4,6,10,9,7,8)]); levels(tumorVol$eventName)

#----- plot individual patient timelines - prediction/prob over visits-----

# pats <- unique(testPredictions[[i]]$pat)
# length(unique(pats))
# checkThese <- unique(testPredictions[[i]]$pat[which(testPredictions[[i]]$correct=='no')]);length(checkThese)
# correct <- pats[!(pats %in% checkThese)]

#------ plot all patients -----------

ggplot(temp_test, aes(x=eventID, y=prediction, color=factor(label), shape=factor(partition))) +
  geom_hline(yintercept=bestPredThreshold) +
  geom_vline(aes(xintercept=patSurv)) +
  geom_vline(aes(xintercept=classLine), color='grey', show.legend = F) +
  geom_point() +
  scale_shape_manual(values=c(16,1)) +
  facet_wrap(~ name, scales='free_x')

#-----  inspect specific patients ---------
# pMany <- ggplot(allPreds[which(allPreds$pat %in% pats[201:256]),], aes(x=eventID, y=prediction, color=factor(label), shape=factor(partition))) +
#   geom_hline(yintercept=testResults[[i]]$cutData['cutoff',1]) +
#   geom_vline(aes(xintercept=patSurv)) +
#   geom_vline(aes(xintercept=classLine), color='grey', show.legend = F) +
#   geom_point() +
#   scale_shape_manual(values=c(16,1)) +
#   # scale_color_manual(values = wes_palette("Darjeeling")[c(2,1)] ) +
#   scale_color_manual(values = c('#F8766D','#00BFC4') ) +
#   facet_wrap(~ name, scales='free_x')
# dev.new();
# print(pMany)

#------ inspect plots -----------
# p <-ggplot(allPreds[which(allPreds$pat %in% c(154,124,275,306,259,268)),], aes(x=eventID, y=prediction, color=factor(label), shape=factor(partition))) + 
# cc <- c(10,122,7,70)
# 
# p <- ggplot(allPreds[which(allPreds$pat %in% c(197)),], aes(x=eventID, y=prediction, color=factor(label), shape=factor(partition))) +
#   geom_hline(yintercept=testResults[[i]]$cutData['cutoff',1], color='grey') +
#   geom_vline(aes(xintercept=patSurv), color=wes_palette("Darjeeling")[1]) +
#   geom_vline(aes(xintercept=classLine), color=wes_palette("Darjeeling")[2], show.legend = F) +
#   geom_point(size=1.5) + geom_text(aes(x=eventID, y=prediction, label=eventID), vjust=-1, hjust=.5, size=3) +
#   scale_shape_manual(values=c(16,1)) +
#   labs(color='truth', shape='dataset') +
#   xlab('Days from Baseline') + ylab('model prediction') +
#   theme_linedraw() +
#   theme(axis.title.y = element_text(angle = 90, vjust=0.5),
#         axis.text = element_text(size=8),
#         axis.title = element_text(size=10),
#         legend.text=element_text(size=9),
#         legend.title = element_text(size=9),
#         legend.margin=margin(l = -0.15, unit='cm'),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank()) +
#   scale_color_manual(values = wes_palette("Darjeeling")[c(2,1)] ) +
#   facet_wrap(~ name, scales='free_x', ncol = 1, strip.position = 'right') +
#   theme(strip.text = element_text(size=8))
# dev.new();
# print(p)
# dev.print(png, file.path(subDir,paste0('example_timeline.','png')), res=600,height=2, width=6, units="in")


iois <- 3829
dataFn <- file.path(outDir, paste0('example_coefficients_pat',iois,'.csv'))

patEvents <- events[events$iois==iois,]
patient <- modelData[modelData$iois==iois,] # features for each visit for one patient in trainData
patPreds <- allPreds[allPreds$pat==iois,]

#
# ----- : check specific days ----
# checkDays <- c(455, 476, 511, 1239, 1234,1449,1491, 1502,1555, 1613, 1628, 1631,1646,1694)
# varsEachDayToPrint <- c( c(1), c(1, 2,3,4,5))
# day <- patient[patient$eventID==1694,] ; day <- day[which(day==1)]; p <- coefs[coefs$variable %in% colnames(day),]; p[order(-p$oddsRatio),]
# 
# 
# 
# #-------- :  print out all patterns for each day for 9-month model --------------
# days <- unique(patient$eventID)
# patPatterns <- lapply(days, function(d) {
#   p <- patient[patient$eventID==d,]
#   if(p$gender=='MALE') p$'genderMALE' <- T
#   if(p$ageDecade=='[30,40)') p$'ageDecade[30,40)' <- T
#   if(p$ageDecade=='[40,50)') p$'ageDecade[40,50)' <- T
#   if(p$ageDecade=='[60,70)') p$'ageDecade[60,70)' <- T
#   if(p$ageDecade=='[70,90]') p$'ageDecade[70,90]' <- T
#   if(p$ethnicity=='WHITE') p$'ethnicityWHITE' <- T
#   # p <- p[which(p==1)]
#   p <- coefs[coefs$variable %in% colnames(p),]
#   p$variable <- renameCoefs(p$variable)
#   p$variable <- gsub('`','', p$variable)
#   p$oddsRatio <- as.numeric(p$oddsRatio)
#   p$variable <- renameCoefs(p$variable)
#   p$variable <- gsub('\\}, \\{',' --> ',p$variable)
#   p$variable <- gsub('<\\{','',p$variable)
#   p$variable <- gsub('\\}>','',p$variable)
#   p <- rbind(p, coefs[coefs$variable=='(Intercept)',])
#   p <- p[order(-p$oddsRatio),]
#   p$daysSinceBaseline <- d
#   p
# })
# patPatterns <- do.call(rbind, patPatterns)
# patPatterns <- patPatterns[,c(3,2,1)]
# write.table(patPatterns, dataFn, sep = ',', row.names = F, col.names = T)

#
#------ : inspect specific patient timeline-----------

# plot prep

patTV <- tumorVol[tumorVol$iois==iois,]
patTV$group <- 1

kps <- patEvents[patEvents$agent=='KPS', ]
kps$eventName <- as.character(kps$eventName)
kps$eventName <- as.factor(kps$eventName)
initialKPS <- kps[kps$eventName=='kpsInitial',]
initialKPS$group <- 1
levels(kps$eventName) <- rev(list('unchanged'='kpsInitial','significantly decreased'='kpsSigDECR', 'decreased' ='kpsDECR', 
                                  'unchanged'='kpsUNCH','increased'='kpsINCR','significantly increased'='kpsSigINCR'))
kps$group <- 1

patChemo <- chemo[chemo$iois==iois,]
patChemo <- patChemo[order(patChemo$eventID),]

chemo_y <- rep(c(0.01, 0.015,0.02,0.25), ceiling(nrow(patChemo)/3))
chemo_y <- chemo_y[1:nrow(patChemo)]
chemo_ylim <- scale_y_continuous(limits = c(0,0.03), expand = c(0, 0))
pred_ylim <-  scale_y_continuous(limits = c(0,1), breaks=seq(0,1,.25))

ylim <- scale_y_continuous(expand = c(0, 0))

# plot constants
survLine <- geom_vline(data=patPreds, aes(xintercept=patSurv), color=wes_palette("Darjeeling1")[1])
predThershold <- geom_hline(yintercept=bestPredThreshold, color='grey')
labelLine <- geom_vline(data=patPreds,aes( xintercept=classLine), color=wes_palette("Darjeeling1")[2], linetype='twodash',show.legend = F)
gTheme <- theme(axis.title.y = element_text(angle = 0, vjust=0.5, size=7),
                axis.text.x = element_text(size=7,margin = margin(0,0,-.25,0,unit = 'cm')),
                axis.text.y = element_text(size=7, margin = margin(0,0,0,-.1,unit='cm')),
                axis.title.x = element_blank(),
                legend.text=element_text(size=7),
                legend.title = element_text(size=7),
                legend.margin=margin(l = -0.3, unit='cm'),
                panel.grid.minor =element_blank())

# xmin <- max( min(patPreds$eventID)-months*30, min(patEvents$endEventID))
xmin <- 0
xmax <- patPreds$patSurv[1] + 10
# xlim <- xlim(c( xmin, xmax))
xlim <- scale_x_continuous(breaks=seq(xmin, xmax, 200), limits = c(xmin,xmax), expand = c(0, 0))

pointsize <- .7
linesize <- .1

#
#------- :: single timeline plots ---------
patPreds$label <- as.character(patPreds$label)
patPreds$label[patPreds$label=='dead'] <- '<= x-months'
patPreds$label[patPreds$label=='alive'] <- '> x-months'
pPred <- ggplot(patPreds, aes(x=eventID, y=prediction, color=factor(label), shape=factor(partition))) + 
  # geom_text(data=patPreds[patPreds$eventID %in% checkDays,],aes(x=eventID, y=prediction, label=eventID), vjust=-1, hjust=.5, size=3) +
  survLine + predThershold + labelLine + theme_linedraw() + gTheme + xlim + pred_ylim +
  theme(panel.grid.major.y= element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.key.size = unit(.8, 'lines')) +
  scale_shape_manual(values=c(19,1)) + 
  labs(color='Truth', shape='Data partition') +
  # xlab('Days from Baseline') + 
  ylab('Model estimate\n(probability of\nresidual survival\n\u2264 9 months') +
  scale_color_manual(values = wes_palette("Darjeeling1")[c(1,2)] , labels=c('\u2264 9 months','> 9 months') )+ 
  geom_point(size=1.5); print(pPred)

pKPS <- ggplot(kps, aes(x=eventID, y=eventName,group=group) ) +
  # geom_text_repel(data=kps[kps$eventID %in% checkDays,],aes(x=eventID, y=eventName, label=eventID), vjust=-1, hjust=.5, size=3) +
  ylab('KPS') + survLine  + labelLine + theme_minimal() + gTheme + xlim +
  geom_point(size=pointsize) +  geom_line(aes(x=eventID, y=eventName),size=linesize) +
  geom_text(data=initialKPS, aes(x=eventID, y='unchanged', label='initial', group=group), vjust=-1, hjust=.5, size=2) +
  scale_y_discrete(drop=FALSE) ; print(pKPS)

pNF <- ggplot(patEvents[patEvents$agent=='neuroFunc', ], aes(x=eventID, y=as.numeric(agent_detail))) +
  # geom_text_repel(data=patEvents[patEvents$eventID %in% checkDays,],aes(x=eventID, y=as.numeric(agent_detail), label=eventID), vjust=-1, hjust=.5, size=3) +
  scale_y_continuous(limits = c(0,4), breaks = seq(0,4,1)) + 
  ylab('Neurologic function') + survLine + labelLine + theme_minimal() + gTheme + xlim +
    geom_point(size=pointsize) +  geom_line(size=linesize); print(pNF)

pON <- ggplot(patEvents[patEvents$agent=='overallNeuro', ], aes(x=eventID, y=as.numeric(agent_detail))) +
  # geom_text(aes(x=eventID, y=as.numeric(agent_detail), label=eventID), vjust=-1, hjust=.5, size=3) +
  ylab('Overall\nneurologic status') + survLine + labelLine + theme_minimal() + gTheme + xlim + 
  scale_y_reverse(limits = c(2,-2), breaks = seq(-2,2,1)) +
    geom_point(size=pointsize) +  geom_line(size=linesize); print(pON)

pMS <- ggplot(patEvents[patEvents$agent=='mentalStatus', ], aes(x=eventID, y=as.numeric(agent_detail))) +
  # geom_text(aes(x=eventID, y=as.numeric(agent_detail), label=eventID), vjust=-1, hjust=.5, size=3) +
  scale_y_continuous(limits = c(0,2), breaks = seq(0,2,1)) + 
  xlab('Days from baseline') +
  ylab('Mental status') + survLine + labelLine + theme_minimal() + xlim +
    theme(axis.title.y = element_text(angle = 0, vjust=0.5, size=7),
                axis.text = element_text(size=7),
                axis.title.x = element_text(size=7),
                legend.text=element_text(size=7),
                legend.title = element_text(size=7),
                legend.margin=margin(l = -0.15, unit='cm'),
                panel.grid.minor.y =element_blank()) +
    geom_point(size=pointsize) +  geom_line(size=linesize); print(pMS)

pInter <- ggplot(patEvents[patEvents$agent=='Surgery' | patEvents$agent=='Radiation', ], aes(x=eventID, y=0)) +
  ylab('Intervention') + survLine  + labelLine + theme_minimal() + gTheme + ylim(c(-1, 1)) +  xlim +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust=0.5),
        axis.ticks.y = element_blank(),
        axis.text.y= element_blank()) +
  geom_point(size=pointsize) + geom_text_repel(aes(label=agent_detail), size=2.5) ; print(pInter)
  # geom_blank()

pTV <- ggplot(patTV, aes(x=eventID, y=eventName, group=group)) +
  # geom_text(aes(x=eventID, y=eventName, label=eventID), vjust=-1, hjust=.5, size=3) +
  ylab("Tumor\nvolume\nrate\nchange\n(mm\u00B3/day)") + survLine  + labelLine + theme_minimal() + gTheme + xlim +
  scale_y_discrete(drop = FALSE) +
  geom_point(size=pointsize) +  geom_line(size=linesize); print(pTV)

pChemo <- ggplot(patChemo) +
  ylab('Drug') + theme_minimal() + gTheme + xlim + survLine + labelLine + 
  geom_rect(aes(xmin=eventID, xmax=endEventID, ymin=rep(c(2,1,0), ceiling(nrow(patChemo)/3))[1:nrow(patChemo)], 
                ymax=rep(c(3,2,1), ceiling(nrow(patChemo)/3))[1:nrow(patChemo)], colour=factor(agent_detail)), size=.2 , fill='transparent') +
  geom_text_repel(aes(x=endEventID, label=agent_detail, y=rep(c(2,1,0), ceiling(nrow(patChemo)/3))[1:nrow(patChemo)], colour=factor(agent_detail)), 
                  size=2.5, vjust=-.7, hjust=.5) +
  # ylim(chemo_ylim) +
  # xlab('Days from baseline') +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x= element_blank(),
        axis.text = element_text(size=7),
        # axis.title.x = element_text(size=7),
        axis.title.y = element_text(angle = 0, vjust=0.5, size=7),
        axis.ticks.y = element_blank(),
        axis.text.y= element_blank()) +
  guides(colour=F); print(pChemo)

#
#--- group all timeline plots ----

# ggarrange(pPred, pChemo, pInter,pTV, pON, pKPS, pNF, pMS, heights=c(1.8,.6,.4,1.6,.7,.7,.7,.3) , ncol=1)
fn <- file.path(outDir,paste0('s8_timeline_expanded_pat',iois,'.','eps'))
cairo_ps(fn, height=6.5, width = 6.5)
exPlot <- ggarrange(pPred, pChemo, pInter,pTV, pON, pKPS, pNF, pMS, heights=c(2,.6,.4,1.6,.7,.7,.7,.3) , ncol=1)
# exPlot <- ggarrange(pPred, pChemo,pTV, pON, pKPS, pNF, pMS, heights=c(2,.6,1.6,.7,.7,.7,.3) , ncol=1)

print(exPlot)
dev.off()


#----- print temporal variables present for a specific day -------
# note that this ignores static variables (get them externally)
day <- patient[patient$eventID==id,] ; day <- day[which(day==576)]; coef9m[coef9m$variable %in% colnames(day),]









#----------- example 1 -----------
pPred <- ggplot(patPreds, aes(x=eventID, y=prediction, color=factor(label), shape=factor(partition))) + 
  geom_text(data=patPreds[patPreds$eventID %in% checkDays,],aes(x=eventID, y=prediction, label=eventID), vjust=-1, hjust=.5, size=3) +
  survLine + predThershold + labelLine + theme_linedraw() + gTheme + xlim + pred_ylim +
  theme(panel.grid.major.y= element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=10, margin=margin(0,15,0,0)),
        axis.text.y = element_text(size=10 ),
        axis.text.x = element_text(size=10),
        legend.key.size = unit(.8, 'lines'),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10),
        legend.key.height=unit(1.1,"line"),
        legend.margin=margin(0,0,0,.02)) +
  labs(color='Truth',shape='Data partition') +
  # xlab('Days from Baseline') + 
  ylab('Model estimate\n(probability of\nresidual survival\n\u2264 9 months') +
  scale_color_manual(values = wes_palette("Darjeeling1")[c(1,2)] , labels=c('\u2264 9 months','> 9 months'))+ 
  geom_point(size=2.8); print(pPred)

pChemo <- ggplot(patChemo) +
  ylab('Drug') + theme_minimal() + gTheme + xlim + survLine + labelLine + 
  geom_rect(aes(xmin=eventID, xmax=endEventID, ymin=rep(c(0,1,2,1), ceiling(nrow(patChemo)/3))[1:nrow(patChemo)], 
                ymax=rep(c(1,2,3,2), ceiling(nrow(patChemo)/3))[1:nrow(patChemo)], colour=factor(agent_detail)), size=.2 , fill='transparent') +
  geom_text_repel(aes(x=endEventID, label=agent_detail, y=rep(c(0,1,2,1), ceiling(nrow(patChemo)/3))[1:nrow(patChemo)], colour=factor(agent_detail)), size=3,vjust=-.7, hjust=.5) +
  # ylim(chemo_ylim) +
  xlab('Days from baseline') + 
  theme(panel.grid.major.y= element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x= element_text(size=10, margin=margin(15,0,0,0)),
        axis.title.y = element_text(size=10,angle = 0, vjust=0.5),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.text.y= element_blank()) +
  guides(colour=F); print(pChemo)

# ggarrange(pPred, pChemo, heights=c(2,.5) , ncol=1)
# dev.print(png, file.path(outDir,paste0('example_timeline1.','png')), res=600,height=4, width=10, units="in")

# fig 3 subplot
cairo_ps(file.path(outDir,paste0('subFig3_timeline_expanded_pat',iois,'.','eps')), height=4.5, width = 10)
sub <- ggarrange(pPred, pChemo, heights=c(2,.45) , ncol=1)
print(sub)
dev.off()


# -------- add custom labels ----------
theseDays <- c( 455,476, 511, 567, 568, 
                1234, 1239, 1281, 1449, 1502, 1631, 1646, 1694, 1736, 1862,1890, 1918)
day <- patient[patient$eventID==1234,] ; day <- day[which(day==1)]; coef9m[coef9m$variable %in% colnames(day),]
patPreds$text <- ''
text <- c('KPS decreased', #455
                'overall neurological status of 1', #476, 511 remembers the last two items
                'mental status of 1', #511
                'decreased volume rate', #567
                '', # 568 remembers overallNeuroNeg1
                'increased volume rate', #1234, only 1, TV went up
                'radiation', #1239, only 1, because TV went up, had radiation
                'decreased volume rate', #1281 (and 1253) remembers TV went up, but TV went down since
                'KPS decreased\nmental status of 1 ', #1449
                'increased volume rate', #1502, only 1
                'surgery\ndecreased volume', #1631, surgery and TV went down
                'KPS significantly decreased\noverall neurological status of -2', #1646
                'increased volume rate\nneurological function of 2\nmental status of 1', # 1694, KPS unchanged from ksp sig decrease, neuro func was 2 and still is 2
                '', #1736 similar patterns to 1694
                '', # 1862 similar patterns to 1694
                'increased volume rate', # 1890, remembers tumor volume rate went down
                '' # 1918 similar patterns to 1694
                )
text <- paste0(theseDays,'\n',text)
ind <- match(theseDays, patPreds$eventID)
patPreds$text[ind] <- text

patPreds$text <- ''
patPreds$text[ind] <- theseDays

# annot$text <- lapply(annot$eventID, function(i) {
#   day <- patient[patient$eventID==i,] 
#   day <- day[which(day==1)]
#   c <- coef9m[coef9m$variable %in% colnames(day),]
#   c <- c[c$oddsRatio > 1,]
#   if (nrow(c)>5) c <- c[1:5,]
#   c$variable <- renameCoefs(c$variable)
#   c$variable <- gsub('\\}, \\{',' --> ',c$variable)
#   c$variable <- gsub('<\\{','',c$variable)
#   c$variable <- gsub('\\}>','',c$variable)
#   text <- paste(c$oddsRatio, c$variable, '\n')
#   t <- paste(text, collapse = '')
# })
# annot$text <- unlist(annot$text)


#----------- example 2 -----------



