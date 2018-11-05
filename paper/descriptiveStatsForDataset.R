#---------------------------------------- Description---------------------------------------- 

# run sequentially

library(gbmSpm)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(ggrepel)
library(egg)

outDir <- '~/Data/Agile/figures'
group <- 'agile' # group name in MySQL configuration file
patTabName <- 'patient'
eventTabName <- 'agent6'


#---- initial data descriptions -----

# get from Agile database
data <- getData(database=group, eventTable=eventTabName, personTable=patTabName) # event data
tumorInfo <- data$event # save tumor location and laterility strings before event cleaning

num_patients <- length(unique(data$person$iois))
num_visits <- nrow(unique(data$event[c('iois','startdate')]))

radiation <- data$event[data$event$agent=='Radiation',]
patientsWithNoRadiation <- unique(data$person$iois[which(!(data$person$iois%in%  radiation$iois ))])

demo <- getDemographics(database=group, personTable=patTabName) # includes MGMT biomarker
demo <- getTumorLocation(demo,tumorInfo) # get first tumor location

initialAge <- lapply(split(data$event, data$event$iois), function(i) min(as.Date(as.character(i$startdate), format='%Y-%m-%d')) )
initialAge <- data.frame(iois=as.integer(names(initialAge)), initialAge=do.call('c', initialAge))

# events to use for spm
events <- data$event
events <- cleanData(events, tType='c') # tType will contain all visits considered in events (all others may have NA at the initial imaging visit)

# attach patient info for each event
# note that only eventName, iois, and eventID columns are used for SPM
events <- merge(events, data$person, by='iois', all.x=T)
events <- prepDemographics(events, demo) # prep for each event, since age does change
events <- prepSurvivalLabels(events) # get survival labels, these also change
events <- getTumorLocation(events,tumorInfo) # get first tumor location

#-------- event stats ----------

cleaned_visits <- nrow(unique(events[c('iois','eventID')]))
totalImageStudies <- nrow(events[events$agent=='measurement',])
num_censored <- length(unique(events$iois[events$status=='ALIVE']))

# did al patients get temodar?
temodar <- chemo[chemo$agent_detail=='TEMODAR Schering-Plough',]
identical(length(which(unique(events$iois) %in% temodar$iois)), length(unique(events$iois)))

# number of follow up months:
followups <- unlist(lapply(unique(events$iois), function(i) {
   start <- min(events[events$iois==i,'startdate']) 
   end <- max(events[events$iois==i,'enddate'])
   as.numeric(end-start)
}))
mean(followups)

# average days between visits
daysBetweenVisits <- unlist(lapply(unique(events$iois), function(i) {
  e <- unique(events$startdate[events$iois==i])
  e <- e[order(e)]
  daysBetween <- unlist(lapply(seq(1,length(e)-1), function(j) e[j+1] - e[j] ))
  mean(daysBetween)
}))
daysBetweenVisits[is.na(daysBetweenVisits)] <- 0
mean(daysBetweenVisits)
median(daysBetweenVisits)


#------ Fig S2 frequency of visits per patient -----
# ignore nvisits, sum visits to get total
visitFreq <- unlist(lapply(unique(events$iois), function(i) length(unique(events$startdate[events$iois==i]) ))) # all data used for spm
mean(visitFreq)
median(visitFreq)
plot(hist(visitFreq, breaks=length(visitFreq)))

freqE <- as.data.frame(visitFreq)
names(freqE) <- c('visits')
freqE$type <- 'overall'
freqE$month <- ' '

tvFreq <- readRDS(file.path('~/Data/Agile/analysis_vol/threshold','tumorVol_threshold_visitFreq.rds')) # freq from tumor vol
freq <- tvFreq$trainPart$freq
freq <- data.frame(visits=rep(freq,length(months)), 
                   month= unlist(lapply(months, function(i) rep(i, length(freq)) ))  )
freq$month <- as.character(freq$month)
freq$type <- 'tumor volumes [train]'
freqT <- freq

spmFreq <- readRDS(file.path('~/Data/Agile/final_model','spm_sampleFreq.rds'))
freq <- lapply(spmFreq, '[[',1)
freq <- lapply(freq, '[[', 3)
freq <- melt(freq)
names(freq) <- c('visits','month')
freq$type <- 'temporal patterns [train]'
freq_all <- freq

freq <- lapply(spmFreq, '[[',2)
freq <- lapply(freq, '[[', 3)
freq <- melt(freq)
names(freq) <- c('visits','month')
freq$type <- 'temporal patterns [test]'
freq_all <- rbind(freq_all,freq, freqT)
# freq_all$month <- factor(freq_all$month, levels=months)

freq_text <- split(freq_all, list(freq_all$month, freq_all$type))
means <- lapply(freq_text, function(i) mean(i$visits))
medians <- lapply(freq_text, function(i) median(i$visits))
freq_text <- data.frame(mean=unlist(means), median=unlist(medians),name=names(means))
freq_text$month <- unlist(lapply(as.character(freq_text$name), function(i) strsplit(i, '\\.')[[1]][1]))
freq_text$type <- unlist(lapply(as.character(freq_text$name), function(i) strsplit(i, '\\.')[[1]][2]))
freq_text$x <- 100
freq_text$x <- as.numeric(freq_text$x)
freq_text$mean <- round(as.numeric(freq_text$mean), digits=1)
freq_text$label <- paste0('mean=',freq_text$mean,', median=', freq_text$median)

freq_text_overall <- data.frame(mean=mean(freqE$visits),median=median(freqE$visits),type='overall',x=135)
freq_text_overall$mean <- round(as.numeric(freq_text_overall$mean), digits=1)
freq_text_overall$label <- paste0('mean=',freq_text_overall$mean,', median=', freq_text_overall$median)

# visits are all thes same across months since same length and gap for three models, so just plot one

pf <- ggplot(freq_all[freq_all$month=='2-month',], aes(x=type, y=visits, color=type)) +
  geom_boxplot() +  coord_flip() + theme(legend.position = 'none') +
  theme_linedraw() +
  geom_text(data=freq_text, aes(x=type,y=x,label=label), vjust=-.8,size=3) +
  theme(axis.text=element_text(size=10, color="#666666"),
        axis.title.y = element_text(margin=margin(0,10,0,0),size=10, angle = 0, vjust=.5),
        axis.title.x = element_text(margin=margin(10,0,0,0),size=10) ) +
  theme(legend.text=element_text(size=9, color='grey30'), 
        legend.title=element_text(size=9),
        legend.position='right',
        legend.key.size=unit(.8,'line'),
        plot.title = element_text(size=15, hjust=0.5),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.margin=margin(l = -0.15, unit='cm')) +
  ylab('Visit Frequency Per Patient') + xlab('Data\nPartition') +
  guides(color=FALSE) +
  scale_y_continuous(breaks=c(seq(0,180,10))) +
  scale_color_manual(values = c(wes_palette("FantasticFox1")[3:5], wes_palette("Darjeeling2")[5]) )
# pf <- pf + facet_grid(month ~ ., scales='free')
print(pf)

overall_pf <- ggplot(freqE, aes(x=type, y=visits, color=type)) +
  geom_boxplot() +  coord_flip() + theme(legend.position = 'none') +
  geom_text(data=freq_text_overall, aes(x=type,y=x,label=label), vjust=-.9,size=3) +
  theme_linedraw() +
  theme(axis.text=element_text(size=10, color="#666666"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin=margin(10,0,0,0),size=10) ) +
  ylab('Visit Frequency Per Patient') + xlab('Data Partition') + 
  guides(color=FALSE) +
  scale_y_continuous(breaks=c(seq(0,600,50))) +
  scale_color_manual(values = wes_palette("Darjeeling2")[5] )
print(overall_pf)

ggarrange(overall_pf, pf, heights=c(.5,3) , ncol=1)

# dev.print(png, file.path(outDir, 'patVisitFreq.png'), res=600,height=4.5, width=8, units="in")
cairo_ps(file.path(outDir,'s2_patVisitFreq.eps'), height=3, width=8)
pp <- ggarrange(overall_pf, pf, heights=c(1,3) , ncol=1)
print(pp)
dev.off()

#--------- tab 1 characteristics of patients with events used for SPM ---------
patients <- data$person[data$person$iois %in% events$iois,]
patients <- merge(patients, demo, by='iois')
patients <- merge(patients, initialAge, by='iois')
nrow(patients)

meanOS <- mean(patients[patients$status=='DECEASED', 'survival']); round(meanOS, digits=1)
sdOS <- sd(patients[patients$status=='DECEASED', 'survival']); round(sdOS, digits=1)

status <- table(patients$status); status; round(100*status/nrow(patients), digits=1)

patients$DOB <- as.Date(as.character(patients$DOB), format='%Y-%m-%d')
patients$initialAge <- getAge(patients$DOB, patients$initialAge)
meanInitialAge <- mean(patients$initialAge); round(meanInitialAge, digits=1)
sdInitialAge <- sd(patients$initialAge); round(sdInitialAge, digits=1)

gender <- table(patients$gender); gender; round(100*gender/nrow(patients), digits=1)

ethnicity <- table(patients$ethnicity); ethnicity; round(100*ethnicity/nrow(patients), digits=1)

mgmt <- table(patients$MSP); mgmt; round(100*mgmt/nrow(patients), digits=1)

idh1 <- table(patients$IDH1); idh1; round(100*idh1/nrow(patients), digits=1)

lateral <- table(patients$tumorLaterality); lateral; round(100*lateral/nrow(patients), digits=1)

locations <- lapply(seq(10,19), function(i) {
  l <- length(which(patients[,i]==T))
  c(l,round(100*l/nrow(patients), digits=1))
})
names(locations) <- colnames(patients)[10:19]
locations

multipleLocations <- tumorLocation[,4:13]
multipleLocations$count <- apply(multipleLocations, 1, function(i) sum(i))
length(which(multipleLocations$count<1))
length(which(multipleLocations$count>2))
length(which(multipleLocations$count==1))
length(which(multipleLocations$count==2))


#---------  fig s3 item frequency for spm ------
events <- data$event
events <- cleanData(events, tType='r') # want rate change 

itemFreq <- as.data.frame(table(events$eventName))
itemFreq <- itemFreq[c(20, 21, 3,4,1,6,2,5, 7:14, 19,18,15:17, 24,23,26,22,25,27,31,30,28,29) ,]
row.names(itemFreq) <- seq(1,nrow(itemFreq))
itemFreq$Var1 <- as.character(itemFreq$Var1)

itemFreq$Var1[22:31] <- lapply(itemFreq$Var1[22:31], function(i) gsub('tumorVolRateChange:','', i))
itemFreq$category <- c( 'radiotherapy', 'surgery', rep('KPS', 6), rep('mental status', 3), rep('neurologic function', 5), 
                        rep('overall neurologic\nstatus', 5), rep('tumor volume\nrate change (mm\u00B3/day)',10) )
itemFreq$item <-c('radiotherapy', 'surgery', 'initial', 'significantly decreased', 'decreased', 'unchanged','increased','significantly increased',
                   paste0(' ',0:2), paste0('  ',0:4), paste0('   ',-2:2),  gsub('tumor vol. rate change:','' , itemFreq$Var1[22:31]))
itemFreq$category <- as.factor(itemFreq$category)
itemFreq$item <- factor(itemFreq$item, levels=rev(unique(itemFreq$item)))

# cairo_ps(file.path(outDir,'itemFreq.eps'), height=3, width=4.5)
ggplot(itemFreq, aes(x=item, y=Freq, fill=category)) +
  geom_bar(stat="identity", position=position_dodge(), width = .8) + coord_flip() + ylim(c(0,4500)) +
  geom_text(aes(label=Freq),hjust=-.05, size=1.75) +
  scale_fill_brewer(type = "qual", palette = 'Paired') +
  theme_bw() +
  theme(axis.title=element_text(size=7), 
        axis.title.y = element_text(angle=0, vjust=.5),
        axis.text=element_text(size=6),
        legend.text=element_text(size=6), 
        legend.title=element_text(size=7),
        legend.margin=margin(l=-0.3, unit='cm'),
        legend.key.size=unit(.8,'line'),
        panel.grid = element_line(size=0.2)) +
  labs(fill='Variable',y='Frequency', x='State') 
dev.print(tiff, file.path(outDir,'s3_itemFreq.tiff'), res=600,height=3, width=5, units="in")
# dev.off()


#---------- number of mined patterns for each final model --------

model2patterns <- readRDS(file.path('~/Data/GBM/Agile/lookBack_allVisits','sup0.1g60l2z4','patterns_rateChange.rds'))
model69patterns <- readRDS(file.path('~/Data/GBM/Agile/lookBack_allVisits','sup0.2g60l3z2','patterns_rateChange.rds'))

length(model2patterns) # num patterns for the 2-month model
length(model69patterns) # for the 6 an 9 month model (had the same spm contraints)