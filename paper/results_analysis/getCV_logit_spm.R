#---------------------- Description ---------------------------------------------------------- 
# find the best model in CV
# use it to get trained model (training ROC over all training partition)
# apply trained model to test data

#' run_logit is called many times, resulting in many log files in their subfolders,
#' one log file for each logit trained using each spm parameter specified. It is expected
#' that the log files were created with the same parameters (e.g., verbosity, training),
#' so that log files are exactly the same format.
#' 
#' This script extracts results based on the log files, which contrain cross validation results.
#' Logit data themselves were not saved since many logits were trained.
#' The best parameters are used to re-run the training (same seed), where logit model/data will be saved,
#' and used for testing.

#' location of logit results logs from cross-validation
#' each subfolder is a cspade parameter group
#' each cspade parameter has multiple log files, one for each tumor volume type
dataDir <- '~/Data/Agile/analysis_spm/logits' 
plotDir <- '~/Data/Agile/analysis_spm'

library(ROCR)
library(gridExtra)
library(grid)
library(wesanderson)
library(viridis)
library(scales) # to truncate geom_bars
library(egg)
library(reshape2) #melt


#------ local functions to format results after reading in -------------------------------------------------------------------------- 

getCVResults <- function(dir, fn, name) {
  cvResultsFile <- file.path(dir,fn)
  cvResults <- read.table(cvResultsFile, header=T, sep=',')
  
  # identify the cspade parameters the logit's features were based on
  temp <- do.call(rbind, lapply(cvResults$model, function(n) regmatches(n, gregexpr('[0-9]+',n))[[1]][2:5] ) )
  colnames(temp) <- c('support', 'gap', 'length', 'size')
  temp <- as.data.frame(temp)
  temp$support <- paste0( '.', as.character(temp$support))
  cvResults <- cbind(cvResults, temp)
  
  cvResults <- cvResults[order(cvResults$model),]
  cvResults$name <- rep(name, nrow(cvResults))
  
  ind <- length(strsplit(as.character(cvResults$model[1]), '\\/')[[1]])
  temp <- strsplit(as.character(cvResults$model),'\\/')
  cvResults$model <- unlist(lapply(temp, '[[', ind))
  cat('duplicates?: ', which(duplicated((cvResults))) )
  return(cvResults)
}

getCVByMonth <- function(cvResults, month, metricColName) {
  cv <- cvResults[cvResults$month==month,]
  cv <- cv[order(-cv[,metricColName]),]
  cv$model <- factor(cv$model, levels=unique(cv$model)) # order levels by decreasing ROC
  cv$name <- as.factor(cv$name)
  cv$group <- paste0(cv$tType, cv$name)
  cv$tType <- as.factor(cv$tType)
  cv$tType <- factor(cv$tType, levels=levels(cv$tType)[c(3,1,2,4)])
  cv
}

rankModels <- function(cvsResults, months, topN, metricColName) {
  r <- lapply(months, function(m) {
    models <- unique(cvsResults$model)
    ranked <- lapply(models, function(m) {
      results <- max(cvsResults[cvsResults$model==m,metricColName])
    })
    names(ranked) <- models
    ranked <- unlist(ranked)
    ranked <- ranked[order(-ranked)]
    names(ranked)[1:topN]
  })
  r <- unlist(r)
  r <- r[!duplicated(r)]
  return(r)
}

#------plotting functions-------------------------------------------------------------------------- 



# https://learnr.wordpress.com/2009/04/29/ggplot2-labelling-data-series-and-adding-a-data-table/
grid.show.layout(Layout)
vplayout <- function(...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}
subplot <- function(x, y) viewport(layout.pos.row = x,
                                   layout.pos.col = y)
mmplot <- function(a, b) {
  vplayout()
  print(a, vp = subplot(1, 1))
  print(b, vp = subplot(2, 1))
}
# end
Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2,0.5), c("null", "null")))
pd <- position_dodge(0.6) # move them .05 to the left and right

# ROC
plotCVResults <- function(cv, month, filename) {
  g <- ggplot(cv, aes(x=model, y=rocAUC, colour=tType, line=tType, group=group) ) +
    geom_point() + geom_line() +
    scale_y_continuous('Average\nArea Under the ROC Curve\n',limits=c(0.7,1),breaks=seq(0.7,1,.05), minor_breaks=seq(0.7,1,0.01)) +
    theme_linedraw() + 
    ggtitle(month) +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.text=element_text(size=15),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,10,0,0)),
          legend.text=element_text(size=12, color='grey30'), 
          legend.title=element_text(size=12),
          # legend.position=c('none'),
          legend.position = c(0.7,0.8),
          legend.key.size=unit(.8,'line'),
          plot.title = element_text(size=15, hjust=0.5)) +  #geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD), width=0.001, position=pd) +
    # geom_line(position=pd) + geom_point(position=pd) + 
    labs(color='tumor volume type') +
    scale_color_manual(values = wes_palette("Darjeeling1")[-4] )
  
  table <- cv[!(duplicated(cv$model)),]
  table <- table[, c('model','size','length','gap','support')]
  # colnames(table) <- c('model', 'events', 'visits', 'days', 'patients')
  table <- melt(table, id='model')
  
  t <- ggplot(table, aes(x=model,y=variable, label=format(value))) +
    geom_text(size = 4.2, color='grey30', hjust=0.6) + theme_bw() + scale_y_discrete('cSPADE\nparameters') +  scale_x_discrete(expand=c(0.04,0)) +
    theme(panel.grid.major = element_blank(), legend.position = "none",
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text=element_text(size=12),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,11,0,5)) ) +
    theme(plot.margin = unit(c(0,.19, 0, 0), "lines")) + xlab(NULL) + ylab(NULL)
  
  dev.new()
  mmplot(g, t)
  dev.print(png, paste0(filename,'.png'), res=800,height=5, width=5.1, units="in")
  # cairo_ps(file.path(paste0(filename,'.eps')), height=5, width = 5.1)
  dev.off()
}

# PR
plotCVResults_PR <- function(cv, month, filename, lp) {
  g <- ggplot(cv, aes(x=model, y=prAUC, colour=tType, line=tType, group=group) ) +
    geom_point() + geom_line() +
    scale_y_continuous('Average\nArea Under the PR Curve\n',limits=c(0.1,0.8),breaks=seq(0.1,0.8,.05), minor_breaks=seq(0.1,0.8,0.01)) +
    theme_linedraw() + 
    ggtitle(month) +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(), 
          axis.text=element_text(size=15),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,13,0,0)),
          legend.text=element_text(size=12, color='grey30'), 
          legend.title=element_text(size=12),
          legend.position=lp,
          # legend.position=c('none'),
          legend.key.size=unit(.8,'line'),
          plot.title = element_text(size=15, hjust=0.5)) +  #geom_errorbar(aes(ymin=ROC-ROCSD, ymax=ROC+ROCSD), width=0.001, position=pd) +
    # geom_line(position=pd) + geom_point(position=pd) + 
    labs(color='tumor volume type') +
    scale_color_manual(values = wes_palette("Darjeeling1")[-4] )
  
  table <- cv[!(duplicated(cv$model)),]
  table <- table[, c('model','size','length','gap','support')]
  # colnames(table) <- c('model', 'events', 'visits', 'days', '% patients')
  table <- melt(table, id='model')
  
  t <- ggplot(table, aes(x=model,y=variable, label=format(value))) +
    geom_text(size = 4.2, color='grey30', hjust=0.6) + theme_bw() + scale_y_discrete('cSPADE\nparameters') +  scale_x_discrete(expand=c(0.04,0)) +
    theme(panel.grid.major = element_blank(), legend.position = "none",
          panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text=element_text(size=12),
          axis.title.y = element_text(angle = 90, vjust = .5, size=12, margin=margin(0,11,0,5)) ) +
    theme(plot.margin = unit(c(0,.19, 0, 0), "lines")) + xlab(NULL) + ylab(NULL)
  
  dev.new()
  mmplot(g, t)
  dev.print(png, paste0(filename,'.png'), res=800,height=5, width=5.1, units="in")
  # cairo_ps(file.path(paste0(filename,'.eps')), height=5, width = 5.1)
  dev.off()
}



# id best model parameters ------------------------------------------------------------

# Parse cross-validation log files
theseDirs <- list.dirs(path=dataDir, full.names=F, recursive=F)

# types and typenames are matched
types <- c('continuous','response','percent','rateChange') # tumor volume types expected in log filename
typenames <- c('mm3','response criteria','percent change','rate change') # corresponding pretty names of tumor volume types
months <- c('2-month','6-month','9-month', '12-month') # month names for classification labels

cvs <- lapply(theseDirs, function(f) {
  cat(f, '\n')
  logFiles <- list.files(path=file.path(dataDir,f), pattern='\\.log$') # locate log files
  
  # extract results for each tumor volume type
  results <- lapply(seq(1,length(types)), function(i) {
    n <- logFiles[grepl(types[i],logFiles)] # read in a log file
    if (length(n) > 0) {
      log <- file.path(dataDir,f,n)
      con <- file(log)
      open(con)
      resultsCols <- read.table(con, skip=46, nrow=1, stringsAsFactors = F)  # expected results colnames
      r2 <- read.table(con, skip=0, nrow=1) # expected location of CV results for each classification
      r6 <- read.table(con, skip=29, nrow=1)
      r9 <- read.table(con, skip=29, nrow=1)
      r12 <- read.table(con, skip=29, nrow=1)
      close(con)
      
      r <- do.call('rbind', list(r2,r6,r9,r12)) # combine and format
      r <- r[,-1]
      colnames(r) <- c(resultsCols[1,])
      
      r$month <- months # attach additional info
      r$tType <- rep(typenames[i], nrow(r))
      r$model <- rep(f, nrow(r))
      temp <- do.call(rbind, lapply(r$model, function(n) regmatches(n, gregexpr('[0-9]+',n))[[1]][2:5] ) )
      colnames(temp) <- c('support', 'gap', 'length', 'size')
      temp <- as.data.frame(temp)
      temp$support <- paste0( '.', as.character(temp$support))
      r <- cbind(r, temp)
      r$name <- 'without'
      r
    } else {
      cat('\t ... file not found for ', typenames[i] , '\n')
    }
  })
  results <- do.call('rbind',results)
})
cvs <- do.call('rbind',cvs) 
cat('duplicates?: ', which(duplicated((cvs))) )
saveRDS(cvs, file=file.path(plotDir, 'spm_logit_cv_results.rds'))


#---  plot spm CV results ----
topN <- 15
top <- rankModels(cvs, months, topN, 'rocAUC') 

topCVS <- cvs[cvs$model %in% top,]
saveRDS(topCVS, file=file.path(plotDir, 'spm_logit_top_cv_results.rds'))

# individual plots
spmROC <- lapply(months, function(m) plotCVResults(getCVByMonth(cvResults=topCVS, m, 'rocAUC'), 
                                         m, 
                                         file.path(plotDir, paste0('spm_CVresults_aucroc_',m) ) ) )

spmPR1 <- plotCVResults_PR(getCVByMonth(cvResults=topCVS, months[1], 'prAUC'),
                           months[1],
                                         file.path(plotDir, paste0('spm_CVresults_aucpr_',months[1]) ),c(0.7,0.8))
spmPR2 <- lapply(months[2:4], function(m) plotCVResults_PR(getCVByMonth(cvResults=topCVS, m, 'prAUC'),
                                            m,
                                            file.path(plotDir, paste0('spm_CVresults_aucpr_',m) ),
                                            c(0.75,0.2)) )



#---  heatmap of cv results ----
# hp <- cvs
# hp$group1 <- paste0('days=',hp$gap,' \nvisits=',hp$length)
# hp$group2 <- paste0('events=', hp$size)
# hp2 <- hp[which(hp$month=='2-month'),]
# hp2$tType <- factor(hp2$tType, levels=c('mm3','percent change','rate change', 'response criteria'),
#                      labels=c('m','pc','rc','r'))
# 
# hp2$support <- round(as.numeric(hp2$support), digits=2)
# p <- ggplot(hp2) +
#   geom_tile() +
#   theme(axis.title.y = element_text(angle=0, vjust=.5)) +
#   theme(legend.text=element_text(size=8, color='grey30'), 
#         legend.title=element_text(size=8),
#         legend.position='right',
#         legend.key.size=unit(.8,'line'),
#         legend.background = element_rect(fill = "transparent", colour = "transparent")) +
#   theme(axis.text=element_text(size=8, color="#666666"),
#         axis.title.y = element_text(margin=margin(0,10,0,0),size=8),
#         axis.title.x = element_text(margin=margin(10,0,0,0),size=8)) +
#   xlab('Tumor Volume Type') + ylab('Support\n(% patients\nwith pattern)') +
#   labs(fill='area under\nthe ROC') +
#   scale_y_continuous(breaks=seq(0,.4,by=0.05)) +
#   scale_fill_viridis() + facet_grid(group1 ~ group2, scales='free', space='free_y') +
#   theme(strip.text.x = element_text(size = 8))
# print(p)
# dev.print(png, file.path(plotDir,paste0('cv_heatmap.','png')), res=600,height=7, width=6, units="in")

