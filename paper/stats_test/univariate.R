library(plyr)

source('renameCoefs.R')
getUnivariateOR <- function(labels, feature) {
  dat <- data.frame(survival=labels, pattern=feature)
  alive <- c(nrow(dat[which(dat$survival=='alive' & dat$pattern==0),]),
             nrow(dat[which(dat$survival=='alive' & dat$pattern==1),]))
  dead <- c(nrow(dat[which(dat$survival=='dead' & dat$pattern==0),]),
            nrow(dat[which(dat$survival=='dead' & dat$pattern==1),]))
  contingency <- rbind(alive,dead)
  
  ft <- fisher.test(contingency)
  
  return(data.frame(univariateOR = ft$estimate,
                    CI95_lower = ft$conf.int[1], 
                    CI95_upper = ft$conf.int[2], 
                    p = ft$p.value))
}

getFisherTests <- function(topcoefs,labels, patterns) {
  results <- lapply(topcoefs$variable, function(p){
    feature <- patterns[,which(colnames(patterns) %in% c(p))]
    print(p)
    if (p=='laterality:  rightTRUE') { # fix naming of this feature for test
      p <- 'laterality:  right'
      feature <- patterns[,which(colnames(patterns) %in% c(p))]
    }
    getUnivariateOR(labels,feature)
    # if (length(feature)==0) {
    #   cat('... error')
    #   NA
    # } else {
    #   getUnivariateOR(labels,feature)
    # }
  })
  results <- do.call('rbind',results)
  topcoefs <- do.call('cbind',list(topcoefs,results))
  topcoefs$roundP <- round(topcoefs$p, digits=5)
  topcoefs$sig <- topcoefs$p < 0.05
  return(topcoefs)
}

# Get data
dataDir <- '~/Data/Agile'
outDir <- '~/Data/Agile/univariate'

spmLogitDir <- file.path(dataDir, 'final_model') # temporal pattern approach
retrainResults <- readRDS(file.path(spmLogitDir, 'retrain_test_results.rds'))  # read in logit data

# get visits which were actually used for model training and testing (patterns required X amount of prior history)
# data holds the binary labels for each pattern that was observed across patients
# data also holds the outcome labels
# together, these two columns in the data are used as input to create contingency tables for the fisher test
temp <- retrainResults$fit2m$fit$trainingData
modelData_2m <- rbind(retrainResults$fit2m$test$data,temp)
colnames(modelData_2m)[which(colnames(modelData_2m)=='.outcome')] <- 'survivalIn60'

temp <- retrainResults$fit6m$fit$trainingData
modelData_6m <- rbind(retrainResults$fit6m$test$data,temp)
colnames(modelData_6m)[which(colnames(modelData_6m)=='.outcome')] <- 'survivalIn180'

temp <- retrainResults$fit9m$fit$trainingData
modelData_9m <- rbind(retrainResults$fit9m$test$data,temp)
colnames(modelData_9m)[which(colnames(modelData_9m)=='.outcome')] <- 'survivalIn270'

rm(temp)

rownames(modelData_2m) <- sub(".*?\\.(.+)", "\\1", rownames(modelData_2m) )
rownames(modelData_6m) <- sub(".*?\\.(.+)", "\\1", rownames(modelData_6m) )
rownames(modelData_9m) <- sub(".*?\\.(.+)", "\\1", rownames(modelData_9m) )


# get the selected features from each model (stored in coefficients text file)
allCoefs <- read.table(file.path(spmLogitDir,'allCoefficients.csv'), sep=',', header=T, stringsAsFactors = F)

vars <- allCoefs$variable # remove extra chars, need to match variable names in train and test data 
vars <- gsub('`','',vars)
allCoefs$variable <- vars

intercept <- which(allCoefs$variable=='(Intercept)') # remove intercept, not doing this test
allCoefs <- allCoefs[-intercept,]

coefs2 <- na.omit(allCoefs[,c(1,2)]) # split coefs by model
coefs6 <- na.omit(allCoefs[,c(1,3)])
coefs9 <- na.omit(allCoefs[,c(1,4)])

# run test
coefs2 <- coefs2[order(coefs2$X2m.rateChange),]
coefs2top <- coefs2[c(1:10,(nrow(coefs2)-9):nrow(coefs2)),] # get only the top patterns
results2m <- getFisherTests(labels=modelData_2m$survivalIn60, topcoefs=coefs2top, patterns=modelData_2m)

coefs6 <- coefs6[order(coefs6$X6m.rateChange),]
coefs6top <- coefs6[c(1:10,(nrow(coefs6)-9):nrow(coefs6)),]
results6m <- getFisherTests(labels=modelData_6m$survivalIn180, topcoefs=coefs6top, patterns=modelData_6m)

coefs9 <- coefs9[order(coefs9$X9m.rateChange),]
coefs9top <- coefs9[c(1:10,(nrow(coefs9)-9):nrow(coefs9)),]
results9m <- getFisherTests(labels=modelData_9m$survivalIn270, topcoefs=coefs9top, patterns=modelData_9m)

write.csv(file=file.path(spmLogitDir,'2-month_top20_univariate_results.csv'),results2m )
write.csv(file=file.path(spmLogitDir,'6-month_top20_univariate_results.csv'),results6m )
write.csv(file=file.path(spmLogitDir,'9-month_top20_univariate_results.csv'),results9m )

# get the top pattern's supports
freqseq2m <- as(readRDS(file=file.path(spmLogitDir,'freqseq_supp0.3_gap60_len3_size3.rds')), 'data.frame') 
freqseq2m$variable <- as.character(freqseq2m$sequence)
coefs2 <- merge(freqseq2m, coefs2)

freqseq6m <- as(readRDS(file=file.path(spmLogitDir,'freqseq_supp0.25_gap60_len3_size4.rds')), 'data.frame') 
freqseq6m$variable <- as.character(freqseq6m$sequence)
coefs6 <- merge(freqseq6m, coefs6)

freqseq9m <- as(readRDS(file=file.path(spmLogitDir,'freqseq_supp0.3_gap60_len3_size4.rds')), 'data.frame') 
freqseq9m$variable <- as.character(freqseq9m$sequence)
coefs9 <- merge(freqseq9m, coefs9)

# save freseq as csv
write.csv(file=file.path(spmLogitDir,'freqseq_supp0.3_gap60_len3_size3.csv'),freqseq2m[, c('sequence','support')] )
write.csv(file=file.path(spmLogitDir,'freqseq_supp0.25_gap60_len3_size4.csv'),freqseq6m[, c('sequence','support')] )
write.csv(file=file.path(spmLogitDir,'freqseq_supp0.3_gap60_len3_size4.csv'),freqseq9m[, c('sequence','support')] )

