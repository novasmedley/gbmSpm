# Find predictive temporal patterns in glioblastoma



`gbmSPM` contains the code used to study temporal patterns in tumor volumes and neurological exams to predict residual survival in a glioblastoma cohort in the paper:

Smedley, Nova F., Benjamin M. Ellingson, Timothy F. Cloughesy, and William Hsu. "Longitudinal Patterns in Clinical and Imaging Measurements Predict Residual Survival in Glioblastoma Patients." *Scientific reports* 8, no. 1 (2018): 14429.

Methods for data preprocessing, sequential pattern mining (SPM) with `arules` and `arulesSequences`, and logit modelling with  `glmnet` and `caret` were converted into an R package named `gbmSPM`.


This `gbmSPM` repository has:

* `R` : R library
* `man` : library documentation
* `vignettes`: vignettes for using the library
* `data` : dummy data
* `paper` : code used in the published paper
* `examples` : Rscripts similar to how the paper used the library

Additional details can be found in Supplemental Materials.


## Getting started

### Install
To just get the R package, in RStudio:

```R
install.packages("devtools")
library(devtools)
install_github("novasmedley/gbmSPM")
```

Now you can run the examples in 'RStudio examples and vignettes.'

To get the full repository, on the terminal, `cd` to your desired directory to put this repository and enter:

```bash
git clone https://github.com/novasmedley/gbmSpm.git
```

Now you can run the examples in 'command line examples.'



## Usage

### gbmSPM R package
As usual, function documentation can be reached by running `?function_name`, in RStudio, e.g., `?getAge`.

#### RStudio examples and vignettes

In RStudio:

```R
library(gbmSpm)

# create output directory
outputDir <- '~/gbm_spm_example'
spmDir <- file.path(outputDir,'spm')
logDir <- file.path(outputDir,'spm_logs')
lapply(c(outputDir, spmDir, logDir), function(i) ifelse(!dir.exists(i), dir.create(i, showWarnings = F, recursive = T), F) )


# load example data
data("fake_data")

# set SPM hyperparameters
minSupp <- 0.4
maxgap <- 60
maxlen <- 2
maxsize <- 2
tType <- 'rate'
suppList <- seq(minSupp, 0.4, .05)

# prep example data
fake_tumorInfo <- fake_data$events
fake_demo <- fake_data$demo
fake_data$events <- cleanData(fake_data$events, tType)
cat('...',nrow(fake_data$events), " events left for SPM after cleaning", '\n')
fake_data <- merge(fake_data$events, fake_data$person, by='iois', all.x=T)
fake_data <- prepDemographics(fake_data, fake_demo)
fake_data <- prepSurvivalLabels(fake_data)
fake_data <- getTumorLocation(fake_data, fake_tumorInfo)

# perform SPM to create data in output directories:
# 1) transactions, 2) patterns, and 3) feature vectors
runSPM(event = fake_data,
       suppList = suppList,
       maxgap = maxgap,
       maxlen = maxlen,
       maxsize = maxsize,
       tType = tType,
       outputDir = spmDir)


# load data back in for inspection
features <- readRDS(file.path(spmDir, 'sup0.4g60l2z2','featureVectors_rateChange.rds'))
```

The transactions file is used by `arules` for generating frequent sequences, see `?createTransactions`. The patterns file is a list of frequent sequences and the days that it was observed for each patient, see `?findPatternDays`.

Feature vectors are clinical visits (delineated by each pair of patient ID and event ID) and whether a frequent sequences was observed or not, and contains other clinical information for that event (e.g., patient age at event ID). These are the inputs to the logit modeling.


For step-by-step examples, see the vignettes: "[Generate sequential patterns](http://htmlpreview.github.io/?https://github.com/novasmedley/gbmSpm/blob/master/vignettes/generate_features.html)"" and "[Predicting residual survival](http://htmlpreview.github.io/?https://github.com/novasmedley/gbmSpm/blob/master/vignettes/survival_prediction.html)", "[Model selection and plotting performance](http://htmlpreview.github.io/?https://github.com/novasmedley/gbmSpm/blob/master/vignettes/model_selection.html)."

#### command line examples
Since the pipeline was repeated while searching for hyperparameters, two components were converted for use via the command line:

1. Generating temporal features:

    To use dummy data, run `example_spm.R ` from the `examples` folder:

    ```
    $ Rscript /PATHTO/example_spm.R --tType 'rate' --minSupp 0.4 --minSuppList 'no' --maxgap 60 --maxlength 2 --maxsize 2--outDir ~/gbm_spm_example
    ```

2. Predicting residual survival

    To use dummy data, run `example_logit_spm.R ` in the `examples` folder:

    ```
    $ Rscript example_logit_spm.R --tType 'rate' --maxgap 60 --maxlength 2 --lmax 0.2 --llength 100 --dataFolder sup0.25g60l2z2 --prefix logits --dir ~/gbm_spm_example_multi --saveLogit yes
    ```

Alternatively, the example scripts can be explored in RStudio. For more details, see in-depth examples in `vignettes`.


### Source code to paper

The `paper` folder contains the R code used in the published paper. Specifically, it has:

* `spm` : creating temporal patterns
* `logit` : logistic regression modeling
* `vol_thresholding` : predicting residual survival using only tumor volume information, aka by thresholding volume
* `results_analysis` : extracting cross-validation results and select best approach
* `stats_test` : testing for differences in ROC curves and univariate analysis of patterns
* `finalModel_retrain.R` : retraining selected approach and apply to test partition
* misc. functions, e.g., creating paper figures

For example, to generate temporal features, `run_spm.R` in the `paper` folder was called in the terminal:

```
$ Rscript run_spm.R --tType 'rate' --minSupp 0.4 --minSupp 0.4 --minSuppList 'no' --maxgap 60 --maxlength 2 --outDir ~/spm_analysis
```

The script first calls a database, performs data preprocessing, mines for patterns, and then creates the feature vectors formatted for logit modeling. The feature vectors included temporal patterns identified by SPM, age, and static variables such as ethnicity.

`run_spm.R` returns nothing, but saves the days for which the patterns occurs, `patternDays` and the corresponding feature vectors, `featureVectors` as .RData objects. These .RData objects are later read by `run_logits_spm.R`.

Thus, bash scripts were used to execute SPM under different conditions, see `run_spm.sh` and `run_spm_all.sh`.


## Citation
If you want to cite this work, please cite the paper:
```
@article{smedley2018longitudinal,
  title={Longitudinal Patterns in Clinical and Imaging Measurements Predict Residual Survival in Glioblastoma Patients},
  author={Smedley, Nova F and Ellingson, Benjamin M and Cloughesy, Timothy F and Hsu, William},
  journal={Scientific reports},
  volume={8},
  number={1},
  pages={14429},
  year={2018},
  publisher={Nature Publishing Group}
}
```
and the repo:
```
@misc{smedleygbmspm,
  title={gbmSPM},
  author={Smedley, Nova F},
  year={2018},
  publisher={GitHub},
  howpublished={\url{https://github.com/novasmedley/gbmSpm}},
}
```
