# Find predictive sequential patterns in imaging and



`gbmSPM` contains the code used to study temporal patterns in tumor volumes and neurological exams to predict residual survival in a glioblastoma cohort in the paper:

Smedley, Nova F., Benjamin M. Ellingson, Timothy F. Cloughesy, and William Hsu. "Longitudinal Patterns in Clinical and Imaging Measurements Predict Residual Survival in Glioblastoma Patients." *Scientific reports* 8, no. 1 (2018): 14429.

Methods for data preprocessing, sequential pattern mining (SPM) with `arules` and `arulesSequences`, and logit modelling with  `glmnet` and `caret` were converted into an R package named `gbmSPM`.


`gbmSPM` has:

* `R` : R library
* `man` : library documentation
* `vignettes`: vignettes for using the library
* `data` : dummy data
* `paper` : code used in the published paper
* `examples` : Rscripts similar to how the paper used the library

Additional details can be found in Supplemental Materials.


## Getting started

#### Install

In RStudio:
```R
install.packages("devtools")
library(devtools)
install_github("novasmedley/gbmSPM")
```

## Usage

### gbmSPM R package
As usual, function documentation can be reached by running `?function_name`, in RStudio, e.g., `?getAge`.

#### vignettes

For step-by-step examples, see the vignettes: "Generate sequential patterns" and "Predicting residual survival", "Model selection and plotting performance."

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
    $ Rscript example_logit_spm.R  XX
    ```

    For more details, see the  vignette.

Alternatively, the example scripts can be explored in RStudio.

### Source code to paper

The `paper` folder contains the R code used in the published paper. Specifically, it has:

* `spm` : creating temporal patterns
* `logit` : logistic regression modeling
* `vol_thresholding` : predicting residual survival using only tumor volume information, aka by thresholding volume
* `results_analysis` : extracting cross-validation results and select best approach
* `stats_test` : testing for differences in ROC curves and univariate analysis of patterns
* `finalModel_retrain.R` : retraining selected approach and apply to test partition
* misc. functions, e.g., creating paper figures

For example, to generate temporal features, `run_spm.R` in the *paper* folder was called in the terminal:

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
