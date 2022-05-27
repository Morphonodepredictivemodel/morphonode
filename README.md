# Morphonode Predictive Model (MPM)
The **Morphonode Predictive Model**, and the related **R** package `morphonode`, is the computational framework for ultrasound signatures detection and malignancy prediction in vulvar cancer.

## Install the latest release from source package

The latest stable release is the morphonode version 1.0.0. You can download this and older versions from the "Releases" panel of this website.
The zip or tar.gz package can be installed (from the home directory) in **R** with:

``` r
install.packages("~/morphonode-1.0.0.tar.gz", repos = NULL)
```

## Installation from source through the R package `devtools`

The latest development version can be installed through **GitHub** from any **R** environment:

``` r
# install.packages("devtools")
devtools::install_github("Morphonodepredictivemodel/morphonode")
```

This will download the source code, build the morphonode tar.gz package, and install it locally.
<br/><br/>

# Using the Morphonode Predictive Model suite

## Morphonode Predictive Model basic usage

The whole MPM suite can be launched in two very simple analysis steps:

- Step 0: &nbsp; library loading.
``` r
library(morphonode)
```
- Step 1: &nbsp; defining the ultrasound profile (in this example we will use a simulated malignant profile).
``` r
x <- new.profile(us.simulate(y = 1))
```
- Step 2: &nbsp; Launching the model!
``` r
mpm <- us.predict(x)
```
The results should look like the following:
```
       Morphonode Predictive Model output
# ------------------------------------------------------------------------------------------- #
|     Prediction (Morphonode-RFC): MALIGNANT (Y = 1)                                          |
|      Estimated prediction error: 0.029     (cutoff: E < 1)                                  |
|  Risk estimate (Morphonode-RBM): 0.974                                                      |
|                      Risk level: HIGH (> 0.29)                                              |
|       Signature (Morphonode-DT): HMR (high metastatic risk)                                 |
|                                                                                             |
|--- Input profile ---------------------------------------------------------------------------|
|                                                                                             |
| [0] w = 1.000 | 10.00 |  7.84 | 1 | 0 | 0 | 1 | 0 | 4 | 2 | 2 | 1 | 1 | 1 | 3 | Y = 1 | HMR |
|                                                                                             |
|--- Top-similar profiles (w: cosine similarity) ---------------------------------------------|
|                                                                                             |
| [1] w = 0.996 | 11.50 |  8.80 | 1 | 0 | 0 | 1 | 0 | 4 | 2 | 1 | 1 | 1 | 1 | 3 | Y = 1 | HMR |
| [2] w = 0.988 | 11.00 |  8.80 | 1 | 0 | 0 | 0 | 0 | 4 | 2 | 2 | 3 | 2 | 1 | 3 | Y = 1 | HMR |
| [3] w = 0.988 |  8.00 |  6.40 | 1 | 0 | 0 | 0 | 0 | 4 | 2 | 1 | 1 | 2 | 1 | 3 | Y = 1 | HMR |
| [4] w = 0.987 |  8.10 |  5.90 | 1 | 0 | 0 | 1 | 0 | 4 | 1 | 3 | 1 | 1 | 1 | 3 | Y = 1 | HMR |
| [5] w = 0.986 | 10.30 |  8.80 | 1 | 0 | 0 | 1 | 1 | 4 | 4 | 2 | 1 | 1 | 1 | 2 | Y = 1 | HMR |
# ------------------------------------------------------------------------------------------- #
```
&nbsp;

The MPM suite is composed by 4 modules:

- **Morphonode-RFC**. Random forest classification (RFC) and prediction error (E) estimate.
  The predicted phenotype can be either malignant (y = 1) or non-malignant (y = 0).
  As a rule of thumb, if E is above or equal to 1, the prediction should be considered as unreliable.
- **Morphonode-RBM**. Malignancy risk estimation by robust binomial modeling (RBM).
  The RBM offers a continuous estimation of y (i.e., the malignancy risk), thus the higher the accordance with the RFC, the higher the prediction
  reliability.
  This module also suggests when the risk reach *moderate* (p > 0.23) or *high* levels (p > 0.29). These cutoffs reflect the optimal risk cutpoint between
  malignant and non-malignant subjects, by maximizing F1 score and Sensitivity/Specificity, respectively.
- **Morphonode-DT**. Decision tree-based metastatic risk signature detection.
  A *high metastatic risk* (HMR) signature is characterized by a high risk of a single metastasis event, whereas
  a *metastatic signature* (MET) is typical of malignancies showing multiple metastasis events.
  Conversely, a *low metastatic risk* (LMR) signature is generally associated with non-malignant phenotypes.
  Finally, a *moderate malignancy risk* (MMR) signature is the group with highest heterogeneity and requires RFC and RBM results to be characterized.
- **Morphonode-SP**. Similarity prolfiling module. The module searches and ranks ultrasound profiles from the given (by default, the simulated) ultrasound
  features dataset. The default function is cosine similarity and the 5 top-similar profiles are shown to screen.
