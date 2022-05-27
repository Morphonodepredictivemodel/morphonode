# Morphonode Predictive Model (MPM)
The **Morphonode Predictive Model**, and the related **R** package `morphonode`, is the computational framework for ultrasound signatures detection and malignancy prediction in vulvar cancer.

## Installation

The latest development version can be installed through **GitHub** from any **R** environment:

``` r
# install.packages("devtools")
devtools::install_github("Morphonodepredictivemodel/morphonode")
```

This will download the source code, build the morphonode tar.gz package, and install it locally.
<br/><br/>

# Using the Morphonode Predictive Model suite

## Morphonode Predictive Model basic usage

The whole MPM suite can be launched in two very simple steps:

- Defining the ultrasound profile (in this example we will use a simulated malignant profile)
``` r
x <- new.profile(us.simulate(y = 1))
```
- Launching the model!
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
<br/><br/>
The MPM suite is composed by 4 modules:

- **Morphonode-RFC**: random forest classification and prediction error (E) estimate.
  The predicted phenotype can be either malignant (y = 1) or non-malignant (y = 0).
  As a rule of thumb, if E is above or equal to 1, the prediction should be considered as unreliable.
- **Morphonode-RBM**: malignancy risk estimation by robust binomial modeling.
  The RBM offers a continuous estimation of y (i.e., the malignancy risk), thus the higher the accordance with the RFC, the higher the prediction
  reliability.
- ...
