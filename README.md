# Morphonode Predictive Model (MPM)
The **Morphonode Predictive Model**, and the related **R** package `morphonode`, is the computational framework for ultrasound signatures detection and malignancy prediction in vulvar cancer.

## Installation

The latest development version can be installed into an **R** environment from **GitHub**:

``` r
# install.packages("devtools")
devtools::install_github("Morphonodepredictivemodel/morphonode")
```

# Using Morphonode Predictive Model suite

## Basic usage

The whole MPM suite can be launched in two very simple steps:

- Defining the ultrasound profile (in this example we will use a simulated malignant profile)
``` r
x <- new.profile(simulate.us(y = 1))
```
- Launching the model!
``` r
mpm <- us.predict(x)
```
