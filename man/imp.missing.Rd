% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPMcore.R
\name{imp.missing}
\alias{imp.missing}
\title{Impute missing values}
\usage{
imp.missing(M, x = NULL, mode = NULL, ...)
}
\arguments{
\item{M}{A matrix or data.frame containing missing values.}

\item{x}{An optional vector that will be attached to M. 
This can be useful if data with missing values can be attached to 
a reference dataset.}

\item{mode}{Either "cat" (cateorical variables) or 
"con" (continuous variables).}

\item{...}{Currently ignored.}
}
\value{
A data.frame with imputed missing values.
}
\description{
Wrapper for the \code{imputeR} function 
\code{\link[imputeR]{impute}}.
}
\examples{

# Sample 30 subjects from the morphonode simulated dataset
data <- mosaic::sample(mpm.us, 30, replace = FALSE, prob = NULL)[, 2:15]

# Entries with missing values
missing <- matrix(c(10.0, 6.3, 1, 0, 0, 0, 0, 1, NA, 2, NA, 2, 3, NA,
                     6.4, 2.1, 1, 0, 0, 0, 0, 1, NA, 2, NA, 1, 1, NA),
                  nrow = 2, byrow = TRUE)
colnames(missing) <- colnames(mpm.us[, 2:15])

# Defining categorical subset
data.cat <- data.frame(apply(mpm.us[, 2:15], 2, factor))

# Imputing missing values
data.cat <- imp.missing(data.cat, x = missing, mode = "cat")
print(tail(data.cat))

}
\author{
Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
}
