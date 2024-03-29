% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MPMutils.R
\name{jaccard}
\alias{jaccard}
\title{Jaccard similarity}
\usage{
jaccard(x, y, ...)
}
\arguments{
\item{x}{A dichotomous vector.}

\item{y}{A dichotomous vector.}

\item{...}{Currently ignored.}
}
\value{
A numeric value corresponding to the Jaccard similarity.
}
\description{
Compute Jaccard similarity between two vectors.
}
\examples{

# Sample two random ultrasound profiles from the default dataset
x <- mosaic::sample(mpm.us, 1, replace = FALSE, prob = NULL)
x <- as.matrix(dichotomize(x[, 2:15]))
y <- mosaic::sample(mpm.us, 1, replace = FALSE, prob = NULL)
y <- as.matrix(dichotomize(y[, 2:15]))

# Compute the euclidean distance
r <- jaccard(x, y)
print(r)

}
\seealso{
\code{\link[morphonode]{similarity}}, 
   \code{\link[morphonode]{euclidean}}
}
\author{
Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
}
