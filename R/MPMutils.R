#  Morphonode Predictive Model (MPM) - The morphonode R package

#  Copyright (C) 2022 Fernando Palluzzi
#  e-mail: <fernando.palluzzi@gmail.com>
#  Bioinformatics facility 
#  Gemelli Science and Technological Park (GSTeP)
#  Fondazione Policlinico Universitario Agostino Gemelli IRCCS,
#  Largo Agostino Gemelli 8, 00168 Roma, Italy

#  morphonode is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  morphonode is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.

# -------------------------------------------------------------------- #


#' @title Build a random forest classifier
#'
#' @description Wrapper for the \code{\link[CMA]{GenerateLearningsets}} 
#'    \code{\link[CMA]{classification}} functions.
#'
#' @param data An (n, m) matrix or data.frame with no outcome attribute.
#' @param status A vector of length n, containing the outcome.
#' @param k number of cross-validation iterations (default = 5).
#' @param method One of the \code{\link[CMA]{GenerateLearningsets}} 
#'    methods, including: "LOOCV", "CV", "MCCV", and "bootstrap" 
#'    (default = "CV").
#' @param ... Currently ignored.
#'
#' @importFrom CMA GenerateLearningsets classification
#' @export
#'
#' @return An object of class A list of objects of class "cloutput" and 
#'    "clvarseloutput", respectively.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{us.predict}} to launch all 
#'    morphonode modules at once.
#'    See also \code{\link[morphonode]{topsim}} for a simple similarity 
#'    search.
#'
#' @references
#' 
#' Slawski M, Daumer M, Boulesteix AL. CMA - a comprehensive Bioconductor 
#' package for supervised classification with high dimensional data. 
#' BMC Bioinformatics 9, 439 (2008).
#' <https://doi.org/10.1186/1471-2105-9-439>
#'
#' @examples
#' 
#' # Extract a subset of 300 subjects and an outcome vector of length 30 
#' # from the default simulated dataset
#' 
#' x <- sample(mpm.us, 300, replace = FALSE, prob = NULL)
#' y <- x$y
#' x <- x[, 2:15]
#' dim(x)
#' length(y)
#' 
#' # Build a 5-fold cross-validation object
#' 
#' CV <- rfc(x, status = y)
#' 
#' # Performances of the first of five predictors
#' 
#' CV1 <- CV[[1]]
#' P <- performance(obs = CV1@y, pred = CV1@yhat)
#' print(P)
#'
rfc <- function(data, status, k = 5, method = "CV") {
	
	L <- GenerateLearningsets(n = nrow(data), status,
							  method = method,
							  fold = k,
							  strat = TRUE)
	
	C <- classification(as.matrix(data), status,
	                    learningsets = L,
	                    classifier = "rfCMA")
	
	return(C)
}

#' @title Compute Brier scores
#'
#' @description Compute Brier scores for a given dataset.
#'
#' @param data An (n, m) matrix or data.frame with no outcome attribute.
#' @param status A vector of length n, containing the outcome.
#' @param k number of cross-validation iterations (default = 5).
#' @param method One of the \code{\link[CMA]{GenerateLearningsets}} 
#'    methods, including: "LOOCV", "CV", "MCCV", and "bootstrap" 
#'    (default = "CV").
#' @param ... Currently ignored.
#'
#' @importFrom CMA GenerateLearningsets classification evaluation
#' @export
#'
#' @return A vector of Brier scores of length n (one value per subject).
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso \code{\link[morphonode]{loss}}
#'
#' @references
#' 
#' Slawski M, Daumer M, Boulesteix AL. CMA - a comprehensive Bioconductor 
#' package for supervised classification with high dimensional data. 
#' BMC Bioinformatics 9, 439 (2008).
#' <https://doi.org/10.1186/1471-2105-9-439>
#' 
#' Brier GW. Verification of forecasts expressed in terms of probability. 
#' Monthly Weather Review. 1950;78(1):1-3. 
#' <https://doi.org/10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2>
#'
#' @examples
#' 
#' # Extract a subset of 300 subjects and an outcome vector of length 30 
#' # from the default simulated dataset
#' 
#' x <- sample(mpm.us, 300, replace = FALSE, prob = NULL)
#' y <- x$y
#' x <- x[, 2:15]
#' print(dim(x))
#' print(length(y))
#' 
#' # Compute brier scores
#' E <- brier(x, y)
#' print(quantile(E))
#'
brier <- function(data, status, k = 5, method = "CV") {
	C <- rfc(data = data, status = status, k = k, method = method)
	E <- evaluation(C, scheme = "observationwise",
	                   measure = "brier score")@score
	return(E)
}

#' @title Loss function
#'
#' @description Compute the loss function needed for prediction error 
#'    estimation.
#'
#' @param p A vector of predicted outcome values.
#' @param y A vector of observed outcome values.
#' @param method Loss function definition. One between "log" (default) 
#'    and "sqerror".
#' @param ... Currently ignored.
#'
#' @details The subject-level prediction error, calculated through the 
#'    \code{\link[morphonode]{topsim}} function, is strongly correlated 
#'    with the loss function values. However, while the former can be 
#'    only calculated for the entire training set, the latter can be 
#'    computed for the new input profile(s). Therefore, the prediction 
#'    error (E) for the input is calculated as: 
#'    E = b0 + b*L, where L is the loss function value.
#'    L can be currently defined as either "log": 
#'    L = -1*(y*log(p) + (1 - y)*log(1 - p)) (default), or 
#'    "sqerror": L = (p - y)^2.
#'    Additionally, the cost is calculated as either average loss sum(L)/n 
#'    or root mean squared error sqrt(sum(L)/n), for "log" and "sqerror", 
#'    respectively.
#'
#' @import randomForest
#' @importFrom stats cor
#' @importFrom lsa cosine
#' @importFrom imputeR impute
#' @export
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item "loss", loss function values;
#' \item "cost", cost function value.
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso \code{\link[morphonode]{us.predict}}, 
#'    \code{\link[morphonode]{brier}}
#'
#' @examples
#' 
#' # RBM prediction vs. reality over the simulated dataset
#' 
#' p <- predict(mpm.rbm$fit, dichotomize(mpm.us[2:15], asFactor = TRUE),
#'              type = "response")
#' 
#' L <- loss(p, mpm.us$y)
#' 
#' print(quantile(L$loss))
#' 
#' print(L$cost)
#' 
#' # Overall RBM performances
#' 
#' y.hat <- ifelse(p > 0.5, 1, 0)
#' P <- performance(obs = mpm.us$y, pred = y.hat)
#' 
#' print(P)
#'
loss <- function(p, y, method = "log") {
	if (method == "log") {
		if (p == 1) {
			p <- 0.9999
		} else if (p == 0) {
			p <- 0.0001
		}
		L = -1*(y*log(p) + (1 - y)*log(1 - p))
		C = sum(L)/length(p)
	} else if (method == "sqerror") {
		L = (p - y)^2
		C = sqrt(sum(L)/length(p))
	} else {
		stop("invalid method")
	}
	return(list(loss = L, cost = C))
}

#' @title Build a randomForest object
#'
#' @description Build a Random Forest Classifier (RFC) of class 
#'    \code{randomForest}.
#'
#' @param model An R formula of shape y ~ x1 + x2 + ... + xn.
#'    You can use \code{model <- formula("y ~ .")} to include all the 
#'    variables of the input dataset.
#' @param data A matrix or data.frame with subjects as rows and 
#'    ultrasound features as columns. The outcome should be named "y".
#' @param n Number of trees to be generated by bootstrap.
#' @param m Number of randomly sampled variables per tree branching.
#' @param vset An optional data.frame of the same shape as the input data. 
#'    This will be used as external validation set. 
#'    A validation set can be also extracted from the input dataset, 
#'    using the \code{\link[morphonode]{vpart}} function.
#' @param ... Currently ignored.
#'
#' @details This is the class of RFC used in the 
#'    \code{\link[morphonode]{us.predict}} and 
#'    \code{\link[morphonode]{rfc.predict}} fuctions. 
#'    The default randomForest object in the \code{morphonode} package 
#'    is an ensemble of 5 RFCs trained over a simulated dataset of 
#'    948 subjects (508 non-malignant and 440 malignant profiles), using 
#'    a nested 5-fold cross-validation scheme. The training/validation 
#'    details and performances of this dataset are stored in the 
#'    \code{\link[morphonode]{mpm.rfc}} object.
#'
#' @import randomForest
#' @export
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item "RFC", an object of class \code{randomForest};
#' \item "performance", a list containing RFC performances.
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso \code{\link[morphonode]{us.predict}}, 
#'    \code{\link[morphonode]{rfc.predict}}, 
#'    \code{\link[morphonode]{vpart}}
#'
#' @references
#' 
#' Fragomeni SM, Moro F, Palluzzi F, Mascilini F, Rufini V, Collarino A, 
#' Inzani F, Giacò L, Scambia G, Testa AC, Garganese G (2022). 
#' Evaluating the risk of inguinal lymph node metastases before surgery 
#' using the Morphonode Predictive Model: a prospective diagnostic study. 
#' Ultrasound xx Xxxxxxxxxx xxx Xxxxxxxxxx. 00(0):000-000.
#' <https://doi.org/00.0000/00000000000000000000>
#' 
#' Liaw A, Wiener M. Classification and Regression by randomForest (2002). 
#' R News, 2(3):18-22. <https://doi.org/10.1023/A:1010933404324>
#'
#' @examples
#' 
#' # Extract a subset of 500 subjects and an outcome vector of length 30 
#' # from the default simulated dataset
#' 
#' x <- sample(mpm.us, 500, replace = FALSE, prob = NULL)
#' x <- x[, 2:16]
#' model <- formula("y ~ .")
#' 
#' # Data partitioning (75% training set, 25% validation set)
#' 
#' x <- vpart(x, p = 0.75)
#' 
#' # RFC building (1000 bootstrapped trees, 3 random variables per split)
#' 
#' rfc <- buildPredictor(model, x$training.set, n = 1000,
#'                       vset = x$validation.set)
#' print(rfc$performance)
#' 
#' \donttest{
#' 
#' # RFC building (10000 bootstrapped trees, 3 random variables per split)
#' 
#' rfc <- buildPredictor(model, x$training.set, vset = x$validation.set)
#' print(rfc$performance)
#' 
#' }
#'
buildPredictor <- function(model, data, n = 10000, m = 3, vset = NULL) {
	
	RFC <- randomForest(model, data = data,
						ntree = n,
						mtry = m,
						keep.forest = TRUE,
						proximity = TRUE,
						importance = TRUE)
	print(RFC)
	if (is.null(vset)) {
		return(RFC)
	} else {
		C <- predict(RFC, vset)
		P <- performance(vset$y, C)
		return(list(RFC = RFC, performance = P))
	}
}

#' @title Data partitioning utility
#'
#' @description Extract a random partition from an input dataset.
#' @param data An input data.frame.
#' @param p Proportion of rows to be extracted from the input dataset.
#' @param shuffle A logical value. If TRUE, the input rows are randomly 
#'    shuffled before data partitioning.
#'
#' @return A list of 2 data.frames:
#' \enumerate{
#' \item "training.set", the portion of the input data defined by p;
#' \item "validation.set", the portion of the input data defined by 1-p.
#' }
#'
#' @export
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso \code{\link[morphonode]{buildPredictor}}, 
#'
#' @examples
#' 
#' # Extract a subset of 300 subjects and an outcome vector of length 30 
#' # from the default simulated dataset
#' 
#' x <- sample(mpm.us, 300, replace = FALSE, prob = NULL)
#' 
#' # Data partitioning
#' 
#' x <- vpart(x, p = 0.75)
#' print(dim(x$training.set))
#' print(dim(x$validation.set))
#'
vpart <- function(data, p, shuffle = FALSE) {
	if (shuffle) {
		data <- data[sample(1:nrow(data), replace = FALSE),]
	}
	j <- sample(nrow(data), size = p*nrow(data), replace = FALSE)
	X <- data[j,]
	V <- data[-j,]
	return(list(training.set = X, validation.set = V))
}

#' @title Predictor performance calculation
#'
#' @description Compute sensitivity, speciticity, precision, F1 score, 
#'    and accuracy for a set of predictions.
#'
#' @param obs A vector of observed values.
#' @param pred A vector of predicted values.
#' @param C 2x2 contingency table (alternative to \code{obs} and 
#'    \code{pred}). If \code{obs} and \code{pred} are given, the 
#'    contingency table is automatically computed. The contingency 
#'    table has the observed values as rows and the predicted ones as 
#'    columns. By default, true negatives are located at position [1, 1], 
#'    while true positives at [2, 2] (see below).
#' @param y Contingenty table orientation. If \code{y = "0,1"} (default), 
#'    true negatives are located at position [1, 1], while true positives 
#'    at position [2, 2]. If \code{y = "1,0"}, these positions are inverted.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A list of 6 objects:
#' \enumerate{
#' \item "ctable", 2x2 contingency table of predicted vs. observed values;
#' \item "sensitivity", Se = TP/(TP + FN);
#' \item "specificity", Sp = TN/(TN + FP);
#' \item "precision", PPV = TP/(TP + FP), 
#'    also called "positive predictive value";
#' \item "NPV", NPV = TN/(TN + FN), 
#'    "negative predictive value";
#' \item "F1", F1 = 2*Se*PPV/(Se + PPV);
#' \item "accuracy", (TP + TN)/(TP + TN + FP + FN).
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso \code{\link[morphonode]{us.predict}}
#'
#' @examples
#' 
#' # RBM performance
#' p <- predict(mpm.rbm$fit, dichotomize(mpm.us[2:15], asFactor = TRUE),
#'              type = "response")
#' y.hat <- ifelse(p > 0.5, 1, 0)
#' rbm <- performance(mpm.us$y, y.hat)
#' 
#' # Compare RBM and RFC F1 scores
#' print(rbm$F1)
#' print(mpm.rfc$performance$F1)
#'
performance <- function(obs = NULL, pred = NULL, C = NULL, y = "0,1") {
	# Confusion Table
	if (class(C)[1] != "matrix" & class(C)[1] != "table") {
		C <- table(obs, pred)
	}
	
	#print(CT)
	#cat("\n")
	
	if (y == "0,1") {
		TP <- C[2, 2]
		FN <- C[2, 1]
		FP <- C[1, 2]
		TN <- C[1, 1]
	} else if (y == "1,0") {
		TP <- C[1, 1]
		FN <- C[1, 2]
		FP <- C[2, 1]
		TN <- C[2, 2]
	}
	
	C <- matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE)
	
	Se <- TP/(TP + FN)                # Sensitivity
	Sp <- TN/(TN + FP)                # Specificity
	A <- (TP + TN)/sum(C)             # Accuracy
	P <- TP/(TP + FP)                 # Precision (PPV)
	N <- TN/(TN + FN)                 # Negative predictive value
	F1 <- (2*Se*P)/(Se + P)           # F1-score
	
	W <- list(ctable = C, sensitivity = Se, specificity = Sp,
	          precision = P, NPV = N, F1 = F1, accuracy = A)
	return(W)
}

#' @title Euclidean distance
#'
#' @description Compute the euclidean distance between two vectors.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A numeric value corresponding to the euclidean distrance.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso \code{\link[morphonode]{similarity}}, 
#'    \code{\link[morphonode]{jaccard}}
#'
#' @examples
#' 
#' # Sample two random ultrasound profiles from the default dataset
#' x <- sample(mpm.us, 1, replace = FALSE, prob = NULL)
#' x <- as.matrix(x[, 2:15])
#' y <- sample(mpm.us, 1, replace = FALSE, prob = NULL)
#' y <- as.matrix(y[, 2:15])
#' 
#' # Compute the euclidean distance
#' d <- euclidean(x, y)
#' print(d)
#'
euclidean <- function(x, y) {
	return(sqrt(sum((y - x)^2)))
}

#' @title Jaccard similarity
#'
#' @description Compute Jaccard similarity between two vectors.
#'
#' @param x A dichotomous vector.
#' @param y A dichotomous vector.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A numeric value corresponding to the Jaccard similarity.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso \code{\link[morphonode]{similarity}}, 
#'    \code{\link[morphonode]{euclidean}}
#'
#' @examples
#' 
#' # Sample two random ultrasound profiles from the default dataset
#' x <- sample(mpm.us, 1, replace = FALSE, prob = NULL)
#' x <- as.matrix(dichotomize(x[, 2:15]))
#' y <- sample(mpm.us, 1, replace = FALSE, prob = NULL)
#' y <- as.matrix(dichotomize(y[, 2:15]))
#' 
#' # Compute the euclidean distance
#' r <- jaccard(x, y)
#' print(r)
#'
jaccard <- function(x, y) {
	return(sum(1*(y == x))/length(x))
}

#' @title Compute different similarity coefficients
#'
#' @description Compute the similarity between two vectors.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param f Similarity function. One among "cosine" (default), 
#'    "jaccard" (for dichotomous vectors only), "pearson", "spearman",
#'    "kendall", or "euclidean".
#' @param ceiling If \code{f = "euclidean"}, the similarity is 
#'    computed as 1/distance. This argument limits the similarity to a 
#'    very high value, in case the euclidean distance is equal to 0.
#'    The default value is 1E+07.
#' @param ... Currently ignored.
#'
#' @importFrom lsa cosine
#' @importFrom stats cor
#' @export
#'
#' @return A numeric value corresponding to the similarity coefficient.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @references
#' 
#' Leydesdorff L (2005). Similarity Measures, Author Cocitation 
#' Analysis,and Information Theory. In: JASIST 56(7), pp.769-772.
#' <https://doi.org/10.48550/arXiv.0911.4292>
#'
#' @seealso \code{\link[morphonode]{jaccard}}, 
#'    \code{\link[morphonode]{euclidean}}
#'
#' @examples
#' 
#' # Sample two random ultrasound profiles from the default dataset
#' x <- sample(mpm.us, 1, replace = FALSE, prob = NULL)
#' x <- as.matrix(x[, 2:15])
#' y <- sample(mpm.us, 1, replace = FALSE, prob = NULL)
#' y <- as.matrix(y[, 2:15])
#' 
#' # Compute the cosine similarity
#' r <- similarity(x, y)
#' print(r)
#'
similarity <- function(x, y, f = "cosine", ceiling = 1E+07) {
	if (f == "cosine") {
		r <- lsa::cosine(x, y)
	} else if (f == "jaccard") {
		r <- jaccard(x, y)
	} else if (f == "euclidean") {
		r <- euclidean(x, y)
		r <- ifelse(r == 0, ceiling, 1/r)
	} else {
		r <- cor(x, y, method = f)
	}
	return(r)
}

#' @title Compute bootstrap standard errors
#'
#' @description Compute bootstrap standard errors (SE) for a given 
#'    (generalized) linear model.
#'
#' @param fit An object of class \code{glm} or \code{lm}.
#' @param boot A bootstrapped model object (see the examples section).
#' @param probs Bootstrap intervals (default = c(0.025, 0.975)).
#' @param z0 z-score under null hypothesis (default = 1.96).
#' @param b0 Effect size under null hypothesis (default = 0).
#' @param ... Currently ignored.
#'
#' @import mosaic
#' @export
#'
#' @return A data.frame of 9 columns:
#' \enumerate{
#' \item "Variable", variable name;
#' \item "Estimate", parameter (effect size) estimate;
#' \item "se.boot", bootstrap standard error;
#' \item "lower", confidence interval lower bound;
#' \item "upper", confidence interval upper bound;
#' \item "conf.level", confidence level;
#' \item "method", estimation method;
#' \item "z", z-score = (estimate - b0)/SE;
#' \item "P", 2-sided p-value; i.e., 2*pnorm(-abs(z)).
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @references
#' 
#' Pruim R, Kaplan DT, Horton NJ (2017). The mosaic Package: Helping 
#' Students to 'Think with Data' Using R. The R Journal, 9(1), 77–102. 
#' <https://journal.r-project.org/archive/2017/RJ-2017-024/index.html>
#'
#' @seealso See \code{\link[morphonode]{p.boot}} for performance indices 
#'    bootstrap confidence intervals.
#'    See also \code{\link[mosaic]{do}} for resampling.
#'
#' @examples
#' 
#' # Dichotomized dataset creation
#' x <- dichotomize(mpm.us, asFactor = TRUE)
#' 
#' # Model specification
#' model <- formula(paste0(c("y ~ shortAxis + cortical + hilum + ",
#'                           "inflammatoryStroma + extracapsularSpread + ",
#'                           "ecostructure + FID + VFL + corticalThickening + ",
#'                           "vascularPattern + CMID + shape + grouping + ",
#'                           "colorScore"), collapse = ""))
#' 
#' # Binomial model fitting (MLE)
#' fit <- glm(model, data = x, family = "binomial")
#' 
#' # Binomial model fitting
#' n.reps <- 100
#' boot <- do(n.reps) * coef(glm(model, data = resample(x), family = "binomial"))
#' 
#' # Bootstrap SE calculation
#' SE <- boot.se(fit, boot)
#' print(SE)
#'
boot.se <- function(fit, boot, probs = c(0.025, 0.975), z0 = 1.96, b0 = 0) {
	conf.level <- probs[2] - probs[1]
	b0 <- summary(fit)$coefficients
	coef <- NULL
	for (i in 1:nrow(b0)) {
		se.boot <- sd(boot[, i])
		#hist(z)
		q <- qnorm(probs)
		lower <- b0[i, 1] - q[2]*se.boot
		upper <- b0[i, 1] - q[1]*se.boot
		Variable <- colnames(boot)[i]
		Estimate <- b0[i, 1]
		coef <- rbind(coef, data.frame(Variable, Estimate, se.boot,
		              lower, upper, conf.level, method = "Bootstrap SE",
		              z = Estimate/se.boot,
		              P = 2*pnorm(-abs(Estimate/se.boot))))
	}
	return(coef)
}

sensitivity <- function(tp, fn) {
	return(tp/(tp + fn))
}

specificity <- function(tn, fp) {
	return(tn/(tn + fp))
}

ppv <- function(tp, fp) {
	return(tp/(tp + fp))
}

npv <- function(tn, fn) {
	return(tn/(tn + fn))
}

accuracy <- function(tp, fp, fn, tn) {
	return((tp + tn)/(tp + fp + fn + tn))
}

f1 <- function(tp, fp, fn, tn) {
	Se <- tp/(tp + fn)
	Sp <- tn/(tn + fp)
	P <- tp/(tp + fp)
	F1 <- (2*Se*P)/(Se + P)
	return(F1)
}

fpr <- function(fp, tn) {
	return(fp/(fp + tn))
}

fdr <- function(tp, fp) {
	return(fp/(fp + tp))
}

fnr <- function(tp, fn) {
	return(fn/(fn + tp))
}

fnc <- function(tp, tn, fn) {
	return(fn/(tp + tn))
}

p.f1 <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fn <- C[2, 1]
	fp <- C[1, 2]
	tn <- C[1, 1]
	p <- f1(tp, fp, fn, tn)
	return(p)
}

p.mode <- function(x) {
	f <- table(x)
	return(names(f)[which(f == max(f))])
}

p.accuracy <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fn <- C[2, 1]
	fp <- C[1, 2]
	tn <- C[1, 1]
	p <- accuracy(tp, fp, fn, tn)
	return(p)
}

p.sensitivity <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fn <- C[2, 1]
	p <- sensitivity(tp, fn)
	return(p)
}

p.specificity <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tn <- C[1, 1]
	fp <- C[1, 2]
	p <- specificity(tn, fp)
	return(p)
}

p.plr <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fp <- C[1, 2]
	fn <- C[2, 1]
	tn <- C[1, 1]
	p <- sensitivity(tp, fn)
	q <- specificity(tn, fp)
	plr <- p/(1 - q)
	return(plr)
}

p.nlr <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fp <- C[1, 2]
	fn <- C[2, 1]
	tn <- C[1, 1]
	p <- sensitivity(tp, fn)
	q <- specificity(tn, fp)
	nlr <- (1 - p)/q
	return(nlr)
}

p.ppv <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fp <- C[1, 2]
	p <- ppv(tp, fp)
	return(p)
}

p.npv <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	fn <- C[2, 1]
	tn <- C[1, 1]
	p <- npv(tn, fn)
	return(p)
}

p.fpr <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	fp <- C[1, 2]
	tn <- C[1, 1]
	p <- fpr(fp, tn)
	return(p)
}

p.fdr <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fp <- C[1, 2]
	p <- fdr(tp, fp)
	return(p)
}

p.fnr <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fn <- C[2, 1]
	p <- fnr(tp, fn)
	return(p)
}

p.fnc <- function(x, j, real = NULL) {
	data <- x[j,]
	if (is.null(real)) {
		real <- data[, 2]
		data <- data[, 1]
	}
	C <- table(real, data)
	tp <- C[2, 2]
	fn <- C[2, 1]
	tn <- C[1, 1]
	p <- fnc(tp, tn, fn)
	return(p)
}

#' @title Bootstrap confidence interfals of classifier performances
#'
#' @description Compute bootstrap confidence intervals for various 
#'    classifier performance indices.
#'
#' @param x An data.frame of n rows (subjects) and 2 columns: the first 
#'    column contains the predicted values for an outcome y (i.e., y.hat), 
#'    and the second column contains the observed ("true") values for y.
#' @param rep Number of bootstrap iterations (default = 5000). A high 
#'    number of iterations is needed for reliable estimations. 
#'    If \code{rep} < 5000, estimation errors might occur.
#' @param ci.method Method used for bootstrap confidence interval estimation 
#'    (default = "bca"; i.e., adjusted bootstrap percentile).
#' @param formula Performance index. One among: "f1" (default), "accuracy", 
#'    "sensitivity", "specificity", "ppv" (or "precision"; i.e., positive 
#'    predictive value), "npv" (negative predictive value), "plr" (positive 
#'    likelihood ratio), "nlr" (hegative likelihood ratio), "fpr" (false 
#'    positive rate), "fdr" (false discovery rate), "fnr" (flase negative 
#'    rate), "fnc" (false negative cost = TP/(FN + FP), see 
#'    Fragomeni et al. 2022).
#' @param ... Currently ignored.
#'
#' @import boot
#' @export
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item "boot", an object of class \code{boot}, containing bootstrap 
#'    replicates of the given statistic;
#' \item "ci", an object of class \code{bootci}, containing bootstrap 
#'    confidence interval estimations.
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @references
#' 
#' Fragomeni SM, Moro F, Palluzzi F, Mascilini F, Rufini V, Collarino A, 
#' Inzani F, Giacò L, Scambia G, Testa AC, Garganese G (2022). 
#' Evaluating the risk of inguinal lymph node metastases before surgery 
#' using the Morphonode Predictive Model: a prospective diagnostic study. 
#' Ultrasound xx Xxxxxxxxxx xxx Xxxxxxxxxx. 00(0):000-000.
#' <https://doi.org/00.0000/00000000000000000000>
#' 
#' Davison AC and Hinkley DV (1997). Bootstrap Methods and Their Application, 
#' Chapter 5. Cambridge University Press.
#' 
#' DiCiccio TJ and Efron B (1996). Bootstrap confidence intervals 
#' (with Discussion)._Statistical Science_, *11*, 189-228.
#' 
#' Efron B (1987). Better bootstrap confidence intervals (with 
#' Discussion). Journal of the American Statistical Association, 
#' *82*, 171-200. <https://doi.org/10.2307/2289144>
#' 
#' @seealso See \code{\link[boot]{boot.ci}} for further details.
#' 
#' @examples
#' 
#' \donttest{
#' 
#' # Predicted ultrasound phenotype from the default morphonode RFC1
#' y.hat <- predict(mpm.rfc$rfc$RFC1, mpm.rfc$validation$V1)
#'
#' # Actual ultrasound phenotype values
#' y <- mpm.rfc$validation$V1$y
#'
#' # Input preparation
#' Y <- data.frame(y.hat, y)
#'
#' # F1 score bootstrap confidence intervals
#' F1 <- p.boot(Y)
#' print(F1$boot$t0)              # F1 score observed value
#' print(F1$ci$bca[4:5])          # F1 score bca confidence interval
#' 
#' # Accuracy bootstrap confidence intervals
#' A <- p.boot(Y, formula = "accuracy")
#' print(A$boot$t0)               # Accuracy observed value
#' print(A$ci$bca[4:5])           # Accuracy bca confidence interval
#' 
#' }
#'
p.boot <- function(x, rep = 5000, ci.method = "bca", formula = "f1") {
	if (formula == "f1") {
		R <- boot(data = x, statistic = p.f1, R = rep)
	} else if (formula == "accuracy") {
		R <- boot(data = x, statistic = p.accuracy, R = rep)
	} else if (formula == "sensitivity") {
		R <- boot(data = x, statistic = p.sensitivity, R = rep)
	} else if (formula == "specificity") {
		R <- boot(data = x, statistic = p.specificity, R = rep)
	} else if (formula %in% c("ppv", "precision")) {
		R <- boot(data = x, statistic = p.ppv, R = rep)
	} else if (formula == "npv") {
		R <- boot(data = x, statistic = p.npv, R = rep)
	} else if (formula == "plr") {
		R <- boot(data = x, statistic = p.plr, R = rep)
	} else if (formula == "nlr") {
		R <- boot(data = x, statistic = p.nlr, R = rep)
	} else if (formula == "fpr") {
		R <- boot(data = x, statistic = p.fpr, R = rep)
	} else if (formula == "fdr") {
		R <- boot(data = x, statistic = p.fdr, R = rep)
	} else if (formula == "fnr") {
		R <- boot(data = x, statistic = p.fnr, R = rep)
	} else if (formula == "fnc") {
		R <- boot(data = x, statistic = p.fnc, R = rep)
	} else {
		stop("invalid formula")
	}
	ci <- boot.ci(R, type = ci.method)
	return(list(boot = R, ci = ci))
}

new.entry <- function(x, m) {
	v <- vector()
	n <- nrow(x)
	j <- sample(1:n, m, replace = FALSE)
	for (k in 1:m) {
		v <- c(v, x[j[k], k])
	}
	return(v)
}

#' @title Ultrasound profile simulation
#'
#' @description Ultrasound profile simulation based on a reference dataset.
#'
#' @param reps Number of simulated vectors to be generated (default = 1).
#' @param x Reference ultrasound features data.frame. By default, a 
#'    dataset of 948 simulated profiles (object \code{mpm.us}) is used.
#' @param y An optional numeric value defining the phenotype to be 
#'    simulated: 0 for non-malignant and 1 for malignant (default = NULL).
#' @param signature An optional string defining the metastatic risk 
#'    signature to be simulated. One among: "LMR" (low metastatic risk), 
#'    "MMR" (moderate metastatic risk), "HMR" (high metastatic risk), 
#'    "MET" (metastatic signature). Defaults to NULL.
#' @param features A vector of indices defining the ultasound features 
#'     within the reference data.frame (default = 2:15).
#' @param header A vector of new names for the ultrasound features vector. 
#'     If NULL (default), feature names will be imported from the reference.
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A simulated ultrasound features vector.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{new.profile}} for ultrasound 
#'    profile creation and \code{\link[morphonode]{us.predict}} for 
#'    subject classification.
#' 
#' @examples
#' 
#' # Simulate a generic ultrasound profile
#' u <- us.simulate()
#' print(u)
#' 
#' # Simulate a non-malignant ultrasound profile
#' u0 <- us.simulate(y = 0)
#' print(u0)
#' 
#' # Simulate a malignant ultrasound profile
#' u1 <- us.simulate(y = 1)
#' print(u1)
#' 
#' # Simulate a high metastatic risk profile (single metastasis event)
#' u.hmr <- us.simulate(signature = "HMR")
#' print(u.hmr)
#' 
#' # Simulate a metastatic ultrasound profile (multiple metastasis events)
#' u.met <- us.simulate(signature = "MET")
#' print(u.met)
#'
us.simulate <- function(reps = 1, x = mpm.us, y = NULL, signature = NULL,
                        features = 2:15, header = NULL) {
	if (!is.null(x)) {
		if (!is.null(y)) x <- x[x$y == y,]
		if (!is.null(signature)) x <- x[x$signature == signature,]
		x <- x[, features]
		v <- vector()
		#n <- nrow(x)
		m <- ncol(x)
		while (reps > 0) {
			v <- rbind(v, new.entry(x, m))
			reps <- reps - 1
		}
		if (is.null(header)) colnames(v) <- colnames(x)
		return(v)
	} else {
		stop("Input dataset not found!")
	}
}
