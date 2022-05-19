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

rfc <- function(data, status, k = 5, method = "CV", brier = 0) {
	
	L <- GenerateLearningsets(n = nrow(data), status,
							  method = method,
							  fold = k,
							  strat = TRUE)
	
	C <- classification(as.matrix(data), status,
	                    learningsets = L,
	                    classifier = "rfCMA")
	
	return(C)
}

p.mode <- function(x) {
	f <- table(x)
	return(names(f)[which(f == max(f))])
}

brier <- function(data, status, k = 5, method = "CV") {
	C <- rfc(data = data, status = status, k = k, method = method)
	E <- evaluation(C, scheme = "observationwise",
	                   measure = "brier score")@score
	return(E)
}

loss <- function(p, y, method = "logistic") {
	if (method == "logistic") {
		L = -1*(y*log(p) + (1 - y)*log(1 - p))
		E = sum(L)/length(p)
	} else if (method == "rmse") {
		L = (p - y)^2
		E = sqrt(sum(L)/length(p))
	} else {
		stop("invalid method")
	}
	return(list(L = L, E = E))
}

buildPredictor <- function(model, data, vset = NULL, y = NULL, n = 10000, m = 3, asFactor = TRUE) {
	RFC <- randomForest(model, data = data,
						ntree = n,
						mtry = m,
						keep.forest = TRUE,
						proximity = TRUE,
						importance = TRUE)

	if (is.null(vset)) {
		return(RFC)
	} else {
		C <- predict(RFC, vset)
		P <- performance(vset[, y], C)
		return(list(RFC = RFC, performance = P))
	}
}

vpart <- function(data, p, shuffle = FALSE) {
	if (shuffle) {
		data <- data[sample(1:nrow(data), replace = FALSE),]
	}
	j <- sample(nrow(data), size = p*nrow(data), replace = FALSE)
	X <- data[j,]
	V <- data[-j,]
	return(list(training.set = X, validation.set = V))
}

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
	P <- TP/(TP + FP)                 # Precision
	F1 <- (2*Se*P)/(Se + P)           # F1-score
	
	W <- list(ctable = C, sensitivity = Se, specificity = Sp,
	          precision = P, F1 = F1, accuracy = A)
	return(W)
}

euclidean <- function(x, y) {
	return(sqrt(sum((y - x)^2)))
}

jaccard <- function(x, y, dichotomize.y = FALSE) {
	return(sum(1*(y == x))/length(x))
}

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

boot.se <- function(fit, boot, probs = c(0.025, 0.975)) {
	conf.level <- probs[2] - probs[1]
	b0 <- summary(fit)$coefficients
	coef <- NULL
	for (i in 1:nrow(b0)) {
		se.boot <- sd(boot[, i])
		#hist(z)
		#Sys.sleep(5)
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

glmBoot <- function(x, model = NULL, y = NULL, n = NULL, family = "binomial", link = "", nrep = 5000) {
	if (is.null(model)) {
		if (is.null(y)) {
			stop("Please, define either 'model' or 'y'.")
		}
		j <- which(colnames(x) == y)
		if (!is.null(n)) {
			k <- which(colnames(x) == n)
			model <- paste0(c(model, " + offset(log(", n, "))"), collapse = "")
			link <- "log"
			model <- paste0(colnames(x[, c(-j, -k)]), collapse = " + ")
		} else {
			model <- paste0(colnames(x[, -j]), collapse = " + ")
		}
		model <- paste0(c(y, "~", model), collapse = " ")
	}
	fam.check <- FALSE
	if (family %in% c("gaussian", "binomial", "poisson",
	                  "quasibinomial", "quasipoisson")) {
		fam.check <- TRUE
	}
	if (link != "" & fam.check) {
		link <- paste0(c("(link = ", link, ")"), collapse = "")
	}
	
	if (family == "nb") {
		if (nrep > 1) {
			boot <- do(nrep) * glm.nb(model, data = resample(x), link = "log")
		}
		fit <- glm.nb(model, data = x, link = nb.link)
	} else if (fam.check) {
		family <- parse(text = paste0(c(family, link), collapse = ""))
		if (nrep > 1) {
			boot <- do(nrep) * glm(model, data = resample(x), family = eval(family))
		}
		fit <- glm(model, data = x, family = eval(family))
	} else {
		stop("An unavailable family was chosen.")
	}
	
	if (nrep > 1) {
		coef <- boot.se(fit, boot)
	} else {
		coef <- summary(fit)$coefficients
		#rownames(coef)[1] <- "Intercept"
		coef <- data.frame(Variable = rownames(coef),
						   Estimate = coef[, 1],
						   SE = coef[, 2],
						   z = coef[, 3],
						   P = coef[, 4])
		#rownames(coef)[1] <- NULL
	}
	
	return(list(coef = coef, model = model, fit = fit))
}

glm.exp <- function(fit, digits = 3) {
	
	if (fit$family$family == "binomial") {
		what <- "OR"
	} else if (fit$family$family == "poisson") {
		what <- "Rate"
    } else {
		stop("Only available for logistic or poisson models.")
    }
    
    coef <- stats::coef(fit)
    CI <- stats::confint(fit)
    est <- cbind(coef = coef, CI)
    est.exp <- round(exp(est), digits)
    colnames(est.exp)[1] <- what
    return(est.exp)
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
	return(list(R = R, ci = ci))
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

simulate.us <- function(reps = 1, x = mn2.us, y = NULL, signature = NULL,
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

