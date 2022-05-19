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

show.mn2env <- function(x = NULL, brief = TRUE) {
	if (is.null(x)) {
		cat("\n### Morphonode 2 - Ultrasound Predictors ###\n\n")
		print(env.us()[, c(1, 4, 5)])
		cat("\n")
	} else {
		
	}
}

new.profile <- function(u = NULL) {
	v <- env.us()
	if (is.null(u)) {
		w <- vector()
		
		for (j in 1:nrow(v)) {
			if (v[j, 4] == "numeric") {
				line <- "is.na(suppressWarnings(as.numeric(u)))"
			} else {
				vals <- c(strsplit(v[j, 4], "")[[1]], "")
				line <- "!(u %in% vals)"
			}
			
			u <- readline(prompt = paste0(c("\n", v[j, 1], " [", v[j, 2], "]:\n"),
			                              collapse = ""))
			
			while (eval(parse(text = line))) {
				message("\nInvalid entry.\nPlease, retry or use Ctrl+C to exit.")
				u <- readline(prompt = paste0(c("\n", v[j, 1],
				                              " [", v[j, 2], "]:\n"),
				                              collapse = ""))
			}
			w <- c(w, u)
		}
		w[w == ""] <- "-1"
		w <- as.numeric(w)
		
	} else {
		w <- as.numeric(u)
	}
	m <- which(w == -1)
	names(w) <- c("shortAxis", "cortical", "hilum", "inflammatoryStroma",
	              "extracapsularSpread", "ecostructure", "FID", "VFL",
	              "corticalThickening", "vascularPattern", "CMID",
	              "shape", "grouping", "colorScore")
	#names(w) <- v$descriptions
	u <- list(ultrasound = w, missing = m)
	return(u)
}

env.us <- function() {
	descriptions <- c("Short axis", "Cortical thickness",
	                  "Nodal core sign",
	                  "Perinodal hyperechogenic ring",
	                  "Cortical interruption",
	                  "Echogenicity", "Focal intranodal deposit",
	                  "Vascular flow localization",
	                  "Cortical thickening",
	                  "Vascular flow architecture pattern",
	                  "Cortical-Medullar interface distortion",
	                  "Shape", "Grouping", "Color score")
	                  #, "Medulla thickness", "Long axis"
	
	units <- c("mm", "mm", "0: absent, 1: present", "0: absent, 1: present",
	           "0: absent, 1: present", "0: homogeneous, 1: inhomogeneous",
	           "0: absent, 1: hyperechoic, 2: hypoechoic, 3: both",
	           "0: non vascularized, 1: central (hilum), 2: peripheral, 3: extranodal, 4: combined",
	           "0: not evaluable, 1: absent, 2: focal, 3: concentric, 4: eccentric",
	           "0: non vascularized, 1: longitudinal axis, 2: scattered, 3: branched, 4: chaotic",
	           "1: absent, 2: focal, 3: diffused, 4: not visible",
	           "1: elliptical, 2: rounded, 3: asymmetrical",
	           "1: assente, 2: intermediate, 3: complete",
	           "1, 2, 3, 4, 5")
	           #, "mm", "mm"
	
	variables <- c("shortAxis", "cortical", "hilum",
	               "inflammatoryStroma", "extracapsularSpread",
	               "ecostructure", "FID", "VFL", "corticalThickening",
	               "vascularPattern", "CMID", "shape", "grouping",
	               "colorScore")
	
	values <- c("numeric", "numeric", "01", "01", "01", "01",
	            "0123", "01234", "01234", "01234", "1234",
	            "123", "123", "12345")
	
	basal <- c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
	
	env <- data.frame(descriptions, units, variables, values, basal)
	return(env)
}

f2n <- function(u, i = 1:2, j = 3:14, adj = -1) {
	n <- as.numeric(u[i])
	f <- as.numeric(u[j]) + adj
	return(c(n, f))
}

f2n.df <- function(x) {
	rows <- rownames(x)
	x <- as.data.frame(t(apply(x, 2, as.numeric)))
	rownames(x) <- rows
	return(x)
}

set.rfcdata <- function(u, levels = NULL, ref = NULL, E = NULL) {
	U <- data.frame(u$ultrasound)
	U$shortAxis <- as.numeric(U$shortAxis)
	U$cortical <- as.numeric(U$cortical)
	U$hilum <- ifelse(U$hilum == 1, 0, 1)
	names(U)[3] <- "hilumAbsence"
	for (j in 3:14) {
		U[, j] <- factor(U[, j], levels = levels[[j]])
	}
	return(list(ultrasound = U, missing = u$missing))
}

set.rbmdata <- function(u, levels, short = 8, cortical = 2,
                        ist = 1, ecs = 1, hab = 0, eco = 1,
                        vp = c(1, 2, 3), vfl = c(2, 3, 4),
                        ct = c(2), fid = c(1, 2, 3), cmid = c(2, 3, 4),
                        shape = c(3), cs = 3, grouping = c(2, 3),
                        ref = NULL, asFactor = FALSE) {
	U <- data.frame(u$ultrasound)
	U$shortAxis <- as.numeric(U$shortAxis)
	U$cortical <- as.numeric(U$cortical)
	U$colorScore <- as.numeric(U$colorScore)
	U <- dichotomize(U, short = short, cortical = cortical,
                     ist = ist, ecs = ecs, hab = hab, eco = eco,
                     vp = vp, vfl = vfl, ct = ct, fid = fid,
                     cmid = cmid, shape = shape, cs = cs,
                     grouping = grouping, asFactor = asFactor)
    return(list(ultrasound = U, missing = u$missing))
}

imp.missing <- function(M, x = NULL, mode = NULL) {
	if (!is.null(x)) M <- rbind(M, x)
	if (mode == "cat") {
		M <- suppressMessages(imputeR::impute(as.matrix(M),
		                                      cFun = "rpartC",
		                                      ini = "majority"))
	} else if (mode == "con") {
		M <- suppressMessages(imputeR::impute(as.matrix(M),
		                                      lmFun = "lassoR",
		                                      ini = "mean"))
	}
	M <- as.data.frame(M$imp)
	return(M)
}

set.missing <- function(v, levels = NULL, ref = NULL, con = 1:2,
                        cat = 3:14, missing = -1, na = NA,
                        asNumeric = TRUE) {
	if (length(v$missing) > 0) {
		if (!is.null(ref)) {
			v$ultrasound[v$ultrasound == missing] <- na
			ref <- rbind(ref, v$ultrasound)
			if (!is.null(levels)) {
				for (j in cat) {
					ref[, j] <- factor(ref[, j], levels = levels[[j]])
				}
			}
			if (v$missing[1] < con[length(con)]) {
				W <- imp.missing(ref[, con], x = v$ultrasound[con],
				                 mode = "con")
				W <- W[nrow(ref) + 1,]
			} else {
				W <- ref[nrow(ref), con]
			}
			if (v$missing[length(v$missing)] > 2) {
				Z <- imp.missing(ref[, cat], x = v$ultrasound[cat],
				                 mode = "cat")
				Z <- Z[nrow(ref) + 1,]
			} else {
				Z <- ref[nrow(ref), cat]
			}
			v$ultrasound <- cbind(W, Z)
			v$missing <- numeric()
		} else if (!is.null(levels)) {
			v$ultrasound[v$missing] <- unlist(lapply(levels[v$missing],
			                                         function(x) x[1]))
		} else {
			v[v == missing] <- na
			if (asNumeric) return(as.numeric(v))
		}
		return(v)
		
	} else {
		warning("no missing values found.")
		return(v)
	}
}

metProbs <- function(is, ecs, h) {
	w <- is + ecs + h
	if (w > 1) {
		return(list(p = 1.00, ci95 = "1.00, 1.00"))
	} else if (is == 1) {
		return(list(p = 0.90, ci95 = "0.69, 0.97"))
	} else if (ecs == 1) {
		return(list(p = 0.89, ci95 = "0.58, 0.95"))
	} else if (h == 0) {
		return(list(p = 0.86, ci95 = "0.62, 0.93"))
	}
}

mmrProb <- function(vp, vfl, eco, ct, fid, cmid, shape, cs, grp, dct = 2) {
	w <- sum(ifelse(vp + vfl > 1, 1, 0), eco, ct, fid, cmid, shape, cs, grp)
	if (w >= dct) {
		return(list(k = "MMR1", p = 0.55, ci95 = "0.46, 0.64", y = 1))
	} else {
		return(list(k = "MMR0", p = 0.15, ci95 = "0.06, 0.25", y = 0))
	}
}

dichotomize <- function(x, short = 8, cortical = 2,
                        ist = 1, ecs = 1, hab = 0, eco = 1,
                        vp = c(1, 2, 3), vfl = c(2, 3, 4),
                        ct = c(2), fid = c(1, 2, 3), cmid = c(2, 3, 4),
                        shape = c(3), cs = 3, grouping = c(2, 3),
                        asFactor = FALSE) {
	
	if (asFactor) {
		a <- "factor("
		b <- ", levels = c(0, 1))"
	} else {
		a <- b <- ""
	}
	v <- list(c("shortAxis", "short", " >= "),
	          c("cortical", "cortical", " >= "),
	          c("inflammatoryStroma", "ist", " == "),
	          c("extracapsularSpread", "ecs", " == "),
	          c("hilum", "hab", " != "),
	          #c("hilumAbsence", "hab", " == "),
	          c("ecostructure", "eco", " == "),
	          c("vascularPattern", "vp", " %in% "),
	          c("VFL", "vfl", " %in% "),
	          c("corticalThickening", "ct", " %in% "),
	          c("FID", "fid", " %in% "),
	          c("CMID", "cmid", " %in% "),
	          c("shape", "shape", " %in% "),
	          c("colorScore", "cs", " >= "),
	          c("grouping", "grouping", " %in% "))
	
	for (j in 1:length(v)) {
		if (!is.null(eval(parse(text = paste0(c("x$", v[[j]][1]),
		                        collapse = ""))))) {
			command <- paste0(c("x$", v[[j]][1], " <- ", a,
								"ifelse(x$", v[[j]][1], v[[j]][3],
								v[[j]][2], ", 1, 0)", b),
								collapse = "")
			#print(command)
			eval(parse(text = command))
		}
	}
	
	return(x)
}

uss <- function(x, dichotomous = FALSE, dct = 2,
                short = 8, cortical = 2, ist = 1, ecs = 1, hab = 0,
                eco = 1, vp = c(1, 2, 3), vfl = c(2, 3, 4),
                ct = c(2), fid = c(1, 2, 3), cmid = c(2, 3, 4),
                shape = c(3), cs = 3, grouping = c(2, 3)) {
	
	if (!dichotomous) x <- dichotomize(x, short, cortical, ist, ecs,
	                                   hab, eco, vp, vfl, ct, fid,
	                                   cmid, shape, cs, grouping)
	signature <- vector()
	prob <- vector()
	CI95 <- vector()
	y.dtr <- vector()
	
	for (j in 1:(nrow(x))) {
		if (x$inflammatoryStroma[j] == 1 | x$extracapsularSpread[j] == 1 | x$hilum[j] == 0) {
			k <- "MET"
			p <- metProbs(x$inflammatoryStroma[j],
			              x$extracapsularSpread[j],
			              x$hilum[j])
			ci95 <- p$ci95
			p <- p$p
			y <- 1
		} else {
			if (x$cortical[j] == 0) {
				k <- "LMR"
				p <- 0.04
				ci95 <- "0.00, 0.10"
				y <- 0
			} else if (x$cortical[j] == 1 & x$shortAxis[j] == 0) {
				mmr <- mmrProb(x$vascularPattern[j], x$VFL[j],
				               x$ecostructure[j], x$corticalThickening[j],
				               x$FID[j], x$CMID[j], x$shape[j],
				               x$colorScore[j], x$grouping[j],
				               dct = dct)
				k <- mmr$k
				p <- mmr$p
				ci95 <- mmr$ci95
				y <- mmr$y
			} else if (x$cortical[j] == 1 & x$shortAxis[j] == 1 & x$vascularPattern[j] == 0) {
				mmr <- mmrProb(x$vascularPattern[j], x$VFL[j],
				               x$ecostructure[j], x$corticalThickening[j],
				               x$FID[j], x$CMID[j], x$shape[j],
				               x$colorScore[j], x$grouping[j],
				               dct = dct)
				k <- mmr$k
				p <- mmr$p
				ci95 <- mmr$ci95
				y <- mmr$y
			} else if (x$cortical[j] == 1 & x$shortAxis[j] == 1 & x$vascularPattern[j] == 1 & x$corticalThickening[j] == 0 & x$VFL[j] == 0) {
				mmr <- mmrProb(x$vascularPattern[j], x$VFL[j],
				               x$ecostructure[j], x$corticalThickening[j],
				               x$FID[j], x$CMID[j], x$shape[j],
				               x$colorScore[j], x$grouping[j],
				               dct = dct)
				k <- mmr$k
				p <- mmr$p
				ci95 <- mmr$ci95
				y <- mmr$y
			} else if (x$cortical[j] == 1 & x$shortAxis[j] == 1 & x$vascularPattern[j] == 1 & x$corticalThickening[j] == 0 & x$VFL[j] == 1) {
				k <- "HMR"
				p <- 0.81
				ci95 <- "0.52, 0.90"
				y <- 1
			} else if (x$cortical[j] == 1 & x$shortAxis[j] == 1 & x$vascularPattern[j] == 1 & x$corticalThickening[j] == 1) {
				k <- "HMR"
				p <- 0.81
				ci95 <- "0.52, 0.90"
				y <- 1
			} else {
				k <- "None"
				p <- -1
				ci95 <- "-1, -1"
				y <- -1
			}
		}
		signature <- c(signature, k)
		prob <- c(prob, p)
		CI95 <- c(CI95, ci95)
		y.dtr <- c(y.dtr, y)
	}
	return(list(signature = signature, p = prob, ci95 = CI95, y = y.dtr))
}

topsim <- function(v, x = mn2.us, f = "cosine", k = 5, p = 0.7, features = 2:15, check.data = TRUE) {
	if (check.data) v <- f2n(v)
	x$R <- apply(x[, features], 1, function(w) similarity(v, w, f = f))
	x <- x[order(x$R, decreasing = TRUE),]
	if (f != "euclidean" & p > 0 & p <= 1) x <- x[x$R >= p,]
	if (k > 0) x <- x[1:k,]
	return(x)
}

ranksim <- function(u, v = NULL, x = mn2.us, k = 5, p = 0.7,
                    j = 2:15, d = c(2:6, 9, 10, 11),
                    signature = NULL, check.data = TRUE,
                    orderbyDistance = FALSE) {
	if (!is.null(signature)) {
		x <- x[x$signature == signature,]
	}
	if (check.data) u <- f2n(u)
	if (is.null(v)) {
		x$R <- apply(x[, j], 1, function(w) similarity(u[j - 1], w,
		                                               f = "cosine"))
	} else {
		z <- dichotomize(x)
		x$R <- apply(z[, j], 1, function(w) similarity(v[j - 1], w,
		                                               f = "jaccard"))
	}
	x <- x[x$R >= p,]
	x$D <- apply(x[, d], 1, function(w) euclidean(u[d - 1], w))
	if (orderbyDistance) {
		x <- x[order(x$D),]
	} else {
		x <- x[order(x$R, decreasing = TRUE),]
	}
	x <- x[1:k,]
	return(x)
}

rfc.predict <- function(u, rfc, recover = TRUE) {
	p <- factor(, levels = c(0, 1))
	for (j in 1:length(rfc)) {
		p <- c(p, predict(rfc[[j]], u))
	}
	y.hat <- p.mode(p)
	#print(y.hat)
	oob <- unlist(lapply(rfc,
		          function(x) x$err.rate[length(x$err.rate[, 1]), 1]))
	if (recover & length(y.hat) > 1) {
		E <- order(oob, decreasing = TRUE)
		i <- c(1)
		while (length(y.hat) > 1 & length(p[E[-i]]) > 1) {
			y.hat <- p.mode(p[E[-i]])
			#print(i)
			#print(E[-i])
			i <- c(i, i[length(i)] + 1)
		}
		oob <- E[-i]
		p <- p[oob]
		L <- length(oob) + 1
		if (L > 1) {
			warning(paste0("Prediction achieved using the best ", L,
			               " classifiers."))
		} else {
			warning(paste0("Prediction achieved using the best classifier."))
		}
	}
	names(p) <- NULL
	return(list(y.hat = as.numeric(y.hat), decisions = p, oob.err = oob))
}

blanks <- function(x, blank = " ", chars = 5, digits = -1) {
	if (digits > 0) {
		x <- format(round(as.numeric(x), digits = digits), nsmall = digits)
	}
	chars <- chars - nchar(x)
	chars <- ifelse(chars >= 0, chars, 0)
	return(paste0(c(rep(blank, chars), x), collapse = ""))
}

printline <- function(x, j, w, y, m, c = 5, d = 2) {
	if (j < 10) {
		line = "| ["
	} else {
		line = "|["
	}
	line <- paste0(c(line, as.character(j), "] w = ", w, " | ",
	               blanks(x[1], chars = c, digits = d), " | ",
	               blanks(x[2], chars = c, digits = d), " | ",
	               x[3], " | ", x[4], " | ", x[5], " | ", x[6], " | ",
	               x[7], " | ", x[8], " | ", x[9], " | ", x[10], " | ",
	               x[11], " | ", x[12], " | ", x[13], " | ", x[14], " | ",
	               "Y = ", y, " | ", m, " |\n"),
	               collapse = "")
	return(line)
}

printlines <- function(x, c = 5, d = 2, score.digits = 3, features = 2:15) {
	w <- c()
	for (j in 1:nrow(x)) {
		if (j < 10) {
			line = "| ["
		} else {
			line = "|["
		}
		score <- format(round(as.numeric(x$R[j]), digits = score.digits),
		                nsmall = score.digits)
		w <- c(w, printline(x = x[j, features], j = j, w = score,
		                    y = x$y[j], m = x$signature[j],
		                    c = c, d = d))
	}
	return(paste0(w, collapse = ""))
}

printout <- function(u, x, y, E, p, mrs, k = 5, wmax = 1,
                     sim = "cosine", features = 2:15) {
	x <- x[order(x$R, decreasing = TRUE),]
	E <- format(round(as.numeric(E), digits = 3), nsmall = 3)
	p <- format(round(as.numeric(p), digits = 3), nsmall = 3)
	wmax <- format(round(as.numeric(wmax), digits = 3), nsmall = 3)
	if (mrs %in% c("MMR0", "MMR1")) mrs <- "MMR"
	
	if (y == 1) {
		msg.y <- "MALIGNANT (Y = 1)"
		spc.y <- 4
	} else {
		msg.y <- "NOT MALIGNANT (Y = 0)"
		spc.y <- 0
	}
	
	if (E < 1) {
		msg.e <- "     (cutoff: E < 1)"
	} else {
		msg.e <- "     Warning: E >= 1"
	}
	
	if (p > 0.29) {
		msg.p <- "HIGH (> 0.29)"
		spc.p <- 7
	} else if (p < 0.23) {
		msg.p <- "LOW (< 0.23)"
		spc.p <- 8
	} else {
		msg.p <- "MODERATE (0.23-0.29)"
		spc.p <- 0
	}
	
	if (mrs == "MET") {
		msg.m <- "MET (metastatic signature)"
		spc.m <- 4
	} else if (mrs == "HMR") {
		msg.m <- "HMR (high metastatic risk)"
		spc.m <- 4
	} else if (mrs %in% c("MMR", "MMR0", "MMR1")) {
		msg.m <- "MMR (moderate metastatic risk)"
		spc.m <- 0
	} else if (mrs == "LMR") {
		msg.m <- "LMR (low metastatic risk)"
		spc.m <- 5
	}
	
	if (sim == "cosine") {
		hyp.s <- 45
	} else {
		hyp.s <- 44
	}
	
	hyphen <- paste0(rep("-", 91), collapse = "")
	
	msg <- paste0(c("\n       Morphonode Predictive Model output",
					"\n# ", hyphen, " #",
					"\n|     Prediction (Morphonode-RFC): ", msg.y,
					rep(" ", 38 + spc.y), "|",
					"\n|      Estimated prediction error: ", E,
					msg.e, rep(" ", 34), "|",
					"\n|  Risk estimate (Morphonode-RBM): ", p,
					rep(" ", 54), "|",
					"\n|                      Risk level: ", msg.p,
					rep(" ", 39 + spc.p), "|",
					"\n|       Signature (Morphonode-DT): ", msg.m,
					rep(" ", 29 + spc.m), "|",
					"\n|", rep(" ", 93), "|",
					"\n|--- Input profile ", rep("-", 75), "|\n",
					"|", rep(" ", 93), "|\n",
					printline(u, j = 0, w = wmax, y = y, m = mrs),
					"|", rep(" ", 93), "|\n",
					"|--- Top-similar profiles (w: ", sim, " similarity) ",
					rep("-", hyp.s), "|\n",
					"|", rep(" ", 93), "|\n",
					printlines(x, features = features),
					"# ", hyphen, " #\n\n"),
					collapse = "")
	cat(msg)
}

us.predict <- function(x, f = "cosine", levels = NULL, ref = NULL,
                       rfc = NULL, rbm = NULL, k = 5, features = 2:15,
                       orderbyDistance = FALSE, uncertainty = "loss",
                       b0 = c(0, 0.028, 0.013), b = 0.055, rho = 0.9,
                       verbose = TRUE, wmax = 1) {
	
	if (is.null(levels)) levels <- RFC5$rfc$RFC1$forest$xlevels
	if (is.null(ref)) ref <- mn2.us[, 2:15]
	if (is.null(rfc)) rfc <- RFC5$rfc
	if (is.null(rbm)) rbm <- mn2.rbm$fit
	
	if (length(x$missing) > 0) {
		suppressWarnings(x <- set.missing(x, levels = levels, ref = ref))
	} else {
		x$ultrasound <- data.frame(t(x$ultrasound))
	}
	suppressWarnings(u <- set.rfcdata(x, levels = levels, ref = ref))
	suppressWarnings(v <- set.rbmdata(x, levels = levels, ref = ref,
	                                  asFactor = TRUE))
	p <- predict(rbm, v$ultrasound, type = "response")
	mrs <- uss(f2n.df(v$ultrasound), dichotomous = TRUE)$signature
	
	if (mrs == "MET") {
		i <- 3
		if (f == "cosine") {
			R <- ranksim(u$ultrasound, k = k, j = features,
			             signature = mrs,
						 orderbyDistance = orderbyDistance)
		} else if (f == "jaccard") {
			R <- ranksim(u$ultrasound, v$ultrasound, k = k,
			             j = features, signature = mrs,
			             orderbyDistance = orderbyDistance)
		} else {
			stop("invalid similarity function")
		}
	} else {
		if (mrs == "HMR") i <- 2 else i <- 1
		if (f == "cosine") {
			R <- ranksim(u$ultrasound, k = k, j = features,
						 orderbyDistance = orderbyDistance)
		} else if (f == "jaccard") {
			R <- ranksim(u$ultrasound, v$ultrasound, k = k, j = features,
						 orderbyDistance = orderbyDistance)
		} else {
			stop("invalid similarity function")
		}
	}
	if (uncertainty == "loss") {
		y <- as.numeric(p.mode(R$y[R$R >= rho]))
		u$ultrasound$E <- b0[i] + b*loss(p, y, method = "logistic")$L
	} else if (uncertainty == "similarity") {
		if (nrow(R) > 3) {
			u$ultrasound$E <- mean(R$E[1:3])
		} else {
			u$ultrasound$E <- mean(R$E)
		}
	} else {
		stop("invalid error estimation method")
	}
	y <- rfc.predict(u$ultrasound, rfc = rfc)
	if (verbose) {
		z <- f2n(u$ultrasound)
		z[3] <- ifelse(z[3] == 1, 0, 1)
		printout(z, R, y$y.hat, u$ultrasound$E, p, mrs, k = k,
		         wmax = wmax, sim = f, features = features)
	}
	return(list(prediction = y, E = u$ultrasound$E, p = p,
	            signature = mrs, profiles = R))
}






