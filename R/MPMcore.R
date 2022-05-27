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

#' @title Create a new ultrasound profile
#'
#' @description Creates a new ultrasound profile, used as input for the 
#' \code{\link[morphonode]{us.predict}} function.
#' 
#' @param u A vector of length 14, corresponding to the input ultrasound 
#'    profile. The order of the elements in u must be: 
#'    short axis [numeric], cortical thickness [numeric], 
#'    nodal core sign (hilum) [dichotomous], 
#'    perinodal hyperechogenic ring [dichotomous], 
#'    cortical interruption [dichotomous], echogenicity [dichotomous], 
#'    focal intranodal deposit [categorical], 
#'    vascular flow localization [categorical], 
#'    cortical thickening [categorical], 
#'    vascular flow architecture pattern [categorical], 
#'    cortical-medullar interface distortion [categorical], 
#'    shape [categorical], grouping [categorical], color score [ordinal].
#'    Value -1 can be entered for missing values.
#'    If \code{u = NULL}, a guided interactive prompt will be launched.
#'    In the interactive mode, typing return without entering any value 
#'    will introduce a missing value (this is not possible for short axis 
#'    and cortical thickness, since they are necessary for a reliable 
#'    prediction).
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item "ultrasound", ultrasound features vector;
#' \item "missing", missing values vector (empty, if no missing values 
#'                  are found).
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{show.mpmenv}} for ultrasound 
#' variable description.
#' See also \code{\link[morphonode]{us.predict}} to launch the full 
#'    Morphonode Predictive Model suite.
#'
#' @examples
#'
#' ### Profile with missing data
#'
#' u <- newProfile(c(10.0, 6.3, 1, 0, 0, 0, -1, 1, 2, 2, 3, -1, -1, -1))
#' print(u)
#'
#' ### New profile from simulated data
#' 
#' # High metastatic risk profile
#' u.hmr <- newProfile(us.simulate(signature = "HMR"))
#' print(u.hmr)
#'
#' # Low metastatic risk profile
#' u.lmr <- newProfile(us.simulate(signature = "LMR"))
#' print(u.lmr)
#'
#' # Malignant profile
#' u1 <- newProfile(us.simulate(y = 1))
#' print(u1)
#'
#' # Non-malignant profile
#' u0 <- newProfile(us.simulate(y = 0))
#' print(u0)
#'
newProfile <- function(u = NULL, ...) {
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

#' @title Morphonode ultrasound variables description
#'
#' @description Shows a description of the utrasound variables used 
#' by the **morphonode** package.
#'
#' @param brief Logical value enabling a compact view (default = TRUE).
#' @param ... Currently ignored.
#'
#' @export
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @examples
#' show.mpmenv()
#'
show.mpmenv <- function(brief = TRUE, ...) {
	if (brief) {
		cat("\n### Morphonode 2 - Ultrasound Predictors ###\n\n")
		print(env.us()[, c(1, 4, 5)])
		cat("\n")
	} else {
		print(env.us())
	}
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

#' @title Impute missing values
#'
#' @description Wrapper for the \code{imputeR} function 
#' \code{\link[imputeR]{impute}}.
#' 
#' @param M A matrix or data.frame containing missing values.
#' @param x An optional vector that will be attached to M. 
#'    This can be useful if data with missing values can be attached to 
#'    a reference dataset.
#' @param mode Either "cat" (cateorical variables) or 
#'    "con" (continuous variables).
#' @param ... Currently ignored.
#'
#' @importFrom imputeR impute
#' @export
#'
#' @return A data.frame with imputed missing values.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @examples
#' 
#' # Sample 30 subjects from the morphonode simulated dataset
#' data <- sample(mpm.us, 30, replace = FALSE, prob = NULL)[, 2:15]
#' 
#' # Entries with missing values
#' missing <- matrix(c(10.0, 6.3, 1, 0, 0, 0, 0, 1, NA, 2, NA, 2, 3, NA,
#'                      6.4, 2.1, 1, 0, 0, 0, 0, 1, NA, 2, NA, 1, 1, NA),
#'                   nrow = 2, byrow = TRUE)
#' colnames(missing) <- colnames(data)
#' 
#' # Defining categorical subset
#' data.cat <- data.frame(apply(data[, 3:14], 2, factor))
#' 
#' # Imputing missing values
#' data.cat <- imp.missing(data.cat, x = missing[, 3:14], mode = "cat")
#' print(tail(data.cat))
#'
imp.missing <- function(M, x = NULL, mode = NULL, ...) {
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

#' @title Cope with possible missing data in the ultrasound profile
#'
#' @description This function imputes missing data in the ultrasound 
#' profile, creating a new profile with imputed missing values.
#' If no missing values are found, it will simply send a message and 
#' return the input profile.
#' 
#' @param v An ultrasound profile generated by 
#'    \code{\link[morphonode]{newProfile}}, i.e., a list of two objects, 
#'    containing a vector of length 14 (\code{list$ultasound}), 
#'    corresponding to the input ultrasound profile (14 ultrasound 
#'    variables or "features"), and a vector of missing value indices 
#'    (\code{list$missing}).
#' @param levels A list of length 14, corresponding to the levels of 
#'    each ultrasound variable. Needed for categorical variables (factors); 
#'    for continuous variables, it should assume the nominal value of 0. 
#'    The default levels variable \code{mpm.levels} can be used.
#' @param ref A data.frame representing the reference dataset. The 
#'    ultrasound profile will be attached to the reference dataset 
#'    before the imputation. This argument is required to impute missing 
#'    features.
#' @param con Vector of the indices corresponding to continuous variables 
#'    in the \code{list$ultasound} vector (default = 1:2).
#' @param cat Vector of the indices corresponding to categorical variables 
#'    in the \code{list$ultasound} vector (default = 3:14).
#' @param missing Value used to mark missing data (default = -1).
#' @param na Value used for "not available" data (default = NA). This 
#'    will be used to dubstitute \code{missing} within the ultrasound 
#'    vector before the imputation.
#' @param asNumeric Logical value used to convert every value in the 
#'    ultrasound vector to class "numeric". This argument is used only 
#'    if \code{ref = NULL} and \code{levels = NULL}.
#' @param ... Currently ignored.
#'
#' @details Automatic imputation is necessary to improve RFC-based 
#'    (malignancy prediction) and RBM-based (metastatic risk evaluation) 
#'    estimations. Imputation is currently forbidden for short axis and 
#'    cortical thickness (i.e., the first two ultrasound features), since 
#'    they have a critical role in the prediction, estimation and 
#'    signature detection processes. Hence, their actual value must be 
#'    entered for a reliable prediction. Although permitted, the imputation 
#'    is discouraged for the following three features: 
#'    nodal core sign (i.e., hilum presence), perinodal hyperechogenic 
#'    ring (i.e., the presence of inflammatory stroma), and 
#'    cortical interruption (i.e., extracapsular spread).
#'    These features define a strongly metastatic profile with possible 
#'    multiple metastases (i.e., the "MET" signature) that are hardly 
#'    imputable from the other ultrasound variables.
#'
#' @importFrom imputeR impute
#' @export
#'
#' @return An ultrasound profile with imputed missing values.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{set.rfcdata}} for random forest 
#' (RFC) data preparation and \code{\link[morphonode]{set.rbmdata}} for 
#' robust binomial model (RBM) data preparation.
#'
#' @examples
#' 
#' # Create an ultrasound profile with missing values
#' 
#' u <- newProfile(c(10.0, 6.3, 1, 0, 0, 0, -1, 1, 2, 2, 3, 1, -1, -1))
#' print(u)
#' 
#' # Fix missing values with the default simulated dataset as reference
#' # (ultrasound features only: \code{mpm.us} attributes 2 to 15).
#' # Default levels are provided by the \code{mpm.levels} object.
#' 
#' v <- set.missing(u, ref = mpm.us[, 2:15], levels = mpm.levels)
#' print(v)
#'
set.missing <- function(v, ref = NULL, levels = NULL, con = 1:2,
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
		message("  No missing values found.")
		return(v)
	}
}

#' @title Ultrasound profile preparation for random forest classification
#'
#' @description Prepare a new ultrasound profile for RFC prediction. 
#'    This function includes missing values check and fix (see 
#'    \code{\link[morphonode]{set.missing}}).
#'
#' @param u New ultrasound profile generated by 
#'    \code{\link[morphonode]{newProfile}}.
#' @param levels A list of length 14, corresponding to the levels of 
#'    each ultrasound variable. Needed for categorical variables (factors); 
#'    for continuous variables, it should assume the nominal value of 0. 
#'    The default levels variable \code{\link[morphonode]{mpm.levels}} 
#'    can be used.
#' @param ref Reference ultrasound features dataset as a (n, 14) 
#'    data.frame, with n being the number of subjects (rows). 
#'    The default simulated dataset \code{\link[morphonode]{mpm.us}} 
#'    can be used.
#' @param ... Currently ignored.
#'
#' @import randomForest
#' @export
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item "ultrasound", ultrasound features vector;
#' \item "missing", indices of missing values (empty if no missing 
#'    values are found).
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{set.rbmdata}} for robust binomial 
#'    model data preparation. 
#'    See also \code{\link[morphonode]{us.predict}} to launch all 
#'    morphonode modules at once.
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
#' @examples
#' 
#' # Prepare a simulated malignant ultrasound profile
#' x <- newProfile(us.simulate(y = 1))
#' 
#' # Set the new profile for RFC prediction
#' u <- set.rfcdata(x, ref = mpm.us[, 2:15], levels = mpm.levels)
#' print(u)
#'
set.rfcdata <- function(u, levels = NULL, ref = NULL) {
	if (length(u$missing) > 0) {
		suppressMessages(u <- set.missing(u, levels = levels, ref = ref))
	} else {
		u$ultrasound <- data.frame(t(u$ultrasound))
	}
	U <- data.frame(u$ultrasound)
	U$shortAxis <- as.numeric(U$shortAxis)
	U$cortical <- as.numeric(U$cortical)
	#U$hilum <- ifelse(U$hilum == 1, 0, 1)
	#names(U)[3] <- "hilumAbsence"
	for (j in 3:14) {
		U[, j] <- factor(U[, j], levels = levels[[j]])
	}
	return(list(ultrasound = U, missing = u$missing))
}

#' @title Ultrasound profile preparation for metastatic risk estimation
#'
#' @description Prepare a new ultrasound profile for metastatic risk 
#'    estimation using robust binomial modeling.
#'    This function includes missing values check and fix (see 
#'    \code{\link[morphonode]{set.missing}}).
#'
#' @param u New ultrasound profile generated by 
#'    \code{\link[morphonode]{newProfile}}.
#' @param levels A list of length 14, corresponding to the levels of 
#'    each ultrasound variable. Needed for categorical variables (factors); 
#'    for continuous variables, it should assume the nominal value of 0. 
#'    The default levels variable \code{\link[morphonode]{mpm.levels}} 
#'    can be used.
#' @param short Numeric value corresponding to the short axis cutoff 
#'    in millimeters (default = 8).
#' @param cortical Numeric value corresponding to the cortical thickness 
#'    cutoff in millimeters (default = 2).
#' @param ist Dichotomous value {0, 1} for the presence of inflammatory 
#'    stroma (perinodal hyperechogenic ring; default = 1).
#' @param ecs Dichotomous value {0, 1} for the presence of extracapsular 
#'    spread (cortical interruption; default = 1).
#' @param hab Dichotomous value {0, 1} for the absence of the hilum 
#'    (nodal core sign; default = 0).
#' @param eco Dichotomous value {0, 1} for heterogeneous echogenicity 
#'    (echostructure; default = 1).
#' @param vp Categorical value (integers from 0 to 4) associated to a 
#'    high-risk vascular flow architecture pattern (default = c(1, 2, 3)).
#' @param vfl Categorical value (integers from 0 to 4) associated to a 
#'    high-risk vascular flow localization (default = c(2, 3, 4)).
#' @param ct Categorical value (integers from 0 to 4) associated to a 
#'    high-risk cortical thickening (default = 2).
#' @param fid Categorical value (integers from 0 to 3) associated to a 
#'    high-risk focal intranodal deposit (default = c(1, 2, 3)).
#' @param cmid Categorical value (integers from 0 to 4) associated to a 
#'    high-risk cortical-medullar interface distortion (default = c(2, 3, 4)).
#' @param shape Categorical value (integers from 1 to 3) associated to a 
#'    high-risk shape (default = 3).
#' @param cs Ordinal value (integers from 1 to 5) associated to a 
#'    high-risk color score (default = 3).
#' @param grouping Categorical value (integers from 1 to 3) associated 
#'    to a high-risk grouping (default = c(2, 3)).
#' @param ref Reference ultrasound features dataset as a (n, 14) 
#'    data.frame, with n being the number of subjects (rows). 
#'    The default simulated dataset \code{\link[morphonode]{mpm.us}} 
#'    can be used.
#' @param asFactor Logical value. If TRUE, data.frame columns are converted 
#'    to factors (default = FALSE).
#' @param ... Currently ignored.
#'
#' @export
#'
#' @return A list of 2 objects:
#' \enumerate{
#' \item "ultrasound", dichotomized ultrasound features vector;
#' \item "missing", indices of missing values (empty if no missing 
#'    values are found).
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{set.rfcdata}} for random forest  
#'    classifier data preparation. 
#'    See also \code{\link[morphonode]{us.predict}} to launch all 
#'    morphonode modules at once.
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
#' @examples
#' 
#' # Prepare a simulated malignant ultrasound profile
#' x <- newProfile(us.simulate(y = 1))
#' print(x)
#' 
#' # Set the new profile for RFC prediction
#' v <- set.rbmdata(x, ref = mpm.us[, 2:15], levels = mpm.levels)
#' print(v)
#'
set.rbmdata <- function(u, levels, short = 8, cortical = 2,
                        ist = 1, ecs = 1, hab = 0, eco = 1,
                        vp = c(1, 2, 3), vfl = c(2, 3, 4),
                        ct = c(2), fid = c(1, 2, 3), cmid = c(2, 3, 4),
                        shape = c(3), cs = 3, grouping = c(2, 3),
                        ref = NULL, asFactor = FALSE) {
	if (length(u$missing) > 0) {
		suppressMessages(u <- set.missing(u, levels = levels, ref = ref))
	} else {
		u$ultrasound <- data.frame(t(u$ultrasound))
	}
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

#' @title Generates a dichotomous ultrasound feature data.frame
#'
#' @description Convert an ultrasound feature data.frame (each row is 
#'    an ultrasound vector) into a dichotomous data.frame. 
#'    Defalut dichotomization cutoffs are computed as described in 
#'    Fragomeni et al. (2022).
#' 
#' @param x An (n, 14) ultrasound features data.frame, where n is the 
#'    number of subjects.
#' @param short Numeric value corresponding to the short axis cutoff 
#'    in millimeters (default = 8).
#' @param cortical Numeric value corresponding to the cortical thickness 
#'    cutoff in millimeters (default = 2).
#' @param ist Dichotomous value {0, 1} for the presence of inflammatory 
#'    stroma (perinodal hyperechogenic ring; default = 1).
#' @param ecs Dichotomous value {0, 1} for the presence of extracapsular 
#'    spread (cortical interruption; default = 1).
#' @param hab Dichotomous value {0, 1} for the absence of the hilum 
#'    (nodal core sign; default = 0).
#' @param eco Dichotomous value {0, 1} for heterogeneous echogenicity 
#'    (echostructure; default = 1).
#' @param vp Categorical value (integers from 0 to 4) associated to a 
#'    high-risk vascular flow architecture pattern (default = c(1, 2, 3)).
#' @param vfl Categorical value (integers from 0 to 4) associated to a 
#'    high-risk vascular flow localization (default = c(2, 3, 4)).
#' @param ct Categorical value (integers from 0 to 4) associated to a 
#'    high-risk cortical thickening (default = 2).
#' @param fid Categorical value (integers from 0 to 3) associated to a 
#'    high-risk focal intranodal deposit (default = c(1, 2, 3)).
#' @param cmid Categorical value (integers from 0 to 4) associated to a 
#'    high-risk cortical-medullar interface distortion (default = c(2, 3, 4)).
#' @param shape Categorical value (integers from 1 to 3) associated to a 
#'    high-risk shape (default = 3).
#' @param cs Ordinal value (integers from 1 to 5) associated to a 
#'    high-risk color score (default = 3).
#' @param grouping Categorical value (integers from 1 to 3) associated 
#'    to a high-risk grouping (default = c(2, 3)).
#' @param asFactor Logical value. If TRUE, data.frame columns are converted 
#'    to factors (default = FALSE).
#' @param ... Currently ignored.
#'
#' @details This function is internally used to estimate the malignancy 
#'    risk through the robust binomial model (Morphonode-RBM).
#'    Dichotomization is performed to avoid the extremely low frequancy 
#'    of some levels in categorical variables.
#'
#' @export
#'
#' @return An ultrasound profile with imputed missing values.
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{uss}} for metastatic risk 
#'    signature detection (Morphonode-DT module).
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
#' @examples
#' 
#' # Create a dichotomous version of subjects with metastatic signature, 
#' # from the default simulated dataset.
#' 
#' M <- dichotomize(mpm.us[mpm.us$signature == "MET", 2:15])
#' print(head(M))
#'
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

#' @title Ultrasound signatures (uss) detection
#'
#' @description This function implements the decision tree described in 
#'    Fragomeni et al. 2022 (Morphonode-DT module) to detect metastatic 
#'    risk signatures (MRSs). 
#'    The first split detects subjects with a "MET" (metastatic) signature. 
#'    These individuals show at least one among the three metastatic 
#'    markers (see details) and a high risk of malignancy (86-100%), 
#'    usually coming with multiple metastatic lymph nodes. 
#'    The other decision tree branches are defined on the base of five 
#'    key ultrasound features (see details), with a malignancy risk (MR)
#'    signature ranging from low ("LMR", 0-10%) to moderate ("MMR", 6-25%) 
#'    to high ("HMR", 52-90%), implying a single metastatic lymph node 
#'    in most malignancies.
#'    The main goal of MRSs is to predict single-metastatic event 
#'    malignancies (HMR) and multiple-metastatic event malignancies (MET). 
#' 
#' @param x An (n, 14) ultrasound features data.frame, where n is the 
#'    number of subjects.
#' @param dichotomous Logical value. The first step of this function is 
#'    to dichotomize the input data.frame. Set \code{dichotomous = TRUE} 
#'    if the input data.frame is already dichotomous (default = FALSE).
#' @param dct Numeric value. Moderate risk signature predicts a baseline 
#'    risk of malignancy equal to 0.16 (CI95%: 0.06-0.25). According to 
#'    Fragomeni et al. (2022), if at least 2 other ultrasound features 
#'    than those used in the decision tree (referred to as diagnostic 
#'    covariates) are above their optimal threshold, the malignancy risk 
#'    increases to 0.55 (CI95%: 0.46-0.64). This signature will be marked 
#'    as MMR1, in contrast to the basal MMR signature (MMR0).
#'    The argument \code{dct} controls the number of diagnostic covariates 
#'    needed to switch from MMR0 to MMR1 (default = 2).
#' @param short Numeric value corresponding to the short axis cutoff 
#'    in millimeters (default = 8).
#' @param cortical Numeric value corresponding to the cortical thickness 
#'    cutoff in millimeters (default = 2).
#' @param ist Dichotomous value {0, 1} for the presence of inflammatory 
#'    stroma (perinodal hyperechogenic ring; default = 1).
#' @param ecs Dichotomous value {0, 1} for the presence of extracapsular 
#'    spread (cortical interruption; default = 1).
#' @param hab Dichotomous value {0, 1} for the absence of the hilum 
#'    (nodal core sign; default = 0).
#' @param eco Dichotomous value {0, 1} for heterogeneous echogenicity 
#'    (echostructure; default = 1).
#' @param vp Categorical value (integers from 0 to 4) associated to a 
#'    high-risk vascular flow architecture pattern (default = c(1, 2, 3)).
#' @param vfl Categorical value (integers from 0 to 4) associated to a 
#'    high-risk vascular flow localization (default = c(2, 3, 4)).
#' @param ct Categorical value (integers from 0 to 4) associated to a 
#'    high-risk cortical thickening (default = 2).
#' @param fid Categorical value (integers from 0 to 3) associated to a 
#'    high-risk focal intranodal deposit (default = c(1, 2, 3)).
#' @param cmid Categorical value (integers from 0 to 4) associated to a 
#'    high-risk cortical-medullar interface distortion (default = c(2, 3, 4)).
#' @param shape Categorical value (integers from 1 to 3) associated to a 
#'    high-risk shape (default = 3).
#' @param cs Ordinal value (integers from 1 to 5) associated to a 
#'    high-risk color score (default = 3).
#' @param grouping Categorical value (integers from 1 to 3) associated 
#'    to a high-risk grouping (default = c(2, 3)).
#' @param ... Currently ignored.
#'
#' @details The core method of the Morphonode-DT model is implemented 
#'    in this function. A series of binary branching points define the 
#'    metastatic risk signature (MRS) of the subject.
#'    The first branch point is based on the evaluation of three metastatic 
#'    markers: the absence of the nodal core sign (hilum), the presence 
#'    of the perinodal hyperechogenic ring, and the presence of cortical 
#'    interruption. If at least one of these conditions are true, the 
#'    ultrasound profile has a high malignancy risk (86-100%) with 
#'    multiple metastatic lymph nodes (Fragomeni et al. 2022). This is 
#'    referred to as the metastatic (MET) signature.
#'    A cortical thickness below 2 mm defines a low metastatic risk (LMR) 
#'    signature (0.04, CI95%: 0.00-0.10). LMR subjects are mostly mot 
#'    malignant. 
#'    Based on the values of short axis, vascular flow architecture pattern, 
#'    cortical thickening, and vascular flow localization, two more 
#'    signatures are defined: (i) moderate metastatic risk (MMR; 0.16, 
#'    CI95%: 0.06-0.25), and (ii) high metastatic risk (HMR; 0.81, 
#'    CI95%: 0.52-0.90). Similarly to the MET signature, HMR subjects are 
#'    mostly malignant, but chraracterized by a single metastatic event.
#'    MRSs should be always compared to the output of the random forest 
#'    classifier (Morphonode-RFC module) and robust binomial model 
#'    (Morphonode-RBM module). 
#'    The main advantage of MRSs is the prediction of either multiple 
#'    (MET signature) or single (HMR signature) metastasis events.
#'
#' @export
#'
#' @return A list of 4 objects:
#' \enumerate{
#' \item "signature", metastatic risk signature (MRS);
#' \item "p", MRS-associated malignancy risk (evaluated as positive 
#'    predictive value, according to Fragomeni et al. 2022);
#' \item "ci95", 95% confidence intervals of p;
#' \item "y.uss", naive guess of the outcome (0: non-malignant, 
#'    1: malignant) based on the MRS (this will be less accurate than 
#'    the RFC-based prediction).
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{us.predict}} to launch all 
#'    morphonode modules at once.
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
#' @examples
#' 
#' # Extract 5 random subjects from the default simulated dataset
#' x <- sample(mpm.us[, 2:15], 5, replace = FALSE, prob = NULL)
#' print(x)
#' 
#' # Assign a metastatic risk signature to each subject in the dataset
#' mrs <- uss(x)
#' x$signature <- mrs$signature
#' print(x)
#'
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
	y.uss <- vector()
	
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
		y.uss <- c(y.uss, y)
	}
	return(list(signature = signature, p = prob, ci95 = CI95, y.uss = y.uss))
}

#' @title Similarity search
#'
#' @description Vector similarity search across a reference dataset.
#'
#' @param v An input vector of length equal to the number of columns 
#'    of the reference data.frame (see below).
#' @param x Reference data.frame. If no reference is specified, 
#'    the default simulated dataset (object \code{\link[morphonode]{mpm.us}}) 
#'    will be used.
#' @param f Similarity function. Available functions: "cosine", "jaccard", 
#'    "euclidean", "pearson", "spearman", "kendall" (default = "cosine").
#' @param k Numeric value defining the number of top-k profiles to 
#'    return after similarity ranking (default = 5).
#' @param p Continuous value from 0 to 1 representing the minimum similarity 
#'    value for an ultrasound profile to be included in the output 
#'    (default = 0.7).
#' @param features Indices of the features (columns) in \code{x} used to 
#'    compute similarity (default = 2:15).
#' @param check.data Logical value. If TRUE (default), the input data 
#'    type is checked.
#' @param ... Currently ignored.
#'
#' @importFrom lsa cosine
#' @importFrom stats cor
#' @export
#'
#' @return A list of 4 objects:
#' \enumerate{
#' \item "signature", metastatic risk signature (MRS);
#' \item "p", MRS-associated malignancy risk (evaluated as positive 
#'    predictive value, according to Fragomeni et al. 2022);
#' \item "ci95", 95% confidence intervals of p;
#' \item "y.uss", naive guess of the outcome (0: non-malignant, 
#'    1: malignant), based on the MRS (this will be less accurate than 
#'    the RFC-based prediction).
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{us.predict}} to launch all 
#'    morphonode modules at once.
#'    See also \code{\link[morphonode]{ranksim}} for ultrasound profile 
#'    similarity ranking.
#'
#' @references
#' 
#' Leydesdorff L (2005). Similarity Measures, Author Cocitation 
#' Analysis,and Information Theory. In: JASIST 56(7), pp.769-772.
#' <https://doi.org/10.48550/arXiv.0911.4292>
#'
#' @examples
#' 
#' # Prepare a simulated malignant ultrasound profile
#' x <- newProfile(us.simulate(y = 1))
#' u <- set.rfcdata(x, ref = mpm.us[, 2:15], levels = mpm.levels)
#' 
#' # Top-similar profiles
#' sim <- topsim(u$ultrasound)
#' print(sim)
#'
topsim <- function(v, x = mpm.us, f = "cosine", k = 5, p = 0.7,
                   features = 2:15, check.data = TRUE) {
	if (check.data) v <- f2n(v)
	x$R <- apply(x[, features], 1, function(w) similarity(v, w, f = f))
	x <- x[order(x$R, decreasing = TRUE),]
	if (f != "euclidean" & p > 0 & p <= 1) x <- x[x$R >= p,]
	if (k > 0) x <- x[1:k,]
	return(x)
}

#' @title Rank ultrasound profiles by similarity
#'
#' @description Filter, rank, and return the k top-similar ultrasound 
#'    profiles respect to the input one, by searchin across a reference 
#'    dataset. This function implements the similarity profiling module 
#'    (Moprhonode-SP).
#'
#' @param u An ultrasound vector generated by 
#'    \code{\link[morphonode]{set.rfcdata}}.
#' @param v An ultrasound profile generated by 
#'    \code{\link[morphonode]{set.rbmdata}} (default = NULL). 
#'    If this argument is not NULL, the core similarity function will 
#'    be switched from "cosine" to "jaccard".
#' @param x Reference ultrasound data.frame. If no reference is specified, 
#'    the default simulated dataset (object \code{\link[morphonode]{mpm.us}}) 
#'    will be used.
#' @param k Numeric value defining the number of top-k profiles to 
#'    return after similarity ranking (default = 5).
#' @param p Continuous value from 0 to 1 representing the minimum similarity 
#'    value for an ultrasound profile to be included in the output 
#'    (default = 0.7).
#' @param j Indices of the features (columns) in \code{x} used to compute 
#'    profile similarity (default = 2:15).
#' @param d Indices of the features (columns) in \code{x} used to compute 
#'    euclidean distance (default = c(2:6, 9, 10, 11)).
#' @param signature One among "LMR", "MMR", "HMR", "MET" (default = NULL). 
#'    Resctrict the similarity search to a given metastatic risk signature 
#'    (if NULL, no signature restriction is applied).
#' @param check.data Logical value. If TRUE (default), the input data 
#'    type is checked.
#' @param orderbyDistance Logical value. If TRUE, the k top-similar 
#'    profiles are finally ordered by increasing euclidean distance 
#'    (default = FALSE).
#' @param ... Currently ignored.
#'
#' @details The input ultrasound profile is compared to each entry in 
#'    the reference dataset by computing pairwise similarity. By default, 
#'    cosine similarity is used, while jaccard similarity is enabled if 
#'    a binary vector is given. The hits are then filtered by minimum 
#'    similarity (by default, > 0.7) and pairwise euclidean distance 
#'    between them is computed. Results are then ranked by either 
#'    decreasing similarity (default) or increasing distance.
#'
#' @importFrom lsa cosine
#' @importFrom stats cor
#' @export
#'
#' @return A list of 4 objects:
#' \enumerate{
#' \item "signature", metastatic risk signature (MRS);
#' \item "p", MRS-associated malignancy risk (evaluated as positive 
#'    predictive value, according to Fragomeni et al. 2022);
#' \item "ci95", 95% confidence intervals of p;
#' \item "y.uss", naive guess of the outcome (0: non-malignant, 
#'    1: malignant), based on the MRS (this will be less accurate than 
#'    the RFC-based prediction).
#' }
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
#' Fragomeni SM, Moro F, Palluzzi F, Mascilini F, Rufini V, Collarino A, 
#' Inzani F, Giacò L, Scambia G, Testa AC, Garganese G (2022). 
#' Evaluating the risk of inguinal lymph node metastases before surgery 
#' using the Morphonode Predictive Model: a prospective diagnostic study. 
#' Ultrasound xx Xxxxxxxxxx xxx Xxxxxxxxxx. 00(0):000-000.
#' <https://doi.org/00.0000/00000000000000000000>
#' 
#' Leydesdorff L (2005). Similarity Measures, Author Cocitation 
#' Analysis,and Information Theory. In: JASIST 56(7), pp.769-772.
#' <https://doi.org/10.48550/arXiv.0911.4292>
#'
#' @examples
#' 
#' # Prepare a simulated malignant ultrasound profile
#' x <- newProfile(us.simulate(y = 1))
#' u <- set.rfcdata(x, ref = mpm.us[, 2:15], levels = mpm.levels)
#' v <- set.rbmdata(x, ref = mpm.us[, 2:15], levels = mpm.levels)
#' print(u)
#' print(v)
#' 
#' # Rank using cosine similarity and the default simulated reference
#' Rc <- ranksim(u$ultrasound)
#'
#' # Rank using jaccard similarity and the default simulated reference
#' Rj <- ranksim(u$ultrasound, v$ultrasound)
#'
#' # Compare the two rankings
#' print(Rc)
#' print(Rj)
#'
ranksim <- function(u, v = NULL, x = mpm.us, k = 5, p = 0.7,
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

#' @title Random Forest Classifier-based prediction
#'
#' @description Core function of the random forest classifier for 
#'    malignancy prediction (Morphonode-RFC module). 
#'    Given an ultrasound vector generated by 
#'    \code{\link[morphonode]{set.rfcdata}}, this function yields a 
#'    prediction of malignancy (y = 1) or non-malignancy (y = 0).
#'
#' @param u An ultrasound vector generated by 
#'    \code{\link[morphonode]{set.rfcdata}}.
#' @param rfc Random forest classifier as an object of class 
#'    \code{randomForest}. The default classifier \code{mpm.rfc$rfc} 
#'    can be used (see details).
#' @param recover Logical value. If TRUE (default) the predictors with 
#'    least out-of-bag error get the highest priority, in case of a tie.
#' @param ... Currently ignored.
#'
#' @details The default classifier (\code{rfc = mpm.rfc$rfc}) is a set 
#'    of 5 RFCs are used to predict subject's phenotype (0: non-malignant, 
#'    1: malignant). Each RFC is trained through a 5-fold nested 
#'    cross-validation procedure over 10000 random trees, with 3/14 
#'    randomly chosen variables per tree branching. Each of the 5 RFCs 
#'    achieves and independent prediction and the majority wins.
#'    The input is the default simulated dataset (object 
#'    \code{\link[morphonode]{mpm.us}}), of 948 subjects 
#'    (508 non-malignant and 440 malignant profiles) and 14 ultrasound 
#'    features. The dataset includes also the expected phenotype (y), 
#'    the related metastatic risk signature (signature), and the Brier 
#'    score (E) calculated during the cross-validation procedure.
#'
#' @import randomForest
#' @export
#'
#' @return A list of 3 objects:
#' \enumerate{
#' \item "y.hat", final prediction;
#' \item "decisions", prediction of each single RFC;
#' \item "oob.error", out-of-bag error of each classifier in the ensemble.
#' }
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
#' # Prepare a simulated malignant ultrasound profile
#' x <- newProfile(us.simulate(y = 1))
#' u <- set.rfcdata(x, ref = mpm.us[, 2:15], levels = mpm.levels)
#' print(u)
#' 
#' # Predict subject's phenotype
#' P <- rfc.predict(u$ultrasound, rfc = mpm.rfc$rfc)
#' print(P)
#'
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
	if (E < 0.001) {
		E <- "< 0.001"
	} else {
		E <- format(round(as.numeric(E), digits = 3), nsmall = 3)
	}
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
	
	if (E == "< 0.001") {
		msg.e <- "   (cutoff: E < 1)"
	} else if (E < 1) {
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

#' @title Morhonode Predictive Model (MPM) launcher
#'
#' @description The \code{us.predict} function launches the 4 MPM modules: 
#'    (i) malignancy prediction (Morphonode-RFC), (ii) malignancy risk 
#'    estimation (Morphonode-RBM), (iii) malignancy risk signature 
#'    detection (Morphonode-DT), and (iv) similarity profiling 
#'    (Morphonode-SP). 
#'    The MPM structure is described in Fragomeni et al. (2022). 
#'    See also the details section.
#'
#' @param x Ultrasound profile generated by the function 
#'    \code{\link[morphonode]{newProfile}}).
#' @param f Similarity profiling core function: one between "cosine" 
#'    (default) and "jaccard". The former directly compares ultrasound 
#'    profiles, while the latter uses dichotomized versions of them 
#'    (see also \code{\link[morphonode]{dichotomize}}).
#' @param levels A list of length 14, corresponding to the levels of 
#'    each ultrasound variable. Needed for categorical variables (factors); 
#'    for continuous variables, it should assume the nominal value of 0. 
#'    If NULL (default), the internal \code{\link[morphonode]{mpm.levels}} 
#'    variable will be used.
#' @param ref Reference ultrasound features dataset as a (n, 14) 
#'    data.frame, with n being the number of subjects (rows). 
#'    If NULL (default), the internal \code{\link[morphonode]{mpm.us}} 
#'    variable will be used.
#' @param rfc Random forest classifier as an object of class 
#'    \code{randomForest}. If NULL (default), the internal 
#'    \code{\link[morphonode]{mpm.rfc}} variable will be used.
#' @param rbm Robust binomial model as a fitted model object of class 
#'    \code{glm}. If NULL (default), the internal 
#'    \code{\link[morphonode]{mpm.rbm}} variable will be used.
#' @param k Numeric value defining the number of top-k profiles to 
#'    return after similarity ranking (default = 5).
#' @param features Indices of the features to be used for the similarity 
#'    profiling (default = 2:15).
#' @param orderbyDistance Logical value enabling sorting by increasing 
#'    euclidean distance of the top-similar ultrasound profiles 
#'    (default = FALSE).
#' @param uncertainty Function used to compute the RFC prediction error. 
#'     It can be one among "loss" (default) and "similarity". The former 
#'     estimates the error on a new prediction based on a parametric linear 
#'     relationship between the loss function and the observed (reference 
#'     dataset) error. The latter estimates the error as the average of 
#'     the observed errors of the top-3 similar profilesfrom the reference 
#'     dataset (non-parametric).
#' @param b0 Baseline uncertainty. A vector of three values representing 
#'     the intercept parameter when \code{uncertainty = "loss"}. 
#'     This value depends on the metastatic risk signature and it is 
#'     fixed to 0 for the LMR one. The b0 values for MMR, HMR, and MET 
#'     signatures are the first, second, and third element of the b0 
#'     argument, respectively (default = c(0, 0.028, 0.013)).
#' @param b Uncertainty coefficient. A numeric value indicating the 
#'     linear coefficient for the loss-to-error conversion equation.
#' @param rho Numeric value between 0 and 1 denoting the minimum required 
#'     similarity coefficient to make a naive estimation of the true 
#'     outcome y = {0, 1}. This guess is needed as a reference value for 
#'     uncertainty calculation (default = 0.9).
#' @param wmax Nominal value for the self-correlation coefficient (for 
#'     visualization purposes only; default = 1).
#' @param verbose If TRUE (default), a user-frienly summary of the 
#'     prediction and estimation results is printed to screen.
#' @param ... Currently ignored.
#'
#' @details The MPM classifier (Morphonode-RFC module) is based on an 
#'    ensemble of 5 RFCs. Each RFC is trained over 10000 random trees, 
#'    with 3/14 randomly chosen variables per tree branching. The 5 RFCs 
#'    yield independent predictions and the majority wins. This module 
#'    provides a dichotomous phenotype classification in malignant (y = 1) 
#'    and non-malignant (y = 0), and an estimation of the prediction 
#'    error (E). Similarly, the Morphonode-RBM module provides a 
#'    continuous estimate of malignancy (i.e., p = malignancy risk), 
#'    through a binomial model with robust bootstrap standard error 
#'    estimation (5000 bootstrap iterations). Optimal cutpoint estimations 
#'    define two thresholds for p (three risk intervals): low risk (p < 0.23), 
#'    moderate risk (0.23 <= p <= 0.29), and high risk (p > 0.29). 
#'    In addition, Morphonode-DT model defines four metastatic risk 
#'    signatures, strongly associated with the metastasis rate in the 
#'    corresponding subjects. LMR (low metastatic risk) and MMR (moderate 
#'    metastatic risk) signatures are associated with none-to-low metastasis 
#'    rates. Conversely, HMR (high metastatic risk) and MET (metastatic) 
#'    signatures are associated with a high risk of single and multiple 
#'    metastatic events (lymph nodes), respectively. 
#'    Finally, the Morphonode-SP module ranks ultrasound profiles from 
#'    the reference dataset (by default, the internal simulated dataset) 
#'    by similarity with respect to the input profile. This provides a 
#'    supplementary support to the classification process, having only a 
#'    secondary role compared to the other three modules. Generally, 
#'    the majority of similar profiles should have the same outcome (y) 
#'    as the input one.
#'
#' @import randomForest
#' @importFrom stats cor
#' @importFrom lsa cosine
#' @importFrom imputeR impute
#' @export
#'
#' @return A list of 5 objects:
#' \enumerate{
#' \item "prediction" (Morphonode-RFC module), a list including:
#'    \itemize{
#'        \item \code{y.hat}: the final malignancy prediction,
#'        \item \code{decisions}: the predictions of each RFC in the ensemble,
#'        \item \code{oob.err}: out-of-bag errors of each RFC in the ensemble;
#'     }
#' \item "E", estimated overall prediction error (Morphonode-RFC module);
#' \item "p", estimated malignancy risk(Morphonode-RBM module);
#' \item "signature", metastatic risk signature (Morphonode-DT module);
#' \item "profiles", data.frame containing the top-k similar profiles 
#'    sorted by similarity (Morphonode-SP module).
#'    This data.frame includes:
#'    \itemize{
#'        \item ID: numeric value identifying a subject,
#'        \item the 14 ultrasound features characterizing each subject,
#'        \item y: the observed outcome,
#'        \item E: subject-level estimated prediction error (Brier score),
#'        \item R: similarity coefficient with the input profile,
#'        \item D: euclidean distance from the input profile.
#'     }
#' }
#'
#' @author Fernando Palluzzi \email{fernando.palluzzi@gmail.com}
#'
#' @seealso See \code{\link[morphonode]{newProfile}} to create a new 
#'    ultrasound profile. 
#'    See also \code{\link[morphonode]{us.simulate}} for ultrasound 
#'    data simulation.
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
#' # Create a simulated malignant ultrasound profile
#' x <- newProfile(us.simulate(y = 1))
#' 
#' # Lauch the Morhonode Predictive Model
#' u <- us.predict(x)
#'
us.predict <- function(x, f = "cosine", levels = NULL, ref = NULL,
                       rfc = NULL, rbm = NULL, k = 5, features = 2:15,
                       orderbyDistance = FALSE, uncertainty = "loss",
                       b0 = c(0, 0.028, 0.013), b = 0.055, rho = 0.9,
                       wmax = 1, verbose = TRUE) {
	
	if (is.null(levels)) levels <- mpm.levels
	if (is.null(ref)) ref <- mpm.us[, 2:15]
	if (is.null(rfc)) rfc <- mpm.rfc$rfc
	if (is.null(rbm)) rbm <- mpm.rbm$fit
	
	suppressMessages(u <- set.rfcdata(x, levels = levels, ref = ref))
	suppressMessages(v <- set.rbmdata(x, levels = levels, ref = ref,
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
		u$ultrasound$E <- b0[i] + b*loss(p, y, method = "log")$loss
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
		#z[3] <- ifelse(z[3] == 1, 0, 1)
		printout(z, R, y$y.hat, u$ultrasound$E, p, mrs, k = k,
		         wmax = wmax, sim = f, features = features)
	}
	return(list(prediction = y, E = u$ultrasound$E, p = p,
	            signature = mrs, profiles = R))
}
