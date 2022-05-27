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

#' @title Morphonode simulated ultrasound feature dataset
#'
#' @description A simulated dataset of 948 ultrasound profiles 
#'    (440 malignant and 508 non-malignant) used to build the default 
#'    ensemble random forest classifier (RFC) and the robust binomial 
#'    model (RBM). This dataset was generated as a 4-fold expansion of 
#'    the original ultrasound feature dataset of 237 groin samples (75 
#'    malignant and 162 non-malignant) from Fragomeni et al. (2022), 
#'    using the morphonode function \code{\link[morphonode]{simulate.us}}.
#' @name mpm.us
#' @usage mpm.us
#' @docType data
#' @format
#' "mpm.us" is a data.frame of 948 rows (simulated ultrasound profiles) 
#' and 18 columns, including: a progressive number used as unique profile 
#' identifier (ID), 14 ultrasound features used for RFC and RBM building, 
#' expected simulation phenotype used as ground truth (y = {0: non-malignant, 
#' 1: malignant}), metastatic risk signature associated to each simulated 
#' ultrasound profile (signature), subject-level prediction error calculated 
#' as Brier score (E).
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
#' # Default simulated dataset phenotype frequencies
#' print(table(mpm.us$y))
#' 
#' # Default simulated dataset signatures frequencies
#' print(table(mpm.us$signature))
#' 
#' # Default simulated dataset prediction error quartiles
#' print(quantile(mpm.us$E))
#'

NULL
