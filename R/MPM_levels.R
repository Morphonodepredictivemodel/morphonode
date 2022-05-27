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

#' @title Ultrasound features levels
#'
#' @description A list of 14 elements containing the levels of each 
#'    ultrasound feature.
#' @name mpm.levels
#' @usage mpm.levels
#' @docType data
#' @format
#' "mpm.levels" is a list of 14 vectors, each containing the levels of 
#' an ultrasound feature. Continuous features (i.e., short axis and cortical 
#' thickness) have the nominal value of 0. Levels are used by predictive 
#' functions to switch from numeric to factor classes.
#'
#' @examples
#' 
#' # Generate a simulated malignant ultrasound profile
#' x <- new.profile(simulate.us(y = 1))
#' 
#' # Prepare the profile for RFC malignancy prediction
#' u <- set.rfcdata(x, ref = mpm.us[, 2:15], levels = mpm.levels)
#' print(u)
#' 
#' # Predict subject's phenotype
#' P <- rfc.predict(u$ultrasound, rfc = mpm.rfc$rfc)
#' print(P)
#'

NULL
