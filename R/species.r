#' Creates species
#'
#' @details Creates a species, which is defined by its niche parameters and its
#' colonisation/extinction functions.
#'
#' Note that for gaussian functions in multiple variables, `scale` is always a scalar, but
#' `location` must be a location vector of length `nx`, and `width` must be a variance covariance
#' matrix with `dim = c(nx, nx)`.
#'
#' The resource use function is set automatically by default, but can be set to an arbitrary
#' function; see [resource use functions][ruf()] for details on how this function should behave.
#'
#' @param col Named list of colonisation parameters, must include `fun`, 
#'		the colonisation function (e.g., [ce_gaussian()])
#' @param ext Named list of colonisation parameters, must include `fun`, 
#'		the extinction function (e.g., [ce_constant()])
#' @param alpha Active dispersal ability
#' @param beta Passive dispersal ability
#' @param r_use Resource use scaling parameter, one per resource axis
#' @param r_trans A transformation function for converting resources into niche dimensions
#' @return An S3 object of class 'species', which contains the following named elements:
#'   * `col`: The colonisation function, takes a state matrix R and returns a vector of
#'		colonisation rates
#'   * `ext`: The extinction function
#'   * `par_c`: colonisation niche parameters
#'   * `par_e`: extinction niche parameters
#'   * `alpha`: Active dispersal ability
#'   * `beta`: Passive dispersal ability
#'   * `r_use`: Resource use rate per niche axis
#'	 * `r_trans`: Resource-to-niche transformation function
#'
#' Additionally, the following attributes:
#'    * `niche_max`: the maximum possible value of the fundamental niche
#' @examples NULL
#' @export
species = function(col, ext, alpha, beta, r_use, r_trans = identity) {
	x = structure(list(), class = "species")
	x$par_c = col[names(col) != 'fun']
	x$par_e = ext[names(ext) != 'fun']
	x$alpha = alpha
	x$beta = beta
	if(length(r_use) == 1)
		r_use = rep(r_use, length(x$par_c$location))
	x$r_use = r_use
	.check_species_params(x)
	x$r_trans = r_trans
	x$col = do.call(col[['fun']], x$par_c)
	x$ext = do.call(ext[['fun']], x$par_e)
	
	attr(x, "niche_max") = f_niche(x, N = x$par_c$location)
	return(x)
}
