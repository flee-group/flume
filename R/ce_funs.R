#' Colonisation and extinction functions
#'
#' @name ce_funs
#' @details Colonisation/extinction functions take a number of forms. Each form has its own
#' set of parameters that are provided when the function is created (e.g., when creating a species).
#' These functions represent the niche of a species, so the parameters are fixed.
#'
#' Regardless of their form, the functions all have the same signature.
#'
#' The parameter lists contain different named items depending on the type of function
#'
#' * 'constant'
#'      - `scale`: the height of the function
#' * 'linear'
#'      - `a`: the y-intercept
#'      - `b`: the slope
#' * 'gaussian'
#'      - `scale`: the maximum value of the curve
#'      - `mean`: the centre of the curve, x-value where y == `scale`
#'      - `sd`: the width of the curve
#'
#' @param x state variable matrix; input for c/e functions; one row per site, one column per variable
#' @param parm Parameter list, a named list containing parameters for the c/e functions; see 'details'
#'
#' @return For c/e functions, a vector of colonisation/extinction rates with length == `nrow(x)`
#' For all others, a c/e function of the desired form
#' @examples
#' fun = ce_linear(parm=list(a=0, b = 1))
#' R = seq(0,1,0.1)
#' fun(R)
ce_linear = function(parm) {
	if(!'a' %in% names(parm) || !'b' %in% names(parm))
		stop("linear functions require named parameters 'a' and 'b'")
	function(x) return(parm$a + parm$b*x)
}

ce_constant = function(parm) {
	if(!'scale' %in% names(parm))
		stop("constant functions require named parameter 'scale'")
	function(x) return(parm$scale)
}

ce_gaussian = function(parm) {
	if(!'scale' %in% names(parm) || !'mean' %in% names(parm) || !'sd' %in% names(parm))
		stop("gaussian functions require named parameters 'scale', 'mean', and 'sd'")
	function(x) parm$scale * exp(-((x - parm$mean)^2)/(2 * parm$sd^2))
}









#' Generate colonisation and extinction functions
#'
#' @param
#' @return A colonization or extinction function. These functions take two arguments:
#'     * `x`: the state variable matrix, one column per dimension given by `nx`
#'     * `param`: a parameter list; can be optional if the c/e function needs no additional parameters
#' @examples
# generate_colfun = function(type = c('linear', 'gaussian'), n_species = 2) {
# 	type = match.arg(type)
#
# 	if(type == 'linear' && (n_species > 3 || nx > 1))
# 		stop("linear col/ext functions only support maximum 2 species and one variable")
# }
#
#' Create a species pool with corresponding colonisation and extinction functions
#'
#' @param n_species The number of species
#'
#' @return An S3 obect of class `speciespool`, which is a list of [species][create_species()] objects.
# generate_species_pool = function(n_species = 2, nx = 1,
# 				c_type =  c('linear', 'gaussian', 'constant'),
# 				e_type = c('constant', 'linear', 'quadratic')) {
# 	c_type = match.arg(c_type)
#
# 	if(c_type == 'linear') {
#
# 	}
# }
