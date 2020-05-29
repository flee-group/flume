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

#' @rdname ce_funs
ce_constant = function(parm) {
	if(!'scale' %in% names(parm))
		stop("constant functions require named parameter 'scale'")
	function(x) return(parm$scale)
}

#' @rdname ce_funs
ce_gaussian = function(parm) {
	if(!'scale' %in% names(parm) || !'mean' %in% names(parm) || !'sd' %in% names(parm))
		stop("gaussian functions require named parameters 'scale', 'mean', and 'sd'")
	function(x) parm$scale * exp(-((x - parm$mean)^2)/(2 * parm$sd^2))
}










