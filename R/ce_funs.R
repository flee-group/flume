#' Colonisation and extinction functions
#'
#' @name ce_funs
#' @details Colonisation/extinction functions take a number of forms. Each form has its own
#' set of parameters that are provided when the function is created (e.g., when creating a species).
#' These functions represent the niche of a species, so the parameters are fixed.
#' @param x state variable matrix; input for c/e functions; one row per site, one column per variable.
#' @param location Location optimum for functions with an optimum
#' @param breadth Breadth of the function (e.g., standard deviation or vcv matrix for gaussian)
#' @param scale Height of the col/ext function
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
ce_constant = function(scale) {
	function(x) {
		if(!is.matrix(x))
		   x = matrix(x, ncol = 1)
		return(rep(scale, nrow(x)))
	}
}

#' @rdname ce_funs
ce_gaussian = function(location, breadth, scale) {

	if(length(location) == 1) {
		return(function(x) scale * exp(-((x - location)^2)/(2 * breadth^2)))
	} else {
		if(!requireNamespace("mvtnorm", quietly = TRUE))
			stop("Multidimensional niches require the 'mvtnorm' and 'cubature' packages")
		return(function(x) scale * mvtnorm::dmvnorm(x, location, breadth) /
			   	mvtnorm::dmvnorm(location, location, breadth))
	}

}










