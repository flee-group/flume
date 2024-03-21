#' Colonisation and extinction functions
#'
#' @name ce_funs
#' @details Colonisation/extinction functions take a number of forms. Each form has its own
#' set of parameters that are provided when the function is created (e.g., when creating a species).
#' These functions represent the niche of a species, so the parameters are fixed.
#' 
#' Note that these functions expect an input matrix in terms of niche axes not resources. If
#' (e.g. for ratio niches) transformations are used, they must be applied before calling these
#' functions. For a higher-level function that accepts resource state as input, see [f_niche()].
#' 
#' For `ce_gaussian`, the `scale` parameter can be negative. This is used, e.g., when one wants a
#' gaussian extinction niche. This will cause the niche the be 'flipped', with the usual peak
#' of the gaussian curve at the lowest point on the curve. The scale in this case determines the 
#' height of this curve above zero and can be interpreted as a background extinction rate when
#' conditions are ideal.
#' 
#' Displacement is an additive transformation applied to the curve, raising or lowering the entire
#' curve by a constant. For curves with a negative scale, it is 2x the scale, so that the curve minimum 
#' is equal to the scale and the flat tails have a height of 2x scale.
#' @param x Niche axis matrix; input for c/e functions; one row per site, one column per axis
#' @param location Location optimum for functions with an optimum
#' @param breadth Breadth of the function (e.g., standard deviation or vcv matrix for gaussian)
#' @param scale Height of the col/ext function, see 'details'.
#' @param displacement Zero-adjustment of col/ext function, see 'details.'
#' @param nr Number of resource dimensions
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
ce_constant = function(scale, nr) {
	if(any(scale < 0)) {
		scale = abs(scale)
		warning("Scale cannot be negative for constant niches, using abs(scale)")
	}
	function(x) {
		if(is.data.frame(x))
			x = as.matrix(x)
		if(!is.matrix(x))
		   x = matrix(x, ncol = nr)
		return(rep(scale, nrow(x)))
	}
}

#' @rdname ce_funs
ce_gaussian = function(location, breadth, scale, displacement) {
	if(length(scale) != 1)
		stop("Scale must be a single value")
	
	if(missing(displacement))
		displacement = ifelse(scale < 0, 2 * abs(scale), 0)

	if(length(displacement) != 1)
		stop("displacement must be a single value")
	
	if(length(location) == 1) {
		f = function(x) displacement + scale * exp(-(x - location)^2/(2 * breadth^2))
	} else {
		if(!requireNamespace("mvtnorm", quietly = TRUE))
			stop("Multidimensional niches require the 'mvtnorm' and 'cubature' packages")
		if(!is(breadth, "matrix") || 
		   (length(location) != ncol(breadth) && length(location) != nrow(breadth))) {
			stop("Non conforming arguments: location must be a mean vector and breadth a variance-covariance matrix")
		}
		f = function(x) {
			displacement + scale * mvtnorm::dmvnorm(x, location, breadth) /
					mvtnorm::dmvnorm(location, location, breadth)
		}
	}
	return(f)
}










