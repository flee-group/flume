
#' Creates species
#'
#' @details Creates a species, which is defined by its niche parameters and its colonisation/extinction functions
#'
#' Note that for gaussian functions in multiple variables, `scale` is always a scalar, but
#' `mean` must be a location vector of length `nx`, and `sd` must be a variance covariance matrix
#' with `dim=c(nx,nx)`.
#'
#' @param c_type Form of [colonisation function][ce_funs]
#' @param e_type Form of [extinction function][ce_funs]
#' @param c_par Named parameter list; constant parameters to include in the c/e functions, see [colonisation functions][ce_funs]
#' @param e_par Named parameter list; constant parameters to include in the c/e functions, see [extinction functions][ce_funs]
#' @return An S3 object of class 'species', which contains the following named elements:
#'       `col`: The colonisation function, takes a state matrix R and returns a vector of colonisation rates
#'       `ext`: The extinction function
#'       `c_par`: colonisation niche parameters
#'       `e_par`: extinction niche parameters
#' @examples
#' create_species('linear', 'constant', list(a=0, b=1), e_par = list(scale = 0.2))
#' @export
create_species = function(c_type, e_type, c_par, e_par) {
	x = list()
	x$col = switch(c_type,
		"constant" = ce_constant(c_par),
		"linear" = ce_linear(c_par),
		"gaussian" = ce_gaussian(c_par),
		stop("unknown colonisation function type:", c_type)
	)

	x$ext = switch(e_type,
		"constant" = ce_constant(e_par),
		"linear" = ce_linear(e_par),
		"gaussian" = stop('gaussian extinction functions not recommended'),
		stop("unknown extinction function type:", e_type)
	)

	x$c_par = c_par
	x$e_par = e_par

	class(x) = c('species', class(x))
	return(x)
}

#' Plot species niches
#' @param x A species
#' @param R a resource state matrix
#' @param axis Which resource axis (i.e., column in R) to plot along
#'@export
plot.species = function(x, R, axis = 1, ...) {
	yc = x$col(R)
	ye = x$ext(R)
	yl = range(c(yc, ye))
	plot(R[,axis], yc, xlab = "Resource concentration", ylab = "Colonisation/extinction rate", ylim=yl, ...)
}
