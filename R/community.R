#' Generates a species pool with all possible species in a metacommunity
#'
#' In general species in a species pool all have the same niche shapes (c_type and e_type). Species with
#' different shapes in possible in principle, however; in this case a species pool can be created with
#' a set of default species and then individual species added or replaced using [create_species()].
#'
#' @param n_species The number of species
#' @param nx The number of niche axes
#' @param c_type The type of [colonisation function][ce_funs]
#' @param e_type The type of [extinction function][ce_funs]
#' @param scale A vector, length 2, with the scales (maximum values) for the colonisation and extinction functions
#' @param xmin Vector, length=nx, the minimum values for each resource for scaling the fitness functions
#' @param xmax Vector, length=nx, the maximum values for each resource for scaling the fitness functions
#'
#' @return An S3 obect of class `speciespool`, which is a list of [species][create_species()] objects.
create_species_pool = function(n_species = 2, nx = 1, c_type = 'linear', e_type = 'constant', scale = c(1, 0.2),
			xmin = rep(0, nx), xmax=rep(1, nx)) {
	if((c_type == 'linear' || e_type == 'linear') && (n_species > 2 || nx > 1))
		stop("Linear functions are only supported for <=2 species and 1 variable")

	c_pars = switch(c_type,
		'constant' = lapply(1:n_species, function(x) list(scale = scale[1])),
		'linear' = {
			x = c(xmin[1],  xmax[1])
			y1 = c(0, scale[1])
			y2 = c(scale[1], 0)
			b = c((y1[2] - y1[1])/(x[2] - x[1]), (y2[2] - y2[1])/(x[2] - x[1]))
			a = c(y1[1] - b[1] * x[1], y2[1] - b[2] * x[1])
			mapply(function(a,b) list(a=a, b=b), a,b)},
		'gaussian' = {
			stop("Gaussian not implemented yet")
			scale = rep(scale[1], n_species)
			## mean = this is difficult to do generally for multiple variables; simple for one
			## maybe to start we only consider a single, or we at least consider variables independently
		},
	)
}

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
