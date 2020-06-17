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

	c_pars = .generate_ce_params(c_type, n_species, scale[1], xmin, xmax)
	e_pars = .generate_ce_params(e_type, n_species, scale[2], xmin, xmax)

	comm = mapply(create_species, c_type = c_type, e_type = e_type, c_par = c_pars, e_par = e_pars, SIMPLIFY = FALSE)
	class(comm) = c("speciespool", class(comm))
	return(comm)
}


#' Create parameter sets for CE functions
#' @keywords internal
.generate_ce_params = function(type, n, scale, xmin, xmax) {
	switch(type,
		'constant' = lapply(1:n, function(x) list(scale = scale)),
		'linear' = {
		   	x = c(xmin,  xmax)
		   	y1 = c(0, scale)
		   	y2 = c(scale, 0)
		   	b = c((y1[2] - y1[1])/(x[2] - x[1]), (y2[2] - y2[1])/(x[2] - x[1]))
		   	a = c(y1[1] - b[1] * x[1], y2[1] - b[2] * x[1])
		   	mapply(function(a,b) list(a=a, b=b), a,b, SIMPLIFY=FALSE)},
		'gaussian' = {
		   	stop("Gaussian not implemented yet")
		   	scale = rep(scale[1], n_species)
		   	## mean = this is difficult to do generally for multiple variables; simple for one
		   	## maybe to start we only consider a single, or we at least consider variables independently
	})
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


#' Plot species independent stable envelopes
#' @param x A species pool
#' @param R an optional resource state matrix
#' @param axis Which resource axis (i.e., column in R) to plot along
#' @export
plot.speciespool = function(x, R, axis = 1, ...) {
	if(missing(R)) R = matrix(seq(0, 1, length.out=50), ncol=1)
	ypl = lapply(x, function(sp) sp$col(R) - sp$ext(R))
	yl = c(0, max(unlist(ypl)))

	# colours from colorbrewer
	cols = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")
	if(length(cols) > length(x))
		cols = cols[1:length(x)]

	args = .default_plot_pool_options()
	args$x = 0
	args$y = 0
	args$xlim = range(R[,axis])
	args$ylim = yl
	par(mar = c(5,4,4,6))
	do.call(plot, args)

	.make_line = function(y, col, args) {
		args$x = R[,axis]
		args$y = y
		args$col = col
		args$type = 'l'
		do.call(lines, args)
	}
	invisible(mapply(.make_line, ypl, cols, MoreArgs = list(args = args)))

	xpd = par()$xpd
	par(xpd = TRUE)
	legend("topright", legend = 1:length(x), col = cols, title = "species",
		   lty = args$lty, lwd = args$lwd, bty = "n", inset=c(-0.15,0), cex=0.7)
	par(xpd = xpd)

}


#' Plot species niches
#' @param x A species
#' @param R a resource state matrix
#' @param axis Which resource axis (i.e., column in R) to plot along
#'@export
plot.species = function(x, R, axis = 1, ...) {
	if(missing(R)) R = matrix(seq(0, 1, length.out=50), ncol=1)
	yc = x$col(R)
	ye = x$ext(R)
	yl = range(c(yc, ye))

	args = .default_plot_species_options(...)
	cols = args$col
	args$x = R[,axis]
	args$y = yc
	if(!"ylim" %in% names(args)) args$ylim = yl
	args$col = cols[1]
	par(mar = c(5,4,4,6))
	do.call(plot.default, args)

	args$y = ye
	args$col = cols[2]
	do.call(lines, args)

	xpd = par()$xpd
	par(xpd = TRUE)
	legend("topright", legend = c("colonisation", "extinction"), col = cols,
		   lty = args$lty, lwd = args$lwd, bty = "n", inset=c(-0.15,0), cex=0.7)
	par(xpd = xpd)
}


#' Set default plot options when not user-specified
#' @keywords internal
.default_plot_pool_options = function(...) {
	dots = 	.default_plot_species_options(...)
	dots$ylab = "Dominant eigenvalue"
	dots$type = "n"
	dots$col = NULL
	return(dots)
}

#' Set default plot options when not user-specified
#' @keywords internal
.default_plot_species_options = function(...) {
	dots = list(...)
	nms = names(dots)
	if(!"xlab" %in% nms) dots$xlab = "Resource concentration"
	if(!"ylab" %in% nms) dots$ylab = "Colonisation/extinction rate"
	if(!"type" %in% nms) dots$type = "l"
	if(!"bty" %in% nms) dots$bty = "n"
	if(!"lwd" %in% nms) dots$lwd = 1.5
	if(!"col" %in% nms) {
		dots$col = c("#1f78b4", "#e31a1c")
	} else {
		if(length(dots$col) == 1)
			dots$col = rep(dots$col, 2)
	}
	return(dots)
}

