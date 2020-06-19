#' Greate a community for simulations
#' @param n_species The number of species
#' @param nx The number of niche axes
#' @param xmin Vector, length=nx, the minimum values for each resource for scaling the fitness functions
#' @param xmax Vector, length=nx, the maximum values for each resource for scaling the fitness functions
#' @param ... Additional parameters to pass to [species_pool()]
#' @return An object of class `community`, with the following items:
#' * `species` a [species pool][species_pool()]
#' * `competition` a [competition()] matrix
#'
#' Additionally, the following attributes:
#' * `xmin` the minimum value of the environmental gradient
#' * `xmax` the maximum value of the environmental gradient
#' @examples
#' comm = community()
#' plot(comm)
#' @export
metacommunity = function(n_species = 2, nx = 1, xmin = rep(0, nx), xmax=rep(1, nx), ...) {
	comm = structure(list(), class = "metacommunity")

	comm$species = species_pool(n_species = n_species, nx = nx, xmin=xmin, xmax=xmax, ...)
	comm$competition = competition(comm$species, xmin, xmax)

	attr(comm, "xmin") = xmin
	attr(comm, "xmax") = xmax
	return(comm)
}




#' Generates a species pool with all possible species in a metacommunity
#'
#' In general species in a species pool all have the same niche shapes (c_type and e_type). Species with
#' different shapes in possible in principle, however; in this case a species pool can be created with
#' a set of default species and then individual species added or replaced using [species()].
#'
#' @param n_species The number of species
#' @param nx The number of niche axes
#' @param c_type The type of [colonisation function][ce_funs]
#' @param e_type The type of [extinction function][ce_funs]
#' @param scale A vector, length 2, with the scales (maximum values) for the colonisation and extinction functions
#' @param xmin Vector, length=nx, the minimum values for each resource for scaling the fitness functions
#' @param xmax Vector, length=nx, the maximum values for each resource for scaling the fitness functions
#'
#' @return A list of [species()]
#' @examples
#' ## create a default 2-species pool
#' spp = flume:::species_pool()
species_pool = function(n_species = 2, nx = 1, c_type = 'linear', e_type = 'constant', scale = c(1, 0.2),
						xmin = rep(0, nx), xmax=rep(1, nx)) {
	if((c_type == 'linear' || e_type == 'linear') && (n_species > 2 || nx > 1))
		stop("Linear functions are only supported for <=2 species and 1 variable")

	c_pars = .generate_ce_params(c_type, n_species, scale[1], xmin, xmax)
	e_pars = .generate_ce_params(e_type, n_species, scale[2], xmin, xmax)

	mapply(species, c_type = c_type, e_type = e_type, c_par = c_pars, e_par = e_pars, SIMPLIFY = FALSE)
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
#' sp = flume:::species('linear', 'constant', list(a=0, b=1), e_par = list(scale = 0.2))
#' plot(sp)
#' @export
species = function(c_type, e_type, c_par, e_par) {
	x = structure(list(), class = "species")
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

	return(x)
}




#' Produce a default competition matrix
#'
#' Competition is defined by computing the overlap of species niches, which are in turn defined as c(R) - m(R).
#' This function computes the niche for each species, then integrates for every pairwise pair
#' @param sp A [species_pool()]
#' @param xmin The minima for the resources
#' @param xmax The maxima for resources
#'
#' @examples
#' spp = flume:::species_pool()
#' competition(spp, 0, 1)
#'
#' @return a matrix giving pairwise competition coefficients
competition = function(sp, xmin, xmax) {
	comp = matrix(0., nrow=length(sp), ncol=length(sp))
	for(i in 1:(length(sp) - 1)) {
		si = sp[[i]]
		comp[i,i] = integrate(.pairwise_comp(si, si), xmin, xmax)$value
		for(j in (i+1):length(sp)) {
			sj = sp[[j]]
			comp[i,j] = integrate(.pairwise_comp(si, sj), xmin, xmax)$value
			comp[j,i] = comp[i,j]
		}
	}
	comp[j,j] = integrate(.pairwise_comp(sj, sj), xmin, xmax)$value
	return(comp)
}


#' Generate a function to integrate to determine competition between two species
#' @param sp1 First species
#' @param sp2 Second species
#' @keywords internal
.pairwise_comp = function(s1, s2) {
	function(x) {
		l1 = s1$col(x) - s1$ext(x)
		l2 = s2$col(x) - s2$ext(x)
		ifelse(l1 > l2, l2, l1)
	}
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




