#' Greate a community for simulations
#' @details The `r_scale` parameter must either be a single value (in which case it is repeated for all species/resources),
#' a vector of length equal to `nx` (allowing all species to have identical effects, but that differ by resource),
#' or a matrix with nrow == n_species and ncol == nx
#' @param n_species The number of species
#' @param nx The number of niche axes
#' @param xmin Vector, length=nx, the minimum values for each resource for scaling the fitness functions
#' @param xmax Vector, length=nx, the maximum values for each resource for scaling the fitness functions
#' @param r_scale Scales for [resource use functions][ruf()], see 'details'
#' @param alpha Active dispersal ability for each species
#' @param beta Passive dispersal ability for each species
#' @param ... Additional parameters to pass to [species_pool()]
#' @return An object of class `community`, with the following items:
#' * `species` a [species pool][species_pool()]
#' * `competition` a [competition()] matrix
#'
#' Additionally, the following attributes:
#' * `niche_max` the maximum value of the niche along the x interval given
#' * `xmin` the minimum value of the environmental gradient
#' * `xmax` the maximum value of the environmental gradient
#' @examples
#' comm = community()
#' plot(comm)
#' @export
metacommunity = function(n_species = 2, nx = 1, xmin = rep(0, nx), xmax=rep(1, nx), r_scale = 0.05, alpha = 0,
					beta = 0.5, ...) {
	comm = structure(list(), class = "metacommunity")

	stopifnot(length(r_scale) == 1 || length(r_scale) == nx || dim(r_scale) == c(n_species, nx))
	if(!is(r_scale, "matrix")) {
		r_scale = matrix(r_scale, nrow=n_species, ncol=nx, byrow = TRUE)
	}

	comm$species = species_pool(n_species = n_species, nx = nx, xmin=xmin, xmax=xmax, alpha = alpha, beta = beta, ...)
	comm$competition = .compute_comp_matrix(comm$species, xmin, xmax)
	comm$r_scale = r_scale

	# default boundary condition is to return zero for all sites
	# n is the number of sites
	comm$boundary = function(n=1) matrix(0, nrow=n, ncol=n_species)

	attr(comm, "niche_max") =
		sapply(comm$species, function(sp) {
			optimise(function(r, sp) {
				r = matrix(r, nrow=1, ncol=1)
				f_niche(sp, r)[1,1]
			}, c(xmin, xmax), sp = sp, maximum = TRUE)$objective
		})
	attr(comm, "xmin") = xmin
	attr(comm, "xmax") = xmax

	return(comm)
}

#' Generate a random community
#' @details This function will ensure that all sites have at least one species. If many species have low prevalence, it
#' may be that the initial random draw produces many zero-richness sites, in which case one species will be selected at
#' random for each of these sites. This can have the side effect that the prevalences in the end do not match the desired
#' prevalence very closely.
#' @param rn A [river_network()]
#' @param mc A [metacommunity()]
#' @param prevalence A vector of length 1 or `length(mc)`, what proportion of sites should be occupied on average by each species.
#' @return A site by species matrix, with sites taken from `rc` and species from `mc`
#' @examples
#' Q = rep(1, 4)
#' adj = matrix(0, nrow = 4, ncol = 4)
#' adj[1,2] = adj[2,3] = adj[4,3] = 1
#' rn = river_network(adj, Q)
#' comm = community()
#' site_by_species(rn) = random_community(rn, mc)
#' plot(rn, variable = "site_by_species")
#' @export
random_community = function(rn, mc, prevalence = 0.25) {
	i = nrow(adjacency(rn))
	j = length(mc$species)
	com = do.call(cbind, mapply(function(jj, pr)
		sample(0:1, i, replace = TRUE, prob = c(1-pr, pr)), 1:j, prevalence, SIMPLIFY = FALSE))

	if(all(prevalence < 1/i))
		warning("All species have low prevalence; random communities might not match desired prevalences")

	r0 = which(rowSums(com) == 0)
	if(length(r0) > 0) {
		for(k in r0)
			com[k, sample(1:length(j))] = 1
	}
	return(com)
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
#' @param alpha Active dispersal ability for each species
#' @param beta Passive dispersal ability for each species
#'
#' @return A list of [species()]
#' @examples
#' ## create a default 2-species pool
#' spp = flume:::species_pool()
species_pool = function(n_species = 2, nx = 1, c_type = 'linear', e_type = 'constant', scale = c(1, 0.2),
						xmin = rep(0, nx), xmax=rep(1, nx), alpha, beta) {
	if((c_type == 'linear' || e_type == 'linear') && (n_species > 2 || nx > 1))
		stop("Linear functions are only supported for <=2 species and 1 variable")

	c_pars = .generate_ce_params(c_type, n_species, scale[1], xmin, xmax)
	e_pars = .generate_ce_params(e_type, n_species, scale[2], xmin, xmax)

	mapply(species, c_type = c_type, e_type = e_type, c_par = c_pars, e_par = e_pars, alpha = alpha, beta = beta,
		   SIMPLIFY = FALSE, USE.NAMES = FALSE)
}


#' Returns dispersal parameters for a metacommunity
#' @param x A [metacommunity()]
#' @return Named list with two elements, `alpha` and `beta`, each one a vector of dispersal parameters for each species
#' @export
#' @examples
#' comm = community()
#' dispersal_params(comm)
dispersal_params = function(x) {
	alpha = sapply(x$species, function(y) y$alpha)
	beta = sapply(x$species, function(y) y$beta)
	list(alpha=alpha, beta=beta)
}

#' Creates species
#'
#' @details Creates a species, which is defined by its niche parameters and its colonisation/extinction functions
#'
#' Note that for gaussian functions in multiple variables, `scale` is always a scalar, but
#' `mean` must be a location vector of length `nx`, and `sd` must be a variance covariance matrix
#' with `dim=c(nx,nx)`.
#'
#' The resource use function is set automatically by default, but can be set to an arbitrary function; see
#' [resource use functions][ruf()] for details on how this function should behave.
#'
#' @param c_type Form of [colonisation function][ce_funs]
#' @param e_type Form of [extinction function][ce_funs]
#' @param c_par Named parameter list; constant parameters to include in the c/e functions, see [colonisation functions][ce_funs]
#' @param e_par Named parameter list; constant parameters to include in the c/e functions, see [extinction functions][ce_funs]
#' @param alpha Active dispersal ability
#' @param beta Passive dispersal ability
#' @return An S3 object of class 'species', which contains the following named elements:
#'       `col`: The colonisation function, takes a state matrix R and returns a vector of colonisation rates
#'       `ext`: The extinction function
#'       `c_par`: colonisation niche parameters
#'       `e_par`: extinction niche parameters
#'       `alpha`: Active dispersal ability
#'       `beta`: Passive dispersal ability
#' @examples
#' sp = flume:::species('linear', 'constant', list(a=0, b=1), e_par = list(scale = 0.2))
#' plot(sp)
#' @export
species = function(c_type, e_type, c_par, e_par, alpha=0, beta=0) {
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
	x$alpha = alpha
	x$beta = beta

	return(x)
}

#' Make a competition matrix
#' @param sp A [species_pool()]
#' @param xmin The minimum for integration
#' @param xmax The maximum for integration
#' @keywords internal
#' @return a matrix giving pairwise competition coefficients
.compute_comp_matrix = function(sp, xmin, xmax) {
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


#' @export
f_niche = function(x, ...)
	UseMethod("f_niche", x)


#' Functions for determining the niche of species
#' @name niche
#' @param x A [species()] or [metacommunity()]
#' @param R A site by resource matrix
#' @examples
#' comm = metacommunity()
#' st = matrix(seq(0, 1, length.out = 4), ncol = 1, dimnames = list(NULL, 'R'))
#'
#' # for a single species
#' f_niche(comm$species[[1]], st)
#'
#' # for multiple species
#' f_niche(comm, st)
#' @return A vector (for a single species) or matrix (for a metacommunity) giving the value(s)
#' of the niche at each site
#' @export
f_niche.species = function(x, R) {
	x$col(R) - x$ext(R)
}

#' @rdname niche
#' @export
f_niche.metacommunity = function(x, R) {
	do.call(cbind, lapply(x$species, f_niche, R=R))
}
