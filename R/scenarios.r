#' @rdname niches
#' @name niches
#' @title Generate scenarios for species niches
#' @param nsp Number of species
#' @param nr Number of resources
#' @param location Location parameter for the niche, see 'details'.
#' @param breadth Niche breadth, see 'details'.
#' @param scale_c Niche height for colonisation, see 'details'.
#' @param scale_e Niche height for extinction, see 'details'.
#' @param r_use Resource use scaling parameter, see 'details
#' @param r_lim Minimum and maximum possible values for each niche dimension; can be a vector
#'		in the case of a single resource, otherwise a matrix with 2 columns and nr rows.
#' @param static Vector, the indices of resources that are static, `NA` if none. Any resources
#'		that are not static are by definition dynamic
#' @param ratio A 2-column matrix; each row is the indices of a pair of resources. The ratio of
#'		the first to the second will define a single niche dimension. `NA` if none.
#' @param ... Additional named niche parameters to override the defaults in `niches_custom()`.
#' @details
#' Scenarios handle the complexity of generating niches for you. For custom niches, there is some
#' in allowable parameter values and dimensions. Details follow.
#'
#' `location`: Required for custom niches; the location of the niche optimum for each
#' species-resource combination.
#' Must be a matrix with one row per species (`nsp`), one column per niche axis (`nr`). For a
#' single niche axis, a vector with one entry per species is also accepted. The location is the
#' mean of the Gaussian colonisation function.
#'
#' `breadth`: Niche breadth, or the standard deviation of the Gaussian colonisation function;
#' larger values indicate species can occur in a wider variety of environments. For a single
#' niche dimension, either a scalar (all species have the same breadth) or a vector of length `nsp`.
#' For multivariate niches, the following are possible:
#'
#'   * A **scalar**: all species have the same breadth in each niche dimension
#'   * A **vector** of length `nsp`: each species has a different breadth, but breadth is the same
#'   	for each niche axis.
#'   * A **vector** of length `nr`: niche breadth varies by axis, but all species have the same
#'   	breadth for each niche axis. If `nr == nsp`, then this option is not possible, instead the
#'   	former applies.
#'   * A **matrix** of `nsp` rows and `nr` columns: each species-resource combo has a unique
#'   	breadth, but all resources are orthogonal. In other words, performance for a given resource
#'   	cannot depend on the concentration of any other resource.
#'   * A **list** of length `nsp`, each element is a square symmetric `nr` by `nr` matrix. Similar
#'   	to the matrix above, but a full variance-covariance matrix is supplied for each species,
#'   	describing the breadth of the niche along each axis but also how the axes covary.
#'
#' `scale_c`, `scale_e`: The (relative) height of the Gaussian colonisation function or
#' the constant extinction function, must be a positive real number.
#' Can be supplied as a scalar (all species have the same scale) or a vector of length `nsp`.
#' If missing a default value of 0.5 (for colonisation) or 0.2 (for extinction) will be used.
#'
#' `r_use`: The scale of resource consumption; larger values indicate
#' faster consumption of resources, see [ruf()] for more details. This can be a scalar (all species
#' consume all resources at the same rate relative to niche position), a vector of length `nsp`
#' (species behave differently, but all resources within species are consumed identically),
#' a vector of length `nr` (all species have identical behaviour, but each resource is consumed at
#' a different rate), or a matrix with `nsp` rows and `nr` columns. Note that static resources
#' ignore this parameter.
#' @return A list of niche parameters, suitable for passing on to [species()].
#' @examples
#' niches_uniform(nsp = 4)
#' niches_custom(nsp = 2, nr = 2, location = matrix(c(1,2,3,4), nrow=2))
#' @export
niches_custom = function(nsp, nr, location, breadth = 1, scale_c = 0.5, scale_e = 0.2,
		r_use = 0.05, r_lim = matrix(rep(c(0, 1), each = nr), ncol = 2), static = NA, ratio = NA) {

	# TODO
	# handle static and ratio

	if(!is.matrix(location)) {
		if(nr != 1)
			stop("If nr > 1, location must be a matrix")
		location = matrix(location, ncol = 1)
	}

	location = lapply(1:nsp, function(i) location[i, ])
	breadth = .check_breadth(breadth, nsp, nr)
	scale_c = .check_scale(scale_c, nsp)
	scale_e = .check_scale(scale_e, nsp)
	r_use = .check_r_use(r_use, nsp, nr)

	if(!is.matrix(r_lim))
		r_lim = matrix(r_lim, ncol = 2)

	if(nrow(r_lim) != nr)
		stop("r_lim must be a matrix with 1 row per resource and 2 columns")

	list(location = location, breadth = breadth, scale_c = scale_c, scale_e = scale_e,
		r_use = r_use, r_lim = r_lim)
}

#' @rdname niches
#' @name niches
#' @export
niches_uniform = function(nsp = 2, nr = 1, r_lim = c(0, 1), ...) {
	if(nr > 1)
		stop("niches_uniform does not yet support multiple resources, use niches_random instead.")

	## for small numbers of species, we create some buffer around the edges
	if(nsp < 6) {
		location = seq(r_lim[1], r_lim[2], length.out = nsp + 2)
		location = location[2:(nsp + 1)]
	} else {
		location = seq(r_lim[1], r_lim[2], length.out = nsp)
	}
	pars = list(...)
	if(! "breadth" %in% names(pars) & nsp > 1)
		pars$breadth = location[2] - location[1]
	pars$nsp = nsp
	pars$nr = nr
	pars$location = location
	pars$r_lim = r_lim
	do.call(niches_custom, pars)
}

#' @rdname niches
#' @name niches
#' @export
niches_random = function(nsp = 2, nr = 1, r_lim = c(0, 1), ...) {
	if(!is.matrix(r_lim) && length(r_lim == 2))
		r_lim = cbind(rep(r_lim[1], nr), rep(r_lim[2], nr))

	location = do.call(rbind, lapply(1:nsp, function(i) runif(nr, r_lim[, 1], r_lim[, 2])))
	pars = list(...)
	pars$location = location
	pars$r_lim = r_lim
	pars$nsp = nsp
	pars$nr = nr
	if(! "breadth" %in% names(pars)) {
		sd_avg = (r_lim[, 2] - r_lim[, 1]) / nsp
		sd_sd = sd_avg / 4
		if(nr == 1) {
			pars$breadth = rnorm(nsp, sd_avg, sd_sd)
		} else {
			pars$breadth = lapply(1:nsp, function(i) {
				y = matrix(0, nrow = nr, ncol = nr)
				diag(y) = rnorm(nr, sd_avg, sd_sd)
				y
			})
		}
	}
	do.call(niches_custom, pars)
}

#' @rdname dispersal
#' @name dispersal
#' @title Generate scenarios for species' dispersal abilities
#' @param nsp Number of species
#' @param alpha A scalar or a vector, active dispersal ability for each species
#' @param beta A scalar or a vector, passive dispersal ability for each species
#' @param ... Additional named dispersal parameters to override the defaults in
#'		`dispersal_custom()`.
#' @return A list of dispersal parameters
#' @examples
#' dispersal_custom(nsp = 3)
#' @export
dispersal_custom = function(nsp, alpha = 0.05, beta = 0.5) {
	if(length(alpha) == 1)
		alpha = rep(alpha, nsp)
	if(length(beta) == 1)
		beta = rep(beta, nsp)
	if(length(alpha) != nsp)
		stop("length(alpha) must be 1 or the number of species")
	if(length(beta) != nsp)
		stop("length(beta) must be 1 or the number of species")
	list(alpha = alpha, beta = beta)
}


#' @rdname community_scenarios
#' @name community_scenarios
#' @title Generate scenarios for site by species matrices
#' @details community_random will ensure that all sites have at least one species.
#' If many species have low prevalence, it may be that the initial random draw produces
#' many zero-richness sites, in which case one species will be selected at random for
#' each of these sites. This can have the side effect that the prevalences in the end do
#' not match the desired prevalence very closely.
#'
#' community_equilibrium simply places species at all sites where they are expected to be present
#' at equilibrium, in the absence of competition (i.e., the fundamental niche is positive).
#' @param rn A [river_network()]
#' @param mc A [metacommunity()]
#' @param prevalence A vector of length 1 or `length(mc[['species']])`, what proportion of sites
#' should be occupied on average by each species.
#' @return A site by species matrix, with sites taken from `rc` and species from `mc`
#' @examples
#' Q = rep(1, 4)
#' adj = matrix(0, nrow = 4, ncol = 4)
#' adj[1,2] = adj[2,3] = adj[4,3] = 1
#' rn = river_network(adj, Q)
#' mc = metacommunity()
#' site_by_species(rn) = community_random(rn, mc)
#' plot(rn, variable = "site_by_species")
#' @export
community_random = function(rn, mc, prevalence = 0.25) {
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

#' @rdname community_scenarios
#' @export
community_equilibrium = function(rn, mc) {
	niches = f_niche(mc, state(rn))
	(f_niche(mc, state(rn)) > 0) * 1
}

