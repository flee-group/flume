#' @rdname niches
#' @name niches
#' @title Generate scenarios for species niches
#' @param nsp Number of species
#' @param nr Number of resources
#' @param nn Number of niche dimensions
#' @param location Location parameter for the niche, see 'details'.
#' @param breadth Niche breadth, see 'details'.
#' @param scale_c Niche height for colonisation, see 'details'.
#' @param scale_e Niche height for extinction, see 'details'.
#' @param dry_c Factor to multiply colonisation rate by when a reach is dry
#' @param dry_e Factor to multiply extinction rate by when a reach is dry
#' @param r_use Resource use scaling parameter, see 'details
#' @param static Vector of resource indices that are static, default is empty vector
#' @param alpha A scalar or a vector, active dispersal ability for each species
#' @param beta A scalar or a vector, passive dispersal ability for each species
#' @param r_trans A function for transforming niche dimensions, default is identity
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
#' For ratio niche dimensions, note that `nr` will be larger than the number of niche dimensions
#' specified in `location` and `breadth`.  Ratio dimensions will always be inserted at the end
#' of the resource matrix, regardless of the position of the original resources, and it is necessary
#' that the `location` and `breadth` parameters reflect this. For example, if the niche dimensions
#' for three resources are the r1:r2 ratio and r3, the `location` parameter should be of length 2,
#' with the first entry being the location for r3, and the second for r1:r2. To avoid confusion,
#' it is recommended (but not required) that users only construct ratio niches from resources that
#' are specified at the end.
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
niches_custom = function(nsp, nr, nn = nr, location, breadth = 1, scale_c = 0.5, scale_e = 0.2, 
		dry_c = 0, dry_e = Inf, r_use = 0.05, static = integer(), alpha = 0.05, beta = 0.5, 
		r_trans = identity) {

	if(!is.matrix(location)) {
		if(nn != 1)
			stop("If nr > 1, location must be a matrix")
		location = matrix(location, ncol = 1)
	}

	if(length(alpha) == 1)
		alpha = rep(alpha, nsp)
	if(length(beta) == 1)
		beta = rep(beta, nsp)
	if(length(dry_c) == 1)
		dry_c = rep(dry_c, nsp)
	if(length(dry_e) == 1)
		dry_e = rep(dry_e, nsp)
	if(length(alpha) != nsp)
		stop("length(alpha) must be 1 or the number of species")
	if(length(beta) != nsp)
		stop("length(beta) must be 1 or the number of species")
	if(length(dry_c) != nsp)
		stop("length(dry_c) must be 1 or the number of species")
	if(length(dry_e) != nsp)
		stop("length(dry_e) must be 1 or the number of species")

	location = lapply(1:nsp, function(i) location[i, ])
	breadth = .check_breadth(breadth, nsp, nn)
	scale_c = .check_scale(scale_c, nsp)
	scale_e = .check_scale(scale_e, nsp)
	r_use = .check_r_use(r_use, nsp, nr, static)


	splist = list()
	for(i in 1:nsp) {
		c_par = list(location = location[[i]], breadth = breadth[[i]], scale = scale_c[[i]], 
			dry = dry_c[[i]])
		e_par = list(scale = scale_e[[i]], dry = dry_e[[i]])
		splist[[i]] = species(c_par, e_par, alpha[i], beta[i], r_use[[i]], r_trans)
	}
	splist
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
	do.call(niches_custom, pars)
}

#' @rdname niches
#' @name niches
#' @export
niches_random = function(nsp = 2, nr = 1, r_lim = c(0, 1), ...) {
	if(!is.matrix(r_lim) && length(r_lim == 2))
		r_lim = cbind(rep(r_lim[1], nr), rep(r_lim[2], nr))


	pars = list(...)
	location = do.call(rbind, lapply(1:nsp, function(i) runif(nr, r_lim[, 1], r_lim[, 2])))
	pars$location = location
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
#' @param allow_empty_sites Logical, if FALSE then empty sites will be filled with a single species
#' at random; this may result in slight departures from the desired prevalence.
#' @return A site by species matrix, with sites taken from `rc` and species from `mc`
#' @examples
#' Q = rep(1, 4)
#' adj = matrix(0, nrow = 4, ncol = 4)
#' adj[1,2] = adj[2,3] = adj[4,3] = 1
#' rn = river_network(adj, Q)
#' st0=rep(0.4,nrow(adj))
#' state(rn, "resources") = st0
#' mc = metacommunity()
#' state(rn, "species") = community_random(rn, mc)
#' plot(rn, variable = "species")
#' state(rn, "species") = community_equilibrium(rn, mc)
#' plot(rn, variable = "species")
#' @export
community_random = function(rn, mc, prevalence = 0.25, allow_empty_sites = TRUE) {
	i = nrow(adjacency(rn))
	j = length(mc$species)
	com = do.call(cbind, mapply(function(jj, pr)
		sample(0:1, i, replace = TRUE, prob = c(1-pr, pr)), 1:j, prevalence, SIMPLIFY = FALSE))

	if(all(prevalence < 1/i))
		warning("All species have low prevalence; random communities might not match desired prevalences")
	r0 = which(rowSums(com) == 0)
	if(length(r0) > 0 & !allow_empty_sites) {
		for(k in r0)
			com[k, sample(1:j)] = 1
	}
	return(com)
}

#' @rdname community_scenarios
#' @export
community_equilibrium = function(rn, mc) {
	niches = f_niche(mc, state(rn, "resources"))
	val = (f_niche(mc, state(rn, "resources")) > 0) * 1
	colnames(val) = attr(mc, "sp_names")
	val
}


#' @rdname r_transform
#' @name r_transform
#' @title Transform resources into niche dimensions
#'
#' @keywords internal
ratio_transform = function(ratio) {
	# take the ratio of the designated columns, and append them to the end, dropping the originals
	return(function(x) {
		cbind(x[,-as.vector(ratio)], 
			apply(ratio, 1, function(i) {
				val = x[,i[1]] / x[, i[2]]
				val[x[, i[2]] == 0] = Inf
				val
			}))
	})
}