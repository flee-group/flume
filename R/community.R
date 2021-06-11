#' Create a metacommunity
#' @param location Location parameter for the niche, see 'details'.
#' @param breadth Niche breadth, see 'details'.
#' @param scale_c Niche height for colonisation, see 'details'.
#' @param scale_e Niche height for extinction, see 'details'.
#' @param alpha A scalar or a vector, active dispersal ability for each species
#' @param beta A scalar or a vector, passive dispersal ability for each species
#' @param r_scale Resource use scaling parameter, see 'details
#' @param niche_lim Minimum and maximum possible values for each niche dimension; can be a vector
#'		in the case of a single resource, otherwise a matrix with 2 columns and nrx rows.
#' @details The metacommunity describes the possible biological space for a `flume` model. It
#' consists of a few named elements:
#' * `species`: a list of [species()] objects, describing the niche and dispersal parameters
#' for each possible species
#' * `competition`: A square species interaction matrix describing the effect of each species
#' on each other species; generated using the overlap between species'
#' [fundamental niches](f_niche.species()).
#' * `boundary`: A function that returns a site by species matrix giving the colonisation flux
#' from outside the river network; by default returns zero for all sites/species.
#'
#' Each species' fundamental niche is defined as the difference between a Gaussian colonisation
#' function and a constant extinction function; a species will be present when c - e > 0.
#' Niches can be automatically generated (see [generate_niches()]; alternatively, to specify
#' the parameters manually, the following can be specified.
#'
#' `location`: Required; the location of the niche optimum for each species-resource combination.
#' Must be a matrix with one row per species (`nsp`), one column per niche axis (`nrx`). For a
#' single niche axis, a vector with one entry per species is also accepted. The location is the
#' mean of the Gaussian colonisation function.
#'
#' `scale_c`, `scale_e`: The (relative) height of the Gaussian colonisation function or
#' the constant extinction function, must be a positive real number.
#' Can be supplied as a scalar (all species have the same scale) or a vector of length `nsp`.
#' If missing a default value of 0.5 (for colonisation) or 0.2 (for extinction) will be used.
#'
#' `breadth`: Optional, niche breadth; standard deviation of the Gaussian colonisation function;
#' larger values indicate species can occur in a wider variety of environments. For a single
#' niche dimension, either a scalar (all species have the same breadth) or a vector of length `nsp`.
#' For multivariate niches, the following are possible:
#'
#'   * A **scalar**: all species have the same breadth in each niche dimension
#'   * A **vector** of length `nsp`: each species has a different breadth, but breadth is the same
#'   	for each niche axis.
#'   * A **vector** of length `nrx`: niche breadth varies by axis, but all species have the same
#'   	breadth for each niche axis. If `nrx == nsp`, then this option is not possible, instead the
#'   	former applies.
#'   * A **matrix** of `nsp` rows and `nrx` columns: each species-resource combo has a unique
#'   	breadth, but all resources are orthogonal. In other words, performance for a given resource
#'   	cannot depend on the concentration of any other resource.
#'   * A **list** of length `nsp`, each element is a square symmetric `nrx` by `nrx` matrix. Similar
#'   	to the matrix above, but a full variance-covariance matrix is supplied for each species,
#'   	describing the breadth of the niche along each axis but also how the axes covary.
#'
#' Dispersal parameters are `alpha` and `beta`, for specifying active (i.e., either up- or
#' downstream) or passive (downstream only) dispersal, respectively. These can either be a scalar,
#' in which case all species have the same dispersal ability, or a vector of length `nsp`.
#'
#' Resource use by each species can be scaled using the `r_scale` parameter; larger values indicate
#' faster consumption of resources, see [ruf()] for more details. This can be a scalar (all species
#' consume all resources at the same rate relative to niche position), a vector of length `nsp`
#' (species behave differently, but all resources within species are consumed identically),
#' a vector of length `nrx` (all species have identical behaviour, but each resource is consumed at
#' a different rate), or a matrix with `nsp` rows and `nrx` columns. Note that immutable niche axes
#' (see [river_network()]) ignore this parameter.
#'
#' @return A metacommunity object
#' @examples
#' \donttest{
#'		# defaults
#'		mc = metacommunity()
#'
#'		# multiple species
#'		mc = metacommunity(location = c(1, 2))
#'
#'		# multiple resources
#'		nlim = matrix(c(0, 0, 1, 1), ncol = 2) # we need to specify the limits manually
#'		mc = metacommunity(location = matrix(c(1, 2), ncol = 2), niche_lim = nlim)
#'
#'		# multiple species and resources with variance covariance matrices
#'		loc = matrix(c(1, 2, 1, 3, 2, 4), ncol = 2) # 3 species, 2 resources
#'		br_sp1 = matrix(c(1, 0.2, 0.2, 1), ncol = 2) # vcv matrix for one species
#'		breadth = list(br_sp1, br_sp1, br_sp1) # assume the same matrix for all species
#'		mc = metacommunity(location = loc, breadth = bre, niche_lim = nlim)
#' }
#' @export
metacommunity = function(location, breadth = 1, scale_c = 0.5, scale_e = 0.2, alpha = 0.05,
						 beta = 0.5, r_scale = 0.05, niche_lim = c(0, 1)) {
	comm = structure(list(), class = "metacommunity")

	## niche parameters
	if(!is.matrix(location))
		location = matrix(location, ncol = 1, dimnames = list(names(location), "r1"))
	nsp = nrow(location)
	nrx = ncol(location)
	if(is.null(rownames(location))) {
		spnames = paste0("sp", 1:nsp)
	} else {
		spnames = rownames(location)
	}

	if(is.null(colnames(location))) {
		rnames = paste0("r", 1:nrx)
	} else {
		rnames = colnames(location)
	}

	location = lapply(1:nsp, function(i) location[i, ])
	breadth = .check_breadth(breadth, nsp, nrx)
	scale_c = .check_scale(scale_c, nsp)
	scale_e = .check_scale(scale_e, nsp)
	r_scale = .check_r_scale(r_scale, nsp, nrx)

	if(!is.matrix(niche_lim))
		niche_lim = matrix(niche_lim, ncol = 2)

	if(nrow(niche_lim) != nrx)
		stop("niche_lim must be a matrix with 1 row per resource and 2 columns")

	if(length(alpha) == 1)
		alpha = rep(alpha, nsp)
	if(length(beta) == 1)
		beta = rep(beta, nsp)
	if(length(alpha) != nsp)
		stop("length(alpha) must be 1 or the number of species")
	if(length(beta) != nsp)
		stop("length(beta) must be 1 or the number of species")

	## create a list of species
	comm[["species"]] = mapply(species, location = location, breadth = breadth, scale_c = scale_c,
		   scale_e = scale_e, alpha = alpha, beta = beta,
		   r_scale = lapply(seq_len(nrow(r_scale)), function(i) r_scale[i, ]), SIMPLIFY = FALSE,
		   USE.NAMES = FALSE)
	comm[["competition"]] = .compute_comp_matrix(comm$species, niche_lim[, 1], niche_lim[, 2])
	comm[["boundary"]] = function(n=1) matrix(0, nrow = n, ncol = n_species)
	attr(comm, "niche_lim") = niche_lim
	attr(comm, "spnames") = spnames
	attr(comm, "rnames") = rnames
	rownames(comm[["competition"]]) = colnames(comm[["competition"]]) = attr(comm, "spnames")
	comm
}



#' @name check_par
#' @title Check parameter dimensions
#' @rdname check_par
#' @keywords internal
.check_scale = function(p, nsp) {
	if(length(p) == 1)
		p = rep(p, nsp)
	if(length(p) != nsp)
		stop("Scale parameters must be missing, a single value, or 1 value per species")
	p
}

#' @rdname check_par
#' @keywords internal
.check_breadth = function(p, nsp, nrx) {
	if(length(p) == 1)
		p = rep(p, nsp)

	if(!is.list(p) && length(p) == nsp) {
		p = matrix(p, nrow = nsp, ncol = nrx)
	} else if(!is.list(p) && length(p) == nrx) {
		p = matrix(p, nrow = nsp, ncol = nrx, byrow = TRUE)
	}

	## univariate niches
	if(nrx == 1) {
		if(!(is.matrix(p) && nrow(p) == nsp && ncol(p) == 1))
			stop("Invalid niche breadth: must be a scalar, vector, or one-column matrix")
		## convert back to a vector so it plays nicely with mapply
		p = p[, 1]
	} else {
		## for multivariate niches
		p = .check_breadth_multi(p, nsp, nrx)
	}
	p
}

#' @rdname check_par
#' @keywords internal
.check_breadth_multi = function(p, nsp, nrx) {
	if(is.matrix(p)) {
		if(nrow(p) != nsp || ncol(p) != nrx)
			stop("Invalid niche breadth; see ?metacommunity")
		p = lapply(1:nsp, function(i) {
			mat = matrix(0, nrow = nrx, ncol = nrx)
			diag(mat) = p[i, ]
			mat
		})
	}

	## test proper dimensions of the breadth param for multivariate niches
	## must be a list, one entry per species
	## each entry a square matrix with dimension nrx
	if(!(is.list(p) &&
		length(p) == nsp &&
		all(sapply(p, is.matrix)) &&
		all((sapply(p, nrow) - nrx) == 0) &&
		all((sapply(p, ncol) - nrx) == 0))) {
		stop("Invalid niche breadth; see ?metacommunity")
	}
	p
}

#' @rdname check_par
#' @keywords internal
.check_r_scale = function(p, nsp, nrx) {
	if(length(p) == 1)
		p = rep(p, nsp)

	if(!is.matrix(p) && length(p) == nsp) {
		p = matrix(p, nrow = nsp, ncol = nrx)
	} else if(!is.matrix(p) && length(p) == nrx) {
		p = matrix(p, nrow = nsp, ncol = nrx, byrow = TRUE)
	}

	if(!(is.matrix(p) && nrow(p) == nsp && ncol(p) == nrx))
		stop("r_scale must be a matrix with one row per species and one column per niche axis")
	p
}

#' @rdname check_par
#' @keywords internal
.check_species_params = function(x) {
	loc = x$par_c$location
	bre = x$par_c$breadth
	csc = x$par_c$scale
	rsc = x$r_scale
	esc = x$par_e$scale
	alpha = x$alpha
	beta = x$beta

	if(length(csc) != 1 | length(esc) != 1 | length(alpha) != 1 | length(beta) != 1)
		stop("The following parameters must be single values: ",
			"c_scale, e_scale, alpha, beta")

	if(any(bre < 0) | csc < 0 | esc < 0 | alpha < 0 | beta < 0)
		stop("The following parameters must not be negative: ",
			"breadth, scale_c, scale_e, alpha, beta")

	if(length(loc) == 1 & length(bre) != 1)
		stop("dimension mismatch: breadth and r_scale must be length 1 with a single niche",
				"dimension")
	if(length(loc) > 1 & (!is.matrix(bre) | !identical(dim(bre), rep(length(loc), 2))))
		stop("dimension mismatch: breadth must be a square matrix with one row/column per ",
			"niche dimension")
	if(length(rsc) != length(loc))
		stop("dimension mismatch: r_scale must have the same length as location")
}




#' Returns dispersal parameters for a metacommunity
#' @param x A [metacommunity()]
#' @return Named list with two elements, `alpha` and `beta`, each one a vector of dispersal
#'		parameters for each species
#' @export
#' @examples
#' comm = community()
#' dispersal_params(comm)
dispersal_params = function(x) {
	alpha = sapply(x$species, function(y) y$alpha)
	beta = sapply(x$species, function(y) y$beta)
	list(alpha = alpha, beta = beta)
}

#' Creates species
#'
#' @details Creates a species, which is defined by its niche parameters and its
#' colonisation/extinction functions.
#'
#' Note that for gaussian functions in multiple variables, `scale` is always a scalar, but
#' `location` must be a location vector of length `nx`, and `width` must be a variance covariance
#' matrix with `dim = c(nx, nx)`.
#'
#' The resource use function is set automatically by default, but can be set to an arbitrary
#' function; see [resource use functions][ruf()] for details on how this function should behave.
#'
#' @param location Location of the niche optimum; one per niche axis
#' @param breadth Niche breadth, or variance-covariance matrix for multivariate
#' @param scale_c Scale of the colonisation portion of the niche
#' @param scale_e Scale of the extinction portion of the niche
#' @param alpha Active dispersal ability
#' @param beta Passive dispersal ability
#' @param r_scale Resource use scaling parameter, one per niche axis
#' @return An S3 object of class 'species', which contains the following named elements:
#'   * `col`: The colonisation function, takes a state matrix R and returns a vector of
#'		colonisation rates
#'   * `ext`: The extinction function
#'   * `par_c`: colonisation niche parameters
#'   * `par_e`: extinction niche parameters
#'   * `alpha`: Active dispersal ability
#'   * `beta`: Passive dispersal ability
#'   * `r_use`: Resource use rate per niche axis
#'
#' Additionally, the following attributes:
#'    * `niche_max`: the maximum possible value of the fundamental niche
#' @examples NULL
#' @export
species = function(location, breadth, scale_c, scale_e, alpha, beta, r_scale) {
	x = structure(list(), class = "species")
	x$par_c = list(location = location, breadth = breadth, scale = scale_c)
	x$par_e = list(scale = scale_e)
	x$alpha = alpha
	x$beta = beta
	if(length(r_scale) == 1)
		r_scale = rep(r_scale, length(location))
	x$r_scale = r_scale
	.check_species_params(x)
	x$col = ce_gaussian(location, breadth, scale_c)
	x$ext = ce_constant(scale_e)

	attr(x, "niche_max") = f_niche(x, location)
	return(x)
}

#' Make a competition matrix
#' @param sp A list of [species()]
#' @param xmin The minimum for integration
#' @param xmax The maximum for integration
#' @keywords internal
#' @return a matrix giving pairwise competition coefficients
.compute_comp_matrix = function(sp, xmin, xmax) {
	if(length(xmin) == 1) {
		integration_fun = integrate
		integral_name = "value"
	} else {
		if(!requireNamespace("cubature", quietly = TRUE))
			stop("Multidimensional niches require the 'mvtnorm' and 'cubature' packages")
		integral_name = "integral"
		if(length(xmin) < 4) {
			integration_fun = cubature::pcubature
		} else {
			integration_fun = cubature::hcubature
		}

	}
	## guard against case when there is only one species
	if(length(sp) == 1) {
		comp = matrix(integration_fun(.pairwise_comp(sp[[1]], sp[[1]]), xmin, xmax)[[integral_name]],
			nrow = 1, ncol = 1)
	} else {
		comp = matrix(0., nrow = length(sp), ncol = length(sp))
		for(i in seq_len(length(sp) - 1)) {
			si = sp[[i]]
			comp[i, i] = integration_fun(.pairwise_comp(si, si), xmin, xmax)[[integral_name]]
			for(j in (i + 1):length(sp)) {
				sj = sp[[j]]
				comp[i, j] = integration_fun(.pairwise_comp(si, sj), xmin, xmax)[[integral_name]]
				comp[j, i] = comp[i, j]
			}
		}
		comp[j, j] = integration_fun(.pairwise_comp(sj, sj), xmin, xmax)[[integral_name]]
	}
	return(comp)
}


#' Generate a function to integrate to determine competition between two species
#' @param sp1 First species
#' @param sp2 Second species
#' @keywords internal
.pairwise_comp = function(s1, s2) {
	function(x) {
		l1 = f_niche(s1, x)
		l2 = f_niche(s2, x)
		l1[l1 < 0] = 0
		l2[l2 < 0] = 0
		ifelse(l1 > l2, l2, l1)
	}
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
	do.call(cbind, lapply(x$species, f_niche, R = R))
}
