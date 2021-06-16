#' Create a metacommunity
#' @param niches The [niche scenario](niches) to use.
#' @param dispersal The [dispersal scenario](dispersal) to use
#' @param sp_names A vector of species names
#' @param r_names A vector of resource names
#' @param niche_args A list of named arguments to be passed to the niche function
#' @param dispersal_args A list of named arguments to be passed to the dispersal function

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
metacommunity = function(nsp = 2, nr = 1, niches = niches_uniform, dispersal = dispersal_custom,
			sp_names = paste0("sp", 1:nsp), r_names = paste0("r", 1:nr), niche_args = list(),
			dispersal_args = list()) {
	comm = structure(list(), class = "metacommunity")

	## niche parameters
	niche_args$nsp = nsp
	niche_args$nr = nr
	n_params = do.call(niches, niche_args)

	## dispersal params
	dispersal_args$nsp = nsp
	d_params = do.call(dispersal, dispersal_args)
	## create a list of species
	comm[["species"]] = mapply(species, location = n_params$location, breadth = n_params$breadth,
		scale_c = n_params$scale_c, scale_e = n_params$scale_e, r_use = n_params$r_use,
		alpha = d_params$alpha, beta = d_params$beta, SIMPLIFY = FALSE)
	attr(comm, "sp_names") = sp_names
	attr(comm, "r_names") = r_names
	attr(comm, "r_lim") = n_params$r_lim
	attr(comm, "n_species") = nsp
	attr(comm, "n_resources") = nr
#	attr(comm, "n_niche") = nn ##TODO
	comm[["competition"]] = .compute_comp_matrix(comm)

	## for now no immigration from outside the metacommunity
	comm[["boundary"]] = function(n=1) matrix(0, nrow = n, ncol = nsp)
	comm
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
#' @param r_use Resource use scaling parameter, one per niche axis
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
species = function(location, breadth, scale_c, scale_e, alpha, beta, r_use) {
	x = structure(list(), class = "species")
	x$par_c = list(location = location, breadth = breadth, scale = scale_c)
	x$par_e = list(scale = scale_e)
	x$alpha = alpha
	x$beta = beta
	if(length(r_use) == 1)
		r_use = rep(r_use, length(location))
	x$r_use = r_use
	.check_species_params(x)
	x$col = ce_gaussian(location, breadth, scale_c)
	x$ext = ce_constant(scale_e)

	attr(x, "niche_max") = f_niche(x, location)
	return(x)
}

#' Make a competition matrix
#' @param x A [metacommunity()]
#' @keywords internal
#' @return a matrix giving pairwise competition coefficients
.compute_comp_matrix = function(x) {
	sp = x[["species"]]
	xmin = attr(x, "r_lim")[, 1]
	xmax = attr(x, "r_lim")[, 2]
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
	rownames(comp) = colnames(comp) = attr(x, "spnames")
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
