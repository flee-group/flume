#' Create a metacommunity
#' @param niches The [niche scenario](niches) to use.
#' @param dispersal The [dispersal scenario](dispersal) to use
#' @param sp_names A vector of species names
#' @param r_names A vector of resource names
#' @param niche_args A list of named arguments to be passed to the niche function
#' @param dispersal_args A list of named arguments to be passed to the dispersal function
#' @param comp_scale Scale for the strength of competition, see 'details'. 
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
#' The `comp_scale` is treated as a multiplier, changing the value of the competition matrix.
#' By default the `comp_scale` parameter is 1/nsp, meaning pairwise competition gets less important
#' with the number of species (but overall competition may still increase, because pairwise
#' interactions are summed). This parameter must be:
#'
#' * a single value, interpreted as the overall strength of competition in the metacommunity,
#' * a vector of length `nsp`, interpreted as the competitive strength of each species, or
#' * a matrix with `nsp` rows and `nsp` columns, if asymmetric competition or fine control over
#' the strength of each interaction is desired. For asymmetric competition, changing the value of
#' `comp_scale[i,j]` will change the effect of species `i` on species `j`, but not the other way
#' around.
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
			dispersal_args = list(), comp_scale = 1/nsp) {
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
		alpha = d_params$alpha, beta = d_params$beta, MoreArgs = list(r_trans = n_params$r_trans), 
		SIMPLIFY = FALSE)
	attr(comm, "sp_names") = sp_names
	attr(comm, "r_names") = r_names
	attr(comm, "r_lim") = n_params$r_lim
	attr(comm, "r_types") = n_params$r_types
	attr(comm, "r_ratio") = n_params$ratio
	attr(comm, "n_species") = nsp
	attr(comm, "n_resources") = nr

	if("ratio" %in% attr(comm, "r_types")) {
		nn = nr - nrow(n_params$ratio)
		attr(comm, "niche_types") = "normal"
		i_r = (nn - nrow(n_params$ratio) + 1):nn  # ratio indices always go on the end
		attr(comm, "niche_types")[i_r] = "ratio"
		n_names_ratio = apply(matrix(r_names[n_params$ratio], ncol=2), 1, 
			function(x) paste(x, collapse=":"))
		if(length(i_r) < nn) {
			i_nr = 1:(nn - nrow(n_params$ratio)) # non-ratio incides
			attr(comm, "niche_names") = c(r_names[i_nr], n_names_ratio)
		} else
			attr(comm, "niche_names") = n_names_ratio
	} else {
		attr(comm, "niche_names") = r_names
		attr(comm, "niche_types") = attr(comm, "r_types")
	}

	## compute reasonable niche limits for plotting and integration
	## this is done rather simply by looking for +/- 2 sds beyond any niche location
	nlocs = niche_par(comm, "location")
	nsds = niche_par(comm, "sd")
	nmins = nlocs - 2 * nsds
	nmaxes = nlocs + 2 * nsds

	# niches/resources have a natural limit of zero
	# note that this means that all habitat variables MUST be on a ratio scale
	# for example, do not use temperature in C, use K
	nmins[nmins < 0] = 0
	attr(comm, "niche_lim") = cbind(apply(nmins, 2, min), apply(nmaxes, 2, max))
	attr(comm, "niche_lims") = list(min = nmins, max = nmaxes)

	comm[["competition"]] = .compute_comp_matrix(comm, comp_scale)

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
#' @param r_use Resource use scaling parameter, one per resource axis
#' @param r_trans A transformation function for converting resources into niche dimensions
#' @return An S3 object of class 'species', which contains the following named elements:
#'   * `col`: The colonisation function, takes a state matrix R and returns a vector of
#'		colonisation rates
#'   * `ext`: The extinction function
#'   * `par_c`: colonisation niche parameters
#'   * `par_e`: extinction niche parameters
#'   * `alpha`: Active dispersal ability
#'   * `beta`: Passive dispersal ability
#'   * `r_use`: Resource use rate per niche axis
#'	 * `r_trans`: Resource-to-niche transformation function
#'
#' Additionally, the following attributes:
#'    * `niche_max`: the maximum possible value of the fundamental niche
#' @examples NULL
#' @export
species = function(location, breadth, scale_c, scale_e, alpha, beta, r_use, r_trans = identity) {
	x = structure(list(), class = "species")
	x$par_c = list(location = location, breadth = breadth, scale = scale_c)
	x$par_e = list(scale = scale_e)
	x$alpha = alpha
	x$beta = beta
	if(length(r_use) == 1)
		r_use = rep(r_use, length(location))
	x$r_use = r_use
	.check_species_params(x)
	x$r_trans = r_trans
	x$col = ce_gaussian(location, breadth, scale_c)
	x$ext = ce_constant(scale_e, length(location))

	attr(x, "niche_max") = f_niche(x, N = location)
	return(x)
}

#' Make a competition matrix
#' @param x A [metacommunity()]
#' @param comp_scale The competition scale
#' @keywords internal
#' @return a matrix giving pairwise competition coefficients
.compute_comp_matrix = function(x, comp_scale) {
	nsp = attr(x, "n_species")

	if(length(comp_scale) == 1) {
		comp_scale = matrix(comp_scale, nrow=nsp, ncol = nsp)
	}
	if(length(comp_scale) == nsp) {
		comp_scale = matrix(comp_scale, nrow = nsp, ncol = nsp)
	}
	if(!is.matrix(comp_scale) || nrow(comp_scale) != nsp || ncol(comp_scale) != nsp)
		stop("comp_scale must be missing, a single value, a vector of length nsp,", 
			"or an nsp by nsp matrix")

	sp = x[["species"]]
	xmin = attr(x, "niche_lims")$min
	xmax = attr(x, "niche_lims")$max
	if(ncol(xmin) == 1) {
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
		comp = matrix(integration_fun(.pairwise_comp(sp[[1]], sp[[1]]), 
			xmin[1,], xmax[1,])[[integral_name]], nrow = 1, ncol = 1)
	} else {
		comp = matrix(0., nrow = length(sp), ncol = length(sp))
		for(i in seq_len(length(sp) - 1)) {
			si = sp[[i]]
			comp[i, i] = integration_fun(.pairwise_comp(si, si), xmin[i,], xmax[i,])[[integral_name]]
			for(j in (i + 1):length(sp)) {
				# need to help the integrator with reasonable limits
				# the only part that matters is the largest niche min and the smallest niche max
				xmn = pmax(xmin[i,], xmin[j,])
				xma = pmin(xmax[i,], xmax[j,])
				sj = sp[[j]]
				comp[i, j] = integration_fun(.pairwise_comp(si, sj), xmn, xma)[[integral_name]]
				comp[j, i] = comp[i, j]
			}
		}
		comp[j, j] = integration_fun(.pairwise_comp(sj, sj), xmin[j,], xmax[j,])[[integral_name]]
	}
	rownames(comp) = colnames(comp) = attr(x, "spnames")
	return(comp * comp_scale)
}


#' Generate a function to integrate to determine competition between two species
#' @param sp1 First species
#' @param sp2 Second species
#' @keywords internal
.pairwise_comp = function(s1, s2) {
	function(x) {
		l1 = f_niche(s1, N = x)
		l2 = f_niche(s2, N = x)
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
#' @param R A site by resource matrix, not needed if `N` is supplied
#' @param N A niche axis matrix; if not supplied, will be computed by transforming `R`
#' @param component One of "lambda", "col", or "ext" specifying whether the total niche (C - E),
#' 		only C, or only E should be returned
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
f_niche.species = function(x, R, N, component = c("lambda", "col", "ext")) {
	component = match.arg(component)
	if(missing(N))
		N = x$r_trans(R)
	if(component == "lambda") {
		val = x$col(N) - x$ext(N)
	} else if(component == "col") {
		val = x$col(N)
	} else {
		val = x$ext(N)
	}
	val
}

#' @rdname niche
#' @export
f_niche.metacommunity = function(x, R, N, component = c("lambda", "col", "ext")) {
	component = match.arg(component)
	res = do.call(cbind, lapply(x$species, f_niche, R = R, N = N, component = component))
	colnames(res) = attr(x, "sp_names")
	res
}


#' @export
niche_par = function(x, ...)
	UseMethod("niche_par", x)

#' Functions for accessing species niche parameters
#' @name nparams
#' @param x A species or metacommunity
#' @export
niche_par.species = function(x, par = c("location", "breadth", "sd", "scale")) {
	par = match.arg(par)
	
	if(par == "sd") {
		val = x$par_c[["breadth"]]
		if(is.matrix(val))
			val = diag(val)
	} else if(par == "scale"){
		val = c(c = x$par_c$scale, e = x$par_e$scale)
	} else {
		val = x$par_c[[par]]
	}
	val
}

#' @name nparams
#' @export
niche_par.metacommunity = function(x, par = c("location", "breadth", "sd", "scale")) {
	par = match.arg(par)
	val = lapply(x$species, niche_par.species, par)
	if(par == "scale") {
		val = do.call(cbind, val)
		rownames(val) = c("c", "e")
		colnames(val) = attr(x, "sp_names")
	} else if(par != "breadth") {
		val = do.call(rbind, val)
		rownames(val) = attr(x, "sp_names")
		colnames(val) = attr(x, "niche_names")
	}
	val
}