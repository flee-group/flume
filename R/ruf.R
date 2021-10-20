#' Generate a resource use function
#'
#' @details Resource use functions describe the impact a species has on the state of a location.
#' This impact will commonly be zero (for static habitat features, for example), but for many
#' resources (e.g., nutrient concentrations) the result will be negative (i.e., species deplete
#' resources).
#'
#' This function is intended as an example, and can be replaced with a user-defined function.
#' The function must take a site by species matrix and a resource state matrix as its first two
#' arguments, and it must return a matrix with the same dimensions as R giving the instantaneous
#' rate of change in each resource. The units must match whatever is being used for the
#' reaction-transport portion of the model; often something like $g-Resource L^{-1} min^{-1}$.
#'
#' The default behaviour is for species to consume more resources the closer they are to the
#' niche optimum. This is done by computing the value of the niche at the current concentration
#' (`sp$col(R) - sp$ext(R)`) using [f_niche()] and taking the ratio with the `niche_max` attribute
#' for each species. This ratio is then multiplied by the `r_scale` of each species to get the rate
#' of change. This leads to a convenient definition of the r_scale, which is the number of units
#' of resources depleted by a species at its maximum growth, per unit time.
#'
#' @param x A site by species matrix
#' @param R A resource state matrix
#' @param C A [metacommunity()]
#' @return A matrix of the same dimensions as `R` giving the rate of change of each resource
#' @examples
#'	comm = metacommunity()
#'	Q = rep(1, 4)
#'	adj = matrix(0, nrow = 4, ncol = 4)
#'	adj[1,2] = adj[2,3] = adj[4,3] = 1
#'	st = matrix(seq(0, 1, length.out = length(Q)), ncol = 1, dimnames = list(NULL, 'R'))
#'	rn = river_network(adj, Q, state = st)
#'	site_by_species(rn) = matrix(1, nrow = length(Q), ncol = length(comm$species))
ruf = function(x, R, C) {
	niche_max = niche_par(C, "niche_max")
	r_use_scale = niche_par(C, "r_use")
	n_ht = sweep(f_niche(C, R), 2, niche_max, "/")

	## take into account presence-absence
	n_ht = x * n_ht

	# we assume that species that are outside their niche have a negligible effect on resource use
	n_ht[n_ht < 0] = 0

	# this produces the resources consumed, so the rate of change will have the opposite sign
	r_use = n_ht %*% r_use_scale
	-1 * r_use * R
}
