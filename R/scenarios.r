#' Generate a random site by species matrices
#' @details This function will ensure that all sites have at least one species.
#' If many species have low prevalence, it may be that the initial random draw produces
#' many zero-richness sites, in which case one species will be selected at random for
#' each of these sites. This can have the side effect that the prevalences in the end do
#' not match the desired prevalence very closely.
#' @param rn A [river_network()]
#' @param mc A [metacommunity()]
#' @param prevalence A vector of length 1 or `length(mc[['species']])`, what proportion of sites should
#' be occupied on average by each species.
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
random_site_by_species = function(rn, mc, prevalence = 0.25) {
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


#' Create niches for a metacommunity
#' @param nsp The number of species
#' @param nrx The number of niche axes
#' @param xmin The minimum value(s) to use for generating niches
#' @param xmax The maximum value(s) for generating niches
#' @param k The number of niche dimensions
#' @details For a single niche axis, niches will be evenly spaced along the axis. For multiple
#' dimensions, niches locations are randomly generated
#' @return A list with location, breadth, and scale parameters to pass to [metacommunity()].
#' @export
generate_niches = function(nsp, nrx, xmin = rep(0, nrx), xmax = rep(1, nrx)) {

	## for single variable, spread niches evenly
	if(nrx == 1) {
		## for small numbers of species, we create some buffer around the edges
		if(nsp < 6) {
			location = seq(xmin, xmax, length.out = nsp+2)
			location = location[2:(length(location)-1)]
		} else {
			location = seq(xmin, xmax, length.out = nsp)
		}
		breadth = location[2] - location[1]
	} else {
		if(length(xmin) == 1)
			xmin = rep(xmin, nrx)
		if(length(xmax) == 1)
			xmax = rep(xmax, nrx)
		location = do.call(rbind, lapply(1:nsp, function(i) runif(nrx, xmin, xmax)))
		sd_avg = (xmax - xmin) / nsp
		sd_sd = sd_avg / 4
		breadth = lapply(1:nsp, function(i) {
			y = matrix(0, nrow=nrx, ncol=nrx)
			diag(y) = rnorm(nrx, sd_avg, sd_sd)
			y
		})
	}
	list(location = location, breadth = breadth, scale_c = rep(0.5, nsp), scale_e = rep(0.2, nsp))
}

