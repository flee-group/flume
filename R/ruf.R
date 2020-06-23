#' Generate a resource use function
#'
#' @details Resource use functions describe the impact a species has on the state of a location. This impact will commonly
#' be zero (for static habitat features, for example), but for many resources (e.g., nutrient concentrations)
#' the result will be negative (i.e., species deplete resources).
#'
#' The `ruf` in general must take two parameters; the first is a species, the second is a vector of resource values
#'  at a single location; the `ruf` must return a vector of the same length giving the instantaneous rate of change in
#' each resource due to this species. The units must match whatever is being used for the reaction-transport portion of
#' the model; often something like $g-Resource L^{-1} min^{-1}$.
#'
#' The default behaviour is for species to consume more resources the closer they are to the niche optimum. This is done
#' by computing the value of the niche at the current concentration (`sp$col(R) - sp$ext(R)`) and taking the ratio with
#' the `niche_max` attribute (note that for linear niches, this value can exceed one if resource concentrations grow
#' beyond the original definition). This ratio is then multiplied by the `r_scale` to get the rate of change. This
#' leads to a convenient definition of the scale, which is the number of units of resources depleted by a species at its
#' maximum growth, per unit time.
#'
#' @param r_scale The scale of the function, see 'details'
#' @param sp A [species()]
#' @param R A resource state vector
#' @seealso [plot_ruf()]
#' @return A resource use function, suitable for adding to a [species()]. This function takes `sp` and `R` as parameters
#' and returns a vector of the same length as `R` giving the rate of change of each resource.
#' @examples
#' sp = species('linear', 'constant')
ruf = function(r_scale) {

	fu = function(sp, R) {
		n_val = sp$col(R) - sp$
	}

	## hmm, maybe this isn't the way to go... might be better if ruf is just a function
	## that takes species and R and returns the resource use, not a property of the 
	## species
}
