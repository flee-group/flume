#' Colonisation and extinction probability, given the current state of the system
#' @name ceprob
#'
#' @details For debugging and tuning, it can be useful to examine the individual components of the
#' colonisation/extinction rate by setting `components = TRUE`, in which case a named list of site
#' by species matrices is returned, one element per component.
#'
#' @param comm A metacommnunity
#' @param network A river network
#' @param dt Time interval
#' @param components Boolean, normally FALSE, see 'details'
#'
#' @return Normally, a site by species matrix of colonisation/extinction probabilities, but see 'details'
#' @export
col_prob = function(comm, network, dt, components = FALSE) {
	R = state(network)
	nsites = nrow(R)
	nsp = length(comm$species)

	## colonisation rate has two terms, the niche portion and the dispersal portion
	## first the niche portion
	col = do.call(cbind, lapply(comm$species, function(sp) sp$col(R)))


	## here the dispersal portion
	P = prevalence(network)
	A = matrix(dispersal_params(comm)$alpha, nrow=nsites, ncol = nsp, byrow = TRUE) # repeat the alphas across all sites in a matrix
	B = matrix(dispersal_params(comm)$beta, nrow=nsites, ncol = nsp, byrow = TRUE) # same for beta
	Q = matrix(network$discharge, nrow=nsites, ncol = nsp) # repeat discharge by site across all species
	immigration = matrix(as.vector(comm$boundary()), nrow = nsites, ncol = nsp, byrow = TRUE)
	dispersal = P * (A + B*Q + immigration)

	dimnames(col) = dimnames(dispersal) = dimnames(site_by_species(network))
	if(components) {
		return(list(colonisation = col, dispersal = dispersal))
	} else {
		res = 1 - exp(-1 * col * dispersal * dt)
		dimnames(res) = dimnames(site_by_species(network))
		return(res)
	}
}

#' @rdname ceprob
#' @export
ext_prob = function(comm, network, dt, components = FALSE) {
	S = site_by_species(network)
	R = state(network)

	# extinction rate has two components, the stochastic extinction rate and the competition portion
	# stochastic first
	m_i = do.call(cbind, lapply(comm$species, function(sp) sp$ext(R)))

	# competition, which is the sum of the competitive effects of species present in each site
	m_ij = comm$competition ## species by species competition matrix
	comp = S %*% m_ij * S ## this is witchcraft

	dimnames(m_i) = dimnames(comp) = dimnames(S)
	if(components) {
		return(list(stochastic = m_i, competition = comp))
	} else {
		res = 1 - exp(-1 * (m_i + comp) * dt)
		dimnames(res) = dimnames(S)
		return(res)
	}

}

#' Compute time-derivative for resources
#' @details For debugging and tuning, it can be useful to examine the individual components of the
#' derivative rate by setting `components = TRUE`, in which case a named list of site
#' by resource matrices is returned, with separate transport and reaction components.
#' @param comm A metacommnunity
#' @param network A river network
#' @param components Boolean, normally FALSE, see 'details'
#' @return Normally, a site by resource matrix of resource fluxes
#' @export
dRdt = function(comm, network, components = FALSE) {
	S = site_by_species(network)
	R = state(network)
	Q = network$discharge
	Ru = t(adjacency(network)) %*% R
	Qu = t(adjacency(network)) %*% Q
	A = network$area
	l = reach_length(network)

	output = Q * R
	input = apply(Ru, 2, function(x) Qu * x) + network$boundary() * lateral_discharge(network)
	transport = (output - input) / (A*l)

	rxn = ruf(S, R, comm)

	dimnames(transport) = dimnames(rxn) = dimnames(R)
	if(components) {
		return(list(transport = transport, reaction = rxn))
	} else {
		res = rxn - transport
		dimnames(res) = dimnames(R)
		return(res)
	}
}
