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
#' @return Normally, a site by species matrix of colonisation/extinction probabilities, 
#' 	but see 'details'
#' @export
col_prob = function(comm, network, dt, components = FALSE) {
	R = state(network, "resources")
	nsites = nrow(R)
	nsp = length(comm$species)

	## colonisation rate has two terms, the niche portion and the dispersal portion
	## first the niche portion
	col = f_niche(comm, R = R, component = "col")

	## here the dispersal portion
	P = prevalence(network)

	# repeat the alphas and betas across all sites in a matrix
	A = matrix(dispersal_params(comm)$alpha, nrow=nsites, ncol = nsp, byrow = TRUE)
	B = matrix(dispersal_params(comm)$beta, nrow=nsites, ncol = nsp, byrow = TRUE)

	# repeat discharge by site across all species
	Q = matrix(discharge(network), nrow=nsites, ncol = nsp)


	dispersal = P * (A + B*Q) + boundary(network, "species")

	dimnames(col) = dimnames(dispersal) = dimnames(state(network, "species"))
	if(components) {
		return(list(colonisation = col, dispersal = dispersal))
	} else {
		res = 1 - exp(-1 * col * dispersal * dt)
		dimnames(res) = dimnames(state(network, "species"))
		return(res)
	}
}

#' @rdname ceprob
#' @export
ext_prob = function(comm, network, dt, components = FALSE) {
	S = state(network, "species")
	R = state(network, "resources")

	# extinction rate has two components, the stochastic extinction rate and the competition portion
	# stochastic first
	m_i = f_niche(comm, R = R, component = "ext")

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


#' Prepare parameter list for dRdt
#' @param C a metacommunity
#' @param N a river network
#' @return Parameter list, to be passed to dRdt function
#' @keywords internal
.dRdt_params = function(C, N) {
	list(S =  state(N, "species"),
		Q =  state(N, "Q"),
		A =  state(N, "area"),
		l = reach_length(N),
		lQ = boundary(N, "Q"),
		lR = boundary(N, "resources"),
		i_static = which(attr(C, "r_types") == "static"),
		adj = adjacency(N),
		comm = C
	)
}

#' Compute time-derivative for resources
#' @details The required params for the resource derivative are as follows:
#'
#'  * S: a site by species matrix
#'  * Q: a site by discharge vector
#'  * A: a site by cross-sectional area vector
#'  * l: a site by length vector
#'  * lQ: site by lateral input discharge
#'  * lR: site by lateral resource concentration
#'  * i_static: indices of static resources
#'  * adj: the network adjacency matrix
#'  * comm: a [metacommunity()]
#'  * components: logical, optional, default FALSE
#'
#' If `params$components = TRUE`, instead of computing the derivative the function returns
#' the individual components, including the use by species, downstream transport, and input from
#' upstream. The units in this case will be mass/time
#' @param t Time variable
#' @param R A resource state matrix; stored as a vector in column major order
#' @param params Named list; additional parameters, see 'details'
#' @return Normally, a site by resource matrix of resource fluxes, in units of 
#' mass * volume^-1 * time^-1
#' @export
dRdt = function(t, R, params, components = FALSE) {
	S = params$S
	Q = params$Q
	A = params$A
	l = params$l
	lQ = params$lQ
	lR = params$lR
	i_static = params$i_static
	adj = params$adj
	comm = params$comm
	components = ifelse("components" %in% names(params), params$components, FALSE)

	if(!is.matrix(R))
		R = matrix(R, ncol = ncol(lR), nrow=nrow(lR))
	# mass flux out of each site
	output = Q * R

	# when lateral discharge is NEGATIVE (i.e., the stream is shrinking) we export based on the
	# concentration in the stream, not the boundary condition
	# TODO: this will need to be improved for intermittent rivers
	if(any(lQ < 0)) {
		for(i in 1:ncol(lR))
			lR[lQ < 0, i] = R[lQ < 0, i]
	}

	# mass flux into each site
	input = t(adj) %*% output + apply(lR, 2, function(x) lQ * x)

	# convert back to concentration
	transport = (output - input) / (A*l)

	# reaction component, in concentration units
	rxn = ruf(S, R, comm)

	if(length(i_static) > 0) {
		transport[, i_static] = 0
		input[, i_static] = 0
		output[, i_static] = 0
		rxn[, i_static] = 0
	} 

	if(components) {
		# return the components of the derivative, in mass units
		dimnames(output) = dimnames(rxn) = dimnames(input) = dimnames(R)
		return(list(resource_use = rxn*A*l, downstream_transport = output, input_transport = input))
	} else {
		# return the derivative in concentration units
		res = rxn - transport
		dimnames(res) = dimnames(R)
		return(list(res))
	}
}
