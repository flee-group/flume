#' flume: A package for modelling FLUvial Meta-Ecosystems
#'
#' Implementation of a theoretical/mechanistic model for exploring fluvial meta-ecosystems.
#'
#' @section Key functions:
#'
#' * [metacommunity()] Create species pools and associated attributes
#' * [river_network()] Build a river network object
#' * [flume()] Creates a model
#' * [run_simulation()] Runs the simulation
#'
#' @docType package
#' @name flume_package
NULL


#' Create a FLUvial Meta-Ecosystem model
#'
#' @param comm A [metacommunity()]
#' @param network A [river_network()]
#' @param dt The time step for the model, in the same time units as all fluxes
#'
#' @return An S3 object of class 'flume'
#' @export
flume = function(comm, network, dt) {
	structure(list(metacom = comm, network = network, dt = dt), class="flume")
}

#' Run the simulation for a specified number of time steps
#'
#' This function begins with the current state of `x`, runs the model for `nt` steps,
#' then returns a copy of `x` with the updated state and state history.
#'
#' @param x A [flume()] object
#' @param nt The number of time steps
#' @return A modified copy of `x`, with state updated with the results of the simulation
#' @export
run_simulation = function(x, nt) {
	if(nt < 1)
		stop("at least one time step is required")
	MC = x[['metacom']]
	RN = x[['network']]
	dt = x[['dt']]
	R = state(RN)
	S = site_by_species(RN)

	for(tstep in 1:nt) {
		cp = col_prob(MC, RN, dt = dt)
		ep = ext_prob(MC, RN, dt = dt)
		dR = dRdt(MC, RN)

		# no cols where already present, no exts where already absent
		cp[S == 1] = 0
		ep[S == 0] = 0
		cols = matrix(rbinom(length(cp), 1, cp), nrow=nrow(cp))
		exts = matrix(rbinom(length(ep), 1, ep), nrow=nrow(ep))
		S[cols == 1] = 1
		S[exts == 1] = 0

		## euler integration for now
		R = R + dR * dt

		site_by_species(RN) = S
		state(RN) = R
	}
	x[['network']] = RN
	return(x)
}
