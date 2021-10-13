#' Create a FLUvial Meta-Ecosystem model
#'
#' @param comm A [metacommunity()]
#' @param network A [river_network()]
#' @param sp0 Initial site by species matrix
#' @param st0 Initial state of the model
#' @param spb Boundary condition for species
#' @param stb Boundary condition for resources
#' @param dt The time step for the model, in the same time units as all fluxes
#' @details `sp0`, `st0`, `spb`, and `stb` are optional parameters that define the initial
#' states (`*0`) and boundary conditions (`*b`) for the species (`sp*`) and resources (`st*`).
#' Initial states must be a site-by-(species or resource) matrix. Boundary conditions should be
#' as described in [river_network()].
#' 
#' If the species and resource state variables are not defined in the network, they must be
#' specified here. If they are specified here, they will override any state set in the network.
#'
#' Boundary conditions, if not defined, will be set automatically; the boundary condition for
#' resources will be equal to the initial state, and for species will be set to zero globally
#' (i.e., no immigration from outside the network).
#' @return An S3 object of class "flume"
#' @examples
#' data(algae)
#' model = flume(algae$metacommunity, algae$network, algae$sp0, algae$r0)
#' @export
flume = function(comm, network, sp0, st0, spb, stb, dt = 86400) {

	# if(missing(st0)) {
	# 	if(is.null(state(network, "resources")))
	# 		stop("Initial resource state not specified")
	# 	st0 = state(network)
	# }
	# if(ncol(st0) != attr(comm, "n_r"))
	# 	stop("Initial network state doesn't match number of resources in comm")
	# colnames(st0) = attr(comm, "r_names")
	# state(network, "resources") <- NULL
	# state(network, "resources") <- st0
	## initial state for species
	# if(missing(sp0)) {
	# 	if(is.null(state(network, "species")))
	# 		stop("Initial site by species matrix not specified")
	# 	sp0 = state(network, "species")
	# }		
	# 	if(ncol(sp0) != attr(comm, "n_species"))
	# 	stop("Initial site by species doesn't match number of species in comm")
	# colnames(sp0) = attr(comm, "sp_names")
	# state(network, "species") <- NULL
	# state(network, "species") <- sp0

	## define initial states and check dimensionality
	network = .set_flume_initial_state(network, comm, st0, "resources")
	network = .set_flume_initial_state(network, comm, sp0, "species")

	# set resource boundary conditions
	if(missing(stb))
		stb = state(network, "resources")
	i_static = which(attr(comm, "r_types") == "static")
	if(length(i_static) > 0) {
		stb[i_static] = 0
	}
	boundary(network, "resources") = stb

	if(missing(spb))
		spb = state(network, "species") * 0
	boundary(network, "species") = spb
	x = structure(list(metacom = comm, networks = list(network), dt = dt), class = "flume")
	x
}

.set_flume_initial_state = function(rn, comm, init, type = c("resources", "species")) {
	type = match.arg(type)
	n = ifelse(type == "resources", "n_r", "n_species")
	nm = ifelse(type == "resources", "r_names", "sp_names")

	if(missing(init)) {
		if(is.null(state(rn, type)))
			stop("Initial ", type, " state not specified")
		init = state(rn, type)
	}
	if(ncol(init) != attr(comm, n))
		stop("Initial network state doesn't match number of ", type, " in comm")
	colnames(init) = attr(comm, nm)
	state(rn, type) <- NULL
	state(rn, type) <- init
	rn
}

#' Compute post simulation statistics
#' @name flumestats
#' @param x A [flume()]
#' @import data.table
#' @return A [data.table::data.table()] with the desired summary statistics
#' @export
occupancy = function(x, type = c("network", "reach")) {
	type = match.arg(type)
	occ = .reshape_sim(x, variable = "species")
	if(type == "network") {
		keyby = c("species", "time")
	} else {
		keyby = c("species", "reach", "time")
	}
	return(occ[, .(occupancy = sum(occupancy) / .N), keyby = keyby])
}

#' @rdname flumestats
#' @export
resource_summary = function(x) {
	res = .reshape_sim(x, variable = "resources")
	return(res[, .(concentration = mean(concentration)), keyby = .(resource, reach, time)])
}

#' Get a long-format occupancy by time by reach by network dataset
#' @keywords internal
.reshape_sim = function(x, variable = c("species", "resources")) {
	variable = match.arg(variable)
	if(variable == "species") {
		variable.name = "species"
		value.name = "occupancy"
	} else {
		variable.name = "resource"
		value.name = "concentration"
	}
	cores = ifelse(.Platform$OS.type == "unix", parallel::detectCores(), 1)
	data.table::rbindlist(parallel::mclapply(x[["networks"]], function(r) {
		S = state(r, variable, history = TRUE)
		res = data.table::rbindlist(lapply(S, function(s) {
			val = data.table::data.table(s)
			val$reach = seq_len(nrow(val))
			val
		}), idcol = "time")
		res = data.table::melt(res, id.vars = c("reach", "time"), variable.name = variable.name,
			value.name = value.name)
		if(variable.name == "species")
			res[[variable.name]] = sub("V(.+)", "\\1", res[[variable.name]])
		res
	}, mc.cores = cores), idcol = "network")
}


#' Run the simulation for a specified number of time steps
#'
#' This function begins with the current state of `x`, runs the model for `nt` steps,
#' then returns a copy of `x` with the updated state and state history.
#'
#' @param x A [flume()] object
#' @param nt The number of outer time steps, see 'details'
#' @param nt_i The number of inner/integration time steps, see 'details'
#' @param reps The number of replicate simulations to run; by default a single sim is run
#' on a new flume; for continuing simulations, the same number of replicates will be used
#' @details Flume runs at two different time scales to account for the fact that community changes
#' are often slower than resource transport. Thus we have two different parameters here (along
#' with the `dt` parameter from the [flume()] function) for specifying the length of time to run 
#' the model:
#'
#' * `nt` Controls the "outer" time step; this is the number of times the model will simulate
#' colonisations and extinctions, and the number of time steps that will be reported in model 
#' output.
#' * `nt_i` Defines the number of inner time steps for each outer time step; this controls the
#' temporal resolution at which the integration for the transport/reaction model is run.
#' * `dt` (Not specified here, but in the [flume()] function) Determines the length of each outer
#' time step; units are the same units as discharge and other fluxes.
#'
#' For example, if discharge is in m^3/s, then the default `dt = 86400` corresponds to an outer
#' time step of 86400 seconds, or one day. Setting nt = 100 will run the model for 100 days, and 
#' setting nt_i = 144 will integrate resource transport and reaction 144 times each day, 
#' corresponding to a 10-minute temporal resolution.
#' @return A modified copy of `x`, with state updated with the results of the simulation.
#' @export
run_simulation = function(x, nt, nt_i=1, reps, parallel = TRUE, cores = parallel::detectCores()) {

	if(nt < 1 || nt_i < 1)
		stop("at least one time step is required")

	if(missing(reps))
		reps = length(x[["networks"]])

	parallel = parallel && (.Platform$OS.type == "unix")

	# normally flumes are initialized with only a single river network
	# if more reps are desired, duplicate them
	if(reps > 1 && length(x[["networks"]]) == 1)
		x[["networks"]] = lapply(1:reps, function(i) x[["networks"]][[1]])

	if(reps != length(x[["networks"]]))
		warning("'reps' was specified, but it is not equal to the number of existing sims ",
			"in 'x'\n'reps' will be ignored and the number of simulations will be equal to ",
			"length(x$networks) (",length(x$networks), ")")

	if(parallel && reps > 1 && cores > 1) {
		x[["networks"]] = parallel::mclapply(x[["networks"]], .do_sim, comm = x[["metacom"]], 
			dt = x[["dt"]], nt = nt, nt_i = nt_i, mc.cores = cores)
	} else {
		x[["networks"]] = lapply(x[["networks"]], .do_sim, comm = x[["metacom"]], dt = x[["dt"]], 
			nt = nt, nt_i = nt_i)
	}

	return(x)
}


#' Worker function for running simulations in parallel
#' @param network A river network
#' @param comm A metacommunity
#' @param dt The size of the (outer) time step; one "generation" for the metacommunity
#' @param nt The number of outer time steps
#' @param nt_i The number of inner time steps
#' @return A modified copy of `network` with the results of the simulation
#' @keywords internal
.do_sim = function(network, comm, dt, nt, nt_i) {
	R = state(network, "resources")
	S = state(network, "species")

	dt_i = dt/nt_i
	for(t_out in 1:nt) {
		cp = col_prob(comm, network, dt = dt)
		ep = ext_prob(comm, network, dt = dt)

		# transport-reaction terms for each site/resource, in concentration units
		rxn_t = transport_t = 0 * R
		for(t_in in 1:nt_i) {
			dR_comp = dRdt(comm, network, components = TRUE)
			dR = dR_comp$reaction - dR_comp$transport
			rxn_t = rxn_t + dR_comp$reaction * dt_i
			transport_t = transport_t + dR_comp$transport * dt_i
			## euler integration for now
			R = R + dR * dt_i
		}

		# no cols where already present, no exts where already absent
		cp[S == 1] = 0
		ep[S == 0] = 0
		cols = matrix(rbinom(length(cp), 1, cp), nrow=nrow(cp))
		exts = matrix(rbinom(length(ep), 1, ep), nrow=nrow(ep))
		S[cols == 1] = 1
		S[exts == 1] = 0

		# currently we update state only at the end of each outer time step
		state(network, "species") = S
		state(network, "resources") = R
		state(network, "reaction") = rxn_t
		state(network, "transport") = transport_t
	}
	return(network)
}
