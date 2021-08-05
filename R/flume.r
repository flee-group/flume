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
flume = function(comm, network, sp0, st0, spb, stb, dt = 1) {

	## define initial states, and check dimensionality
	if(is.null(state(network)) && missing(st0))
		stop("Initial resource state not specified")

	if(is.null(network[["si_by_sp"]]) && missing(sp0))
		stop("Initial site by species matrix not specified")

	if(missing(st0))
		st0 = state(network)
	if(ncol(st0) != attr(comm, "n_r"))
		stop("Initial network state doesn't match number of resources in comm")
	colnames(st0) = attr(comm, "r_names")
	network = reset_state(network, st0)
	if(missing(stb))
		stb = state(network)

	i_static = which(attr(comm, "r_types") == "static")
	if(length(i_static) > 0) {
		stb[i_static] = 0
	}

	boundary(network) = stb
	if(missing(sp0))
		sp0 = site_by_species(network)
	if(ncol(sp0) != attr(comm, "n_species"))
		stop("Initial site by species doesn't match number of resources in comm")
	colnames(sp0) = attr(comm, "sp_names")
	network = reset_species(network, sp0)
	if(missing(spb))
		spb = matrix(0, nrow = attr(network, 'n_sites'), ncol = attr(comm, 'n_species'), 
			dimnames = list(attr(network, "site_names"), attr(comm, "sp_names")))
	boundary_species(network) = spb
	x = structure(list(metacom = comm, networks = list(network), dt = dt), class = "flume")
	x
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
		fun = site_by_species
		variable.name = "species"
		value.name = "occupancy"
	} else {
		fun = state
		variable.name = "resource"
		value.name = "concentration"
	}
	cores = ifelse(.Platform$OS.type == "unix", parallel::detectCores(), 1)
	data.table::rbindlist(parallel::mclapply(x[["networks"]], function(r) {
		S = fun(r, history = TRUE)
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
#' @param nt The number of time steps
#' @param reps The number of replicate simulations to run; by default a single sim is run
#' on a new flume; for continuing simulations, the same number of replicates will be used
#' @return A modified copy of `x`, with state updated with the results of the simulation.
#' @export
run_simulation = function(x, nt, reps, parallel = TRUE, cores = parallel::detectCores()) {

	if(nt < 1)
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
			dt = x[["dt"]], nt = nt, mc.cores = cores)
	} else {
		x[["networks"]] = lapply(x[["networks"]], .do_sim, comm = x[["metacom"]], dt = x[["dt"]], 
			nt = nt)
	}

	return(x)
}


#' Worker function for running simulations in parallel
#' @param network A river network
#' @param comm A metacommunity
#' @param dt The size of the time step
#' @param nt The number of time steps
#' @return A modified copy of `network` with the results of the simulation
#' @keywords internal
.do_sim = function(network, comm, dt, nt) {
	R = state(network)
	S = site_by_species(network)

	for(tstep in 1:nt) {
		cp = col_prob(comm, network, dt = dt)
		ep = ext_prob(comm, network, dt = dt)
		dR = dRdt(comm, network)

		# no cols where already present, no exts where already absent
		cp[S == 1] = 0
		ep[S == 0] = 0
		cols = matrix(rbinom(length(cp), 1, cp), nrow=nrow(cp))
		exts = matrix(rbinom(length(ep), 1, ep), nrow=nrow(ep))
		S[cols == 1] = 1
		S[exts == 1] = 0

		## euler integration for now
		R = R + dR * dt

		site_by_species(network) = S
		state(network) = R
	}
	return(network)
}
