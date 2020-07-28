#' Create a river network
#'
#' @details The adjacency matrix must be a square matrix, with a row/column for each river network reach. Zero entries
#' indicate no adjacency between reaches. Nonzero entries indicate adjacency, with the row indicating the upstream
#' reach and the column the downstream (i.e., adjacency[i,j] == 1 indicates water flows from i to j). Nonzero elements
#' will be set equal to one.
#' @param adjacency An adjacency matrix; see 'details'
#' @param discharge A vector for discharge values, one per row/column in adjacency
#' @param area Optional vector for cross sectional area values, one per row/column in adjacency
#' @param state Optional, starting state of the network
#' @param length The length of each reach, must be a single value (the model assumes all reaches are the same length)
#' @param skip_checks Logical; if true, no checks for valid topology will be performed
#' @return An S3 object of class 'river_network', with the following attributes:
#' * `adjacency` The adjacency matrix
#' * `discharge` A discharge vector
#' * `area` Cross sectional area at each node
#' * `.state` The state history of the network; access with [state][state.river_network()].
#' @examples
#' Q = rep(1, 4)
#' adj = matrix(0, nrow = 4, ncol = 4)
#' adj[1,2] = adj[2,3] = adj[4,3] = 1
#' rn = river_network(adj, Q)
#' plot(rn)
#' @export
river_network = function(adjacency, discharge, area, state, length = 1, skip_checks = FALSE) {
	if(!skip_checks) {
		if(!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
			stop("Packages igraph and Matrix are required for topology checks; please install to proceed. ",
				"If no topology checks are required, disable this error with skip_checks = TRUE.")
		}
		if(!.validate_adjacency(adjacency))
			stop("Adjacency matrix failed validation; see '?river_network' for the requirements.")
		stopifnot(length(discharge) == nrow(adjacency))
	}
	if(any(! adjacency %in% c(0,1))) {
		adjacency[adjacency != 0] = 1
	}

	if(missing(area)) {
		if(!requireNamespace("WatershedTools", quietly = TRUE))
			stop("WatershedTools is required to estimate area from discharge; ",
			"use devtools::install_github('mtalluto/WatershedTools') to install")

		area = WatershedTools::hydraulic_geometry(discharge)
		area = area$depth * area$width
	}

	rn = structure(list(.adjacency = adjacency, discharge = discharge, area = area,
						.state = list(), .length = length), class = "river_network")

	## by default, boundary condition is set to equal the starting state
	if(!missing(state)) {
		state(rn) = state
		rn$boundary = function() return(state)
	} else {
		rn$boundary = function() return(rep(0, nrow(rn[['.adjacency']])))
	}
	return(rn)
}

#' River network adjacency matrix
#'
#' @param x A river network
#' @param weighted boolean, default to FALSE, if true returns the weighted adjacency matrix
#' @return A matrix Y, such that non-zero values in Y[i,j] indicate flow from i to j
#' @export
adjacency = function(x, weighted = FALSE) {
	adj = x[['.adjacency']]
	# if(!weighted)  ## removed weights from adjacency matrix
	# 	adj = (adj > 0) * 1
	return(adj)
}


#' Length of every reach in the river network
#'
#' Currently, all reaches are assumed to have the same length
#' @param x A river network
#' @return a vector of reach lengths
#' @export
reach_length = function(x) {
	return(rep(x[['.length']], nrow(adjacency(x))))
}

#' Perform checks on a river network adjacency matrix
#' @param adjacency An adjacency matrix
#' @importFrom Matrix rowSums colSums
#' @keywords internal
#' @return Logical, TRUE if the validation passes
.validate_adjacency = function(adjacency) {
	if(!(is.matrix(adjacency) || is(adjacency, "Matrix")) ||   ## must be standard or Matrix package matrix
	   any(adjacency < 0) ||                 ## negative distances not allowed
	   any(is.infinite(adjacency)) ||        ## nor infinite/NA
	   any(rowSums(adjacency) > 1) ||        ## a node may be upstream of at most one other node
	   any(rowSums(adjacency) + colSums(adjacency) == 0) ## nodes must not be isolated
	) return(FALSE)

	plt = igraph::graph_from_adjacency_matrix(adjacency, mode = "directed")
	if(!igraph::is_dag(plt))
		return(FALSE)

	return(TRUE)
}


#'@export
state = function(x, ...)
	UseMethod("state", x)
#'@export
'state<-' = function(x, value)
	UseMethod('state<-', x)
#'@export
site_by_species = function(x, ...)
	UseMethod("site_by_species", x)
#'@export
'site_by_species<-' = function(x, value)
	UseMethod('site_by_species<-', x)

#' Setter and getter methods for river network state variables
#' @name state
#' @details By default, setting state will save the current state in the state history, then update current state
#' to `value`.
#' @param x A river network
#' @param history Logical; if TRUE, entire state history is returned
#' @param value The value to update the attribute with
#' @export
state.river_network = function(x, history = FALSE) {
	if(length(x[['.state']]) == 0)
		return(NULL)

	if(history) {
		return(x[['.state']])
	} else {
		return(x[['.state']][[length(x[['.state']])]])
	}
}

#' @rdname state
#' @export
'state<-.river_network' = function(x, value) {
	if(nrow(value) != nrow(adjacency(x)))
		stop("nrow(value) must be equal to the number of nodes in the river network")

	i = length(x[['.state']])
	if(i > 0) {
		## changing state dimensions is not allowed
		stopifnot(ncol(value) == ncol(state(x)))
	}

	x[['.state']][[i + 1]] = value
	return(x)
}

#' @rdname state
#' @export
site_by_species.river_network = function(x, history = FALSE) {
	if(!"si_by_sp" %in% names(x))
		stop("Site by species information is missing for this river network")
	if(history) {
		return(x[["si_by_sp"]])
	} else {
		return(x[['si_by_sp']][[length(x[['si_by_sp']])]])
	}
}

#' @rdname state
#' @export
'site_by_species<-.river_network' = function(x, value) {
	if(!"si_by_sp" %in% names(x))
		x[['si_by_sp']] = list()

	i = length(x[['si_by_sp']])
	if(i > 0) {
		## changing state dimensions is not allowed
		stopifnot(identical(dim(value), dim(site_by_species(x))))
	}
	x[['si_by_sp']][[i + 1]] = value
	return(x)
}

#' Compute current prevalence for all species in a river network
#' @param x A river network
#' @return A site by species matrix; values indicate the prevalence in surrounding sites
prevalence = function(x) {
	## for prevalence, adjacency is upstream, downstream or self
	adj = adjacency(x) + t(adjacency(x))
	diag(adj) = 1

	adj %*% site_by_species(x)
}

