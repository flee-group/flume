#' Create a river network
#'
#' @details The adjacency matrix must be a square matrix, with a row/column for each river network reach. Zero entries
#' indicate no adjacency between reaches. Nonzero entries indicate adjacency, with the value equal to the distance between
#' reaches.
#' @param adjacency An adjacency matrix; see 'details'
#' @param discharge A vector for discharge values, one per row/column in adjacency
#' @param state Optional, starting state of the network
#' @param skip_checks Logical; if true, no checks for valid topology will be performed
#' @return An S3 object of class 'river_network', with the following attributes:
#' * `adjacency` The adjacency matrix
#' * `discharge` A discharge vector
#' * `.state` The state history of the network; access with [state][state.river_network()].
#' @examples
#' Q = rep(1, 4)
#' adj = matrix(0, nrow = 4, ncol = 4)
#' adj[1,2] = adj[2,3] = adj[4,3] = 1
#' rn = river_network(adj, Q)
#' plot(rn)
#' @export
river_network = function(adjacency, discharge, state, skip_checks = FALSE) {
	if(!skip_checks) {
		if(!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
			stop("Packages igraph and Matrix are required for topology checks; please install to proceed. ",
				"If no topology checks are required, disable this error with skip_checks = TRUE.")
		}
		if(!.validate_adjacency(adjacency))
			stop("Adjacency matrix failed validation; see '?river_network' for the requirements.")
		stopifnot(length(discharge) == nrow(adjacency))
	}

	rn = structure(list(adjacency = adjacency, discharge = discharge,
						.state = list()), class = "river_network")
	if(!missing(state))
		state(rn) = state
	return(rn)
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
	if(nrow(value) != nrow(x[['adjacency']]))
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
