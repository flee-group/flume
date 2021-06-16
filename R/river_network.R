#' Create a river network
#'
#' @details The adjacency matrix must be a square matrix, with a row/column for each river network
#' reach. Zero entries indicate no adjacency between reaches. Nonzero entries indicate adjacency,
#' with the row indicating the upstream reach and the column the downstream 
#' (i.e., adjacency[i,j] == 1 indicates water flows from i to j). Nonzero elements will be set
#' equal to one.
#' @param adjacency An adjacency matrix; see 'details'
#' @param discharge Optional discharge values, see [discharge.river_network()]
#' @param area Optional vector for cross sectional area values, one per reach
#' @param length The length of each reach, must be a single value (the model assumes all reaches
#' are the same length)
#' @param layout Optional, matrix of x-y coordinates used for plotting the network
#' @param site_names Optional vector of site names
#' @param skip_checks Logical; if true, no checks for valid topology will be performed
#' @return An S3 object of class 'river_network', with the following members:
#' * `adjacency` The adjacency matrix
#' * `discharge` A discharge vector
#' * `area` Cross sectional area at each node
#' * `.state` The state history of the network; access with [state][state.river_network()].
#'
#' Additionally, the `layout` attribute, if used, is a matrix of coordinates for plotting
#' @examples
#' Q = rep(1, 4)
#' adj = matrix(0, nrow = 4, ncol = 4)
#' adj[1,2] = adj[2,3] = adj[4,3] = 1
#' rn = river_network(adj, Q)
#' plot(rn)
#' @export
river_network = function(adjacency, discharge, area, length = 1, layout, site_names,
	skip_checks = FALSE) {

	if(any(! as(adjacency, "vector") %in% c(0,1))) {
		adjacency[adjacency != 0] = 1
	}

	if(!skip_checks) {
		if(!requireNamespace("igraph", quietly = TRUE) || 
					!requireNamespace("Matrix", quietly = TRUE)) {
			stop("Packages igraph and Matrix are required for topology checks; please install to ",
				"proceed. If no topology checks are required, disable this error with ",
				"skip_checks = TRUE.")
		}
		if(!.validate_adjacency(adjacency))
			stop("Adjacency matrix failed validation; see '?river_network' for the requirements.")
	}

	if(missing(site_names))
		site_names = paste0("n", 1:nrow(adjacency))

	rn = structure(list(.adjacency = adjacency, .state = list(), .length = length),
		class = "river_network")

	if(!missing(discharge))
		discharge(rn) = discharge
	if(!missing(area))
		cs_area(rn) = area

	if(!missing(layout))
		attr(rn, "layout") = layout
	attr(rn, "n_sites") = nrow(adjacency)
	attr(rn, "site_names") = site_names

	return(rn)
}

#' River network adjacency matrix
#'
#' @param x A river network
#' @return A matrix Y, such that non-zero values in Y[i,j] indicate flow from i to j
#' @export
adjacency = function(x) {
	return(x[['.adjacency']])
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
	if(!(is.matrix(adjacency) || is(adjacency, "Matrix")) ||   ## must be standard matrix/Matrix
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



#' @export
state = function(x, ...)
	UseMethod("state", x)

#' @export
'state<-' = function(x, value)
	UseMethod('state<-', x)

#' @export
boundary = function(x, ...)
	UseMethod("boundary", x)

#' @export
'boundary<-' = function(x, value)
	UseMethod('boundary<-', x)

#' @export
site_by_species = function(x, ...)
	UseMethod("site_by_species", x)

#' @export
'site_by_species<-' = function(x, value)
	UseMethod('site_by_species<-', x)

#' @export
boundary_species = function(x, ...)
	UseMethod("boundary_species", x)

#' @export
'boundary_species<-' = function(x, value)
	UseMethod('boundary_species<-', x)


#' @rdname state
#' @export
reset_state = function(x, value) {
	if(!is.matrix(value))
		value = matrix(value, nrow=nrow(adjacency(x)), ncol=1)
	attr(x, "r_names") = colnames(value)
	x[['.state']] = list()
	state(x) = value
	x
}

#' @rdname state
#' @export
reset_species = function(x, value) {
	if(!is.matrix(value))
		value = matrix(value, nrow=nrow(adjacency(x)), ncol=1)
	attr(x, "sp_names") = colnames(value)
	x[['si_by_sp']] = list()
	site_by_species(x) = value
	x
}



#' Setter and getter methods for river network state variables
#' @name state
#' @rdname state
#' @title Get/set the (resource) state of a river network
#' @details By default, setting state will save the current state in the state history, then update
#' current state to `value`. `reset_state` erases the state variable and sets a new one, and does
#' not check that the dimensions make sense. Other methods update state and do dimension checking.
#'
#' `boundary()` returns the boundary condition of the river network, as a site by resource matrix.
#' `boundary() <-` sets the boundary; it must be a site by resoure matrix in the same format as
#' state.
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
boundary.river_network = function(x) {
	return(x[['.boundary']])
}

#' Check state matrix for errors
#' @param x The state matrix to check
#' @param nn The number of nodes in the network
#' @param nr The number of resources, if NA, this is not checked
#' @keywords internal
.check_state = function(x, nn, nr = NA) {
	if(!is.matrix(x))
		stop("Can only assign a matrix to state")
	if(nrow(x) != nn)
		stop("State matrix must have one row per reach in the network")
	if(!is.na(nr) && ncol(x) != nr)
		stop("State matrix must have one column per resource")
}

#' @rdname state
#' @export
'state<-.river_network' = function(x, value) {
	i = length(x[['.state']])
	.check_state(value, nrow(adjacency(x)), ifelse(i == 0, NA, ncol(state(x))))
	rownames(value) = attr(x, "site_names")
	colnames(value) = attr(x, "r_names")
	x[['.state']][[i + 1]] = value
	return(x)
}

#' @rdname state
#' @export
'boundary<-.river_network' = function(x, value) {
	.check_state(value, nrow(adjacency(x)), ncol(state(x)))
	rownames(value) = attr(x, "site_names")
	colnames(value) = attr(x, "r_names")
	x[['.boundary']] = value
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
boundary_species.river_network = function(x) {
	return(x[['.boundary_sp']])
}

#' @rdname state
#' @export
'boundary_species<-.river_network' = function(x, value) {
	.check_state(value, nrow(adjacency(x)), ncol(site_by_species(x)))
	rownames(value) = attr(x, "site_names")
	colnames(value) = attr(x, "sp_names")
	x[['.boundary_sp']] = value
	return(x)
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
	rownames(value) = attr(x, "site_names")
	colnames(value) = attr(x, "sp_names")
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


#' Number of nodes in a river network
#' @param x A river network
#' @return The number of nodes in the network
#' @export
length.river_network = function(x) {
	return(nrow(adjacency(x)))
}
