#' Create a river network
#'
#' @details The adjacency matrix must be a square matrix, with a row/column for each river network
#' reach. Zero entries indicate no adjacency between reaches. Nonzero entries indicate adjacency,
#' with the row indicating the upstream reach and the column the downstream 
#' (i.e., adjacency[i,j] == 1 indicates water flows from i to j). Nonzero elements will be set
#' equal to one.
#' @param adjacency An adjacency matrix; see 'details'
#' @param discharge Optional discharge values, see [discharge.river_network()]
#' @param area Optional cross sectional area values, see [discharge.river_network()]
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
river_network = function(adjacency, discharge, area, length = 1, layout, 
	site_names = paste0("si", 1:nrow(adjacency)), skip_checks = FALSE) {

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

	rn = structure(list(.adjacency = adjacency, .length = length), class = "river_network")

	if(!missing(discharge))
		discharge(rn) = discharge
	if(!missing(area))
		area(rn) = area

	if(!missing(layout))
		attr(rn, "layout") = layout
	attr(rn, "n_sites") = nrow(adjacency)
	attr(rn, "names_sites") = site_names

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
'state<-' = function(x, ..., value)
	UseMethod('state<-', x)

#' @export
boundary = function(x, ...)
	UseMethod("boundary", x)

#' @export
'boundary<-' = function(x, ..., value)
	UseMethod('boundary<-', x)



#' Setter and getter methods for river network state variables
#' @name state
#' @rdname state
#' @title Get/set the (resource) state of a river network
#' @details By default, setting state will save the current state in the state history, then update
#' current state to `value`. Setting state to NULL erases the state variable.
#'
#' `boundary()` returns the boundary condition of the river network, as a site by resources
#' or site by species matrix.
#' `boundary() <-` sets the boundary; it must be a site by var matrix in the same format as
#' state.
#' @param x A river network
#' @param var The state variable to extract/modify; must be one of `resources`, `
#' @param history Logical; if TRUE, entire state history is returned
#' @param value The value to update the variable with
#' @export
state.river_network = function(x, var, history = FALSE) {
	if(missing(var)) {
		warning("Calling state(x) with no variable name is deprecated; use state(x, 'resources')")
		var = "resources"
	}
	allowed_vars = c("resources", "species", "Q", "area", "reaction", "transport")
	if(length(var) != 1 || !var %in% allowed_vars)
		stop("var must be one of the following: ", paste(allowed_vars, collapse = ' '))
	var = paste0('.', var)

	if(var == ".Q") {
		val = discharge(x, type = ifelse(history, "history", "current"))
	} else if(var == ".area") {
		val = area(x, type = ifelse(history, "history", "current"))
	} else {
		if(length(x[[var]]) == 0)
			return(NULL)
		
		if(history) {
			val = x[[var]]
		} else {
			val = x[[var]][[length(x[[var]])]]
		}
	}
	val
}


#' @rdname state
#' @export
'state<-.river_network' = function(x, var, value) {
	if(missing(var)) {
		warning("Calling state(x) with no variable name is deprecated; use state(x, 'resources')")
		var = "resources"
	}
	if(var == "Q")
		stop("Setting Q with state<- is not permitted, use 'discharge(x) <- value' instead.")
	if(var == "area")
		stop("Setting area with state<- is not permitted, use 'area(x) <- value' instead.")
	allowed_vars = c("resources", "species", "reaction", "transport")
	if(length(var) != 1 || !var %in% allowed_vars)
		stop("var must be one of the following:", paste(allowed_vars))
	cnames = paste0("names_", var)
	var_st = paste0('.', var)

	if(is.null(value)) {
		x[[var_st]] = list()
		return(x)
	}

	i = length(x[[var_st]])
	if(i == 0) {
		x[[var_st]] = list()
		if(!is.matrix(value))
			value = matrix(value, nrow=nrow(adjacency(x)), ncol=1)
		attr(x, cnames) = colnames(value)
	} else {
		.check_state(value, nrow(adjacency(x)), ncol(state(x, var)))
	}

	rownames(value) = attr(x, "names_sites")
	colnames(value) = attr(x, cnames)
	x[[var_st]][[i + 1]] = value
	return(x)
}


#' @rdname state
#' @export
boundary.river_network = function(x, var) {
	if(missing(var)) {
		warning("Calling boundary(x) with no variable name is deprecated;",
			"use boundary(x, 'resources')")
		var = "resources"
	}
	allowed_vars = c("resources", "species", "Q")
	if(length(var) != 1 || !var %in% allowed_vars)
		stop("var must be one of the following:", paste(allowed_vars))

	if(var == "Q") {
		val = lateral_discharge(x)
	} else {
		var = paste0('.b_', var)
		val = x[[var]]
	}
	return(val)
}

#' @rdname state
#' @export
'boundary<-.river_network' = function(x, var, value) {
	if(missing(var)) {
		warning("Calling boundary(x) with no variable name is deprecated;",
			"use boundary(x, 'resources')")
		var = "resources"
	}
	allowed_vars = c("resources", "species")
	if(length(var) != 1 || !var %in% allowed_vars)
		stop("var must be one of the following:", paste(allowed_vars))
	cnames = paste0("names_", var)
	var_b = paste0('.b_', var)

	.check_state(value, nrow(adjacency(x)), ncol(state(x, var)))
	rownames(value) = attr(x, "names_sites")
	colnames(value) = attr(x, cnames)
	x[[var_b]] = value
	return(x)
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
		stop("State matrix must have one column per state dimension")
}


#' Compute current prevalence for all species in a river network
#' @param x A river network
#' @return A site by species matrix; values indicate the prevalence in surrounding sites
prevalence = function(x) {
	## for prevalence, adjacency is upstream, downstream or self
	adj = adjacency(x) + t(adjacency(x))
	diag(adj) = 1

	adj %*% state(x, "species")
}


#' Number of nodes in a river network
#' @param x A river network
#' @return The number of nodes in the network
#' @export
length.river_network = function(x) {
	return(nrow(adjacency(x)))
}








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
boundary_species.river_network = function(x) {
	.Deprecated("boundary(x, 'species')")
	return(x[['.boundary_sp']])
}

#' @rdname state
#' @export
'boundary_species<-.river_network' = function(x, value) {
	.Deprecated("boundary(x, 'species') <- value")
	.check_state(value, nrow(adjacency(x)), ncol(site_by_species(x)))
	rownames(value) = attr(x, "site_names")
	colnames(value) = attr(x, "sp_names")
	x[['.boundary_sp']] = value
	return(x)
}

#' @rdname state
#' @export
'site_by_species<-.river_network' = function(x, value) {
	.Deprecated("state(x, 'species') <- value")
	if(!".species" %in% names(x))
		x[['.species']] = list()

	i = length(x[['.species']])
	if(i > 0) {
		## changing state dimensions is not allowed
		stopifnot(identical(dim(value), dim(site_by_species(x))))
	}
	rownames(value) = attr(x, "site_names")
	colnames(value) = attr(x, "sp_names")
	x[['.species']][[i + 1]] = value
	return(x)
}


#' @rdname state
#' @export
site_by_species.river_network = function(x, history = FALSE) {
	.Deprecated("state(x, 'species')")
	if(!".species" %in% names(x))
		stop("Site by species information is missing for this river network")
	if(history) {
		return(x[[".species"]])
	} else {
		return(x[['.species']][[length(x[['.species']])]])
	}
}

#' @rdname state
#' @export
reset_state = function(x, value) {
	.Deprecated("state(x, 'resources') <- NULL; state(x, 'resources') <- value")
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
	.Deprecated("state(x, 'species') <- NULL; state(x, 'species') <- value")
	if(!is.matrix(value))
		value = matrix(value, nrow=nrow(adjacency(x)), ncol=1)
	attr(x, "sp_names") = colnames(value)
	x[['.species']] = list()
	site_by_species(x) = value
	x
}

