#' Create a river network
#'
#' @details The adjacency matrix must be a square matrix, with a row/column for each river network reach. Zero entries
#' indicate no adjacency between reaches. Nonzero entries indicate adjacency, with the value equal to the distance between
#' reaches.
#' @param adjacency An adjacency matrix; see 'details'
#' @param discharge A vector for discharge values, one per row/column in adjacency
#' @param skip_checks Logical; if true, no checks for valid topology will be performed
river_network = function(adjacency, discharge, skip_checks = FALSE) {
	if(!skip_checks) {
		if(!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("Matrix", quietly = TRUE)) {
			stop("Packages igraph and Matrix are required for topology checks; please install to proceed. ",
				"If no topology checks are required, disable this error with skip_checks = TRUE.")
		}
		if(!.validate_adjacency(adjacency))
			stop("Adjacency matrix failed validation; see '?river_network' for the requirements.")
		stopifnot(length(discharge) == nrow(adjacency))
	}

	rn = structure(list(adjacency = adjacency, discharge = discharge), class = "river_network")
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
