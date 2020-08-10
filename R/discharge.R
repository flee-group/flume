#' @export
discharge = function(x, ...)
	UseMethod("discharge", x)

#' @export
'discharge<-' = function(x, value)
	UseMethod('discharge<-', x)

#' @export
cs_area = function(x, ...)
	UseMethod("cs_area", x)

#' @export
'cs_area<-' = function(x, value)
	UseMethod('cs_area<-', x)


#' Get and set discharge and cross-sectional area of a river network
#'
#' @details There are three possible models for discharge/cross-sectional area, depending on the
#' kind of input provided.
#'
#' If discharge is provided as a vector or a single-column matrix, then discharge will be treated as constant
#' throughout any simulations. In this case, then length of the input must equal the number of nodes in the
#' river network.
#'
#' If discharge is provided as a matrix with ncol > 1, then a variable discharge model will be used. The input in this
#' case is a hydrograph, where the rows are nodes in the river network and the columns are time steps. Time
#' steps will be recycled, so if a model is run for a more steps than is present in the discharge matrix, then discharge
#' will restart at the first columns of the matrix.
#'
#' The final possibility is to provide a function that takes a single parameter, the time step, as input and returns
#' a discharge vector with one element per node in the river network.
#'
#' Cross sectional area can either be missing or provided in the same format as discharge. If missing, then cross
#' sectional area will be computed using hydrological scaling relationships; see [geometry()] for details.
#'
#' @param x A [river_network()]
#' @param value Discharge/cross-sectional area (see 'details')
#' @name discharge
#' @return A discharge vector
#' @export
discharge.river_network = function(x) {
	Q = .get_hydro_attr(x, ".Q")
}

#' @rdname discharge
#' @export
'discharge<-.river_network' = function(x, value) {
	if(is.vector(value))
		value = matrix(value, ncol=1)

	if(is.function(value)) {
		attr(x, "discharge_model") = "function"
	} else if(is.matrix(value)) {
		if(any(value < 0))
			stop("Negative discharge is not allowed")
		if(!nrow(value) == length(x))
			stop("nrow(discharge) must equal length(x)")
		if(ncol(value) > 1) {
			attr(x, "discharge_model") = "variable"
		} else {
			attr(x, "discharge_model") = "constant"
		}
	} else {
		stop("Unknown input for discharge")
	}
	x[['.Q']] = value
	return(x)
}

#' @rdname discharge
#' @return Cross-sectional area vector
#' @export
cs_area.river_network = function(x) {
	if('.area' %in% ls(x)) {
		A = .get_hydro_attr(x, '.area')
	} else {
		A = geometry(discharge(x))
		A = A$width * A$depth
	}
	return(A)
}

#' @rdname discharge
#' @export
'cs_area<-.river_network' = function(x, value) {
	if(is.vector(value))
		value = matrix(value, ncol = 1)

	if(!is.function(value) && any(value < 0))
		stop("Negative areas not allowed")

	if(attr(x, "discharge_model") == "constant") {
		if(!is.matrix(value) || ! ncol(value) == 1 || !nrow(value) == length(x))
			stop("Areas must be a 1-column matrix or a vector with length == length(x)")
	} else if(attr(x, "discharge_model") == "variable") {
		if(! is.matrix(value) || !nrow(value) == length(x))
			stop("Areas must be a matrix with nrow == length(x)")
	} else {
		if(! is.function(value))
			stop("Areas must be a function")
	}
	x[['.area']] = value
	return(x)
}

#' Compute hydraulic geometry from discharge
#' @param Q A vector of discharge values.
#' @references Raymond, PA et al. 2012. Scaling the gas transfer velocity and hydraulic
#' 		geometry in streams and small rivers. *Limnology and Oceanography: Fluids and
#' 		Environments*. **2**:41-53.
#' @return A data frame with discharge, velocity, depth, and width
#' @export
geometry = function(Q) {
	va <- -1.64
	vb <- 0.285
	da <- -0.895
	db <- 0.294
	wa <- 2.56
	wb <- 0.423

	velocity = exp(va + vb * log(discharge))
	depth = exp(da + db * log(discharge))
	width = exp(wa + wb * log(discharge))
	ind = which(discharge == 0)
	velocity[ind] = depth[ind] = width[ind] = 0
	data.frame(discharge, velocity, depth, width)
}

#' Pull out discharge/cross sectional area
#' @param x A river network
#' @param attr Which attribute to get
#' @keywords internal
.get_hydro_attr = function(x, attr = c('.Q', '.area')) {
	attr = match.arg(attr)
	tm = length(x[['.state']]) ## determine what time step we are at from how many states we have saved
	if(attr(x, "discharge_model") == "constant") {
		val = x[[attr]]
	} else if(attr(x, "discharge_model") == "variable") {
		if(tm > ncol(x[[attr]])) {
			tm = tm %% ncol(x[[attr]])
			if(tm == 0)
				tm = ncol(x[[attr]])
		}
		val = x[[attr]][,tm]
	} else {
		val = x[[attr]](tm)
	}
	return(val)
}


