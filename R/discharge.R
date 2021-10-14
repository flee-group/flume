#' @export
discharge = function(x, ...)
	UseMethod("discharge", x)

#' @export
'discharge<-' = function(x, value)
	UseMethod('discharge<-', x)

#' @export
area = function(x, ...)
	UseMethod("area", x)

#' @export
'area<-' = function(x, value)
	UseMethod('area<-', x)


#' Get and set discharge and cross-sectional area of a river network
#'
#' @details There are three possible models for discharge/cross-sectional area, depending on the
#' kind of input provided.
#'
#' If discharge is provided as a vector or a single-column matrix, then discharge will be treated 
#' as constant throughout any simulations. In this case, then length of the input must equal the 
#' number of nodes in the river network.
#'
#' If discharge is provided as a matrix with ncol > 1, then a variable discharge model will be used.
#' The input in this case is a hydrograph, where the rows are nodes in the river network and the 
#' columns are time steps. Time steps will be recycled, so if a model is run for a more steps than
#' is present in the discharge matrix, then discharge will restart at the first columns of the 
#' matrix.
#'
#' The final possibility is to provide a function that takes a single parameter, the time step, as 
#' input and returns a discharge vector with one element per node in the river network.
#'
#' Cross sectional area can either be missing or provided in the same format as discharge. If 
#' missing, then cross sectional area will be computed using hydrological scaling relationships; 
#' see [geometry()] for details.
#'
#' Note that if cross-sectional area is provided, then the discharge model is changed (e.g., from 
#' constant to variable), it is strongly recommended (and not checked!) that the user ALSO update
#' cross-sectional area to match the new model, otherwise undefined behavior will occur.
#'
#' `area(x) <- NULL` will delete any provided cross sectional area and revert to the default 
#' behavior, computing area using [geometry()]
#'
#' @param x A [river_network()]
#' @param type The type of record to get; either the raw data, the current state, or the whole 
#'		history
#' @param value Discharge/cross-sectional area (see 'details')
#' @name discharge
#' @return A discharge vector
#' @export
discharge.river_network = function(x, type = c("raw", "current", "history")) {
	type = match.arg(type)
	if(type == "raw") {
		val = x[[".Q"]]
	} else {
		val = .get_hydro_attr(x, ".Q", type)
	}
	val
}

#' @rdname discharge
#' @export
'discharge<-.river_network' = function(x, value) {
	if(is.vector(value) && is.atomic(value))
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
area.river_network = function(x, type = c("raw", "current", "history")) {
	type = match.arg(type)
	missing_area = !('.area' %in% ls(x, all.names=TRUE))

	if(type == "raw") {
		val = x[[".area"]]
	} else if(missing_area) {
		if(type == "current") {
			val = .compute_area(state(x, "Q"))
		} else {
			val = apply(state(x, "Q", history = TRUE), 2, .compute_area)
		}
	} else {
		val = .get_hydro_attr(x, '.area', type)
	}
	val
}

.compute_area = function(Q) {
	geom = geometry(Q)
	geom$depth * geom$width
}

#' @rdname discharge
#' @export
'area<-.river_network' = function(x, value) {
	if(is.null(value)) {
		x[['.area']] = NULL
		return(x)
	}

	if(is.vector(value) && is.atomic(value))
		value = matrix(value, ncol = 1)

	if(!is.function(value) && any(value < 0))
		stop("Negative areas not allowed")

	if(attr(x, "discharge_model") == "constant") {
		if(!is.matrix(value) || ! ncol(value) == 1 || !nrow(value) == length(x))
			stop("Areas must be a 1-column matrix or a vector with length == length(x)")
	} else if(attr(x, "discharge_model") == "variable") {
		if(!(is.matrix(value) && identical(dim(value), dim(x[[".Q"]]))))
			stop("Areas must be a matrix with dimensions identical to dim(discharge(x))")
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

	velocity = exp(va + vb * log(Q))
	depth = exp(da + db * log(Q))
	width = exp(wa + wb * log(Q))
	ind = which(Q == 0)
	velocity[ind] = depth[ind] = width[ind] = 0
	data.frame(Q, velocity, depth, width)
}

#' Pull out discharge/cross sectional area
#' @param x A river network
#' @param attr Which attribute to get
#' @keywords internal
.get_hydro_attr = function(x, attr = c('.Q', '.area'), type) {
	attr = match.arg(attr)
	## determine what time step we are at from how many states we have saved
	tm = length(state(x, "resources", history = TRUE))
	tm = ifelse(tm == 0, 1, tm)
	mod = attr(x, "discharge_model")

	if(mod == "constant") {
		val = x[[attr]][,1]
		if(type == "history") {
			val = matrix(val, ncol = tm)
		}
	} else if(mod == "variable") {
		val = x[[attr]]
		nc = ncol(val)
		if(type == "current") {
			i = ifelse(tm <= nc, tm, tm %% nc)
			i = ifelse(i == 0, nc, i)
		} else {
			i = if(tm <= nc) seq_len(tm) else c(rep(seq_lem(nc), tm %/% nc), seq_len(tm %% nc))
		}
		val = val[,i]
	} else {
		fun = x[[attr]]
		if(type == "current") {
			val = fun(tm)
		} else {
			val = sapply(1:tm, fun)
		}
	}
	return(val)
}


#' Compute lateral discharge
#'
#' Assumed to be the difference between Q of a site and the sum of Q of upstream sites
#' @param x A river network
#' @param tm A time step
#' @return A vector of discharge values
#' @export
lateral_discharge = function(x) {
	Q = state(x, "Q")
	Qu = t(adjacency(x)) %*% Q
	as.vector(Q - Qu)
}
