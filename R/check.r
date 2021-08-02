#' @name check_par
#' @title Check parameter dimensions
#' @rdname check_par
#' @keywords internal
.check_scale = function(p, nsp) {
	if(length(p) == 1)
		p = rep(p, nsp)
	if(length(p) != nsp)
		stop("Scale parameters must be missing, a single value, or 1 value per species")
	p
}

#' @rdname check_par
#' @keywords internal
.check_ratio = function(r, nr) {
	if(!is.matrix(r))
		r = matrix(r, nrow=1)
	if(any(r %% 1 != 0) )
		stop("ratio niches must be indicated by a matrix of integer indices")
	if(any(r > nr))
		stop("All entries in ratio must be <= nr")
	if(ncol(r) != 2)
		stop("ratio must be a 2-column matrix")
	r
}

#' @rdname check_par
#' @keywords internal
.check_breadth = function(p, nsp, nrx) {
	if(length(p) == 1)
		p = rep(p, nsp)

	if(!is.list(p) && length(p) == nsp) {
		p = matrix(p, nrow = nsp, ncol = nrx)
	} else if(!is.list(p) && length(p) == nrx) {
		p = matrix(p, nrow = nsp, ncol = nrx, byrow = TRUE)
	}

	## univariate niches
	if(nrx == 1) {
		if(!(is.matrix(p) && nrow(p) == nsp && ncol(p) == 1))
			stop("Invalid niche breadth: must be a scalar, vector, or one-column matrix")
		## convert back to a vector so it plays nicely with mapply
		p = p[, 1]
	} else {
		## for multivariate niches
		p = .check_breadth_multi(p, nsp, nrx)
	}
	p
}

#' @rdname check_par
#' @keywords internal
.check_breadth_multi = function(p, nsp, nrx) {
	if(is.matrix(p)) {
		if(nrow(p) != nsp || ncol(p) != nrx)
			stop("Invalid niche breadth; see ?niches")
		p = lapply(1:nsp, function(i) {
			mat = matrix(0, nrow = nrx, ncol = nrx)
			diag(mat) = p[i, ]
			mat
		})
	}

	## test proper dimensions of the breadth param for multivariate niches
	## must be a list, one entry per species
	## each entry a square matrix with dimension nrx
	if(!(is.list(p) &&
		length(p) == nsp &&
		all(sapply(p, is.matrix)) &&
		all((sapply(p, nrow) - nrx) == 0) &&
		all((sapply(p, ncol) - nrx) == 0))) {
		stop("Invalid niche breadth; see ?metacommunity")
	}
	p
}

#' @rdname check_par
#' @keywords internal
.check_r_use = function(p, nsp, nrx) {
	if(length(p) == 1)
		p = rep(p, nsp)

	if(!is.matrix(p) && length(p) == nsp) {
		p = matrix(p, nrow = nsp, ncol = nrx)
	} else if(!is.matrix(p) && length(p) == nrx) {
		p = matrix(p, nrow = nsp, ncol = nrx, byrow = TRUE)
	}

	if(!(is.matrix(p) && nrow(p) == nsp && ncol(p) == nrx))
		stop("r_use must be a matrix with one row per species and one column per niche axis")

	# convert to a list, one entry per row
	apply(p, 1, function(x) x, simplify = FALSE)
}

#' @rdname check_par
#' @keywords internal
.check_species_params = function(x) {
	loc = x$par_c$location
	bre = x$par_c$breadth
	csc = x$par_c$scale
	rsc = x$r_use
	esc = x$par_e$scale
	alpha = x$alpha
	beta = x$beta

	if(length(csc) != 1 | length(esc) != 1 | length(alpha) != 1 | length(beta) != 1)
		stop("The following parameters must be single values: ",
			"c_scale, e_scale, alpha, beta")

	if(any(bre < 0) | csc < 0 | esc < 0 | alpha < 0 | beta < 0)
		stop("The following parameters must not be negative: ",
			"breadth, scale_c, scale_e, alpha, beta")

	if(length(loc) == 1 & length(bre) != 1)
		stop("dimension mismatch: breadth and r_scale must be length 1 with a single niche",
				"dimension")
	if(length(loc) > 1 & (!is.matrix(bre) | !identical(dim(bre), rep(length(loc), 2))))
		stop("dimension mismatch: breadth must be a square matrix with one row/column per ",
			"niche dimension")
}
