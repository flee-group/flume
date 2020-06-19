#' Plot a river network
#' @param x A [river_network()]
#' @param ... Additional arguments to [igraph::plot.igraph()]
#' @export
plot.river_network = function(x, ...) {
	if(!requireNamespace("igraph", quietly = TRUE)) {
		stop("Package 'igraph' is required for plotting, please install it and try again")
	}

	args = .default_river_plot_options(...)
	args$x = graph_from_adjacency_matrix(x$adjacency, mode = "directed")
	wt = x$discharge[2:nrow(x$adjacency)]
	args$edge.width = (wt / max(wt)) * args$edge.width
	do.call(plot, args)
}

#' Plot species independent stable envelopes
#' @param x A [metacommunity()]
#' @param R an optional resource state matrix
#' @param axis Which resource axis (i.e., column in R) to plot along
#' @export
plot.metacommunity = function(x, R, axis = 1, ...) {
	x = x$species
	if(missing(R)) R = matrix(seq(0, 1, length.out=50), ncol=1)
	ypl = lapply(x, function(sp) sp$col(R) - sp$ext(R))
	yl = c(0, max(unlist(ypl)))

	# colours from colorbrewer
	cols = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")
	if(length(cols) > length(x))
		cols = cols[1:length(x)]

	args = .default_plot_pool_options()
	args$x = 0
	args$y = 0
	args$xlim = range(R[,axis])
	args$ylim = yl
	par(mar = c(5,4,4,6))
	do.call(plot, args)

	.make_line = function(y, col, args) {
		args$x = R[,axis]
		args$y = y
		args$col = col
		args$type = 'l'
		do.call(lines, args)
	}
	invisible(mapply(.make_line, ypl, cols, MoreArgs = list(args = args)))

	xpd = par()$xpd
	par(xpd = TRUE)
	legend("topright", legend = 1:length(x), col = cols, title = "species",
		   lty = args$lty, lwd = args$lwd, bty = "n", inset=c(-0.15,0), cex=0.7)
	par(xpd = xpd)

}


#' Plot species niches
#' @param x A species
#' @param R a resource state matrix
#' @param axis Which resource axis (i.e., column in R) to plot along
#' @export
plot.species = function(x, R, axis = 1, ...) {
	if(missing(R)) R = matrix(seq(0, 1, length.out=50), ncol=1)
	yc = x$col(R)
	ye = x$ext(R)
	yl = range(c(yc, ye))

	args = .default_plot_species_options(...)
	cols = args$col
	args$x = R[,axis]
	args$y = yc
	if(!"ylim" %in% names(args)) args$ylim = yl
	args$col = cols[1]
	par(mar = c(5,4,4,6))
	do.call(plot.default, args)

	args$y = ye
	args$col = cols[2]
	do.call(lines, args)

	xpd = par()$xpd
	par(xpd = TRUE)
	legend("topright", legend = c("colonisation", "extinction"), col = cols,
		   lty = args$lty, lwd = args$lwd, bty = "n", inset=c(-0.15,0), cex=0.7)
	par(xpd = xpd)
}


#' Set default plot options when not user-specified
#' @keywords internal
.default_river_plot_options = function(...) {
	dots = 	list(...)
	nms = names(dots)
	if(!"edge.width" %in% nms) dots$edge.width = 10
	if(!"vertex.color" %in% nms) dots$vertex.color = "#7BA08C"
	if(!"edge.color" %in% nms) dots$edge.color = "#a6bddb"
	if(!"edge.arrow.size" %in% nms) dots$edge.arrow.size = 0.2 * dots$edge.width
	return(dots)
}


#' Set default plot options when not user-specified
#' @keywords internal
.default_plot_pool_options = function(...) {
	dots = 	.default_plot_species_options(...)
	dots$ylab = "Dominant eigenvalue"
	dots$type = "n"
	dots$col = NULL
	return(dots)
}

#' Set default plot options when not user-specified
#' @keywords internal
.default_plot_species_options = function(...) {
	dots = list(...)
	nms = names(dots)
	if(!"xlab" %in% nms) dots$xlab = "Resource concentration"
	if(!"ylab" %in% nms) dots$ylab = "Colonisation/extinction rate"
	if(!"type" %in% nms) dots$type = "l"
	if(!"bty" %in% nms) dots$bty = "n"
	if(!"lwd" %in% nms) dots$lwd = 1.5
	if(!"col" %in% nms) {
		dots$col = c("#1f78b4", "#e31a1c")
	} else {
		if(length(dots$col) == 1)
			dots$col = rep(dots$col, 2)
	}
	return(dots)
}
