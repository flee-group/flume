#' Plot the results of a flume simulation
#' @param x A [flume()]
#' @param variable Variable to plot; 'occupancy' plots occupancy by species over time, 'resources' plots
#' resource concentration over time
#' @param type How to compute the plot; either shows the average for the entire network, or plots by reach.
#' @return A ggplot2 object
#' @export
plot.flume = function(x, variable = c("occupancy", "resources"), type = c("network", "reach")) {
	if(!requireNamespace("ggplot2", quietly=TRUE))
		stop("Package ggplot2 is required for this functionality")

	variable = match.arg(variable)
	type = match.arg(type)

	if(variable == "occupancy") {
		.occupancy_plot(x, type=type)
	} else {
		.resource_plot(x)
	}
}

.occupancy_plot = function(x, type) {
	occ = occupancy(x, type=type)
	if(type == 'network') {
		pl = ggplot2::ggplot(occ, ggplot2::aes(x = time, y = occupancy, colour = species))
	} else {
		pl = ggplot2::ggplot(occ, ggplot2::aes(x = time, y = occupancy, colour = as.factor(reach))) +
			ggplot2::facet_grid(.~species)
	}
	pl + ggplot2::geom_line() + ggplot2::theme_minimal() + ggplot2::ylim(0,1)
}

.resource_plot = function(x) {
	res = resource_summary(x)
	ggplot2::ggplot(res, ggplot2::aes(x = time, y = concentration, colour = as.factor(reach))) +
		ggplot2::geom_line() + ggplot2::theme_minimal() + ggplot2::facet_wrap(.~resource, scales = 'free')
}


#' Plot a river network
#' @details The argument 'variable' can either be a column number from the state variable matrix, a column name from
#' the state variable matrix, or the special name "site_by_species", which produces a plot of the network with
#' species presence-absence.
#' @param x A [river_network()]
#' @param variable If state is defined, the column to use for plotting; see 'details'
#' @param t Optional, time step to print; if missing, defaults to most recent
#' @param zlim Optional, z limits for colour scales when plotting a state variable
#' @param ... Additional arguments to [igraph::plot.igraph()]
#' @export
plot.river_network = function(x, variable = 1, t, zlim, ...) {
	if(!requireNamespace("igraph", quietly = TRUE)) {
		stop("Package 'igraph' is required for plotting, please install it and try again")
	}

	args = .default_river_plot_options(...)
	args$x = igraph::graph_from_adjacency_matrix(adjacency(x), mode = "directed")

	## weight the plot by discharge for the UPSTREAM node
	up_nodes = as.integer(igraph::ends(args$x, igraph::E(args$x))[,1])
	wt = discharge(x)[up_nodes]

	args$edge.width = (wt / max(wt)) * args$edge.width
	if(!("layout" %in% args) && "layout" %in% names(attributes(x)))
		args$layout = attr(x, "layout")

	## colour scale, if available and desired
	if(variable == "site_by_species") {
		if(missing(t)) {
			S = site_by_species(x)
		} else {
			S = site_by_species(x, TRUE)[[t]]
		}
		nsp = ncol(S)
		# colours from colorbrewer
		cols = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")
		args$vertex.label.color = '#444444'

		par(mfrow = .set_mfrow(nsp))
		for(i in 1:nsp) {
			args$vertex.color = rep('white', nrow(S))
			args$vertex.color[S[,i] == 1] = cols[(i %% length(cols))+1]
			do.call(plot, args)
		}
	}
	else {
		if(!is.null(state(x)) && !("vertex.color" %in% names(args)) && requireNamespace("scales", quietly = TRUE)) {
			if(missing(t)) {
				R = state(x)[,variable]
			} else {
				R = state(x, TRUE)[[t]][,variable]
			}
			if(missing(zlim))
				zlim = range(R)
			args$vertex.color = scales::col_numeric("PuBu", zlim)(R)
			args$vertex.label.color = rev(scales::col_numeric("YlOrBr", zlim)(R))
		} else {
			message("Install the scales package for plotting the state variable with a colour scale")
			args$vertex.color = "#7BA08C"
		}
		do.call(plot, args)
	}
}

#' Plot species independent stable envelopes
#' @param x A [metacommunity()]
#' @param R an optional resource state matrix
#' @param axis Which resource axis (i.e., column in R) to plot along
#' @param res Plotting resolution, how many points along the x-axis to plot; more produces smoother
#'		lines.
#' @export
plot.metacommunity = function(x, R, axis = 1, res = 100, ...) {
	if(missing(R)) {
		R = matrix(0, nrow = res, ncol = attr(x, "n_r"))
		R[, axis] = seq(attr(x, "r_lim")[axis, 1], attr(x, "r_lim")[axis, 2],
			length.out = res)
	}
	pdat = lapply(x$species, function(sp)
		data.table::data.table(R = R[, axis], lam = f_niche(sp, R)))
	names(pdat) = attr(x, "spnames")
	pdat = data.table::rbindlist(pdat, idcol = "species")

	p1 = ggplot2::ggplot(pdat, ggplot2::aes(x = R, y = lam.V1, colour = factor(species))) +
		ggplot2::geom_line() + ggplot2::scale_colour_brewer(type = "qual", palette = "Set2") +
		ggplot2::theme_minimal() + ggplot2::xlab(attr(x, "rnames")[axis]) +
		ggplot2::labs(colour = "Species") + 
		ggplot2::ylab("Dominant Eigenvalue") + ggplot2::geom_hline(ggplot2::aes(yintercept = 0))

	comp = x$competition
	comp[upper.tri(comp)] = NA
	comp = reshape2::melt(comp)
	comp = comp[complete.cases(comp), ]
	colnames(comp) = c("sp1", "sp2", "competition")

	p2 = ggplot2::ggplot(comp) +
		ggplot2::geom_tile(ggplot2::aes(x = sp1, y = sp2, fill = competition)) +
		ggplot2::scale_fill_viridis_c(option = "plasma") + ggplot2::theme_minimal() +
		ggplot2::xlab("") + ggplot2::ylab("")

	gridExtra::grid.arrange(p1, p2, nrow = 1)
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
	if(!"edge.color" %in% nms) dots$edge.color = "#a6bddb"
	if(!"edge.arrow.size" %in% nms) dots$edge.arrow.size = 0.1 * dots$edge.width
	return(dots)
}


# #' Set default plot options when not user-specified
# #' @keywords internal
# .default_plot_pool_options = function(...) {
# 	dots_orig = list(...)
# 	nms = names(dots_orig)
# 	dots = 	.default_plot_species_options(...)
# 	if(!"ylab" %in% nms) dots$ylab = "Dominant eigenvalue"
# 	if(!"type" %in% nms) dots$type = "n"
# 	dots$col = NULL
# 	# print(dots_orig)
# 	# print(dots)
# 	return(dots)
# }

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


#' Choose some sensible number of plots for a river network
.set_mfrow = function(n) {
	if(n <= 3) {
		nr = 1
	} else if(n <= 6) {
		nr = 2
	} else if(n <= 12) {
		nr = 3
	} else if(n <= 24) {
		nr = 4
	} else {
		ratio = 10/16
		nr = ceiling(sqrt(ratio * n))
	}
	nc = ceiling(n / nr)
	return(c(nr, nc))
}
