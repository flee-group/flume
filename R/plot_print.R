#' Plot the results of a flume simulation
#' @param x A [flume()]
#' @param variable Variable to plot; 'occupancy' plots occupancy by species over time, 'resources' plots
#' resource concentration over time
#' @param type How to compute the plot; either shows the average for the entire network, or plots by reach.
#' @return A ggplot2 object
#' @export
plot.flume = function(x, variable = c("occupancy", "resources"), type = c("network", "reach")) {
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
#' @details The argument 'variable' can either be a column number from the state variable matrix, 
#' a column name from the state variable matrix, or the special name "species", which produces a 
#' plot of the network with species presence-absence.
#' @param x A [river_network()]
#' @param variable If state is defined, the column to use for plotting; see 'details'
#' @param t Optional, time step to print; if missing, defaults to most recent
#' @param zlim Optional, z limits for colour scales when plotting a state variable
#' @param ... Additional arguments to [igraph::plot.igraph()]
#' @examples
#' Q = rep(1, 4)
#' adj = matrix(0, nrow = 4, ncol = 4)
#' adj[1,2] = adj[2,3] = adj[4,3] = 1
#' rn = river_network(adj, Q)
#' plot(rn)
#' plot.river_network(rn)
#' @export
plot.river_network = function(x, variable = 1, t, zlim, ...) {
	if(!requireNamespace("igraph", quietly = TRUE)) {
		stop("Package 'igraph' is required for plotting, please install it and try again")
	}

	args = .default_river_plot_options(...)
		args$x = igraph::graph_from_adjacency_matrix(adjacency(x), mode = "directed")

	## weight the plot by discharge for the UPSTREAM node
	if(is.null(rownames(adjacency(x)))) {
		node_names = as.character(1:nrow(adjacency(x)))
	} else {
		node_names = rownames(adjacency(x))
	}
	up_nodes = igraph::ends(args$x, igraph::E(args$x))[,1]
	up_nodes = match(up_nodes, node_names)
	wt = discharge(x)[up_nodes]

	args$edge.width = (wt / max(wt)) * args$edge.width
	if(!("layout" %in% args) && "layout" %in% names(attributes(x)))
		args$layout = attr(x, "layout")

	## colour scale, if available and desired
	if(variable == "site_by_species") {
		warning("variable = 'site_by_species' is deprecated, use variable = 'species' instead.")
		variable = "species"
	}
	if(variable == "species") {
		if(missing(t)) {
			S = state(x, "species")
		} else {
			S = state(x, "species", TRUE)[[t]]
		}
		nsp = ncol(S)
		# colours from colorbrewer
		cols = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")

		par(mfrow = .set_mfrow(nsp))
		for(i in 1:nsp) {
			args$vertex.color = rep('white', nrow(S))
			args$vertex.color[S[,i] == 1] = cols[(i %% length(cols))+1]
			do.call(plot, args)
		}
	}
	else {
		if(!is.null(state(x, "resources")) && 
					!("vertex.color" %in% names(args)) && 
					requireNamespace("scales", quietly = TRUE)) {
			if(missing(t)) {
				R = state(x, "resources")[,variable]
			} else {
				R = state(x, "resources", history = TRUE)[[t]][,variable]
			}
			if(missing(zlim))
				zlim = range(R)
			args$vertex.color = scales::col_numeric("PuBu", zlim)(R)
			
		} else {
			if(!requireNamespace("scales", quietly = TRUE))
				message("Install the scales package for plotting the state variable with a colour scale")
			args$vertex.color = "#7BA08C"
		}
		do.call(plot, args)
	}
}

#' Plot species independent stable envelopes
#' @param x A [metacommunity()]
#' @param axis Which niche axis to plot along
#' @param res Plotting resolution, how many points along the x-axis to plot; more produces smoother
#'		lines.
#' @param default_r Default values for niche dimensions not being plotted along. If not specified
#'		the default is to use the midpoint of possible values for each niche axis.
#' @param xlim A vector of two, optional axis limits for plotting
#' @param ylim A vector of two, optional axis limits for plotting
#' @export
plot.metacommunity = function(x, axis = 1, res = 100, default_r, lwd = 1, xlim, ylim) {
	# loc = niche_par(x, "location")
	# sc = niche_par(x, "sd")
	if(missing(xlim))
		xlim = attr(x, "niche_lim")[axis, ]

	if(missing(default_r))
		default_r = rowMeans(attr(x, "niche_lim"))

	R = matrix(rep(default_r, each = res), nrow = res, ncol = length(default_r))
	R[,axis] = seq(xlim[1], xlim[2], length.out = res)

	niches = f_niche(x, N = R)
	colnames(niches) = attr(x, "sp_names")
	pldat = as.data.frame(niches)
	pldat$r = R[,axis]
	pldat = reshape2::melt(pldat, id.vars = "r", variable.name = "species")

	if(missing(ylim)) {
		ylim = c(-0.3*max(pldat$value), max(pldat$value))
	}

	i = which(pldat$value < ylim[1] | pldat$value > ylim[2])
	if(length(i) > 0)
		pldat = pldat[-i,]

	p1 = ggplot2::ggplot(pldat) + 
		ggplot2::geom_line(ggplot2::aes(x = r, y = value, colour = species), size = lwd) +
		ggplot2::scale_colour_brewer(type = "qual", palette = "Set2") +
		ggplot2::theme_minimal() + ggplot2::xlab(attr(x, "niche_names")[axis]) + 
		ggplot2::labs(colour = "Species") + ggplot2::ylim(ylim[1], ylim[2]) +
		ggplot2::ylab("Dominant Eigenvalue") + 
		ggplot2::geom_hline(ggplot2::aes(yintercept = 0), size = 1.2)

	comp = x$competition
	comp[upper.tri(comp)] = NA
	comp = reshape2::melt(comp)
	comp = comp[complete.cases(comp), ]
	colnames(comp) = c("sp1", "sp2", "competition")
	comp$sp1 = factor(attr(x, "sp_names")[comp$sp1], levels = attr(x, "sp_names"))
	comp$sp2 = factor(attr(x, "sp_names")[comp$sp2], levels = attr(x, "sp_names"))

	p2 = ggplot2::ggplot(comp) +
		ggplot2::geom_tile(ggplot2::aes(x = sp1, y = sp2, fill = competition)) +
		ggplot2::scale_fill_viridis_c(option = "plasma") + ggplot2::theme_minimal() +
		ggplot2::xlab("") + ggplot2::ylab("")

	gridExtra::grid.arrange(p1, p2, nrow = 1)
}



#' Plot species niches
#' @param x A species
#' @param axis Which niche axis (i.e., column in R) to plot along
#' @param res Plotting resolution, how many points along the x-axis to plot; more produces smoother
#'		lines.
#' @return A ggplot2 object; additional modifications to the plot can be made using usual ggplot2
#' 		syntax.
#' @export
plot.species = function(x, R, axis = 1, res = 100, lwd = 1) {
	loc = niche_par(x, "location")
	sc = niche_par(x, "sd")

	xlim = c(-2, 2) * sc[axis] + loc[axis]
	R = sapply(loc, function(X) rep(X, res))
	R[, axis] = seq(xlim[1], xlim[2], length.out = res)

	pl = data.frame(R = rep(R[,axis], 2), ce_rate = c(x$col(R), x$ext(R)), 
		type = rep(c("colonisation", "extinction"), each = nrow(R)))

	ggplot2::ggplot(pl) + 
		ggplot2::geom_line(ggplot2::aes(x = R, y = ce_rate, colour = type), size = lwd) +
		ggplot2::scale_colour_discrete(c("#386cb0", "#bf5b17")) + ggplot2::xlim(xlim[1], xlim[2]) +
		ggplot2::xlab("Resource concentration") + ggplot2::theme_minimal() + 
		ggplot2::ylab("Colonisation/extinction rate") +
		ggplot2::theme(legend.title = ggplot2::element_blank())
}


#' Set default plot options when not user-specified
#' @keywords internal
.default_river_plot_options = function(...) {
	dots = 	list(...)
	nms = names(dots)
	if(!"edge.width" %in% nms) dots$edge.width = 10
	if(!"edge.color" %in% nms) dots$edge.color = "#a6bddb"
	if(!"edge.arrow.size" %in% nms) dots$edge.arrow.size = 0.1 * dots$edge.width
	if(!"vertex.label.color" %in% nms) dots$vertex.label.color ="#444444"
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
