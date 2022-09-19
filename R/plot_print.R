#' Plot the results of a flume simulation
#'
#' @param x A [flume()]
#' @param type Which plot to produce; see 'details'.
#' 
#' @details `flume` comes with a collection of default plots. All of them are built around the
#' [summarise()] function, which can be used to extract the data from a fitted `flume`. The
#' following plot `type`s are recognized (partial matches are allowed):
#'
#' 	* **occupancy**: (plots the proportion of sites occupied by each in the network against time,  
#'		with quantiles if replciates are present)
#' 
#' All plots return a `ggplot`, so further customisation is easy.
#' @return A ggplot2 object
#' @export
plot.flume = function(x, type = "occupancy") {
	pl_types = list(occupancy = .pl_occupancy, ef_time = .pl_ef_time, bef = .pl_bef)
	type = match.arg(type, names(pl_types))
	pl_types[[type]](x)
}

#' @keywords internal
.pl_occupancy = function(x) {
	occ = summarise(x, stat = "occupancy", by = c("time", "species"))
	pl = ggplot2::ggplot(occ, ggplot2::aes(x = time, y = occupancy, colour = species))
	if("occupancy_lo" %in% colnames(occ))
		pl = pl + geom_ribbon(aes(ymin = occupancy_lo, ymax = occupancy_hi, fill = species, 
			colour = NULL), alpha = 0.3)
	pl = pl + ggplot2::geom_line() + theme_flume() + ggplot2::ylim(0,1)
	pl
}

#' @keywords internal
.pl_ef_time = function(x) {
	summ = summarise(x, stat = "EF", by = c("time", "resources"))
	summ = data.table::melt(summ, id.vars = c("time"), variable.name = "resource")
	if(any(grepl(".+_lo", summ$resource))) {
		summ$stat = "value"
		pat = ".+_(lo|hi)"
		i = grep(pat, summ$resource)
		summ$stat[i] = sub(pat, "\\1", summ$resource[i])
		summ$resource = sub("ef(_(.*?)_?)(hi|lo)?$", "\\2", summ$resource, perl = TRUE)
		summ = data.table::dcast(summ, time + resource ~ stat)
	}
	col = scales::hue_pal()(1)
	pl = ggplot2::ggplot(summ, ggplot2::aes(x = time, y = value)) 
	if("lo" %in% colnames(summ))
		pl = pl + ggplot2::geom_ribbon(ggplot2::aes(ymin = hi, ymax=lo), fill = col, 
			alpha = 0.4, col = NA)
	pl = pl+ ggplot2::geom_line(colour = col) + theme_flume() + ggplot2::facet_wrap(~resource) + 
		ggplot2::ylab("EF")
	pl
}

.pl_bef = function(x) {
	summ = summarise(x, stat = c("EF", "richness"), 
		by = c("time", "resources", "reach"), quantile = NULL)
	summ = data.table::melt(summ, id.vars = c("network", "reach", "time", "richness"))
	pl = ggplot2::ggplot(summ, ggplot2::aes(x = factor(richness), y = value, fill = variable))
	pl = pl + ggplot2::geom_boxplot() + theme_flume() + ggplot2::ylab("EF")
	pl = pl + ggplot2::xlab("Species Richness") + ggplot2::labs(fill = "Resource")
	pl
}


#' ggplot2 theme for flume
#' @keywords internal
theme_flume = function() ggplot2::theme_minimal()

resource_plot = function(x) {
	res = resource_summary(x)
	ggplot2::ggplot(res, ggplot2::aes(x = time, y = concentration, colour = as.factor(reach))) +
		ggplot2::geom_line() + ggplot2::theme_minimal() + 
		ggplot2::facet_wrap(.~resource, scales = 'free') + ggplot2::labs(colour = "Reach")
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
		## choose a vertex colour scheme if none specified
		if(!("vertex.color" %in% names(args))) {
			if(is.null(state(x, "resources"))) {
				args$vertex.color = "#7BA08C"
			} else {
				R = state(x, "resources", history = TRUE)
				t = if(missing(t)) length(R) else t
				R = R[[t]][,variable]

				if(missing(zlim))
					zlim = range(R)
				args$vertex.color = scales::col_numeric("PuBu", zlim)(R)
			}
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
#' @param plot_comp logical, should the competition matrix be plotted as well?
#' @export
plot.metacommunity = function(x, axis = 1, res = 100, default_r, lwd = 1, xlim, ylim, 
	plot_comp = TRUE) {
	# loc = niche_par(x, "location")
	# sc = niche_par(x, "sd")
	if(missing(xlim))
		xlim = attr(x, "niche_lim")[axis, ]

	if(missing(default_r))
		default_r = rowMeans(attr(x, "niche_lim"))

	R = matrix(rep(default_r, each = res), nrow = res, ncol = length(default_r))
	R[,axis] = seq(xlim[1], xlim[2], length.out = res)

	niches = f_niche(x, N = R)
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

	if(plot_comp) {
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
	} else {
		p1
	}
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
	if(!"edge.width" %in% nms) dots$edge.width = 20
	if(!"edge.color" %in% nms) dots$edge.color = "#a6bddb"
	if(!"edge.arrow.size" %in% nms) dots$edge.arrow.size = 0
	if(!"vertex.label.color" %in% nms) dots$vertex.label.color ="#444444"
	if(!"vertex.label.cex" %in% nms) dots$vertex.label.cex = 0.8
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
