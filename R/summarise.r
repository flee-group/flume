#' Methods for summarising simulation results
#'
#' @param x A fitted flume model 
#' @param by The margin(s) across which to summarise; see 'details'
#' @param stat The statistic(s) to compute; see 'details'
#' @param quantile The quantiles to compute; see 'details'
#' @param network
#'
#' @details The `by` parameter controls how detailed of a summary to create; adding more
#' margins creates a more detailed summary. The default behaviour provides statistics for each
#' time step and each reach, with quantiles computed across replicate networks. Allowable values:
#' 
#'	* time
#'	* reach
#'	* species (incompatible with resources)
#'	* resources (incompatible with species)
#' 
#' The statistics desired interact with the margins selected; certain values for `stat` are only
#' compatible with specific values within `by`. If missing, a default set will be selected given the
#' value in `by`. Allowable values, with compatible margins in parentheses, and required margins
#' in bold:
#'
#'	* occupany (**species**, reach, time)
#'	* EF (resources, reach, time)
#'	* richness (reach, time)
#'
#' If the fitted model includes multiple replicate networks, then summarise by default returns
#' upper and lower quantiles along with the median, summarised across networks. 
#' Only two quantiles may be provided. If quantile is NULL, then all replicates will be returned.
#' @return A data.table with the desired summary statistics
#' @import data.table
#' @export
summarise = function(x, by = c("time", "reach"), stat, quantile = c(0.05, 0.95)) {
	st_choices = list(occupancy = .st_occupancy, EF = .st_ef, richness = .st_richness)
	stat = match.arg(stat, choices = names(st_choices), several.ok = TRUE)
	res = lapply(stat, \(s) st_choices[[s]](x, by, quantile))
	if(length(res) == 1) {
		res = res[[1]]
	} else {
		res = Reduce(function(...) merge(..., all = TRUE), res)
	}
	res
}


# Statistics for flumes
#' @keywords internal
#' @import data.table
.st_richness = function(x, by, quantile) {
	if("species" %in% by)
		stop("Cannot compute richness by species")
	by = by[by != "resources"]
	S = .output_table(x, "species")
	S = lapply(S, function(s) s[, .(richness = sum(occupancy)), keyby = by])
	S = data.table::rbindlist(S, idcol = "network")
	if(length(x[["networks"]]) == 1) {
		res = S[, .(richness = median(richness)), keyby = by]
	} else if(!is.null(quantile)) {
		res = S[, .(richness = median(richness), 
			richness_lo = quantile(richness, quantile[1]), 
			richness_hi = quantile(richness, quantile[2])), keyby = by]
	} else {
		res = S
	}
	res
}

# Statistics for flumes
#' @keywords internal
#' @import data.table
.st_ef = function(x, by, quantile) {
	rn = attr(x[["networks"]][[1]], "names_resources")
	S = .output_table(x, "reaction")
	if("resources" %in% by) {
		ef_n = paste0("ef_", rn)
		by = by[-which(by == "resources")]
		S = lapply(S, function(s) {
			cols = match(rn, colnames(s))
			s = s[, lapply(.SD, sum), .SDcols = cols, keyby = by]
			cols = match(rn, colnames(s))
			colnames(s)[cols] = ef_n
			s
		})
	} else {
		S = lapply(S, function(s) {
			s$ef = 0
			for(r in rn) {
				s[["ef"]] = s[["ef"]] + s[[r]]
				s[[r]] = NULL
			}
			s[, .(EF = sum(ef)), keyby = by]
		})
	}
	S = data.table::rbindlist(S, idcol = "network")
	cols_S = grep("(^ef_)|(EF)", colnames(S))
	if(length(x[["networks"]]) == 1) {
		res = S[, lapply(.SD, median), keyby = by, .SDcols = cols_S]
	} else if(! is.null(quantile)) {
		S_m = S[, lapply(.SD, median), keyby = by, .SDcols = cols_S]
		cols_res = grep("(^ef_)|(EF)", colnames(S_m))
		S_l = S[, lapply(.SD, quantile, probs = quantile[1]), keyby = by, .SDcols = cols_S]
		colnames(S_l)[cols_res] = paste0(colnames(S_l)[cols_res], "_lo")
		S_u = S[, lapply(.SD, quantile, probs = quantile[2]), keyby = by, .SDcols = cols_S]
		colnames(S_u)[cols_res] = paste0(colnames(S_u)[cols_res], "_hi")
		res = cbind(S_m, S_l[,..cols_res], S_u[,..cols_res])
	} else {
		res = S
	}
	res
}




#' @keywords internal
#' @import data.table
.st_occupancy = function(x, by, quantile) {
	S = .output_table(x, "species")
	S = lapply(S, function(s) s[, .(occupancy = sum(occupancy) / .N), keyby = by])
	S = data.table::rbindlist(S, idcol = "network")
	if(length(x[["networks"]]) == 1) {
		res = S[, .(occupancy = median(occupancy)), keyby = by]
	} else if(!is.null(quantile)) {
		res = S[, .(occupancy = median(occupancy), 
			occupancy_lo = quantile(occupancy, quantile[1]), 
			occupancy_hi = quantile(occupancy, quantile[2])), keyby = by]
	} else {
		res = S
	}
	res
}

#' Extract a variable across an entire network or series of networks
#' @param x A [river_network()] or [flume()]
#' @param v The variable to extract
#' @return A data.table
#' @keywords internal
setGeneric(".output_table", function(x, v) {
	standardGeneric(".output_table")
})

## delete when converting these to S4
setOldClass("river_network")
setOldClass("flume")

#' @import data.table
setMethod(".output_table", c(x = "river_network"),
	function(x, v) {
		S = state(x, v, history = TRUE)
		res = rbindlist(lapply(S, function(s) {
				val = data.table(s)
				val$reach = rownames(s)
				val
			}), idcol = "time")
		if(v == "species") {
			res = melt(res, id.vars = c("reach", "time"), 
				variable.name = "species", value.name = "occupancy")
		} else if(v == "resources") {
			res = melt(res, id.vars = c("reach", "time"), 
				variable.name = "resources", value.name = "concentration")			
		}
		res
	}
)

setMethod(".output_table", c(x = "flume"),
	function(x, v) {
		cores = ifelse(.Platform$OS.type == "unix", getOption("mc.cores", 1L), 1L)
		parallel::mclapply(x[["networks"]], .output_table, v = v, 
			mc.cores = cores)
	}
)

