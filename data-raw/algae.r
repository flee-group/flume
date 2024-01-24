#' Algae N:P use, thanks to Thomas Fuß

library(flume)
library(ggplot2)
library(data.table)
data(flume_networks)

# river network
kamp = flume_networks$kamp

# site by ASV matrix
kdat = readRDS("data-raw/algae/kamp_flume.rds")

# site coordinates
sites = read.csv("data-raw/algae/kamp_sites.csv")


# remove the original K03 site, replace with K03plus, (ask Thomas why)
rownames(kdat)[rownames(kdat) == "K03"] = "K03_old"
rownames(kdat)[rownames(kdat) == "K03plus"] = "K03"
kdat$name = sub("K0", "K", rownames(kdat))
kdat = merge(kdat, sites, by = "name")


# change the units to put N and P in the same currency
kdat$SRP = kdat$SRP/1000



# plotting locations of nodes, used to match nodes to sites
pldat = data.frame(attr(kamp, "layout"))
pldat$site = attr(kamp ,"names_sites")


### GRAPHICAL EXPLORATION BELOW, FOLLOWED BY MANUAL MATCHING OF NODES AND SITES

# ggplot(pldat, aes(x = X, y = Y)) + geom_point() + geom_label(aes(label = site)) + 
# 	geom_point(data = kdat, aes(x = x ,y = y), colour = 'blue') + 
# 	geom_text(data = kdat, aes(x = x, y = y, label = name), nudge_x = 1000, colour = "blue", size = 3)

# # zoom to help with points around node 6
# ggplot(pldat, aes(x = X, y = Y)) + geom_point() + geom_label(aes(label = site)) + 
# 	geom_point(data = kdat, aes(x = x ,y = y), colour = 'blue') + 
# 	geom_text(data = kdat, aes(x = x, y = y, label = name), nudge_x = 1000, colour = "blue", size = 3) + 
# 	xlim(4695802, 4710802) + ylim(2841328, 2847328)

# this mapping of sites to flume reaches is just done manually / graphically
# they follow the whole confluence logic, so be sure the connections are correct for groups of 3
# K3: 3
# K2: 5
# K1: 27
# K30 (6) -> K28 (not in network) & K29 (7) #  these need zooming, one will be 7, the other is not in the network
# K6(28) -> K5 (29) & K4 (38)
# K21 (29) -> K19 (33) & K20 (30)
# K9 (41) -> K7 (44) & K8 (43)
# K16 (51) -> K17 (52) & K18 (not in network)
# K15 (47) -> K14 (48) * K13 (50)
# K22 (31) -> K23 (34) & K24 (not in network)
# K33 (8) -> K32 (9) & K31 (not in network)
# K36 (9) -> K35 (10) & K34 (14)
# K42 (18) ->K41 (19) & K40 (23)
# K39 (15) -> K38 (16) & K37 (17)


# map site names to nodes in the kamp network
# sites covering tributaries that are not in the network are NA
name_mapping = data.frame(
	name = c("K1", "K2", "K3", "K4", "K5", "K6", "K7", "K8", "K9", "K13", "K14", "K15", 
		"K16", "K17", "K18", "K19", "K20", "K21", "K22", "K23", "K24",  "K28", "K29", "K30", 
		"K31", "K32", "K33", "K34", "K35", "K36", "K37", "K38", "K39", "K40", "K41", "K42"),
	node = c( 27, 5, 3, 38, 29, 28, 44, 43, 41, 50, 48, 47, 51, 52, NA, 33, 30, 29, 31, 34, NA, 
		   NA, 7, 6, NA, 9, 8, 14, 10, 9, 17, 16, 15, 23, 19, 18))


# add node numbers and Q from the flume network to the algae data
mg_kamp = merge(name_mapping, kdat, by = "name")
mg_kamp = merge(mg_kamp, data.frame(node = 1:52, Q = state(kamp, "Q")), by = "node", all.y = TRUE)


# Some rows are duplicated, beacuse there was more than one algae sampling point
# in a single node from the network
# We aggregate by taking the mean
mg_kamp = data.table(mg_kamp)
mg_kamp$name = NULL
mg_kamp = mg_kamp[, lapply(.SD, mean, na.rm = TRUE), .(node)]
mg_kamp = as.data.frame(mg_kamp)


# Interpolate missing values to all nodes
# a bit of a hacky one off function to compute a discharge weighted mean of the needed columns
wtmean = function(dat, i, cols, wtcol = "Q") {
	wt = dat[[wtcol]][i]/sum(dat[[wtcol]][i])
	res = mapply(\(ii, w) dat[ii, cols] * w, i, wt)
	# mapply returns a weird matrix-list combo, so some modification is needed
	rowSums(matrix(unlist(res), ncol = length(i)))
}

# same, but to interpolate when there are multiple steps
linint = function(dat, i, cols, steps) {
	x = dat[i[1], cols]
	y = dat[i[2], cols]
	mapply(\(xx, yy) seq(xx, yy, length.out = steps + 2)[2:(steps+1)], x, y)
}


## columns to interpolate
j = c(match(c("SRP", "NH4", "NO3", "DIN"), colnames(mg_kamp)), 
		grep("ASV", colnames(mg_kamp))) 


## summary of rules for assigning starting values
## 1. known (measured) values are assigned directly
## 2. unknowns at confluences get a discharge weighted mean
## 3. unknowns directly between two knowns get a discharge weighted mean
## 4. unknowns that are on linear stretches get a linear interpolation of the nearest upstream
##		and downstream known (e.g., known -- unk1 -- unk2 -- known)
## 5. unknowns at the edge of the network (headwaters, the outlet) get a copy of the nearest known






mg_kamp[1:2, j] = mg_kamp[3, j] # nodes 1 and 2 are ds of 3, so just copy state
mg_kamp[25, j] = mg_kamp[27, j] 
mg_kamp[4, j] = wtmean(mg_kamp, c(5, 25), j)
mg_kamp[11:12, j] = mg_kamp[10, j]
mg_kamp[13, j] = wtmean(mg_kamp, c(14, 15), j)
mg_kamp[24, j] = mg_kamp[13, j] 
mg_kamp[20:22, j] = mg_kamp[19, j] 
mg_kamp[26, j] = wtmean(mg_kamp, c(27, 28), j)
mg_kamp[32, j] = wtmean(mg_kamp, c(31, 33), j)
mg_kamp[35:37, j] = mg_kamp[34, j] 
mg_kamp[39:40, j] = linint(mg_kamp, c(38, 41), j, steps = 2)
mg_kamp[42, j] = mg_kamp[43, j] 
mg_kamp[45, j] = wtmean(mg_kamp, c(44, 47), j)
mg_kamp[46, j] = mg_kamp[45, j] 
mg_kamp[49, j] = mg_kamp[48, j] 


# compute niches by ASV
ntop = kdat$DIN / kdat$SRP
j = grep("ASV", colnames(kdat))
asv = kdat[,j]
ntop_mean = sapply(asv, \(x) weighted.mean(ntop, x))
ntop_sd = mapply(\(x, mu, wt) {
	wt = wt/sum(wt)
	sqrt(sum((wt * (x - mu)^2 )))
	}, 
	wt = asv, mu = ntop_mean, MoreArgs = list(x = ntop), SIMPLIFY = TRUE)
# total niche height modified by the total abundance for each asv, with the most abundant species
# normalized to twice the default
niche_scale = 2 * (colSums(asv) / max(colSums(asv)))


niche_params = data.frame(location = ntop_mean, breadth = ntop_sd, scale = niche_scale)
rownames(niche_params) = colnames(asv)



# set up the river network with chemistry & community starting state
r0 = as.matrix(mg_kamp[, c('DIN', "SRP")])
sp0 = as.matrix(mg_kamp[, grep("ASV", colnames(kdat))])
sp0[sp0 > 0] = 1


algae = list(niche_params = niche_params,	r0 = r0, sp0 = sp0)
usethis::use_data(algae, overwrite = TRUE)











