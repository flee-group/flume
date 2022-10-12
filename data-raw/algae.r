#' Algae N:P use, thanks to Thomas Fuß
site_by_species = read.csv(file.path("data-raw", "algae", "spec_cs.csv"), row.names = 1)
site_by_env = read.csv(file.path("data-raw", "algae", "env_cs.csv"), row.names = 1)
catch_area = site_by_env$catch
site_by_species = as.matrix(site_by_species)
adj = read.csv(file.path("data-raw", "algae", "bodingbach.csv"), sep = ";", row.names = 1)
colnames(adj) = rownames(adj)
adj = as.matrix(adj)
nsites = nrow(site_by_env)

# we also define a layout to make plotting nice
layout = matrix(c(0, 0,   2, 3,   1, 3,   -5, 1.5,   -5, 2.5,   -6, 2.5,   -5, 3,   -4, 2,
	-3, 2,   -2, 2.5,   -2, 2,   -1, 2,   1, 2,   -2, 1,   -1, 1.5,   -1, 1,   0, 1),
    byrow = TRUE, nrow = nsites)


## we need the discharge at every node
#for now I choose catchment area, scaled from Raymond paper (see watershed package)
library(units)
# in cubic kilometers per day
Q = exp(-11.9749609 + 0.78 * log(catch_area / (1000^2)))
units(Q) = "km^3 / day"
# convert to m^3/s
Q = set_units(Q, "m^3/s")

z = set_units(exp(-0.895 +  0.294 * log(drop_units(Q))), "m")
w = set_units(exp(2.56 + 0.423 * log(drop_units(Q))), "m")
cs_area = z * w


rn = river_network(adjacency = adj, discharge = drop_units(Q), area = drop_units(cs_area),
	layout = layout, site_names = rownames(adj))

# and we need starting resource concentrations and species distributions
r0 = as.matrix(site_by_env[, c("N_mean", "P_mean")])


# niche_loc is the "optimal" N:P ratio for a species; mean of N:P at all sites where sp is found
ntop = matrix((r0[,1]/r0[,2]), ncol=1)
niche_loc = (t(site_by_species) %*% ntop) / colSums(site_by_species)
colnames(niche_loc) = "N:P"

# niche_breadth is the breadth of the niche; we use sd, could also use the range, CV etc
niche_breadth = sweep(site_by_species, 1, ntop, FUN = `*`)
niche_breadth[niche_breadth == 0] = NA
niche_breadth = apply(niche_breadth, 2, sd, na.rm = TRUE)

niches = data.frame(species = rownames(niche_loc), location = niche_loc[,1], 
	breadth = niche_breadth)


nopts = list(location = niches$location, breadth = niches$breadth)
mc = metacommunity(nsp = nrow(niches), nr = 2, niches = niches_custom, niche_args = nopts,
  sp_names = niches$species, r_names = c("N", "P"), ratio = c(1, 2))

alg_flume = flume(mc, rn, sp0 = site_by_species, st0 = r0)


algae = list(niches = niches, adjacency = adj, discharge = Q, cs_area = cs_area, layout = layout,
	network = rn, r0 = r0, sp0 = site_by_species, metacommunity = mc)




usethis::use_data(algae, overwrite = TRUE)





# TODO: add something like this legend to my own plot functions
# plot(network, variable = "N:P", edge.arrow.size = 0.2, edge.width = 20)
# col = scales::col_numeric("PuBu", c(0,1))(seq(0, 1, 0.1))
# legend("bottomleft", legend=seq(14, 660, 70), fill=col, bty='n', cex=0.7, title="N:P")







