#' Algae N:P use, thanks to Thomas Fuß
site_by_species = read.csv(file.path("data-raw", "algae", "spec_cs.csv"), row.names = 1)
site_by_env = read.csv(file.path("data-raw", "algae", "env_cs.csv"), row.names = 1)
catch_area = site_by_env$catch
site_by_species = as.matrix(site_by_species)
site_by_env = as.matrix(site_by_env$NP_mean)
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


rn = river_network(adjacency = adj, discharge = drop_units(Q)*5, area = drop_units(cs_area),
	layout = layout)

# and we need starting resource concentrations and species distributions
R = matrix(site_by_env$NP_mean, ncol = 1, dimnames = list(rownames(site_by_env), "N:P"))


network = river_network(adjacency = adj, discharge = Q, state = R, layout = layout)
site_by_species(network) = site_by_species

# TODO: decide on boundary conditions; default is same as starting state


# TODO: add something like this legend to my own plot functions
# plot(network, variable = "N:P", edge.arrow.size = 0.2, edge.width = 20)
# col = scales::col_numeric("PuBu", c(0,1))(seq(0, 1, 0.1))
# legend("bottomleft", legend=seq(14, 660, 70), fill=col, bty='n', cex=0.7, title="N:P")


# niche_loc is the "optimal" N:P ratio for a species; mean of N:P at all sites where sp is found
niche_loc = (t(site_by_species) %*% site_by_env) / colSums(site_by_species)
colnames(niche_loc) = "N_to_P"

# niche_breadth is the breadth of the niche; we use sd, could also use the range, CV etc
niche_breadth = sweep(site_by_species, 1, site_by_env, FUN = `*`)
niche_breadth[niche_breadth == 0] = NA
niche_breadth = apply(niche_breadth, 2, sd, na.rm = TRUE)

niches = data.frame(species = rownames(niche_loc), location = niche_loc[,1], 
	breadth = niche_breadth)

algae = list(niches = niches, adjacency = adj, discharge = Q, )
usethis::use_data(algae, overwrite = TRUE)

# mc = metacommunity(niche_loc, niche_breadth, niche_lim = c(0, 650))

# alg_flume = flume(mc, network, dt = 1)
