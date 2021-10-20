test_that("Species flux works", {
	comm = metacommunity()
	Q = c(13,12,8,7,1,4,1,1,2,1)
	nsites = length(Q)
	adj = matrix(0, nrow=nsites, ncol=nsites)
	adj[2,1] = adj[3,2] = adj[9,2] = adj[4,3] = adj[5,4] = adj[6,4] =
		adj[7,6] = adj[8,6] = adj[10,9] = 1

	network = river_network(adjacency = adj, discharge = Q)
	state(network, "resources") = matrix(seq(0, 1, length.out = nsites), ncol=1)
	boundary(network, "resources") = state(network, "resources")
	state(network, "species") = community_random(network, comm, prevalence = c(0.35, 0.65))
	boundary(network, "species") = state(network, "species") * 0
	expect_error(cp <- col_prob(comm, network, dt=1), regex=NA)
	expect_error(ep <- ext_prob(comm, network, dt=1), regex=NA)
	expect_error(Rflux <- dRdt(0, state(network, "resources"), 
		.dRdt_params(comm, network)), regex=NA)
})

test_that("Fluxes with static resources", {
	data(algae)
	loc = apply(algae$r0, 2, function(x) colMeans(x*algae$sp0))
	bre = apply(algae$r0, 2, function(x) apply((x*algae$sp0), 2, sd))
	nopts = list(location = loc, breadth = bre, static = 1, r_lim = t(apply(algae$r0, 2, range)))
	mc = metacommunity(nsp = nrow(algae$niches), nr = 2, niches = niches_custom, niche_args = nopts,
		sp_names = algae$niches$species, r_names = c("N", "P"))
	fl = flume(mc, algae$network, algae$sp0, algae$r0)
	expect_error(dR <- dRdt(0, state(fl$network[[1]], "resources"), 
		.dRdt_params(fl$metacom, fl$network[[1]])), regex = NA)
	expect_true(all(dR[[1]][,1] == 0))
})