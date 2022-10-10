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
	nopts = list(location = loc, breadth = bre)
	mc = metacommunity(nsp = nrow(algae$niches), nr = 2, niches = niches_custom, niche_args = nopts,
		sp_names = algae$niches$species, r_names = c("N", "P"), 
		static = 1, r_lim = t(apply(algae$r0, 2, range)))
	fl = flume(mc, algae$network, algae$sp0, algae$r0)
	expect_error(dR <- dRdt(0, state(fl$network[[1]], "resources"), 
		.dRdt_params(fl$metacom, fl$network[[1]])), regex = NA)
	expect_true(all(dR[[1]][,1] == 0))
})

test_that("Fluxes with intermittency", {
	data(algae)
	loc = apply(algae$r0, 2, function(x) colMeans(x*algae$sp0))[,1]
	bre = apply(algae$r0, 2, function(x) apply((x*algae$sp0), 2, sd))[,1]
	nopts = list(location = loc, breadth = bre, scale_c = 1e-4, scale_e = 1e-6)
	mc = metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, niche_args = nopts,
		sp_names = algae$niches$species, r_names = "N", comp_scale = 0, 
		r_lim = apply(algae$r0, 2, range)[,1])
	i = c(5, 11, 13, 3)
	fl = flume(mc, algae$network, algae$sp0, algae$r0[,1, drop = FALSE])
	cp_before = col_prob(fl$metacom, fl$networks[[1]], fl$dt)
	ep_before = ext_prob(fl$metacom, fl$networks[[1]], fl$dt)
	discharge(fl$networks[[1]])[i] = 0
	cp_after <- col_prob(fl$metacom, fl$networks[[1]], fl$dt)
	ep_after = ext_prob(fl$metacom, fl$networks[[1]], fl$dt)

	expect_true(!all(cp_before[i,] == 0))
	expect_true(!all(ep_before[i,] == 1))
	expect_true(all(cp_after[i,] == 0))
	expect_true(all(ep_after[i,] == 1))
})
