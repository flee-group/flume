test_that("Species flux works", {
	comm = metacommunity()
	Q = c(13,12,8,7,1,4,1,1,2,1)
	nsites = length(Q)
	adj = matrix(0, nrow=nsites, ncol=nsites)
	adj[2,1] = adj[3,2] = adj[9,2] = adj[4,3] = adj[5,4] = adj[6,4] =
		adj[7,6] = adj[8,6] = adj[10,9] = 1

	network = river_network(adjacency = adj, discharge = Q)
	state(network) = matrix(seq(0, 1, length.out = nsites), ncol=1)
	boundary(network) = state(network)
	site_by_species(network) = community_random(network, comm, prevalence = c(0.35, 0.65))
	boundary_species(network) = site_by_species(network) * 0
	expect_error(cp <- col_prob(comm, network, dt=1), regex=NA)
	expect_error(ep <- ext_prob(comm, network, dt=1), regex=NA)

	expect_error(Rflux <- dRdt(comm, network), regex=NA)
})
