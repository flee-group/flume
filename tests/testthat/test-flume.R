test_that("Model creation", {
	sp_pool = metacommunity(n_species = 2, nx = 1, c_type = 'linear', e_type='constant')
	Q = c(13,12,8,7,1,4,1,1,2,1)
	nsites = length(Q)
	adj = matrix(0, nrow=nsites, ncol=nsites)
	adj[2,1] = adj[3,2] = adj[9,2] = adj[4,3] = adj[5,4] = adj[6,4] = adj[7,6] = adj[8,6] = adj[10,9] = 1

	network = river_network(adjacency = adj, discharge = Q)
	network$boundary = function() matrix(seq(0, 1, length.out = nsites), ncol=1, dimnames = list(NULL, 'R'))
	state(network) = network$boundary()
	site_by_species(network) = random_community(network, sp_pool, prevalence = c(0.35, 0.65))

	expect_error(sim <- flume(sp_pool, network, 1), regex=NA)
	expect_error(sim <- run_simulation(sim, 1), regex=NA)

})
