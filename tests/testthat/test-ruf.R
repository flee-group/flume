test_that("Resource use functions", {
	comm = metacommunity(r_scale = matrix(c(0.05, 0.1), nrow=2), c_type='linear')
	Q = rep(1, 4)
	adj = matrix(0, nrow = 4, ncol = 4)
	adj[1,2] = adj[2,3] = adj[4,3] = 1
	st = matrix(seq(0, 1, length.out = length(Q)), ncol = 1, dimnames = list(NULL, 'R'))
	rn = river_network(adj, Q, state = st)
	site_by_species(rn) = random_community(rn, comm)

	expect_error(ru <- ruf(site_by_species(rn), state(rn), comm), regex=NA)

	# overall by default we expect all resource uses to be negative or zero
	expect_true(all(colSums(ru) <= 0))

	# absolute rate of change should be greatest when a greedy species is at it's optimum
	expect_gt(abs(ru[1,1]), abs(ru[4,1]))
})
