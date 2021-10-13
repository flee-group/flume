test_that("Resource use functions", {
	# default 2-species community with uniformly spread niches
	comm = metacommunity()

	# simple 4-node river network
	# sites 2 and 3 are optimal for species 1 and 2
	Q = rep(1, 4)
	adj = matrix(0, nrow = 4, ncol = 4)
	adj[1,2] = adj[2,3] = adj[4,3] = 1
	st = matrix(seq(0, 1, length.out = length(Q)), ncol = 1, dimnames = list(NULL, 'R'))
	rn = river_network(adj)
	state(rn, "resources") = st
	state(rn, "species") = matrix(1, nrow = length(Q), ncol = length(comm$species))
	
	expect_error(ru <- ruf(state(rn, "species"), state(rn, "resources"), comm), regex=NA)

	# overall by default we expect all resource uses to be negative or zero
	expect_true(all(colSums(ru) <= 0))

	# rate of change per unit of resource should be greatest at the niche optimum
	ru_norm = ru / state(rn, "resources")
	expect_gt(abs(ru_norm[2,1]), abs(ru_norm[4,1]))

	# make sure that the presence-absence matrix actually matters
	state(rn, "species") = matrix(c(1,0,0,1,1,1,1,0), nrow = length(Q), ncol = length(comm$species))
	ru2 = ruf(state(rn, "species"), state(rn, "resources"), comm)
	expect_lt(sum(-1*ru2), sum(-1*ru)) #when all species are present, resource use is higher

})
