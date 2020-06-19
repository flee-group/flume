context("Plot methods")

test_that("Species plotting", {
	sp = species('linear', 'constant', list(a = 0, b = 1), list(scale=0.2))
	pl = function() plot(sp)
	vdiffr::expect_doppelganger("Linear/Constant Species Plot", pl)
})

test_that("Community plotting", {
	comm = metacommunity()
	pl = function() plot(comm)
	vdiffr::expect_doppelganger("Default Species Pool Plot", pl)
})

test_that("River Network plotting", {
	adj = matrix(0, nrow = 4, ncol = 4)
	adj[1,2] = adj[2,3] = adj[4,3] = 1
	Q = rep(1, 4)
	rn = river_network(adj, Q)
	layout = matrix(c(-1,1, -0.5, 0.5, 0, 0, 0.5,0.5), ncol=2, byrow=TRUE)
	pl = function() plot(rn, layout = layout)
	vdiffr::expect_doppelganger("River Network Plot", pl)
})
