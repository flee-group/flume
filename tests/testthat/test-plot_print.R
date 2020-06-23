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

adj = matrix(0, nrow = 4, ncol = 4)
adj[1,2] = adj[2,3] = adj[4,3] = 1
Q = rep(1, 4)
rn = river_network(adj, Q)
l_out = matrix(c(-1,1, -0.5, 0.5, 0, 0, 0.5,0.5), ncol=2, byrow=TRUE)

test_that("River Network plotting", {

	pl = function() plot(rn, layout = l_out)
	vdiffr::expect_doppelganger("River Network Plot", pl)

	## add state and plot
	st = matrix(seq(0, 1, length.out = length(Q)), ncol = 1, dimnames = list(NULL, 'R'))
	state(rn) = st
	vdiffr::expect_doppelganger("River Network State Plot", pl)

})

test_that("River Network community plotting", {
	comm <- metacommunity()
	site_by_species(rn) = matrix(c(1,1,1,0,0,0,1,1), ncol=length(comm))
	pl = function() plot(rn, layout = l_out, variable = 'site_by_species')
	vdiffr::expect_doppelganger("River Network Species Plot", pl)
})
