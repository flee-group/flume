context("Plot methods")

test_that("Species plotting", {
	sp <- species('linear', 'constant', list(a = 0, b = 1), list(scale=0.2))
	pl <- function() plot(sp)
	vdiffr::expect_doppelganger("Linear/Constant Species Plot", pl)
})

test_that("Community plotting", {
	comm <- community()
	pl <- function() plot(comm)
	vdiffr::expect_doppelganger("Default Species Pool Plot", pl)
})
