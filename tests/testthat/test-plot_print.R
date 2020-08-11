context("Plotting")
comm = readRDS(system.file("testdata/metacom.rds", package="flume"))
rn = readRDS(system.file("testdata/river_network.rds", package="flume"))

test_that("Species plotting", {
	sp = comm$species[[1]]
	pl = function() plot(sp)
	vdiffr::expect_doppelganger("Linear/Constant Species Plot", pl)
})

test_that("Community plotting", {
	pl = function() plot(comm)
	vdiffr::expect_doppelganger("Default Species Pool Plot", pl)
})


test_that("River Network plotting", {
	pl = function() plot(rn)
	vdiffr::expect_doppelganger("River Network Plot", pl)
})

test_that("River Network community plotting", {
	pl = function() plot(rn, variable = 'site_by_species')
	vdiffr::expect_doppelganger("River Network Species Plot", pl)
})
