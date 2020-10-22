context("Plotting")
comm = readRDS(system.file("testdata/metacom.rds", package="flume"))
rn = readRDS(system.file("testdata/river_network.rds", package="flume"))
sim = readRDS(system.file("testdata/sim.rds", package="flume"))

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

## ggplot figs not working with vdiffr for some reason
# test_that("Occupancy Plotting", {
#  	pl = function() plot(sim, variable = "occupancy")
#  	vdiffr::expect_doppelganger("Sim Occupancy Plot", pl)
#
# })
# test_that("Resource Plotting", {
# 	pl = function() plot(sim, variable = "resources")
# 	vdiffr::expect_doppelganger("Sim Resource Plot", pl)
# })
