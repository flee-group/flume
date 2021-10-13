## workflow:
# test(filter="plot_print")
## if tests fail
# testthat::snapshot_review()

context("Plotting")
comm = readRDS(system.file("testdata/metacom.rds", package="flume"))
rn = readRDS(system.file("testdata/river_network.rds", package="flume"))
sim = readRDS(system.file("testdata/sim.rds", package="flume"))

test_that("Species plotting", {
	sp = comm$species[[1]]
	vdiffr::expect_doppelganger("Linear/Constant Species Plot", plot(sp))
})

test_that("Community plotting", {
	vdiffr::expect_doppelganger("Default Species Pool Plot", plot(comm))
})


test_that("River Network plotting", {
	vdiffr::expect_doppelganger("River Network Plot", plot(rn))
})

test_that("River Network community plotting", {
	vdiffr::expect_doppelganger("River Network Species Plot", 
		plot(rn, variable = 'species'))
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
