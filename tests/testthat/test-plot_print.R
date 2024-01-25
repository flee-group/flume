## workflow:
# test(filter="plot_print")
## if tests fail
# testthat::snapshot_review()

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

test_that("Plot resource concentration over time", {
	pl = function() plot(sim, type = 'resources')
	vdiffr::expect_doppelganger("Sim Resource Plot", pl)
})

test_that("Occupancy Plotting", {
 	pl = function() plot(sim, type = "occupancy")
 	vdiffr::expect_doppelganger("Sim Occupancy Plot", pl)

})

test_that("BEF Plot", {
	pl = function() plot(sim, type = "bef")
	vdiffr::expect_doppelganger("Sim BEF Plot", pl)
	
})
