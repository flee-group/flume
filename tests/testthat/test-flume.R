test_that("Model creation", {
	sp_pool = readRDS(system.file("testdata/metacom.rds", package="flume"))
	network = readRDS(system.file("testdata/river_network.rds", package="flume"))

	sbysp = site_by_species(network)
	R = state(network)

	expect_error(sim <- flume(sp_pool, network, 1), regex=NA)
	expect_error(sim <- run_simulation(sim, 1), regex=NA)

	expect_error(occ_nt <- occupancy(sim), regex=NA)
	expect_true(all(occ_nt$occupancy >= 0 & occ_nt$occupancy <= 1))
	expect_equal(length(unique(occ_nt$species)), 2)
	expect_true(all(range(unique(occ_nt$time)) == c(1,2)))
	expect_identical(occ_nt[time == 1]$occupancy, colSums(sbysp)/nrow(sbysp))

	expect_error(occ_re <- occupancy(sim, "reach"), regex=NA)

	expect_error(res <- resource_summary(sim), regex=NA)
	expect_identical(res[time == 1]$concentration, R[,1])

})



