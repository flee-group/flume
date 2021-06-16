test_that("Model creation", {
	comm = readRDS(system.file("testdata/metacom.rds", package = "flume"))
	network = readRDS(system.file("testdata/river_network.rds", package = "flume"))

	R = state(network)
	sbysp = site_by_species(network)


	# create a model with defaults
	expect_error(sim <- flume(comm, network), regex=NA)

	# create a model specifying states, but with default boundaries
	expect_error(flume(comm, network, st0 = R, sp0 = sbysp), regex=NA)

	# create a model with everything specified
	expect_error(flume(comm, network, st0 = R, sp0 = sbysp, stb = R, spb = sbysp), regex=NA)

	R = state(sim$network[[1]])
	sbysp = site_by_species(sim$network[[1]])

	## errors
	## missing initial state
	nt2 = network
	nt2[['.state']] = NULL
	expect_error(flume(comm, nt2), regex="Initial resource state")
	nt2 = network
	nt2[['si_by_sp']] = NULL
	expect_error(flume(comm, nt2), regex="Site by species")


	expect_error(sim <- run_simulation(sim, 1), regex=NA)

	expect_error(occ_nt <- occupancy(sim), regex=NA)
	expect_true(all(occ_nt$occupancy >= 0 & occ_nt$occupancy <= 1))
	expect_identical(length(unique(occ_nt$species)), length(comm$species))
	expect_true(all(range(unique(occ_nt$time)) == c(1,2)))
	expect_equal(occ_nt[time == 1]$occupancy, colSums(sbysp)/nrow(sbysp), check.names = FALSE)

	expect_error(occ_re <- occupancy(sim, "reach"), regex=NA)

	expect_error(res <- resource_summary(sim), regex=NA)
	expect_equal(res[time == 1]$concentration, R[,1], check.names = FALSE)

})



