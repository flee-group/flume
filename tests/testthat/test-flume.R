test_that("Model creation", {
	comm = readRDS(system.file("testdata/metacom.rds", package = "flume"))
	network = readRDS(system.file("testdata/river_network.rds", package = "flume"))

	R = state(network, "resources")
	S = state(network, "species")


	# create a model with defaults
	expect_error(sim <- flume(comm, network), regex=NA)

	# create a model specifying states, but with default boundaries
	expect_error(flume(comm, network, st0 = R, sp0 = S), regex=NA)

	# create a model with everything specified
	expect_error(flume(comm, network, st0 = R, sp0 = S, stb = R, spb = S), regex=NA)

	R = state(sim$network[[1]], "resources")
	S = state(sim$network[[1]], "species")

	## errors
	## missing initial state
	nt2 = network
	state(nt2, "resources") = NULL
	expect_error(flume(comm, nt2), regex="Initial resources state")
	nt2 = network
	state(nt2, "species") = NULL
	expect_error(flume(comm, nt2), regex="Initial species state")

	expect_error(sim <- run_simulation(sim, 1), regex=NA)

	expect_error(occ_nt <- occupancy(sim), regex=NA)
	expect_true(all(occ_nt$occupancy >= 0 & occ_nt$occupancy <= 1))
	expect_identical(length(unique(occ_nt$species)), length(comm$species))
	expect_true(all(range(unique(occ_nt$time)) == c(1,2)))
	expect_equal(occ_nt[time == 1]$occupancy, colSums(S)/nrow(S), check.names = FALSE)

	expect_error(occ_re <- occupancy(sim, "reach"), regex=NA)

	expect_error(res <- resource_summary(sim), regex=NA)
	expect_equal(res[time == 1]$concentration, R[,1], check.names = FALSE)

})

test_that("Create model with special resource types", {
	data(algae)

	# static
	loc = apply(algae$r0, 2, function(x) colMeans(x*algae$sp0))
	bre = apply(algae$r0, 2, function(x) apply((x*algae$sp0), 2, sd))
	nopts = list(location = loc, breadth = bre, static = 1, r_lim = t(apply(algae$r0, 2, range)))
	mc = metacommunity(nsp = nrow(algae$niches), nr = 2, niches = niches_custom, niche_args = nopts,
		sp_names = algae$niches$species, r_names = c("N", "P"))
	expect_error(fl <- flume(mc, algae$network, algae$sp0, algae$r0), regex=NA)
	expect_error(fl <- run_simulation(fl, 1), regex=NA)
	# first resource is static and shouldn't change, second should change
	expect_equal(state(fl$network[[1]], "resources")[,1], algae$r0[,1], check.attributes = FALSE) 
	expect_false(all(state(fl$network[[1]], "resources")[,2] == algae$r0[,2])) 

	# ratio
	expect_error(fl <- flume(algae$metacommunity, algae$network, algae$sp0, algae$r0), regex=NA)
	expect_error(fl <- run_simulation(fl, 1), regex=NA)
})

