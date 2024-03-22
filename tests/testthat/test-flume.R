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

	# test creation of using abundance as initial values
	expect_warning(sim2 <- flume(comm, network, st0 = R, sp0 = 2*S, stb = R, spb = S), regex="abundance")
	expect_identical(state(sim2$networks[[1]], "species"), S)
	
	## errors
	## missing initial state
	nt2 = network
	state(nt2, "resources") = NULL
	expect_error(flume(comm, nt2), regex="Initial resources state")
	nt2 = network
	state(nt2, "species") = NULL
	expect_error(flume(comm, nt2), regex="Initial species state")

	expect_error(sim <- run_simulation(sim, 1), regex=NA)

})

test_that("Create model with special resource types", {
	algae = readRDS(system.file("testdata/algae-test.rds", package="flume"))
	
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

test_that("Running a model with variable discharge", {
	comm = readRDS(system.file("testdata/metacom.rds", package = "flume"))
	network = readRDS(system.file("testdata/river_network.rds", package = "flume"))

	discharge(network) = cbind(discharge(network), 2*discharge(network), 1.5*discharge(network))
	sim = flume(comm, network)
	expect_error(sim <- run_simulation(sim, nt = 2* ncol(discharge(network))), regex = NA)

	# make sure that state reports the proper history
	### TODO
	## Currently fails, see issue #25
	expect_equal(ncol(state(sim$networks[[1]], "Q", history = TRUE)), 1+2 * ncol(discharge(network)))
})

