comm = readRDS(system.file("testdata/metacom.rds", package = "iflume"))
network = readRDS(system.file("testdata/river_network.rds", package = "iflume"))

## note best to add a boundary condition for species so nothing goes globally extinct
## for this sim, would be good if it was high enough to ensure some colonisation no matter what
bnd = matrix(0.5, ncol = attr(comm,"n_species"), nrow = attr(network,"n_sites"), 
	dimnames = list(attr(network,"names_sites"), attr(comm,"sp_names")))

# run the model for some time steps, with a period of intermittency
nt = 9
Q = matrix(rep(state(network, "Q"), nt), ncol = nt)
inodes = c(4, 6, 7, 10)
itime = 4:6
Q[inodes, itime] = 0



test_that("Sim behaves itself during intermittent phase", {
	sim = flume(comm, network, spb = bnd)
	## run until the first intermittent timestep
	i = 1
	while(all(Q[,i] > 0)) {
		discharge(sim$networks[[1]]) = Q[,i]
		sim = run_simulation(sim, 1)
		i = i + 1
	}
	# and a single intermittent timestep
	discharge(sim$networks[[1]]) = Q[,i]
	sim = run_simulation(sim, 1)
	i = i + 1

	## col prob is zero, ext prob is 1 when discharge = 0, all values are finite
	cp = col_prob(sim$metacom, sim$networks[[1]], sim$dt)
	ep = ext_prob(sim$metacom, sim$networks[[1]], sim$dt)
	expect_true(all(cp[inodes,] == 0))
	expect_true(all(ep[inodes,] == 1))
	expect_true(all(is.finite(cp)))
	expect_true(all(is.finite(ep)))

	# repeat tests after a second intermittent timestep
	discharge(sim$networks[[1]]) = Q[,i]
	sim = run_simulation(sim, 1)
	i = i + 1
	cp = col_prob(sim$metacom, sim$networks[[1]], sim$dt)
	ep = ext_prob(sim$metacom, sim$networks[[1]], sim$dt)
	expect_true(all(cp[inodes,] == 0))
	expect_true(all(ep[inodes,] == 1))
	expect_true(all(is.finite(cp)))
	expect_true(all(is.finite(ep)))

	# finish running to the end to test recovery
	for(j in i:ncol(Q)) {
		discharge(sim$networks[[1]]) = Q[,j]
		sim = run_simulation(sim, 1)		
	}
	cp = col_prob(sim$metacom, sim$networks[[1]], sim$dt)
	ep = ext_prob(sim$metacom, sim$networks[[1]], sim$dt)
	expect_true(all(cp[inodes,] > 0))
	expect_true(all(ep[inodes,] < 1))   ## failing due to bug in test data
	expect_true(all(is.finite(cp)))
	expect_true(all(is.finite(ep)))

	## zero discharge is correctly recorded in history
	expect_equal(discharge(sim$networks[[1]], "history"), Q)    # failing, enhancements needed in main branch

})


test_that("Full sim runs with intermittency and recovery", {

	sim = flume(comm, network, spb = bnd)
	discharge(sim$networks[[1]]) = Q

	## sim runs with intermittency
	expect_warning(sim_r <- run_simulation(sim, nt), regex = NA)

	## all species go extinct when discharge = 0
	# occ = summarise(sim_r, c("time", "occupancy")
	sphist = state(sim_r$networks[[1]], "species", TRUE)
	## itime + 2 because we add one for the initial state (time = 0), add one because
	## the extinctions should happen in the time step after discharge drops
	sptotal_i = Reduce(`+`, sphist[itime + 2])
	expect_equal(sptotal_i[inodes,], sptotal_i[inodes,]*0) # fails, species are present; this test should solve itself after CP and EP are resolved
	# do not expect global extinctions
	expect_gt(sum(sptotal_i), 0)

	## reaction and transport should also go to zero when discharge is zero
	rxn = do.call(cbind, state(sim_r$networks[[1]], "reaction", history = TRUE))
	tr = do.call(cbind, state(sim_r$networks[[1]], "transport", history = TRUE))
	expect_true(all(rxn[inodes, itime + 1] == 0))
	expect_true(all(tr[inodes, itime + 1] == 0))
	expect_true(all(rxn[-inodes, itime + 1] <= 0))
	expect_true(all(tr[-inodes, itime + 1] > 0))
	# and recover when flow resumes
	expect_lt(sum(rxn[inodes, nt + 1]), 0)
	expect_gt(sum(tr[inodes, nt + 1]), 0)

	## resource concentration should be stored as zero (not quite accurate, but prob most useful)
	rhist = do.call(cbind, state(sim_r$networks[[1]], "resources", TRUE))
	expect_equal(sum(rhist[inodes, itime + 1]), 0) # fails, concentrations not zero
	expect_gt(sum(rhist[-inodes, itime + 1]), 0)
	## after resumption of flow we get reasources back in the affected nodes
	expect_true(all(rhist[inodes, ncol(rhist)] > 0))

})




