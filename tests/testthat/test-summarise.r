test_that("EF Summary", {
	sim = readRDS(system.file("testdata/sim.rds", package = "flume"))
	expect_error(summarise(sim, "species", "EF"), "EF by species")
	
	# try to summarise concentration without resources in by, or with species in by
	expect_error(summarise(sim, stat = "concentration"), "must contain 'resources'")
	expect_error(summarise(sim, by = c('species', 'time'), stat = "concentration"), "must not contain 'species'")
	

	expect_error(occ_nt <- summarise(sim, c("species", "time"), "occupancy"), regex=NA)
	expect_true(all(occ_nt$occupancy >= 0 & occ_nt$occupancy <= 1))
	expect_identical(length(unique(occ_nt$species)), length(sim$metacom$species))
	expect_equal(range(unique(occ_nt$time)), c(1,51))

	R = state(sim$network[[1]], "resources")
	S = state(sim$network[[1]], "species")
	maxt = length(state(sim$network[[1]], "species", TRUE))
	expect_equal(occ_nt[time == maxt]$occupancy, unname(colSums(S)/nrow(S)), check.names = FALSE)

	expect_error(occ_re <- summarise(sim, c("species", "reach"), "occupancy"), regex=NA)
	
	# slicing time steps
	expect_error(summarise(sim, stat = "richness", t_steps = 50:55), "simulation bounds")
	t_steps = 40:50
	expect_error(summ <- summarise(sim, stat = "richness", t_steps = t_steps), NA)
	expect_true(all(unique(summ$time) %in% t_steps))
})
