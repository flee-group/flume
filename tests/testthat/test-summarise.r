test_that("EF Summary", {
	sim = readRDS(system.file("testdata/sim.rds", package = "flume"))
	expect_error(summarise(sim, "species", "EF"), "EF by species")
	
	# try to summarise concentration without resources in by, or with species in by
	expect_error(summarise(sim, stat = "concentration"), "must contain 'resources'")
	expect_error(summarise(sim, by = c('species', 'time'), stat = "concentration"), "must not contain 'species'")
})
