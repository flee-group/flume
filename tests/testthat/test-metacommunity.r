## exported functions
# f_niche
# dispersal_params



test_that("Metacommunity creation", {

	# dimension checking
	# to check
	# location
	nargs = list(
		location = c(1, 2),
		breadth = c(1, 1.5),
		scale = list(col = c(0.5, 0.7), ext = c(0.2, 0.3)),
		r_use = c(0.5, 0.4),
		r_lim = c(0, 1)
	)

	dargs = list(
		alpha = c(0.05, 0.07),
		beta = c(0.4, 0.3)
	)

	# single species, single niche axis
	expect_error(metacommunity(nsp = 1, nr = 1, niches = niches_custom,
		niche_args = list(location = 1)), regex = NA)

	# multiple species, single niche axis
	expect_error(metacommunity(nsp = 2, nr = 1, niches = niches_custom,
		niche_args = list(location = nargs$location)), regex = NA)
	expect_error(metacommunity(nsp = 2, nr = 1, niches = niches_custom, niche_args = nargs, 
		dispersal_args = dargs), regex = NA)
})


test_that("Multivariate metacommunities", {
	skip_on_cran()
	skip_if_not(identical(Sys.getenv("SLOW_TESTS"), "TRUE"))


	# breadth
	br_3_2_vcv = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
	br_3_2_l = list(br_3_2_vcv, br_3_2_vcv, br_3_2_vcv)

	nargs_12 = list(
		location = matrix(c(1, 2), ncol = 2),
		breadth = c(1, 1.5),
		r_use =  c(0.5, 0.4),
		r_lim = rbind(c(0, 1), c(0, 1))
	)

	nargs_32 = list(
		location = rbind(nargs_12$location, nargs_12$location + 0.25, nargs_12$location - 0.25),
		breadth = br_3_2_l,
		r_use = rbind(nargs_12$r_use, nargs_12$r_use, nargs_12$r_use)
	)
	# single species, multiple niche axes
	expect_error(metacommunity(nsp = 1, nr = 2, niches = niches_custom, niche_args = nargs_12), 
		regex = NA)

	# multiple species, multiple niche axes
	expect_error(metacommunity(nsp = 3, nr = 2, niches = niches_custom, niche_args = nargs_32),
		regex = NA)
})

test_that("Niches estimated correctly", {
	# data(algae, package = "flume")
	loc = c(0.25, 0.75)
	scale = list(col = 0.5, ext = 0.2)
	comm = metacommunity(niches = niches_custom, 
		niche_args = list(location = loc, scale = scale))
	expect_equal(attr(comm$species[[1]], "niche_max"), scale$col - scale$ext, check.names = FALSE)
	expect_equal(f_niche(comm$species[[1]], loc[1]), attr(comm$species[[1]], "niche_max"))

	# species 2 performs worse at species 1's optimum
	expect_gt(f_niche(comm$species[[1]], loc[1]), f_niche(comm$species[[2]], loc[1]))

	st = matrix(seq(0, 1, 0.05), ncol = 1)
	expect_error(sp_niche <- f_niche(comm$species[[1]], st), regex = NA)
	expect_true(is(sp_niche, "matrix"))
	expect_equal(dim(sp_niche), c(nrow(st), 1))
	expect_identical(sp_niche, comm$species[[1]]$col(st) - comm$species[[1]]$ext(st))

	expect_error(c_niche <- f_niche(comm, st), regex = NA)
	expect_true(is(c_niche, "matrix"))
	expect_equal(dim(c_niche), c(nrow(st), length(comm$species)))
	expect_equal(c_niche[, 1, drop = FALSE], sp_niche, check.attributes = FALSE)
})

test_that("Ratio and static niches", {
	algae = readRDS(system.file("testdata/algae-test.rds", package="flume"))
	nopts = list(location = algae$niches$location, breadth = algae$niches$breadth, static = 1)
	expect_error(
		metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, niche_args = nopts),
		regex = NA)

	nopts = list(location = algae$niches$location, breadth = algae$niches$breadth, ratio = c(1, 2))
	expect_error(
		metacommunity(nsp = nrow(algae$niches), nr = 2, niches = niches_custom, niche_args = nopts),
		regex = NA)

	# various dimension errors
	expect_error(
		metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
			niche_args = nopts), regex = "<= nr")

	nopts2 = nopts
	nopts2$breadth = rep(nopts2$breadth, 2)
	expect_error(
		metacommunity(nsp = nrow(algae$niches), nr = 2, niches = niches_custom, 
			niche_args = nopts2), regex = "Invalid niche breadth")
})

test_that("Asymmetric competition", {
	algae = readRDS(system.file("testdata/algae-test.rds", package="flume"))
	nopts = list(location = algae$niches$location, breadth = algae$niches$breadth)

	# simple scale with all species the same
	expect_error(comm1 <- metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
		niche_args = nopts, comp_scale = 1), regex = NA)
	comm2 <- metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
		niche_args = nopts, comp_scale = 2)
	expect_true(all(comm2$competition >= comm1$competition))
	expect_equal(comm2$competition/2, comm1$competition)

	# change scale per species
	expect_error(comm3 <- metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
		niche_args = nopts, comp_scale = c(1, 2, 1, 1)), regex = "nsp")
	expect_error(comm3 <- metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
		niche_args = nopts, comp_scale = c(1, 2, 1, 1, 1)), regex = NA)
	expect_equal(comm3$competition[2,]/2, comm1$competition[2,])

	# asymmetric pairwise
	comp = matrix(1, nrow = length(nopts$location), ncol = length(nopts$location))
	comp[2,1] = 2
	expect_error(comm4 <- metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
		niche_args = nopts, comp_scale = comp), regex = NA)
	expect_equal(comm4$competition[2,1], comm3$competition[2,1])
	expect_equal(comm4$competition[1,2], comm1$competition[1,2])

	# species 2 facilitates species 1, but not the other way around
	comp[2,1] = -0.2
	expect_error(comm5 <- metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
		niche_args = nopts, comp_scale = comp), regex = NA)
	expect_lt(comm5$competition[2,1], 0)
	expect_gt(comm5$competition[1,2], 0)

})
