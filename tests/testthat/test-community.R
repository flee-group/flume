## exported functions
# f_niche
# dispersal_params

test_that("Species creation", {

	# input validation
	# many parameters must not be negative
	par_c_default = list(location = 0, breadth = 1, scale = 1)
	par_c = par_c_default; par_c$breadth = -1
	par_e_default = list(scale = 1)
	par_e = par_e_default
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = 1), regex = "negative")
	par_c = par_c_default; par_c$scale = -1
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = 1), regex = "negative")
	par_c = par_c_default; par_e$scale = -1
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = 1), regex = "negative")
	par_e = par_e_default
	expect_error(species(par_c, par_e, alpha = -1, beta = 1, r_use = 1), regex = "negative")
	expect_error(species(par_c, par_e, alpha = 1, beta = -1, r_use = 1), regex = "negative")

	# some parameters may be negative or positive
	par_c = par_c_default
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = -1), regex = NA)
	par_c$location = -1
	expect_error(species(par_c, par_e,  alpha = 1, beta = 1, r_use = 1), regex = NA)


	# multivariate niches
	loc = c(1, 2)
	rsc = c(1, 2)
	rsc_wrong = 1:3
	bre = matrix(c(1, 0, 0, 1), nrow = 2)
	bre_wrong = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)
	par_c = list(location = loc, breadth = bre, scale = 1)

	# all dimensions must agree
	par_c$breadth = bre_wrong
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = rsc), 
		regex = "dimension mismatch")
	par_c$breadth = bre; par_c$location = 1
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = rsc), 
		regex = "dimension mismatch")
	par_c$location = loc; par_c$breadth = bre[1,]
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = rsc_wrong), 
		regex = "dimension mismatch")

	# also check that no entries in vcv matrix are negative
	par_c$breadth = -1 * bre
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = rsc), regex = "negative")

	# scale_c, scale_e, alpha and beta must be single values
	par_c = list(location = loc, breadth = bre, scale = c(1, 2))
	expect_error(species(par_c, par_e, alpha = 1, beta = 1, r_use = rsc), regex = "single value")
	par_c$scale = 1
	expect_error(species(par_c, list(scale = c(1, 2)), alpha = 1, beta = 1, r_use = rsc), 
		regex = "single value")
	expect_error(species(par_c, par_e, alpha = c(1, 0), beta = 1, r_use = rsc), 
		regex = "single value")
	expect_error(species(par_c, par_e, alpha = 1, beta = c(1, 0), r_use = rsc), 
		regex = "single value")
})

test_that("Metacommunity creation", {

	# dimension checking
	# to check
	# location
	nargs = list(
		location = c(1, 2),
		breadth = c(1, 1.5),
		scale_c = c(0.5, 0.7),
		scale_e = c(0.2, 0.3),
		r_use = c(0.5, 0.4),
		alpha = c(0.05, 0.07),
		beta = c(0.4, 0.3)
	)
	r_lim = c(0, 1)

	# single species, single niche axis
	expect_error(metacommunity(nsp = 1, nr = 1, niches = niches_custom,
		niche_args = list(location = 1)), regex = NA)

	# multiple species, single niche axis
	expect_error(metacommunity(nsp = 2, nr = 1, niches = niches_custom,
		niche_args = list(location = nargs$location)), regex = NA)
	expect_error(metacommunity(nsp = 2, nr = 1, niches = niches_custom, niche_args = nargs, 
		r_lim = r_lim), regex = NA)
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
		r_use =  c(0.5, 0.4))

	nargs_32 = list(
		location = rbind(nargs_12$location, nargs_12$location + 0.25, nargs_12$location - 0.25),
		breadth = br_3_2_l,
		r_use = rbind(nargs_12$r_use, nargs_12$r_use, nargs_12$r_use)
	)
	# single species, multiple niche axes
	expect_error(metacommunity(nsp = 1, nr = 2, niches = niches_custom, niche_args = nargs_12, 
		r_lim = rbind(c(0, 1), c(0, 1))), regex = NA)

	# multiple species, multiple niche axes
	expect_error(metacommunity(nsp = 3, nr = 2, niches = niches_custom, niche_args = nargs_32),
		regex = NA)
})

test_that("Niches estimated correctly", {
	# data(algae, package = "flume")
	loc = c(0.25, 0.75)
	scale_c = 0.5
	scale_e = 0.2
	comm = metacommunity(niches = niches_custom, 
		niche_args = list(location = loc, scale_c = scale_c, scale_e = scale_e))
	expect_equal(attr(comm$species[[1]], "niche_max"), scale_c - scale_e, check.names = FALSE)
	expect_equal(as.vector(f_niche(comm$species[[1]], loc[1])), attr(comm$species[[1]], "niche_max"))

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
	data(algae)
	nopts = list(location = algae$niches$location, breadth = algae$niches$breadth, static = 1)
	expect_error(
		metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, niche_args = nopts),
		regex = NA)

	nopts = list(location = algae$niches$location, breadth = algae$niches$breadth)
	expect_error(metacommunity(nsp = nrow(algae$niches), nr = 2, niches = niches_custom, 
		niche_args = nopts, ratio = c(1, 2)), regex = NA)

	# various dimension errors
	expect_error(metacommunity(nsp = nrow(algae$niches), nr = 1, niches = niches_custom, 
			ratio = c(1, 2), niche_args = nopts), regex = "<= nr")

	nopts2 = nopts
	nopts2$breadth = rep(nopts2$breadth, 2)
	expect_error(
		metacommunity(nsp = nrow(algae$niches), nr = 2, niches = niches_custom, 
			niche_args = nopts2, ratio = c(1, 2)), regex = "Invalid niche breadth")
})

test_that("Asymmetric competition", {
	data(algae)
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

test_that("Niches/species with intermittency", {
	# try creating a community with invalid drying parameters
	n_params = list(dry_c = c(0, 0), dry_e = c(Inf, Inf, Inf))
	expect_error(mc <- metacommunity(niche_args = n_params), regex = "species")

	# create a community with defaults and test the niches when dry
	mc <- metacommunity()
	R = matrix(rep(0, 2), ncol = 1)
	# inappropriate dry length causes error
	dry = c(TRUE, FALSE, TRUE)
	expect_error(cn <- f_niche(mc, R, component = "col", dry = dry), regex = "length\\(dry\\)")

	dry = c(TRUE, FALSE)
	expect_error(ln <- f_niche(mc, R, dry = dry), regex = NA)
	expect_error(cn <- f_niche(mc, R, component = "col", dry = dry), regex = NA)
	expect_error(en <- f_niche(mc, R, component = "ext", dry = dry), regex = NA)
	expect_true(all(cn[dry,] == 0))
	expect_true(all(cn[!dry,] != 0))
	expect_true(all(is.infinite(en[dry,])))
	expect_true(all(is.finite(en[!dry,])))
	expect_true(all(is.infinite(ln[dry,])))
	expect_true(all(is.finite(ln[!dry,])))


	# create and test community with non defaults
	n_params = list(dry_c = c(0.5, 0.5), dry_e = c(2, 2))
	mc = metacommunity(niche_args = n_params)
	cn = f_niche(mc, R, component = "col", dry = dry)
	en = f_niche(mc, R, component = "ext", dry = dry)
	expect_equal(cn[dry,], cn[!dry,] * n_params$dry_c)
	expect_equal(en[dry,], en[!dry,] * n_params$dry_e)

})