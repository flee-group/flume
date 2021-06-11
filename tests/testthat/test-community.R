## exported functions
# f_niche
# dispersal_params

test_that("Species creation", {

	# input validation
	# many parameters must not be negative
	expect_error(species(location = 0, breadth = -1, scale_c = 1, scale_e = 1, alpha = 1, beta = 1,
		r_scale = 1), regex = "negative")
	expect_error(species(location = 0, breadth = 1, scale_c = -1, scale_e = 1, alpha = 1, beta = 1,
		r_scale = 1), regex = "negative")
	expect_error(species(location = 0, breadth = 1, scale_c = 1, scale_e = -1, alpha = 1, beta = 1,
		r_scale = 1), regex = "negative")
	expect_error(species(location = 0, breadth = 1, scale_c = 1, scale_e = 1, alpha = -1, beta = 1,
		r_scale = 1), regex = "negative")
	expect_error(species(location = 0, breadth = 1, scale_c = 1, scale_e = 1, alpha = 1, beta = -1,
		r_scale = 1), regex = "negative")

	# some parameters may be negative or positive
	expect_error(species(location = 1, breadth = 1, scale_c = 1, scale_e = 1, alpha = 1, beta = 1,
		r_scale = 1), regex = NA)
	expect_error(species(location = -1, breadth = 1, scale_c = 1, scale_e = 1, alpha = 1, beta = 1,
		r_scale = -1), regex = NA)


	# multivariate niches
	loc = c(1, 2)
	rsc = c(1, 2)
	rsc_wrong = 1:3
	bre = matrix(c(1, 0, 0, 1), nrow = 2)
	bre_wrong = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)

	# all dimensions must agree
	expect_error(species(location = loc, breadth = bre_wrong, scale_c = 1, scale_e = 1, alpha = 1,
		beta = 1, r_scale = rsc), regex = "dimension mismatch")
	expect_error(species(location = 1, breadth = bre, scale_c = 1, scale_e = 1, alpha = 1,
		beta = 1, r_scale = rsc), regex = "dimension mismatch")
	expect_error(species(location = loc, breadth = bre, scale_c = 1, scale_e = 1, alpha = 1,
		beta = 1, r_scale = rsc_wrong), regex = "dimension mismatch")
	expect_error(species(location = loc, breadth = bre[1, ], scale_c = 1, scale_e = 1, alpha = 1,
		beta = 1, r_scale = rsc_wrong), regex = "dimension mismatch")

	# also check that no entries in vcv matrix are negative
	expect_error(species(location = loc, breadth = -1 * bre, scale_c = 1, scale_e = 1, alpha = 1,
		beta = 1, r_scale = rsc), regex = "negative")

	# scale_c, scale_e, alpha and beta must be single values
	expect_error(species(location = loc, breadth = bre, scale_c = c(1, 2), scale_e = 1, alpha = 1,
		beta = 1, r_scale = rsc), regex = "single value")
	expect_error(species(location = loc, breadth = bre, scale_c = 1, scale_e = c(1, 2), alpha = 1,
		beta = 1, r_scale = rsc), regex = "single value")
	expect_error(species(location = loc, breadth = bre, scale_c = 1, scale_e = 1, alpha = c(1, 0),
		beta = 1, r_scale = rsc), regex = "single value")
	expect_error(species(location = loc, breadth = bre, scale_c = 1, scale_e = 1, alpha = 1,
		beta = c(1, 0), r_scale = rsc), regex = "single value")
})

test_that("Metacommunity creation", {

	# dimension checking
	# to check
	# location
	loc2 = c(1, 2)

	# breadth
	br2 = c(1, 1.5)

	# scale_c
	scc2 = c(0.5, 0.7)

	# scale_e
	sce2 = c(0.2, 0.3)

	# alpha
	al2 = c(0.05, 0.07)
	# beta
	be2 = c(0.4, 0.3)

	# r_scale
	rsc2 = c(0.5, 0.4)

	# niche_lim
	nl1 = c(0, 1)

	# single species, single niche axis
	expect_error(metacommunity(location = 1), regex = NA)
	expect_error(metacommunity(location = 1, breadth = br2), regex = "niche breadth")
	expect_error(metacommunity(location = 1, scale_c = scc2), regex = "Scale parameters")
	expect_error(metacommunity(location = 1, scale_e = sce2), regex = "Scale parameters")
	expect_error(metacommunity(location = 1, alpha = al2), regex = "alpha")
	expect_error(metacommunity(location = 1, beta = be2), regex = "beta")
	expect_error(metacommunity(location = 1, r_scale = rsc2), regex = "r_scale")
	expect_error(metacommunity(location = 1, niche_lim = rbind(nl1, nl1)), regex = "niche_lim")

	# multiple species, single niche axis
	expect_error(metacommunity(location = loc2), regex = NA)
	expect_error(metacommunity(location = loc2, breadth = br2, scale_c = scc2, scale_e = sce2,
			alpha = al2, beta = be2, r_scale = rsc2, niche_lim = nl1), regex = NA)
	expect_error(metacommunity(location = loc2, breadth = rep(br2, 2), scale_c = scc2,
			scale_e = sce2, alpha = al2, beta = be2, r_scale = rsc2, niche_lim = nl1),
			regex = "niche breadth")
	expect_error(metacommunity(location = loc2, breadth = br2, scale_c = scc2,
			scale_e = rep(sce2, 2), alpha = al2, beta = be2, r_scale = rsc2, niche_lim = nl1),
			regex = "Scale parameters")
	expect_error(metacommunity(location = loc2, breadth = br2, scale_c = scc2,
			scale_e = sce2, alpha = rep(al2, 2), beta = be2, r_scale = rsc2, niche_lim = nl1),
			regex = "alpha")
	expect_error(metacommunity(location = loc2, breadth = br2, scale_c = scc2,
			scale_e = sce2, alpha = al2, beta = rep(be2, 2), r_scale = rsc2, niche_lim = nl1),
			regex = "beta")
	expect_error(metacommunity(location = loc2, breadth = br2, scale_c = scc2,
			scale_e = sce2, alpha = al2, beta = be2, r_scale = rep(rsc2, 2), niche_lim = nl1),
			regex = "r_scale")
	expect_error(metacommunity(location = loc2, breadth = br2, scale_c = scc2,
			scale_e = sce2, alpha = al2, beta = be2, r_scale = rsc2,
			niche_lim = rbind(nl1, nl1)), regex = "niche_lim")
})


test_that("Multivariate metacommunities", {
	skip_on_cran()
	skip_if_not(identical(Sys.getenv("SLOW_TESTS"), "TRUE"))

	# location
	loc_1_2 = matrix(c(1, 2), ncol = 2)
	loc_3_2 = rbind(loc_1_2, loc_1_2 + 0.25, loc_1_2 - 0.25)

	# breadth
	br_1_2 = c(1, 1.5)
	br_3_2_vsp = c(1, 2, 3)
	br_3_2_vr = c(1, 2)
	br_3_2_m = cbind(br_3_2_vsp, br_3_2_vsp)
	br_3_2_vcv = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
	br_3_2_l = list(br_3_2_vcv, br_3_2_vcv, br_3_2_vcv)

	# r_scale
	rsc_1_2 = c(0.5, 0.4)
	rsc_3_2 = rbind(rsc_1_2, rsc_1_2, rsc_1_2)

	# niche_lim
	nl_2 = rbind(c(0, 1), c(0, 1))

	# single species, multiple niche axes
	expect_error(metacommunity(location = loc_1_2, niche_lim = nl_2), regex = NA)
	expect_error(metacommunity(location = loc_1_2, breadth = br_1_2, r_scale = rsc_1_2,
		niche_lim = nl_2), regex = NA)
	expect_error(metacommunity(location = loc_1_2, breadth = rep(br_1_2, 2),
		niche_lim = nl_2), regex = "niche breadth")
	expect_error(metacommunity(location = loc_1_2, r_scale = rep(rsc_1_2, 2),
		niche_lim = nl_2), regex = "r_scale")

	# multiple species, multiple niche axes
	expect_error(metacommunity(location = loc_3_2, niche_lim = nl_2), regex = NA)
	expect_error(metacommunity(location = loc_3_2, breadth = br_3_2_vsp, niche_lim = nl_2),
		regex = NA)
	expect_error(metacommunity(location = loc_3_2, breadth = br_3_2_vr, niche_lim = nl_2),
		regex = NA)
	expect_error(metacommunity(location = loc_3_2, breadth = br_3_2_m, niche_lim = nl_2),
		regex = NA)
	expect_error(metacommunity(location = loc_3_2, breadth = br_3_2_l, niche_lim = nl_2),
		regex = NA)
	expect_error(metacommunity(location = loc_3_2, r_scale = rsc_3_2, niche_lim = nl_2), regex = NA)

	expect_error(metacommunity(location = loc_3_2, breadth = rep(br_3_2_vsp, 2), niche_lim = nl_2),
		regex = "niche breadth")
	expect_error(metacommunity(location = loc_3_2, breadth = rbind(br_3_2_m, br_3_2_m),
		niche_lim = nl_2), regex = "niche breadth")
	expect_error(metacommunity(location = loc_3_2, breadth = c(br_3_2_l, br_3_2_l),
		niche_lim = nl_2), regex = "niche breadth")

	expect_error(metacommunity(location = loc_3_2, r_scale = rsc_3_2[1:2, ], niche_lim = nl_2),
		regex = "r_scale")

})

test_that("Niches estimated correctly", {
	# data(algae, package = "flume")
	loc = c(0.25, 0.75)
	scale_c = 0.5
	scale_e = 0.2
	comm = metacommunity(location = loc, scale_c = scale_c, scale_e = scale_e)
	expect_equal(attr(comm$species[[1]], "niche_max"), scale_c - scale_e, check.names = FALSE)
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
	expect_identical(c_niche[, 1, drop = FALSE], sp_niche)
})
