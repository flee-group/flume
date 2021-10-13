test_that("Custom niches", {

	## exported functions
	# 


	## input validation
	location = c(1, 2)
	breadth = c(1, 1.5)
	scale_c = c(0.5, 0.7)
	scale_e = c(0.2, 0.3)
	r_use = c(0.5, 0.4)
	r_lim = c(0, 1)

	# single species, single resource
	expect_error(niches_custom(nsp = 1, nr = 1, location = location[1]), regex=NA)
	expect_error(niches_custom(nsp = 1, nr = 1, location = location[1], breadth = breadth[1], 
		scale_c = scale_c[1], scale_e = scale_e[1], r_use = r_use[1], r_lim = r_lim), regex=NA)
	expect_error(niches_custom(nsp = 1, nr = 1, location = 1, breadth = breadth),
		regex = "niche breadth")
	expect_error(niches_custom(nsp = 1, nr = 1, location = 1, scale_c = scale_c),
		regex = "Scale parameters")
	expect_error(niches_custom(nsp = 1, nr = 1, location = 1, scale_e = scale_e),
		regex = "Scale parameters")
	expect_error(niches_custom(nsp = 1, nr = 1, location = 1, r_use = r_use),
		regex = "r_use")
	expect_error(niches_custom(nsp = 1, nr = 1, location = 1, r_lim = rbind(r_lim, r_lim)),
		regex = "r_lim")

	# multiple species, single resource
	expect_error(niches_custom(nsp = 2, nr = 1, location = location), regex=NA)
	expect_error(niches_custom(nsp = 2, nr = 1, location = location, breadth = breadth, 
		scale_c = scale_c, scale_e = scale_e, r_use = r_use, r_lim = r_lim), regex=NA)
	expect_error(niches_custom(nsp = 2, nr = 1, location = location, breadth = rep(breadth, 2)),
		regex = "niche breadth")
	expect_error(niches_custom(nsp = 2, nr = 1, location = location, scale_c = rep(scale_c, 2)),
		regex = "Scale parameters")
	expect_error(niches_custom(nsp = 2, nr = 1, location = location, scale_e = rep(scale_e, 2)),
		regex = "Scale parameters")
	expect_error(niches_custom(nsp = 2, nr = 1, location = location, r_use = rep(r_use, 2)),
		regex = "r_use")

	# single species, multiple niche axes
	location = matrix(c(1, 2), ncol = 2)
	breadth = c(1, 1.5)
	r_use =  c(0.5, 0.4)
	r_lim = rbind(c(0, 1), c(0, 1))
	expect_error(niches_custom(nsp = 1, nr = 2, location = location, breadth = breadth,
		r_use = r_use, r_lim = r_lim), regex = NA)
	expect_error(niches_custom(nsp = 1, nr = 2, location = location, breadth = rep(breadth, 2)), 
		regex = "niche breadth")
	expect_error(niches_custom(nsp = 1, nr = 2, location = location, r_use = rep(r_use, 2)), 
		regex = "r_use")

	# multiple species, multiple niche axes
	location = rbind(location, location + 0.25, location - 0.25)
	r_use = rbind(r_use, r_use, r_use)
	br_3_2_vsp = c(1, 2, 3)
	br_3_2_vr = c(1, 2)
	br_3_2_m = cbind(br_3_2_vsp, br_3_2_vsp)
	br_3_2_vcv = matrix(c(1, 0.3, 0.3, 1), nrow = 2)
	br_3_2_l = list(br_3_2_vcv, br_3_2_vcv, br_3_2_vcv)

	expect_error(niches_custom(nsp = 3, nr = 2, location = location, breadth = br_3_2_vsp), 
		regex = NA)
	expect_error(niches_custom(nsp = 3, nr = 2, location = location, breadth = br_3_2_vr),
		regex = NA)
	expect_error(niches_custom(nsp = 3, nr = 2, location = location, breadth = br_3_2_m),
		regex = NA)
	expect_error(niches_custom(nsp = 3, nr = 2, location = location, breadth = br_3_2_l),
		regex = NA)

	expect_error(niches_custom(nsp = 3, nr = 2, location = location, breadth = rep(br_3_2_vsp, 2)),
		regex = "niche breadth")
	expect_error(niches_custom(nsp = 3, nr = 2, location = location, 
		breadth = rbind(br_3_2_m, br_3_2_m)), regex = "niche breadth")
	expect_error(niches_custom(nsp = 3, nr = 2, location = location, 
		breadth = c(br_3_2_l, br_3_2_l)), regex = "niche breadth")
	expect_error(niches_custom(nsp = 3, nr = 2, location = location, breadth = br_3_2_l,
		r_use = r_use[1:2, ]), regex = "r_use")
})

test_that("Static niche dimensions", {
	location = matrix(c(1, 2), ncol = 2)
	breadth = c(1, 1.5)
	r_use =  c(0.5, 0.4)
	r_lim = rbind(c(0, 1), c(0, 1))
	expect_error(ni <- niches_custom(nsp = 1, nr = 2, location = location, breadth = breadth,
		r_use = r_use, r_lim = r_lim, static = 1), regex = NA)
	expect_equal(ni$r_use[[1]][1], 0)
	expect_equal(ni$r_use[[1]][2], r_use[2])

	expect_error(niches_custom(nsp = 1, nr = 2, location = location, breadth = breadth,
		r_use = r_use, r_lim = r_lim, static = -1), regex = "Indices in static")
	expect_error(niches_custom(nsp = 1, nr = 2, location = location, breadth = breadth,
		r_use = r_use, r_lim = r_lim, static = 3), regex = "Indices in static")

})

test_that("Niche scenarios", {
	expect_error(niches_uniform(nsp = 2, nr = 2), regex = "multiple resources")

	# few species
	expect_error(nu <- niches_uniform(nsp = 2, nr = 1), regex = NA)
	expect_equal(unlist(nu$location), c(1/3, 2/3))

	# many species
	expect_error(nu <- niches_uniform(nsp = 10, nr = 1), regex = NA)
	expect_equal(unlist(nu$location), seq(0, 1, length.out = 10))

	# custom niche width
	br = 2.4
	expect_error(nu <- niches_uniform(nsp = 2, nr = 1, breadth = br), regex = NA)
	expect_equal(nu$breadth, rep(br, 2))

	# random niches, single resource
	nsp = 2
	nr = 1
	expect_error(ni_r <- niches_random(nsp = nsp, nr = nr), regex = NA)
	expect_equal(length(unlist(ni_r$location)), nsp * nr)

	# random niches, multiple resources
	nr = 2
	expect_error(ni_r <- niches_random(nsp = nsp, nr = nr), regex = NA)
	expect_equal(length(unlist(ni_r$location)), nsp * nr)

	# random niches, custom breadth
	br = 1
	expect_error(ni_r <- niches_random(nsp = nsp, nr = nr, breadth = br), regex = NA)
	br_o = unlist(lapply(ni_r$breadth, diag))
	expect_equal(br_o, rep(br, length(br_o)))

})

test_that("Dispersal scenarios", {
	alpha = c(0.05, 0.07)
	beta = c(0.4, 0.3)

	expect_error(dispersal_custom(nsp = 1, alpha = alpha[1], beta = beta[1]), regex = NA)
	expect_error(dispersal_custom(nsp = 2, alpha = alpha, beta = beta), regex = NA)
	expect_error(dispersal_custom(nsp = 1, alpha = alpha, beta = beta[1]), regex = "alpha")
	expect_error(dispersal_custom(nsp = 1, alpha = alpha[1], beta = beta), regex = "beta")

})

test_that("Community scenarios", {
	Q = rep(1, 4)
	adj = matrix(0, nrow = 4, ncol = 4)
	adj[1,2] = adj[2,3] = adj[4,3] = 1
	rn = river_network(adj, Q)
	state(rn, "resources") = matrix(seq(0, 1, length.out = 4), ncol=1)
	mc = metacommunity()

	# make sure all sites have species in random communities
	expect_error(rcomm <- community_random(rn, mc, prevalence = c(0.5, 0.25)), 
		regex = NA)
	expect_warning(rcomm <- community_random(rn, mc, prevalence = 0), 
		regex = "low prevalence")

	expect_error(rcomm <- community_random(rn, mc, prevalence = c(0.5, 0.25), 
		allow_empty_sites = FALSE), regex = NA)
	expect_true(all(rowSums(rcomm) > 0))

	expect_error(ecomm <- community_equilibrium(rn, mc), regex = NA)
	expect_identical(ecomm == 1, f_niche(mc, state(rn, "resources")) > 0)
})
