test_that("Species creation", {

	# input validation
	expect_error(species('line', 'constant', NULL, NULL), regex = "line")
	expect_error(species('linear', 'constant', list(scale = 0.2), NULL), regex = "a")
	expect_error(sp <- species('linear', 'gaussian', list(a = 0, b = 1), list(scale=0.2)), regex="gaussian")

	expect_error(sp <- species('linear', 'constant', list(a = 0, b = 1), list(scale=0.2)), regex=NA)
	expect_equal(sp$col(1), 1)
	expect_equal(sp$ext(1), 0.2)

	expect_error(sp <- species('constant', 'constant', list(scale=0.3), list(scale=0.2)), regex=NA)
	expect_equal(sp$col(1), 0.3)

	expect_error(sp <- species('gaussian', 'linear', list(scale=1, mean=0.5, sd=0.2), list(a=1, b=-1)), regex=NA)
	expect_equal(sp$col(0.5), 1)
	expect_equal(sp$ext(1), 0)

	# default, 2 species metacommunity, linear niches
	expect_error(comm <<- metacommunity(c_type="linear"), regex = NA)
	pool = comm$species
	expect_equal(length(pool), 2)
	expect_equal(pool[[1]]$ext(1), pool[[2]]$ext(1))
	R <- matrix(seq(0,1,length.out=20), ncol=1)
	expect_equal(mean(pool[[1]]$col(R)), mean(pool[[2]]$col(R)))
	default_mat = matrix(c(0.30, 0.05, 0.05, 0.30), ncol=2)
	expect_equal(comm$competition, default_mat)

	# gaussian niches, 2 species
	expect_error(comm <- metacommunity(c_type = 'gaussian'), regex = NA)

	# gaussian niches, 3 species
	expect_error(comm <- metacommunity(n_species = 3, c_type = 'gaussian'), regex = NA)

	# gaussian niches, 3 species, 2 variables
	expect_error(comm <- metacommunity(n_species = 3, nx = 2, c_type = 'gaussian'), regex = NA)

	# check that dispersal function works
	expect_error(di <- dispersal_params(comm), regex=NA)
	expect_true(is(di, "list"))
	expect_equal(length(di$alpha), length(di$beta))
	expect_equal(length(di$alpha), length(comm$species))
})

test_that("Niches estimated correctly", {
	comm <- metacommunity(c_type = 'gaussian')
	st = matrix(seq(0, 1, length.out = 4), ncol = 1, dimnames = list(NULL, 'R'))
	expect_error(sp_niche <- f_niche(comm$species[[1]], st), regex=NA)
	expect_true(is(sp_niche, "matrix"))
	expect_equal(dim(sp_niche), c(nrow(st), 1))
	expect_identical(sp_niche, comm$species[[1]]$col(st) - comm$species[[1]]$ext(st))

	expect_error(c_niche <- f_niche(comm, st), regex=NA)
	expect_true(is(c_niche, "matrix"))
	expect_equal(dim(c_niche), c(nrow(st), length(comm$species)))
	expect_identical(c_niche[,1, drop=FALSE], sp_niche)

})

test_that("Random community creation", {
	Q = rep(1, 4)
	adj = matrix(0, nrow = 4, ncol = 4)
	adj[1,2] = adj[2,3] = adj[4,3] = 1
	rn = river_network(adj, Q)

	expect_warning(random_community(rn, comm, 0.01), regex = "low prevalence")
	expect_error(si_sp <- random_community(rn, comm), regex = NA)
	expect_equal(ncol(si_sp), length(comm$species))
	expect_equal(nrow(si_sp), nrow(adj))
})
