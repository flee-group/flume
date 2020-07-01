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

	expect_error(comm <<- metacommunity(), regex = NA)
	pool = comm$species
	expect_equal(length(pool), 2)
	expect_equal(pool[[1]]$ext(1), pool[[2]]$ext(1))
	R <- matrix(seq(0,1,length.out=20), ncol=1)
	expect_equal(mean(pool[[1]]$col(R)), mean(pool[[2]]$col(R)))
	default_mat = matrix(c(0.30, 0.05, 0.05, 0.30), ncol=2)
	expect_equal(comm$competition, default_mat)
})

test_that("Niches estimated correctly", {
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
