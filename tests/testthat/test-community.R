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

	expect_error(comm <- metacommunity(), regex = NA)
	pool = comm$species
	expect_equal(length(pool), 2)
	expect_equal(pool[[1]]$ext(1), pool[[2]]$ext(1))
	R <- matrix(seq(0,1,length.out=20), ncol=1)
	expect_equal(mean(pool[[1]]$col(R)), mean(pool[[2]]$col(R)))
	default_mat = matrix(c(0.30, 0.05, 0.05, 0.30), ncol=2)
	expect_equal(comm$competition, default_mat)
})
