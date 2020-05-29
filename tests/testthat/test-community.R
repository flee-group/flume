test_that("Species creation", {

	# input validation
	expect_error(create_species('line', 'constant', NULL, NULL), regex = "line")
	expect_error(create_species('linear', 'constant', list(scale = 0.2), NULL), regex = "a")
	expect_error(sp <- create_species('linear', 'gaussian', list(a = 0, b = 1), list(scale=0.2)), regex="gaussian")

	expect_error(sp <- create_species('linear', 'constant', list(a = 0, b = 1), list(scale=0.2)), regex=NA)
	expect_equal(sp$col(1), 1)
	expect_equal(sp$ext(1), 0.2)

	expect_error(sp <- create_species('constant', 'constant', list(scale=0.3), list(scale=0.2)), regex=NA)
	expect_equal(sp$col(1), 0.3)

	expect_error(sp <- create_species('gaussian', 'linear', list(scale=1, mean=0.5, sd=0.2), list(a=1, b=-1)), regex=NA)
	expect_equal(sp$col(0.5), 1)
	expect_equal(sp$ext(1), 0)
})
