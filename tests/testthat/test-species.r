test_that("Species creation", {
	col = new("GaussianNicheFun", location = 0, breadth = matrix(1), scale = 1, nr = 1L)
	ext = new("ConstantNicheFun", scale = 1, nr = 1L)
	
	# input validation
	# many parameters must not be negative
	col1 = col
	col1@breadth = -1 * col1@breadth
	expect_error(species(col1, ext, alpha = 1, beta = 1, r_use = 1), regex = "negative")

	col1 = col
	col1@scale = -1 * col1@scale
	expect_error(species(col1, ext, alpha = 1, beta = 1, r_use = 1), regex = "negative")

	ext1 = ext
	ext1@scale = -1 * ext1@scale
	expect_error(species(col, ext1, alpha = 1, beta = 1, r_use = 1), regex = "negative")

	expect_error(species(col, ext, alpha = -1, beta = 1, r_use = 1), regex = "negative")
	expect_error(species(col, ext, alpha = 1, beta = -11, r_use = 1), regex = "negative")
	
	# some parameters may be negative or positive
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = 1), regex = NA)
	col1 = col
	col1@location = -1 * col1@location
	expect_error(species(col1, ext, alpha = 1, beta = 1, r_use = 1), regex = NA)
	

	# multivariate niches
	loc = c(1, 2)
	rsc = c(1, 2)
	rsc_wrong = 1:3
	bre = matrix(c(1, 0, 0, 1), nrow = 2)
	bre_wrong = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)

	# all dimensions must agree
	col@location = loc
	col@breadth = bre_wrong
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = rsc), regex = "dimension mismatch")
	col@breadth = bre
	col@location = 1
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = rsc), regex = "dimension mismatch")
	col@location = loc
	col@breadth = bre[1,, drop = FALSE]
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = rsc_wrong), regex = "dimension mismatch")
	
	# also check that no entries in vcv matrix are negative
	col@breadth = -1*bre
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = rsc), regex = "negative")

	col@breadth = bre
	# scale_c, scale_e, alpha and beta must be single values
	col@scale = rsc
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = rsc), regex = "single value")
	col@scale = 1
	ext@scale = rsc
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = rsc), regex = "single value")

	ext@scale = 1
	expect_error(species(col, ext, alpha = c(1, 0), beta = 1, r_use = rsc), regex = "single value")
	expect_error(species(col, ext, alpha = 1, beta = c(1, 0), r_use = rsc), regex = "single value")
	
	# valid input
	expect_error(species(col, ext, alpha = 1, beta = 1, r_use = rsc), regex = NA)
})

