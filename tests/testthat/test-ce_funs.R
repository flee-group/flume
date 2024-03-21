test_that("Colonisation and extinction functions", {

	# input validation
	expect_error(ce_linear(parm = list(a = 0.5)), regex = "linear")
	expect_error(ce_linear(parm = list(b = 0.5)), regex = "linear")

	# linear function
	expect_error(cfun <- ce_linear(parm = list(a = 0, b = 1)), regex = NA)

	R = c(0, 1)
	expect_error(res <- cfun(R), regex = NA)

	# basic output
	expect_equal(length(res), length(R))
	expect_identical(res, c(0, 1))

	# constant functions
	sc = 0.2
	expect_error(cfun <- ce_constant(scale = sc, nr = 1), regex = NA)
	expect_error(res <- cfun(R), regex = NA)
	expect_identical(res, rep(sc, length(R)))

	# gaussian functions
	expect_error(ce_gaussian(1, 0.2, 0.2), regex = NA)
	
	# multivariate, invalid inputs
	expect_error(ce_gaussian(c(1, 2), 0.2, 0.2), regex = 'variance-covariance matrix')
	expect_error(ce_gaussian(c(1, 2), c(0.1, 0.2), 0.2), regex = 'variance-covariance matrix')
	vcv = matrix(c(0.2, 0, 0, 0.2), ncol = 2)
	expect_error(ce_gaussian(c(1, 2), vcv, c(0.1, 0.2)), regex = 'single value')
	
	# valid input
	expect_error(cfun <- ce_gaussian(c(1, 2), vcv, 0.2), regex = NA)
	
	# negative scales
	sc = 0.2
	loc = 1
	expect_warning(cfun <- ce_constant(scale = -sc, nr = 1), regex = "negative")
	expect_error(cfun <- ce_gaussian(loc, 0.2, sc), regex = NA)
	expect_error(efun <- ce_gaussian(loc, 0.2, -sc), regex = NA)
	expect_equal(efun(Inf), 2 * sc)
	expect_equal(efun(loc), sc)

	loc = c(1,2)
	expect_error(efun <- ce_gaussian(loc, vcv, -sc), regex = NA)
	expect_equal(efun(1e6*loc), 2 * sc) # Inf doesn't work, just try big num
	expect_equal(efun(loc), sc)
	
})
