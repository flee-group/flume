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
	expect_error(cfun <- ce_constant(scale = sc), regex = NA)
	expect_error(res <- cfun(R), regex = NA)
	expect_identical(res, rep(sc, length(R)))

	# gaussian functions
	expect_error(ce_gaussian(1, 0.2, 0.2), regex = NA)
	expect_error(ce_gaussian(c(1, 2), 0.2, 0.2), regex = NA)
	expect_error(ce_gaussian(c(1, 2), c(0.1, 0.2), 0.2), regex = NA)
	expect_error(ce_gaussian(c(1, 2), 0.2, c(0.1, 0.2)), regex = NA)
	expect_error(ce_gaussian(c(1, 2), c(0.1, 0.2), c(0.1, 0.2)), regex = NA)

	expect_error(cfun <- ce_gaussian(0, 1, 0.2), regex = NA)
	expect_equal(cfun(R[1]), 0.2)

})
