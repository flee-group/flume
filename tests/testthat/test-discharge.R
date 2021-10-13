rn = readRDS(system.file("testdata/river_network.rds", package="flume"))

test_that("discharge", {
	## reading works as expected
	expect_error(Q <- discharge(rn), regex=NA)

	## input validation, invalid inputs rejected
	expect_error(discharge(rn) <- -1*Q, "Negative")
	expect_error(discharge(rn) <- Q[1], "nrow")
	expect_error(discharge(rn) <- list(Q), "Unknown")

	## valid input accepted, attributes are set correctly, returns the right value
	expect_error(discharge(rn) <- function(tm) return(Q), regex = NA)
	expect_match(attr(rn, "discharge_model"), "function")
	expect_identical(discharge(rn), Q)

	zeros = rep(0, length(Q))
	expect_error(discharge(rn) <- cbind(Q, zeros), regex = NA)
	expect_match(attr(rn, "discharge_model"), "variable")
	expect_identical(discharge(rn), Q)
	state(rn, "resources") = state(rn, "resources") ## advance the time step of the model
	expect_identical(discharge(rn), zeros)

	expect_error(discharge(rn) <- Q, regex = NA)
	expect_match(attr(rn, "discharge_model"), "constant")
	expect_identical(discharge(rn), Q)


	##
	## valid input accepted
	## read returns the right value
})

test_that("cross sectional area", {
	Q = discharge(rn)
	zeros = rep(0, length(Q))

	## reading works as expected
	expect_error(A <- cs_area(rn), regex=NA)
	expect_error(A2 <- geometry(discharge(rn)), regex=NA)
	A2 = A2$width * A2$depth
	expect_identical(A, A2)
	A[1] = 2*A[1] # change it just to make sure the right vals are returned

	## input validation, invalid inputs rejected, valid inputs accepted for different models
	## constant
	expect_error(cs_area(rn) <- -1*A, regex="Negative")
	expect_error(cs_area(rn) <- cbind(A, A), regex="1-column")
	expect_error(cs_area(rn) <- function(tm) return(A), regex="1-column")
	expect_error(cs_area(rn) <- A, regex=NA)
	expect_identical(A, cs_area(rn))

	## variable
	discharge(rn) = cbind(Q, zeros)
	expect_error(cs_area(rn) <- A, regex = "dimensions")
	expect_error(cs_area(rn) <- cbind(A, zeros), regex = NA)
	expect_identical(cs_area(rn), A)
	state(rn, "resources") = state(rn, "resources") ## advance the time step of the model
	expect_identical(cs_area(rn), zeros)

	## function
	discharge(rn) = function(tm) return(Q)
	expect_error(cs_area(rn) <- A, regex="function")
	expect_error(cs_area(rn) <- function(tm) return(A), regex=NA)
	expect_identical(A, cs_area(rn))

	## resetting
	expect_error(cs_area(rn) <- NULL, regex = NA)
	A = cs_area(rn)
	expect_identical(A, A2)
})
