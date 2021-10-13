network = readRDS(system.file("testdata/river_network.rds", package="flume"))

test_that("simple discharge/area", {
	rn = network
	## reading works as expected
	expect_error(Q <- state(rn, "Q"), regex=NA)

	## input validation, invalid inputs rejected
	expect_error(discharge(rn) <- -1*Q, "Negative")
	expect_error(discharge(rn) <- Q[1], "nrow")
	expect_error(discharge(rn) <- list(Q), "Unknown")

	# reset discharge; setting it via state should fail, other way should work
	expect_error(state(rn, "Q") <- Q*2, regex = "discharge")
	expect_error(discharge(rn) <- Q*2, regex = NA)
	expect_equal(state(rn, "Q"), Q*2)


	### code to test area
	expect_error(A <- state(rn, "area"), regex=NA)
	expect_error(A2 <- geometry(state(rn, "Q")), regex=NA)
	A2 = A2$width * A2$depth
	expect_equal(A, A2)

	expect_error(state(rn, "area") <- A*2, regex = "area")
	expect_error(area(rn) <- -1*A, regex="Negative")
	expect_error(area(rn) <- cbind(A, A), regex="1-column")
	expect_error(area(rn) <- function(time) return(A), regex="1-column")
	expect_error(area(rn) <- A*2, regex=NA)
	expect_identical(A*2, state(rn, "area"))

	## setting null returns us to the default
	expect_error(area(rn) <- NULL, regex=NA)
	expect_true(is.null(area(rn)))
	expect_equal(state(rn, "area"), A2)
})



test_that("Variable discharge model", {
	rn = network
	Q = state(rn, "Q")

	expect_error(discharge(rn) <- cbind(Q, Q*2), regex = NA)
	expect_match(attr(rn, "discharge_model"), "variable")
	expect_true(is.matrix(discharge(rn)))
	expect_equal(state(rn, "Q"), Q)
	state(rn, "resources") = state(rn, "resources") ## advance the time step of the model
	expect_identical(state(rn, "Q"), Q*2)

	## cs_area
	# area() method returns a the raw value (null, it is not set), state the current vector
	expect_true(is.null(area(rn)))
	expect_error(a_curr <- state(rn, "area"), regex=NA)

	# setting with state fails
	expect_error(state(rn, "area") <- a_curr*2, regex = "area")

	# invalid inputs
	expect_error(area(rn) <- a_curr, regex="dimensions")
	expect_error(area(rn) <- cbind(a_curr, a_curr*2), regex=NA)
	expect_equal(cbind(a_curr, a_curr*2), area(rn))

	## make sure recycling happens properly
	state(rn, "resources") = state(rn, "resources") ## advance the time step of the model
	expect_equal(state(rn, "Q"), Q)
	state(rn, "resources") = state(rn, "resources") ## advance the time step of the model
	expect_equal(state(rn, "Q"), Q*2)
})



test_that("Discharge as function", {
	rn = network
	Q = state(rn, "Q")
	A = state(rn, "area")

	expect_error(discharge(rn) <- function(time) return(Q * time), regex = NA)
	expect_match(attr(rn, "discharge_model"), "function")

	# the discharge method should return the actual function so we can use it
	expect_true(is.function(discharge(rn)))

	# the state method returns the values
	# for first time step, equal to Q, second Q*2
	expect_equal(state(rn, "Q"), Q)
	expect_equal(state(rn, "area"), A)
	state(rn, "resources") = state(rn, "resources")
	expect_equal(state(rn, "Q"), Q*2)
	expect_true(all(state(rn, "area") > A))
	expect_error(qhist <- state(rn, "Q", history = TRUE), regex=NA)
	expect_equal(qhist, cbind(Q, Q*2), check.attributes = FALSE) 

	# if we set state to null, then we get time step 1
	state(rn, "resources") = NULL
	expect_equal(state(rn, "Q"), Q)
	expect_equal(state(rn, "area"), A)

})


