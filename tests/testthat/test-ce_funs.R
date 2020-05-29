# library(flume)
#
# test_that("Colonisation and extinction functions", {
#
# 	# input validation
# 	expect_error(generate_colfun(type = 'lienar', n_species = 3), regex = "linear")
# 	expect_error(generate_extfun(type = 'lienar', n_species = 3), regex = "linear")
# 	expect_error(generate_colfun(type = 'lienar', nx = 2), regex = "linear")
# 	expect_error(generate_extfun(type = 'lienar', nx = 3), regex = "linear")
#
# 	R = c(0,0.5,1)
#
# 	# linear function, two species
# 	expect_error(col <- generate_colfun(type='linear', n_species=2), regex = NA)
# 	expect_error(ext <- generate_extfun(type='linear', n_species=2), regex = NA)
# 	# basic output
# 	expect_equal(length(col), 2)
# 	expect_true(all(sapply(col, is.function)))
# 	expect_true(all(sapply(ext, is.function)))
# 	expect_true(all(is.function(col)) && is.function(ext))
#
# 	expect_equal(col(R), -1*(ext(R) - max(ext(R))))
#
# 	# constant functions
#
# 	# gaussian functions
# })
