test_that("River network creation", {

	Q = rep(1, 4)
	# Invalid topologies
	adj = matrix(0, nrow = 4, ncol = 4)
	adj[1,2] = adj[2,3] = adj[4,3] = 1

	# not a matrix
	expect_error(river_network(as.vector(adj), Q), regex = "validation")

	# invalid weights
	adj[2,2] = -1
	expect_error(river_network(adj, Q), regex = "validation")
	adj[2,2] = Inf
	expect_error(river_network(adj, Q), regex = "validation")
	adj[2,2] = 0

	# cycles and loops
	adj[3,1] = 1
	expect_error(river_network(adj, Q), regex = "validation")
	adj[3,1] = 0
	adj[1,1] = 1
	expect_error(river_network(adj, Q), regex = "validation")
	adj[1,1] = 0

	# single node upstream of two downstream nodes
	adj[1,4] = 1
	expect_error(river_network(adj, Q), regex = "validation")
	adj[1,4] = 0

	# isolated nodes
	adj[4,3] = 0
	expect_error(river_network(adj, Q), regex = "validation")
	adj[4,3] = 1

	# invalid discharge
	expect_error(river_network(adj, Q[1:3]), regex = "length\\(discharge\\)")

	# works with regular and Matrix matrices
	expect_error(river_network(adj, Q), regex = NA)
	expect_error(river_network(Matrix::Matrix(adj), Q), regex = NA)
})
