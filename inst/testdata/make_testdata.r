devtools::load_all()
## make a species pool
mc = metacommunity(nsp = 2, nr = 1, niche_args = list(r_use = c(0.3, 0.1)))
saveRDS(mc, "inst/testdata/metacom.rds")


## make a river network
##
Q = c(13,12,8,7,1,4,1,1,2,1)
nsites = length(Q)
adj = matrix(0, nrow=nsites, ncol=nsites)
adj[2,1] = adj[3,2] = adj[9,2] = adj[4,3] = adj[5,4] = adj[6,4] = adj[7,6] = adj[8,6] = 
	adj[10,9] = 1

layout = matrix(c(0,0,0,1,-0.5,2,-0.5,3,-1,4,0,4,-0.5,5,0.5,5,0.5,2,1,3), 
	byrow=TRUE, nrow=nsites)
R = matrix(seq(0, 1, length.out = nsites), ncol=1, dimnames = list(NULL, 'R'))
network = river_network(adjacency = adj, discharge = Q, layout = layout)
state(network, "resources") = R
boundary(network, "resources") = R
state(network, "species") = community_equilibrium(network, mc)
saveRDS(network, "inst/testdata/river_network.rds")

sim = flume(mc, network)
set.seed(45580085)
sim = run_simulation(sim, 50)
saveRDS(sim, "inst/testdata/sim.rds")
