devtools::load_all()
## make a species pool
sp_pool = metacommunity(n_species = 2, nx = 1, c_type = 'linear', e_type='constant')
sp_pool$r_scale = matrix(c(0.3, 0.1), nrow = 2, ncol=1)
saveRDS(sp_pool, "inst/testdata/metacom.rds")


## make a river network
##

Q = c(13,12,8,7,1,4,1,1,2,1)
nsites = length(Q)
adj = matrix(0, nrow=nsites, ncol=nsites)
adj[2,1] = adj[3,2] = adj[9,2] = adj[4,3] = adj[5,4] = adj[6,4] = adj[7,6] = adj[8,6] = adj[10,9] = 1

network = river_network(adjacency = adj, discharge = Q)
layout = matrix(c(0,0,0,1,-0.5,2,-0.5,3,-1,4,0,4,-0.5,5,0.5,5,0.5,2,1,3), byrow=TRUE, nrow=nsites)

R = matrix(seq(0, 1, length.out = nsites), ncol=1, dimnames = list(NULL, 'R'))
network$boundary = function() return(matrix(seq(0, 1, length.out = 10), ncol=1, dimnames = list(NULL, 'R')))
state(network) = R
attr(network, "layout") = layout

sbysp = matrix(0, ncol = 2, nrow=length(Q))
sbysp[c(1,3,5),1] = sbysp[c(2:4,6:10), 2] = 1
site_by_species(network) = sbysp

saveRDS(network, "inst/testdata/river_network.rds")

sim = flume(sp_pool, network, 1)
set.seed(45580085)
sim = run_simulation(sim, 50)
saveRDS(sim, "inst/testdata/sim.rds")
