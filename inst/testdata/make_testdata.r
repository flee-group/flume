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
network$boundary = function() matrix(seq(0, 1, length.out = nsites), ncol=1, dimnames = list(NULL, 'R'))
attr(network, "layout") = layout
state(network) = network$boundary()
site_by_species(network) = random_community(network, sp_pool, prevalence = c(0.35, 0.65))

saveRDS(network, "inst/testdata/river_network.rds")
