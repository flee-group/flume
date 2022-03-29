library(flume)
library(Matrix)
kamp = readRDS("data-raw/networks/kamp_old.rds")
# plot(kamp, vertex.size = 5)

drop = c(2,3,5,7, 63:65, 62, 61, 12, 30:32)
adj = kamp[[".adjacency"]][-drop, -drop]
names = 1:nrow(adj)
colnames(adj) = rownames(adj) = names

kamp_new = river_network(adj, kamp[[".Q"]][-drop], kamp[[".area"]][-drop], kamp[[".length"]], 
	attr(kamp, "layout")[-drop, ], site_names = names) 
plot(kamp_new, vertex.size = 5)
saveRDS("kamp_new", "data_raw/networks/kamp.rds")

mara = readRDS("data-raw/networks/mara_network.rds")
thur = readRDS("data-raw/networks/thur_network.rds")
vjosa = readRDS("data-raw/networks/vjosa_network.rds")
ybbs = readRDS("data-raw/networks/ybbs_network.rds")
flume_networks = list(kamp = kamp_new, mara = mara, thur = thur, vjosa = vjosa, ybbs = ybbs)
usethis::use_data(flume_networks, overwrite = TRUE)
