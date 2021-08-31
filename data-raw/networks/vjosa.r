library(watershed)
library(raster)
library(sf)
library(ggplot2)
library(flume)

# dem1 = raster("data/elevation/eu_dem_v11_E50N10.TIF")
# dem2 = raster("data/elevation/eu_dem_v11_E50N20.TIF")
# dem = merge(dem1, dem2, file = )
# ext = extent(5073000, 5350000, 1874000, 2035000)
# dem = crop(dem, ext, file="data/elevation/vjosa_dem.tif")
thresh = 5e7
dem = raster("data/elevation/vjosa_dem.tif")
vj = delineate(dem, threshold = thresh, reach_len=7000)
Tp = pixel_topology(vj)

## for plotting
# dem_cr = crop(dem, vj)
# dem_df = as.data.frame(rasterToPoints(dem_cr))

# remove small headwater reaches
vj_v = vectorise_stream(vj[["stream"]], Tp=Tp)
# Tr = reach_topology(vj, Tp=Tp)
# hws = which(colSums(Tr != 0) == 0)
# lens = st_length(vj_v)
# hw_drop = hws[which(lens[hws] < units::set_units(3500, "m"))]
# vj[["stream"]][values(vj[["stream"]]) %in% hw_drop] = NA
# Tp = pixel_topology(vj)

# # renumber the reaches
# rids = unique(values(vj$stream))
# rids = sort(rids) ## also gets rid of NA
# vj$stream = match(vj$stream, rids)
# # make into vector to do further revisions in QGIS
# vj_v = vectorise_stream(vj[["stream"]], Tp=Tp)
# vj_v$len = st_length(vj_v)
# st_write(vj_v, "data/vjosa_stream.shp", append = FALSE)


# ggplot() + geom_raster(data=dem_df, aes(x = x, y = y, fill = vjosa_dem)) + 
# 	geom_sf(data = vj_v, aes(colour=catchment)) + scico::scale_fill_scico(palette = "turku")

discharge = read.csv("data/discharge.txt")
discharge = discharge[discharge$watershed == "Vjosa" & discharge$season == "Fall",]
discharge = st_as_sf(discharge, coords = c('x', 'y'), crs=3035)
# ggplot() + geom_raster(data=dem_df, aes(x = x, y = y, fill = vjosa_dem)) + 
	# geom_sf(data = vj_v, aes(colour=catchment)) + scico::scale_fill_scico(palette = "turku") +
	# geom_sf(data = discharge)

discharge = st_snap(discharge, vj_v, tolerance=100)
discharge$ca = catchment(vj, type = "points", y = st_coordinates(discharge))
discharge = discharge[discharge$ca > thresh,]

vj_v$catchment = catchment(vj, type = 'reach', Tp = Tp)
vj_v = cbind(vj_v, hydraulic_geometry(discharge$ca, discharge$discharge, vj_v$catchment))
Tp_r = reach_topology(vj, Tp)



# save intermediate shit in case of a crash
vj_v$CA = units::set_units(vj_v$CA, "km^2")
saveRDS(vj_v, "res/vjosa_stream.rds")
saveRDS(discharge, "res/vjosa_discharge.rds")
writeRaster(vj, "res/vjosa_watershed.grd", overwrite = TRUE)
saveRDS(Tp, "res/vjosa_pix_topo.rds")
saveRDS(Tp_r, "res/vjosa_each_topo.rds")
layout = st_coordinates(lwgeom::st_startpoint(vj_v))
vj_v$Q = units::drop_units(vj_v$Q)
vj_v$z = units::drop_units(vj_v$z)
vj_v$w = units::drop_units(vj_v$w)


# one last attempt to edit the topology manually to get rid of tiny reaches
# Tp_r_bak = Tp_r
# vj_v_bak = vj_v

vj_v = vj_v_bak
Tp_r = Tp_r_bak
layout = st_coordinates(lwgeom::st_startpoint(vj_v))

rdel = list(
	c(2:19,308),
	21:33,
	c(35:66,306:307),
	c(35:66,306:307),
	c(81:89,305),
	91:92,
	c(94:101,303:304),
	c(94:101,303:304),
	c(104:127, 302),
	c(129:139, 299:301),
	141:153,
	141:153,
	c(157:174, 298, 155),
	c(157:174, 298, 155),
	178:187,
	178:187,
	191,
	191,
	273:277,
	273:277,
	287,
	c(193:202, 271),
	c(193:202, 271),
	206,
	206,
	c(210:212, 214)
)

rconn = list(
	c(20, 1),
	c(34,20),
	c(67,34),
	c(80,34),
	c(90,80),
	c(93,90),
	c(102,93),
	c(103,93),
	c(128,103),
	c(140,128),
	c(154,140),
	c(156,140),
	c(176,156),
	c(177,156),
	c(189,177),
	c(190,177),
	c(272,190),
	c(192,190),
	c(279,272),
	c(282,272),
	c(286,283),
	c(203,192),
	c(270,192),
	c(208,207),
	c(209,207),
	c(213,209))

lens = rowSums(Tp_r)
for(i in 1:length(rdel)) {
	dd = rdel[[i]]
	cc = rconn[[i]]
	Tp_r[cc[1], cc[2]] = sum(lens[dd]) + lens[cc[1]]
}

rdel = unique(unlist(rdel))
Tp_r = Tp_r[-rdel, -rdel]
vj_v = vj_v[-rdel,]
layout = st_coordinates(lwgeom::st_startpoint(vj_v))


vjosa = river_network(Tp_r, discharge = vj_v$Q, area = vj_v$z*vj_v$w, length = 7000, 
	layout = layout)
plot(vjosa, edge.arrow.size=0, vertex.size=2, vertex.label=NA)

## notes
# there are many small reaches in the main stem. These need to be combined into simpler reaches
