library(watershed)
library(raster)
library(sf)
library(ggplot2)
library(flume)

dem = raster("data/elevation/mara_dem_102022.tif")
thresh = 1e7
len = 10000
ma = delineate(dem, threshold = thresh, reach_len=len, outlet = c(1040958, -154588))
Tp = pixel_topology(ma)

## for plotting
dem_cr = crop(dem, ma)
dem_df = as.data.frame(rasterToPoints(dem_cr))
colnames(dem_df)[3] = "elev"
ma_v = vectorise_stream(ma[["stream"]], Tp=Tp)

discharge = readRDS("data/mara_sites_final.rds")
discharge = st_snap(discharge, ma_v, tolerance=100)
discharge$ca = catchment(ma, type = "points", y = st_coordinates(discharge))
discharge = discharge[discharge$ca > thresh,]

ggplot() + geom_raster(data=dem_df, aes(x = x, y = y, fill = elev)) + 
	geom_sf(data = ma_v, aes(colour=reach_id)) + scico::scale_fill_scico(palette = "turku") + 
	geom_sf(data = discharge, colour='red', size = 0.75)


ma_v$catchment = catchment(ma, type = 'reach', Tp = Tp)
ma_v = cbind(ma_v, hydraulic_geometry(discharge$ca, discharge$q, ma_v$catchment))
Tp_r = reach_topology(ma, Tp)
layout = st_coordinates(lwgeom::st_startpoint(ma_v))
ma_v$len = st_length(ma_v)
st_write(ma_v, "res/mara.gpkg", append = FALSE)


ma_v$Q = units::drop_units(ma_v$Q)
ma_v$z = units::drop_units(ma_v$z)
ma_v$w = units::drop_units(ma_v$w)
mara = river_network(Tp_r, discharge = ma_v$Q, 
	area = (ma_v$z*ma_v$w), length = len, 
	layout = layout)

	plot(mara, edge.arrow.size=0, vertex.size=2, vertex.label=NA)

saveRDS(mara, "res/mara_network.rds")