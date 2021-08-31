library(watershed)
library(raster)
library(sf)
library(ggplot2)
library(flume)

dem = raster("data/elevation/dhm_at_lamb_10m_2018_epsg3035.tif")
ext = extent(4652000, 4722000, 2740000, 2807000)
dem = crop(dem, ext)
thresh = 2e5
len = 2000
yb = delineate(dem, threshold = thresh, reach_len=len, outlet = c(4690462, 2756217))
Tp = pixel_topology(yb)

## for plotting
dem_cr = crop(dem, yb)
dem_df = as.data.frame(rasterToPoints(dem_cr))
colnames(dem_df)[3] = "elev"
yb_v = vectorise_stream(yb[["stream"]], Tp=Tp)

discharge = read.csv("data/discharge.txt")
discharge = discharge[discharge$watershed == "Ybbs",]
discharge = st_as_sf(discharge, coords = c('x', 'y'), crs=3035)

discharge = st_snap(discharge, yb_v, tolerance=100)
discharge$ca = catchment(yb, type = "points", y = st_coordinates(discharge))
discharge = discharge[discharge$ca > thresh,]

ggplot() + geom_raster(data=dem_df, aes(x = x, y = y, fill = elev)) + 
	geom_sf(data = yb_v, aes(colour=reach_id)) + scico::scale_fill_scico(palette = "turku") + 
	geom_sf(data = discharge, colour='red', size = 0.75)


yb_v$catchment = catchment(yb, type = 'reach', Tp = Tp)
yb_v = cbind(yb_v, hydraulic_geometry(discharge$ca, discharge$discharge, yb_v$catchment))
Tp_r = reach_topology(yb, Tp)
layout = st_coordinates(lwgeom::st_startpoint(yb_v))

yb_v$len = st_length(yb_v)
st_write(yb_v, "res/ybbs.shp", append = FALSE)



# combine reaches
simplify_reach = function(tp, i) {
	hw = i[watershed:::.headwater(tp[i,i])]
	if(length(hw) > 1)
		hw = hw[which.max(hw)]
	out = i[watershed:::.outlet(tp[i,i])]
	conn = watershed:::.upstream(tp, hw)

	# find branches
	all_upstream = function(tp, j) {
		val = unlist(sapply(j, function(k) unname(watershed:::.upstream(tp, k))))
		if(length(val) == 0) {
			return(val)
		} else {
			return(c(val, all_upstream(tp, val)))
		}
	}
	br = unname(unlist(sapply(i, function(j) watershed:::.upstream(tp, j))))
	br = br[which(!br %in% c(i, conn))]
	if(length(br) > 0) {
		br = c(br, all_upstream(tp, br))
		tp[br,] = 0
		tp[,br] = 0
	}

	lens = sapply(conn, function(j) {
		k = watershed:::.downstream(tp, j)
		l = tp[j, k]
		while(k != out) {
			ko = k
			k = watershed:::.downstream(tp, ko)
			l = l + tp[ko, k]
		}
		l
	})

	del = c(i[i != out], br)
	if(length(conn) == 0) {
		connect = data.frame()
	} else {
		connect = data.frame(ds = out, us = conn, length = lens)
	}
	

	list(connect = connect, delete = del)
}

simplify_reaches = function(tp, rs) {
	form = lapply(rs, function(r) simplify_reach(tp, r))
	for(ff in form) {
		if(nrow(ff$connect) > 0) {
			for(i in 1:nrow(ff$connect))
				tp[ff$connect$us[i], ff$connect$ds[i]] = ff$connect$length[i]
		}
	}
	del = unlist(sapply(form, function(ff) ff$delete))
	if(length(del) > 0) {
		tp = tp[-del,]
		tp = tp[,-del]
	}
	list(topology = tp, delete = del)
}

rcomb = list(
	1:3, 5:6, 338:339, 345:346, 371:374, 376:378, 386:387, 347:358, 349:350, 351:352, 368:370, 
	353:356, 361:362, 363:364, 340:343, 389:390, 391:392, 393:394, 10:13, 15:18, 19:20, 21:22, 
	324:325, 327:331, 23:25, 26:28, 29:33, 82:84, 85:87, 35:37, 40:42, 43:44, 45:47, 51:53, 55:57, 
	59:61, 66:70, 77:80, 71:76, 261:279, 280:282, 284:286, 288:291, 292:294, 295:297, 303:304, 
	305:307, 309:313, 92:93, 95:96, 97:99, 101:102, 103:105, 106:108, 109:111, 256:259, 115:116, 
	117:119, 120:121, 123:125, 137:140, 141:146, 147:148, 150:155, 157:159, 161:163, 164:167, 
	249:252, 168:173, 174:177, 178:179, 182:184, 245:246, 185:188, 212:213, 189:190, 192:195, 
	197:199, 201:203, 205:207, 208:210, 215:222, 223:224, 237:238, 240:242, 226:227, 228:229, 
	232:234)

res = simplify_reaches(Tp_r, rcomb)

layout = st_coordinates(lwgeom::st_startpoint(yb_v))
yb_v$Q = units::drop_units(yb_v$Q)
yb_v$z = units::drop_units(yb_v$z)
yb_v$w = units::drop_units(yb_v$w)

ybbs = river_network(res$topology, discharge = yb_v$Q[-res$delete], 
	area = (yb_v$z*yb_v$w)[-res$delete], length = len, 
	layout = layout[-res$delete,])
ybbs_old = river_network(Tp_r, discharge = yb_v$Q, 
	area = (yb_v$z*yb_v$w), length = len, 
	layout = layout)

par(mfrow=c(1,2))
	plot(ybbs, edge.arrow.size=0, vertex.size=2, vertex.label=NA)
	plot(ybbs_old, edge.arrow.size=0, vertex.size=2, vertex.label=NA)

saveRDS(ybbs, "res/ybbs_network.rds")