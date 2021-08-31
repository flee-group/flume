library(watershed)
library(raster)
library(sf)
library(ggplot2)
library(flume)

dem = raster("data/elevation/thur_dem_3035.tif")
thresh = 1e6
len = 5000
th = delineate(dem, threshold = thresh, reach_len=len, outlet = c(4260565, 2706691))
Tp = pixel_topology(th)

## for plotting
dem_cr = crop(dem, th)
dem_df = as.data.frame(rasterToPoints(dem_cr))
colnames(dem_df)[3] = "elev"
th_v = vectorise_stream(th[["stream"]], Tp=Tp)

discharge = read.csv("data/discharge.txt")
discharge = discharge[discharge$watershed == "Thur",]
discharge = st_as_sf(discharge, coords = c('x', 'y'), crs=3035)

discharge = st_snap(discharge, th_v, tolerance=100)
discharge$ca = catchment(th, type = "points", y = st_coordinates(discharge))
discharge = discharge[discharge$ca > thresh,]

ggplot() + geom_raster(data=dem_df, aes(x = x, y = y, fill = elev)) + 
	geom_sf(data = th_v, aes(colour=reach_id)) + scico::scale_fill_scico(palette = "turku") + 
	geom_sf(data = discharge, colour='red', size = 0.75)


th_v$catchment = catchment(th, type = 'reach', Tp = Tp)
th_v = cbind(th_v, hydraulic_geometry(discharge$ca, discharge$discharge, th_v$catchment))
Tp_r = reach_topology(th, Tp)
layout = st_coordinates(lwgeom::st_startpoint(th_v))
th_v$len = st_length(th_v)
st_write(th_v, "res/thur.gpkg", append = FALSE)



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

rcomb = list(2:4, 216:217, 184:186, 214:215, 188:190, 191:195, 209:213, 203:205, 206:207,
	196:197, 200:201, 5:8, 178:180, 9:17, 18:19, 20:22, 25:26, 27:30, 34:35, 147:149, 154:156,
	157:159, 160:162, 166:167, 168:170, 44:45, 145:146, 48:49, 50:53, 55:56, 57:58, 60:61, 63:66, 
	68:69, 141:143, 70:71, 73:76, 78:79, 80:81, 82:83, 84:86, 87:89, 95:98, 134:135, 104:105, 
	131:132, 106:111, 114:115, 119:120, 126:127, 122:124, 117:118)

res = simplify_reaches(Tp_r, rcomb)

th_v$Q = units::drop_units(th_v$Q)
th_v$z = units::drop_units(th_v$z)
th_v$w = units::drop_units(th_v$w)

thur = river_network(res$topology, discharge = th_v$Q[-res$delete], 
	area = (th_v$z*th_v$w)[-res$delete], length = len, 
	layout = layout[-res$delete,])
thur_old = river_network(Tp_r, discharge = th_v$Q, 
	area = (th_v$z*th_v$w), length = len, 
	layout = layout)

par(mfrow=c(1,2))
	plot(thur, edge.arrow.size=0, vertex.size=2, vertex.label=NA)
	plot(thur_old, edge.arrow.size=0, vertex.size=2, vertex.label=NA)

saveRDS(thur, "res/thur_network.rds")